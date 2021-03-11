# Functions 

library(rentrez)
set_entrez_key("e1b887b07de1764a6e68883fce0f9f69d108") # My API key
Sys.getenv("ENTREZ_KEY") 

"%ni%" <- Negate("%in%")

# ---------------------------------------------------------------------------------------------------------------- # 
# Function to assemble taxo and table data

# taxo_otu <- blank_taxo
# table_otu <- blank_table

assemble_data <- function(table_otu, taxo_otu){
  
  # Clean data - adjust to fit all - normalise by amplicon name - add ifelse for control for each category 
  taxo_otu$amplicon <- str_trim(as.character(taxo_otu$definition))
  taxo_otu$sequence <- toupper(taxo_otu$sequence)
  
  # Checks - for sequences: all sequences from the table must be in the csv. (The opposite is not true du to sequence loss while cleaning rapidrun data)
  verif_seq <- table_otu$sequence %in% taxo_otu$sequence
  if( length(verif_seq[verif_seq==FALSE]) > 0 ) stop(paste(" error: some sequences are in table but not in csv file"))
  
  # Checks - for sequences: all MOTUs from the table must be in the csv. (The opposite is not true du to sequence loss while cleaning rapidrun data)
  verif_motu <- table_otu$amplicon %in% taxo_otu$amplicon
  if( length(verif_motu[verif_motu==FALSE]) > 0  ) stop(paste(" error: some MOTU names are in table but not in csv file"))
  
  # Merge data
  sw <- merge(taxo_otu, table_otu, by=c("amplicon", "sequence"))
  
  # Extract all sample names to prepare formatting
  samples2 <- sw %>%
    dplyr::select(starts_with("CPCR"), starts_with("Cext"), starts_with("SPY"), starts_with("CMET"), starts_with("CNEG"), starts_with("Other"), 
                  starts_with("Blank"), starts_with("Blanc"), starts_with("Other")) %>%
    colnames()
  
  # Format 
  sw2 <- sw %>%
    gather(key = "sample_name", value = "count_reads", samples2) %>%
    dplyr::filter(count_reads>0)
  
  # Standardise the % ID column by the same name for automation
  # Rename column database ID
  database_col_name <- grep("best_identity", colnames(sw2), value=T)
  sw2[,"best_identity_database"] <- sw2[[database_col_name]]
  sw2[[database_col_name]] <- NULL
  
  return(sw2)
}

  
# ---------------------------------------------------------------------------------------------------------------- # 
# Function to clean tag-jump

# file_edna <- project_data
# file_other <- other_data
# tag_jump_value = 0.001

clean_tag_jump <- function(file_edna, file_other, tag_jump_value = 0.001){
  
  # Add project type
  file_edna$clean <- "Project"
  file_other$clean <- "Other"
  
  # Bind
  file_all <- rbind(file_edna, file_other)
  
  # Count & clean
  file_edna_clean <- file_all %>%
    group_by(run, amplicon) %>%
    mutate(somme_tot = sum(count_reads), 
           seuil = tag_jump_value * somme_tot) %>%
    ungroup() %>%
    filter(count_reads > seuil) %>%
    filter(clean == "Project") %>% select(-clean)
  
  # discarded
  file_all_discarded <- file_all %>%
    group_by(run, amplicon) %>%
    mutate(somme_tot = sum(count_reads), 
           seuil = tag_jump_value * somme_tot) %>%
    ungroup() %>%
    filter(count_reads < seuil) %>%
    filter(clean == "Project") %>% select(-clean)
  
  return(list(file_edna_clean, file_all_discarded))
}

# ---------------------------------------------------------------------------------------------------------------- # 
# Function to clean index-hoping (inter library tag-jump)

# file_edna <- project_data
# file_other <- other_data
# file_blank <- blank_data
# 
# file_edna <- rbind(project_data, other_data)

clean_index_hoping <- function(file_edna, file_blank){
  
  # Check if `lot` field in metadata_i exists, if not make a dummy considering all runs 
  if(is.null(metadata_i$lot)){
    file_edna[,"lot"] <- "A"
    file_blank[,"lot"] <- "A"
  }
  
  # Format data 
  file_edna$clean <- "other_or_project"
  file_blank$clean <- 'blank'
  
  all <- rbind(file_edna, file_blank)
  
  # Clean tag-jump intra-library first 
  tag_jump_value <- 0.001
  
  #  all2 <- all %>%
  #    group_by(run, sequence) %>%
  #    mutate(somme_tot = sum(count_reads), 
  #           seuil = tag_jump_value * somme_tot) %>%
  #    ungroup() %>%
  #    filter(count_reads > seuil)
  #  
  # For defining the seuil, take only the occurences with more than 10? more than one? 
  table_counts <- all %>%
    # TAG JUMP
    group_by(run, sequence) %>%
    mutate(somme_tot = sum(count_reads), 
           seuil = tag_jump_value * somme_tot) %>%
    ungroup() %>%
    filter(count_reads > seuil) %>%
    # TAG JUMP OVER
    group_by(lot, plaque, sequence) %>%
    summarise(n_reads_tots = sum(count_reads), 
              n_reads_other_project = sum(count_reads[clean == "other_or_project"]), 
              n_reads_blanks = sum(count_reads[clean == "blank"])) %>%
    mutate(seuil_blank = n_reads_blanks/n_reads_tots) %>%
    filter(n_reads_blanks>5) %>% ungroup()
  
  # Now, define a seuil at which the cleaning should happen -- It must be by lot!! TO CORRECT
  # Take only the occurences > 5 >10?? and remove the seuil set at 1? 
  # seuil_blank <- max(table_counts$seuil_blank)
  
  seuil_blank_df <- table_counts %>%
    # A trial??
    filter(seuil_blank != 1) %>%
    group_by(lot) %>%
    summarise(seuil_blank = max(seuil_blank))
  
  # Warning if seuil is too important 
  if(max(seuil_blank_df$seuil_blank) > 0.01) stop(paste(" error: the cleaning threshold using blanks seems too high"))
  
  # Clean the project_data
  project_counts <- all %>%
    # Joint
    left_join(., seuil_blank_df) %>%
    group_by(lot, plaque, sequence) %>%
    mutate(somme_tot_plaque_lot = sum(count_reads), 
           seuil_plaque_lot = seuil_blank*somme_tot_plaque_lot,
           discard = ifelse(count_reads < seuil_plaque_lot, "yes", "no")) %>%
    # Clean for output
    ungroup() %>% 
    # filter(project == project_i) %>% 
    filter(grepl(project_i, project, ignore.case = TRUE)) %>%
    select(-clean)
  
  #
  table(project_counts$discard)
  
  # Filter
  project_clean <- project_counts %>% filter(discard == "no") %>% select(colnames(project_data))
  project_discarded <- project_counts %>% filter(discard == "yes")
  
  # output
  seuil_blank_df$seuil_blank <- round(seuil_blank_df$seuil_blank, 5)
  return(list(project_clean, project_discarded, seuil_blank_df))
}

# ---------------------------------------------------------------------------------------------------------------- # 
# Function to clean taxonomy outputs from ecotag

# Store a file 

clean_taxonomy <- function(file){
  
  # Clean taxo
  file2 <- file %>%
    mutate(rank_ncbi_corrected = case_when(
      !is.na(species_name) & best_identity_database >= 0.99999 ~ "species",
      (!is.na(genus_name) & best_identity_database >= 0.99999 & is.na(species_name)) | (best_identity_database < 0.99999 & best_identity_database >= 0.90 & !is.na(genus_name)) ~ "genus",
      best_identity_database  >= 0.85 & !is.na(family_name) & best_identity_database < 0.90 | (best_identity_database >= 0.90 & is.na(genus_name) & !is.na(family_name)) ~ "family", 
      best_identity_database < 0.85 | is.na(order_name) | is.na(family_name) ~ "higher"
    )) %>% 
    # Set new taxa for ncbi
    mutate(scientific_name_ncbi_corrected = case_when(
      !is.na(species_name) & best_identity_database >= 0.99999 ~ species_name,
      (!is.na(genus_name) & best_identity_database >= 0.99999 & is.na(species_name)) | (best_identity_database < 0.99999 & best_identity_database >= 0.90 & !is.na(genus_name)) ~ genus_name,
      best_identity_database  >= 0.85 & !is.na(family_name) & best_identity_database < 0.96 | best_identity_database >= 0.90 & is.na(genus_name) & !is.na(family_name) ~ family_name, 
      best_identity_database < 0.85 & !is.na(order_name) | best_identity_database >= 0.85 & is.na(family_name) & !is.na(order_name) ~ order_name, 
      best_identity_database < 0.85 & is.na(order_name) | is.na(order_name) ~ scientific_name,
      is.na(order_name) & is.na(family_name) ~ scientific_name
    )) %>%
    # Correct the family name 
    mutate(family_name_corrected = case_when(
      rank_ncbi_corrected %in% c("species", "genus", "family") ~ family_name
    )) %>%
    # Correct the genus name 
    mutate(genus_name_corrected = case_when(
      rank_ncbi_corrected %in% c("species", "genus") ~ genus_name
    )) %>%
    # Correct the genus name 
    mutate(species_name_corrected = case_when(
      rank_ncbi_corrected %in% c("species") ~ species_name
    ))
  
  # Clean the columns a bit for clarity: remove taxids for each rank, stats from swarm, outputs from ecotag
  columns_to_remove <- c("id", "count", "family_name", "genus_name", "rank", "scientific_name", "species_name",
                         grep("taxid", colnames(file2), value=T), 
                         grep("match", colnames(file2), value=T), 
                         grep("species_list", colnames(file2), value=T))
  
  file3 <- file2 %>%
    select(-one_of(columns_to_remove))
  
  # Output
  return(file3)
  
}


# ---------------------------------------------------------------------------------------------------------------- # 
# Function to add a class column

# edna_file <- df_teleo_clean
# Create an archive to avoid running a long time! 

# This outputs a list, the [[1]] is the file, the [[2]] is the updated archive file 
add_class_name_archive <- function(edna_file, archive_file=NULL){
  
  require(taxize) 
  
  # Control if archive file is null or not
  if(is.null(archive_file)){
    archive_file <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("class_name", "scientific_name_ncbi_corrected"))
  }
  
  # Identify the non-fish MOTU
  # All names
  specieslist <- unique(edna_file$scientific_name_ncbi_corrected)
  
  # Are they in archive or not? 
  specieslist_present <- specieslist[specieslist %in% archive_file$scientific_name_ncbi_corrected]
  specieslist_query <- specieslist[specieslist %ni% archive_file$scientific_name_ncbi_corrected]
  
  # the classification
  liste_classification <- classification(specieslist_query, db = "ncbi", rows=1)
  
  # Transform into data frame
  liste_classification <- map_df(liste_classification, ~as.data.frame(.x), .id="initial_name")
  
  # Extract useful data
  liste_classification_class <- liste_classification %>%
    filter(rank == "class") %>%
    rename(class_name = name) %>%
    rename(scientific_name_ncbi_corrected = initial_name) %>%
    select(class_name, scientific_name_ncbi_corrected)
  #   
  #   # Condition - if true, add those new data to the archive 
  #   if(add_new_to_archive == TRUE){
  #     archive_file <- rbind(archive_file, liste_classification_class)
  #   }
  
  archive_file <- rbind(archive_file, liste_classification_class)
  
  # Add the class to the dataframe
  edna_with_class <- edna_file %>%
    left_join(., archive_file, by=c("scientific_name_ncbi_corrected"))
  
  return(list(edna_with_class, archive_file))
  
}

# ---------------------------------------------------------------------------------------------------------------- # 
# Function to clean motus based un user-input thresholds 


# Debug
  #  edna_file <- list_read_step2[[1]]
  #  remove_PCR_blanks = TRUE
  #  min_reads=10
  #  remove_chimera=TRUE
  #  remove_not_fish_taxize = TRUE
  #  min_size_seq = 20
  #  max_size_seq = 150
  #  min_PCR = 1
  #  min_percentage_id = 0.80


clean_motus <- function(edna_file, remove_PCR_blanks = TRUE, min_reads=10, remove_chimera=TRUE, remove_not_fish_taxize = TRUE,
                       min_size_seq = 20, max_size_seq = 150, min_PCR = 1, min_percentage_id = 0.80){
  
  # Fish class kept 
  fish_class <- c("Actinopteri", "Chondrichthyes")
  
  # Verification that the class column exists
  if( remove_not_fish_taxize & !('class_name' %in% colnames(edna_file)) ) stop(paste(" error: the column class_name does not exist. Please use the add_class_name function before or set the remove_not_fish_taxize option to FALSE."))
  
  # Isolate the PCR blanks for the filter
  amplicon_control <- edna_file %>%
    # Base filters
    dplyr::filter(count_reads > min_reads) %>% 
    `if`(remove_chimera, filter(., chimera == "N"), .) %>%
    # Keep only the control samples
    dplyr::filter(!stringr::str_detect(sample_name, "SPY"))
  
  # Verif - condition of non NULL for amplicon_control 
  if(nrow(amplicon_control) == 0){
    amplicon_control <- data.frame(NA)
    amplicon_control$amplicon <- NA
  } else {
    
    # Filer only fish
    fish_controls <- amplicon_control %>% filter(class_name %in% fish_class)
    
    # Print the count of MOTUs in blank
    cat(paste0(
      "\nThere is ", length(unique(amplicon_control$sequence)), " MOTUs in the blanks for the ", unique(edna_file$project), " project",
      "\nincluding ", length(unique(fish_controls$sequence)), " MOTUs belonging to fish taxa"
    ))
  }
  
  # Remove the MOTUs found in the blanks
  edna_file_filtered <- edna_file %>%
    # Remove controls
    filter(str_detect(sample_name, "SPY")) %>% 
    # Read counts
    filter(count_reads > min_reads) %>% 
    # Chimeras
    `if`(remove_chimera, filter(., chimera == "N"), .) %>%
    select(-chimera) %>%
    # Remove MOTUs in blank
    `if`(remove_PCR_blanks, filter(., !sequence %in% amplicon_control$sequence), .) %>%
    # Remove non fish - taxize
    `if`(remove_not_fish_taxize, filter(., class_name %in% fish_class), .) %>%
    # Length sequence
    mutate(length_sequence = nchar(sequence)) %>%
    filter(length_sequence > min_size_seq & length_sequence < max_size_seq) %>%
    # PCR - all dataset -- MAYBE HERE DO IT LATER ON AFTER LIST BINDING??
    group_by(sequence) %>%
    mutate(n_PCR = n_distinct(sample_name)) %>%
    ungroup() %>%
    filter(n_PCR > min_PCR) %>%
    # Select minimum % of ID
    filter(best_identity_database > min_percentage_id) %>%
    # Remove salmo salar or Salmo ID - lab contaminants
    filter(scientific_name_ncbi_corrected %ni% c("Salmo salar", "Salmo"))
  
  return(list(edna_file_filtered, amplicon_control))
}

# ---------------------------------------------------------------------------------------------------------------- # 
# Function to clean motus based on LULU



apply_lulu <- function(edna_file, path_lulu, match_lulu = 84, co_occurence_lulu = 0.95){
  
  require(lulu)
  
  # Create necessary files
  otu_seq <- edna_file %>%
    distinct(definition, sequence) %>%
    mutate(sequence = toupper(sequence))
  
  # Convert to fasta
  fa <- character(2 * nrow(otu_seq))
  fa[c(TRUE, FALSE)] = sprintf(">%s", otu_seq$definition)
  fa[c(FALSE, TRUE)] = as.character(otu_seq$sequence)
  writeLines(fa, paste(path_lulu, "OTU_sequences.fasta", sep=""))
  
  # Launch blastn with bash
  system(paste("makeblastdb -in ", path_lulu, "OTU_sequences.fasta -parse_seqids -dbtype nucl", sep=""))
  
  # Matchlist with blast
  # We can change -task to 1) "blastn" or 2) "megablast"
  system(paste("blastn -db ", path_lulu, "OTU_sequences.fasta -outfmt '6 qseqid sseqid pident' -out  ", path_lulu, "match_list.txt -task megablast -qcov_hsp_perc 60 -perc_identity 60 -query ", path_lulu, "OTU_sequences.fasta", sep=""))
  system(paste("blastn -db ", path_lulu, "OTU_sequences.fasta -outfmt '6 qseqid sseqid pident' -out  ", path_lulu, "match_list_blastn.txt -task blastn -qcov_hsp_perc 60 -perc_identity 60 -query ", path_lulu, "OTU_sequences.fasta", sep=""))
  
  # ------------------------------ # 
  # Analysis
  
  # Format for LULU function
  otutab <- edna_file %>%
    select(definition, count_reads, sample_name) %>% # Here, sample name or sample name no pcr to try 
    spread(sample_name, count_reads) %>%
    replace(., is.na(.), "0") %>%
    as.data.frame()
  
  # Put the otu name in row name
  # Steps because the conversion takes out the rownames
  otutab2 <- otutab
  otutab$definition <- NULL
  otutab <- otutab %>% mutate_if(is.character, as.numeric)
  rownames(otutab) <- otutab2$definition
  
  # Import matchlist
  matchlist <- read.table(paste(path_lulu, "match_list.txt", sep=""), header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
  matchlist_blastn <- read.table(paste(path_lulu, "match_list_blastn.txt", sep=""), header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
  
  # Run LULU 
  # megablast: 11 MOTUs lost, blastn: 25 MOTUs lost 
  lulu_edna <- lulu(otutab, matchlist, minimum_match = match_lulu, minimum_relative_cooccurence = co_occurence_lulu) # 94% is 3 mismatchs
  
  # Clean the logs 
  system("rm lulu.log*")
  
  # Output
  return(lulu_edna)
}

# ---------------------------------------------------------------------------------- # 
#### infos_statistiques ----

infos_statistiques <- function(file){
  Reads <- sum(file$count_reads)
  MOTUs <- file %>% summarise(n = n_distinct(sequence)) %>% pull()
  Species <-  file %>% filter(!is.na(species_name_corrected))%>%summarise(n = n_distinct(species_name_corrected)) %>% pull()
  
  ab <- data.frame(Reads, MOTUs, Species)
  return(ab)
}
