
# Libs 
library(tidyverse)
library(conflicted)
library(taxize) # Make sure you have an API
library(purrr)
library(rentrez)
library(lulu)


# ---------------------------------------------------------------------------------- # 
#### SET THE API TO USE TAXIZE ----

set_entrez_key("e1b887b07de1764a6e68883fce0f9f69d108") # My API key
Sys.getenv("ENTREZ_KEY") 

# ---------------------------------------------------------------------------------- # 
#### READ_DATA ----

# This function reads the data and adds the metadata


read_data <- function(taxo_motu, table_otu,
                      metadata_sampling, metadata_sequencing){
  
  # Clean data - adjust to fit all - normalise by amplicon name - add ifelse for control for each category 
  taxo_motu$amplicon <- as.character(taxo_motu$definition)
  taxo_motu$amplicon <- str_trim(taxo_motu$amplicon)
  taxo_motu$sequence <- toupper(taxo_motu$sequence)
  
  # Checks - for sequences: all sequences from the table must be in the csv. (The opposite is not true due to sequence loss while cleaning rapidrun data)
  verif_seq <- table_otu$sequence %in% taxo_motu$sequence
  if( length(verif_seq[verif_seq==FALSE]) > 0 ) stop(paste(" error: some sequences are in table but not in csv file"))
  
  # Checks - for sequences: all MOTUs from the table must be in the csv. (The opposite is not true due to sequence loss while cleaning rapidrun data)
  verif_motu <- table_otu$amplicon %in% taxo_motu$amplicon
  if( length(verif_motu[verif_motu==FALSE]) > 0  ) stop(paste(" error: some MOTU names are in table but not in csv file"))
  
  # Merge data
  sw <- merge(taxo_motu, table_otu, by=c("amplicon", "sequence"))
  
  # Extract all sample names to prepare formatting
  samples2 <- sw %>%
    dplyr::select(starts_with("CPCR"), starts_with("Cext"), starts_with("SPY"), starts_with("CMET"), starts_with("CNEG")) %>%
    colnames()
  
  # Format 
  sw2 <- sw %>%
    gather(key = "sample_name", value = "count_reads", samples2) %>%
    dplyr::filter(count_reads>0)
  
  # Standardise the % ID column by the same name for automation
  sw2[,"best_identity_database"] <- sw2[,colnames(sw2)[grepl("best_identity", colnames(sw2))]]
  
  # Correct % assignement
  sw3 <- sw2 %>% 
    # Change factor to character
    mutate_if(is.factor, as.character) %>%
    # Set new rank for ncbi
    mutate(new_rank_ncbi = case_when(
      !is.na(species_name) & best_identity_database >= 0.99 ~ "species",
      (!is.na(genus_name) & best_identity_database >= 0.99 & is.na(species_name)) | (best_identity_database < 0.99 & best_identity_database >= 0.90 & !is.na(genus_name)) ~ "genus",
      best_identity_database  >= 0.85 & !is.na(family_name) & best_identity_database < 0.90 | (best_identity_database >= 0.90 & is.na(genus_name) & !is.na(family_name)) ~ "family", 
      best_identity_database < 0.85 | is.na(order_name) | is.na(family_name) ~ "higher"
    )) %>% 
    # Set new taxa for ncbi
    mutate(new_scientific_name_ncbi = case_when(
      !is.na(species_name) & best_identity_database >= 0.99 ~ species_name,
      (!is.na(genus_name) & best_identity_database >= 0.99 & is.na(species_name)) | (best_identity_database < 0.99 & best_identity_database >= 0.90 & !is.na(genus_name)) ~ genus_name,
      best_identity_database  >= 0.85 & !is.na(family_name) & best_identity_database < 0.96 | best_identity_database >= 0.90 & is.na(genus_name) & !is.na(family_name) ~ family_name, 
      best_identity_database < 0.85 & !is.na(order_name) | best_identity_database >= 0.85 & is.na(family_name) & !is.na(order_name) ~ order_name, 
      best_identity_database < 0.85 & is.na(order_name) | is.na(order_name) ~ scientific_name,
      is.na(order_name) & is.na(family_name) ~ scientific_name
    )) %>%
    # Prepare sample name without PCR number
    mutate(sample_name_all_pcr = word(.$sample_name, 1, sep="_")) %>%
    # Count the reads at the filter level (not PCR level) 
    group_by(sample_name_all_pcr, amplicon) %>% 
    mutate(count_reads_all_pcr = sum(count_reads)) %>%
    ungroup() %>%
    # Taxid for NCBI
    mutate(taxid_ncbi = taxid) %>%
    # Add the run metadata
    left_join(., metadata_sequencing, by= "sample_name") %>%
    # Add the metadata
    left_join(., metadata_sampling, by = c("sample_name_all_pcr" = "code_spygen")) %>%
    # Remove temoin / control if needed
    dplyr::filter(!grepl("temoin", station)) %>%
    dplyr::filter(!grepl("control", station)) %>%
    # Correct the family name 
    mutate(new_family_name = case_when(
      new_rank_ncbi %in% c("species", "genus", "family") ~ family_name
    )) %>%
    # Correct the genus name 
    mutate(new_genus_name = case_when(
      new_rank_ncbi %in% c("species", "genus") ~ genus_name
    )) %>%
    # Correct the genus name 
    mutate(new_species_name = case_when(
      new_rank_ncbi %in% c("species") ~ species_name
    ))
  
  # Clean the columns a bit for clarity: remove taxids for each rank, stats from swarm, outputs from ecotag
  columns_to_remove <- c("definition", "id", "taxid", "count", "family", "genus", "order", "rank", "scientific_name",
                         "species", "OTU", "total", "cloud", "length", "abundance", "spread", "quality", "identity", "taxonomy", 
                         "references", "family_name", "genus_name", "species_name", "taxid_ncbi") 
  
  sw4 <- sw3 %>%
    select(-one_of(columns_to_remove)) %>%
    select(-starts_with("taxid_by_db")) %>%
    select(-starts_with("best_identity.db")) %>%
    select(-starts_with("species_list")) %>%
    select(-starts_with("match_count")) %>%
    select(-starts_with("best_identity:"))
  
  
  # Return object
  return(sw4)
}

# ---------------------------------------------------------------------------------- # 
#### ADD_CLASS ----

# This function allows to add the class of each taxa using the package taxize
# It can be time consuming (15 min for 1200 queries) but is more efficient and automatic

add_class_name <- function(edna_file){
  
  # Identify the non-fish MOTU
  # All names
  specieslist <- unique(edna_file$new_scientific_name_ncbi)
  
  # the classification
  liste_classification <- classification(specieslist, db = "ncbi", rows=1)
  
  # Transform into data frame
  liste_classification <- map_df(liste_classification, ~as.data.frame(.x), .id="initial_name")
  
  # Extract useful data
  liste_classification_class <- liste_classification %>%
    filter(rank == "class") %>%
    rename(class_name = name) %>%
    rename(new_scientific_name_ncbi = initial_name) %>%
    select(new_scientific_name_ncbi, class_name)
  
  # Add the class to the dataframe
  edna_with_class <- edna_file %>%
      left_join(., liste_classification_class, by=c("new_scientific_name_ncbi"))

  return(edna_with_class)
  
}


# ---------------------------------------------------------------------------------- # 
#### CLEAN_DATA ----

# ----------------------------------------- # 
# In this function, you choose which filter to apply to the data 
# There is some arguments pre-built to allow user's choice on which filter to apply 

# remove_blanks remove MOTUs found in PCR blanks, default is TRUE
# min_reads sets the minimum number of reads to be filtered, default is 10
# remove_chimera removes sequences identified as chimeras, default is TRUE
# remove_not_fish_manual removes sequences assigned to taxa other than fish manually, you have to fill the none fish orders and taxa and the beggining of this script. Default is FALSE. 
# remove_not_fish_taxize removes sequences assigned to taxa other than fish automatically using taxize. Here, fish comprise both fishes and sharks. You need to have used the add_class_name function before for this option to work. Default is TRUE. 
# ** You can control what taxa are removed at the beginning of the script: there is a dual control for efficiency purpose: first remove all orders not belonging to fish 
# ** Then, as ecotag is not perfect and sometimes assigns NAs to orders to both fish and non-fish sequences, remove all scientitic_name not assign to fish
# min_size_seq is the minimum sequence length to be kept, default is 30
# max_size_seq is the maximum sequence length to be kept, default is 100
# tag_jump is a filter removing all reads within a sequencing run < 0.0001 per sequence, it cleans for tag-jumps, that can wrongly assign some reads to a filter when its just contamination from a more abundant other filter, default is TRUE
# tag_jump_value sets the filter at 0.001: you can change it
# min_PCR set the minimum number of PCRs in the entire dataset necessary to be kept, default is 1
# min_PCR_sample is an alternative PCR filter; where the minimum number of PCR is not calculated on the entire dataset but on a pair of field replicates. Default is 0 (not applied) 
# **It is usually not used in tropical reef eDNA as too much information is lost, but clasfically performed in stream eDNA or for Lengguru as the protocol is different from other studies
# habitat_select is a value allowing to filter samples based on its habitat (provided in metadata). It is a vector. You can visualize the list of habitats using the command: unique(metadata_sampling$habitat). Default is set to keep only marine. 
# min_percentage_id is a value set to remove the sequences with a very low % of ID. There is no consensus on such a value, but some low % sequence assigned to fish could not really correspond to fish species. Default is 80%. 
# delete_gps_col controls the deletion of multiples GPS columns -- to have a small dataset if those data are not needed. Since the protocols are different, there is a lot of ways to store GPS. Default is TRUE
# ----------------------------------------- # 


clean_data <- function(edna_file, remove_blanks = TRUE, min_reads=10, remove_chimera=TRUE, remove_not_fish_manual=FALSE, remove_not_fish_taxize = TRUE,
                       min_size_seq = 30, max_size_seq = 100, tag_jump=TRUE, tag_jump_value = 0.001, min_PCR = 1,
                       min_PCR_sample = 0, habitat_select = c("marine"), min_percentage_id = 0.80, delete_gps_col=TRUE){
  
  # Verify that the class column exists
  if( remove_not_fish_taxize & !('class_name' %in% colnames(edna_file)) ) stop(paste(" error: the column class_name does not exist. Please use the add_class_name function before or set the remove_not_fish_taxize option to FALSE."))
  
  # Isolate the PCR blanks for the filter
  amplicon_control <- edna_file %>%
    # Base filters
    dplyr::filter(count_reads > min_reads) %>% 
    `if`(remove_chimera, filter(., chimera == "N"), .) %>%
    # Threshold 1/1000
    group_by(sample_name_all_pcr, amplicon) %>% 
    mutate(count_reads_all_pcr = sum(count_reads)) %>%
    ungroup() %>%
    group_by(Run, amplicon) %>%
    mutate(somme_tot = sum(count_reads), seuil = tag_jump_value * somme_tot) %>%
    ungroup() %>%
    `if`(tag_jump, filter(., count_reads_all_pcr > seuil), .) %>%
    # Keep only the control samples
    dplyr::filter(!stringr::str_detect(sample_name, "SPY"))
  
  # Verif - condition of non NULL for amplicon_control 
  if(nrow(amplicon_control) == 0){
    amplicon_control <- data.frame(NA)
    amplicon_control$amplicon <- NA
  } else {
    # Print the count of MOTUs in blank
    print(paste0(
      "There is ", length(unique(amplicon_control$amplicon)), " MOTUs in the blanks for the ", unique(edna_file$project_name), " project"
    ))
  }
  
  # Prepare columns to delete -- you can also change it here if necessary
  columns_delete <- c("turbidity", "gps_start", "gps_b", "lat_gps_b", "long_gps_b", "gps_c", "long_gps_c", "lat_gps_d", "gps_half_turn", "longitude_turn", "latitude_end", "longitude_end", 
                      "gps_end", "long_gps_d", "gps_d", "lat_gps_c", "latitude_turn", "data_manager", "gps_owner", "chimera")
  
  # Remove the MOTUs found in the blanks
  edna_file_filtered <- edna_file %>%
    # Remove controls
    filter(str_detect(sample_name, "SPY")) %>% 
    # Read counts
    filter(count_reads > min_reads) %>% 
    # Chimeras
    `if`(remove_chimera, filter(., chimera == "N"), .) %>%
    # Remove MOTUs in blank
    `if`(remove_blanks, filter(., !amplicon %in% amplicon_control$amplicon), .) %>%
    # Remove non fish - manual 
    `if`(remove_not_fish_manual, filter(., !order_name %in% non_fish_orders), .) %>%
    `if`(remove_not_fish_manual, filter(., !new_scientific_name_ncbi %in% non_fish_taxa), .) %>%
    # Remove non fish - taxize
    `if`(remove_not_fish_taxize, filter(., class_name %in% c("Actinopteri", "Chondrichthyes")), .) %>%
    # Length sequence
    mutate(length_sequence = nchar(sequence)) %>%
    filter(length_sequence > min_size_seq & length_sequence < max_size_seq) %>%
    # Tag-jump, 1/1000
    group_by(sample_name_all_pcr, amplicon) %>% 
    mutate(count_reads_all_pcr = sum(count_reads)) %>%
    ungroup() %>%
    group_by(Run, amplicon) %>%
    mutate(somme_tot = sum(count_reads), seuil = tag_jump_value * somme_tot) %>%
    ungroup() %>%
    `if`(tag_jump, filter(., sample_name_all_pcr > seuil), .) %>%
    # PCR - all dataset
    group_by(amplicon) %>%
    mutate(n_PCR = n_distinct(sample_name)) %>%
    ungroup() %>%
    filter(n_PCR > min_PCR) %>%
    # PCR - at the site level
    group_by(amplicon, station) %>% # careful here relative to the standardization -- need to control for data quality to perform the filter correctly
    mutate(n_PCR_samples = n_distinct(sample_name)) %>%
    ungroup() %>%
    filter(n_PCR_samples > min_PCR_sample) %>%
    # Choose habitats
    filter(habitat %in% habitat_select) %>%
    # Select minimum % of ID
    filter(best_identity_database > min_percentage_id) %>%
    # Clean columns
    `if`(delete_gps_col, select(., -one_of(columns_delete)), .)
  
  
  return(edna_file_filtered)
  
}


# ---------------------------------------------------------------------------------- # 
#### PCR_LEVEL ----

# Summarize a dataset at the sample level instead of PCR level

simplify_sample_level <- function(edna_file){
  
  # Filter
  edna_file_simplified <- edna_file %>%
    distinct(amplicon,sample_name_all_pcr, .keep_all = TRUE) %>%
    select(-sample_name, -count_reads)
  
  return(edna_file_simplified)
  
}

# ---------------------------------------------------------------------------------- # 
#### APPLY LULU -----

# To apply this function, you need to have a UNIX OS system (The system function doesn't work on windows)
# And you need to install the blastn tools in your local machine 
# The output is the complete list out of the lulu fonction, see its documention ?lulu

apply_lulu <- function(edna_file, path_lulu, match_lulu = 84, co_occurence_lulu = 0.95){
  
  # Create necessary files
  otu_seq <- edna_file %>%
    distinct(amplicon, sequence) %>%
    mutate(sequence = toupper(sequence))
  
  # Convert to fasta
  fa <- character(2 * nrow(otu_seq))
  fa[c(TRUE, FALSE)] = sprintf(">%s", otu_seq$amplicon)
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
    select(amplicon, count_reads, sample_name) %>%  
    spread(sample_name, count_reads) %>%
    replace(., is.na(.), "0") %>%
    as.data.frame()
  
  # Put the otu name in row name
  # Steps because the conversion takes out the rownames
  otutab2 <- otutab
  otutab$amplicon <- NULL
  otutab <- otutab %>% mutate_if(is.character, as.numeric)
  rownames(otutab) <- otutab2$amplicon
  
  # Import matchlist
  matchlist <- read.table(paste(path_lulu, "match_list.txt", sep=""), header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
  matchlist_blastn <- read.table(paste(path_lulu, "match_list_blastn.txt", sep=""), header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
  
  # Run LULU 
  lulu_edna <- lulu(otutab, matchlist, minimum_match = match_lulu, minimum_relative_cooccurence = co_occurence_lulu) 
  
  # Clean the logs 
  system("rm lulu.log*")
  
  # Output
  return(lulu_edna)
}
