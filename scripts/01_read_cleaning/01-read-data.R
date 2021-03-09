# This script reads eDNA metabarcoding data from the SWARM clustering pipeline 
# It cleans and formats data correctly 


# Lib 
library(tidyverse)
library(data.table)

# Source functions
source("scripts/01_read_cleaning/00_functions.R")

# List the directories 
list_projects_dir <- list.dirs(path = "data/swarm", full.names = TRUE, recursive = F)

# Create outputs list
list_read_step1 <- list()
list_clean_lot_discarded <- list()

directories_multiples_projects <- c("malpelo", "fakarava", "santamarta", "providencia")

# For metadata field 
columns_delete_field_metadata <- c("turbidity", "gps_start", "gps_b", "lat_gps_b", "long_gps_b", "gps_c", "long_gps_c", "lat_gps_d", "gps_half_turn", "longitude_turn", "latitude_end", "longitude_end", 
                    "gps_end", "long_gps_d", "gps_d", "lat_gps_c", "latitude_turn", "data_manager", "gps_owner", "chimera")


metadata_field <- read.csv("metadata/Metadata_eDNA_global_V6.csv")

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- # 
# Step 1 - Assemble & clean
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- # 

# Clean the index-hoping (ie. same plate number, different library)
# Clean the tag-jump (1/1000 threshold for each MOTU within each library)

for(i in 1:length(list_projects_dir)){
  
  # i = 2
  dir_i <- list_projects_dir[[i]]
  project_i <- word(dir_i, sep="/", 3)
  
  # Open files --- DO A WARNING HERE IF 0 FILES 
  files_i <- list.files(path=dir_i, pattern = "(.*)teleo(.*)csv", recursive = F, include.dirs = FALSE)
  
  if(length(files_i)==0){next}
  
  # Metadata file 
  metadata_i <- fread(paste0(dir_i, "/metadata/all_samples.csv"), sep=";", h=F, stringsAsFactors = F) # %>% # There is sometimes a bug where '-' in original metadata ends up being "." after processing - correct it here -- it only happends with read.csv!! (weird) not fread mutate(V3 = gsub("\\-","\\.",V3))
  colnames(metadata_i) <- c("plaque", "run", "sample_name", "project", "marker")
  
  # ----- # Open files 
  
  # For the project 
  project_table <- fread(
    paste0(dir_i, "/", grep(paste0(
      project_i, "(.*)table"), files_i, value=T, ignore.case = TRUE)),
    sep="\t", stringsAsFactors = F, h=T)
  
  project_taxo <- fread(
    paste0(dir_i, "/", grep(paste0(
      project_i, "(.*)ecotag"), files_i, value=T, ignore.case = TRUE)),
    sep="\t", stringsAsFactors = F, h=T)
  
  # -------- # For the other files
  # If there are multiple projects within a directory, the other projects needs to be grouped with the 'Other' files. 
  # Ifelse condition here to filter 
  # TEST ICI: REMOVE THE BLANKS TO BE ADDED: IT IS NOT OTHER FILES 
  
  # CONDITION 1: Other et pas de multiple
  # CONDITION 2: Other et multiple 
  # CONDITION 3: Pas de Other
  # CONDITION 4: Pas de other et multiple
  
  # project_i %in% directories_multiples_projects
  # length(grep("Other(.*)", files_i, value=T)) == 0 # There is no others files
  # length(grep("Other(.*)", files_i, value=T)) > 0 # There is others files
  
  if(project_i %ni% directories_multiples_projects & length(grep("Other(.*)", files_i, value=T)) > 0){
    {
      ##### IF THERE IS A SINGLE PROJECT WITHIN THE DIRECTORY (CLASSIC CASE) & THERE IS OTHERS FILES 
      
      # Other 
      other_table <- fread(paste0(dir_i, "/", grep("Other(.*)table", files_i, value=T)), sep="\t", stringsAsFactors = F, h=T)
      other_taxo <- fread(paste0(dir_i, "/", grep("Other(.*)ecotag", files_i, value=T)), sep="\t", stringsAsFactors = F, h=T)
      
      # Assemble
      other_data <- assemble_data(table_otu = other_table, taxo_otu = other_taxo) %>%
        left_join(., metadata_i)
      
    }
  } else if (project_i %in% directories_multiples_projects & length(grep("Other(.*)", files_i, value=T)) > 0){
    ##### IF THERE IS A MULTIPLE PROJECTS WITHIN THE DIRECTORY (CLASSIC CASE) & THERE IS OTHERS FILES 
    message(paste0("There is multiple projects within the ", project_i, " directory, so other projects will be condensed together"))
    
    # Other - normal
    other_table <- fread(paste0(dir_i, "/", grep("Other(.*)table", files_i, value=T)), sep="\t", stringsAsFactors = F, h=T)
    other_taxo <- fread(paste0(dir_i, "/", grep("Other(.*)ecotag", files_i, value=T)), sep="\t", stringsAsFactors = F, h=T)
    other_data_part1 <- assemble_data(table_otu = other_table, taxo_otu = other_taxo) %>%
      left_join(., metadata_i)
    
    # Other -- the other(s) project(s) also analyzed 
    other_projects <- unique(word(grep(paste0(project_i, "|Other"), files_i, ignore.case = T, value=T, invert=T), 1, sep="_"))
    
    list_other <- lapply(other_projects, function(x){
      
      # x = "Fakarava"
      # Same thing
      other_bis_table <- fread(paste0(dir_i, "/", grep(paste0(x, "(.*)table"), files_i, value=T)), sep="\t", stringsAsFactors = F, h=T)
      other_bis_taxo <- fread(paste0(dir_i, "/", grep(paste0(x, "(.*)ecotag"), files_i, value=T)), sep="\t", stringsAsFactors = F, h=T)
      other_bis_data <- assemble_data(table_otu = other_bis_table, taxo_otu = other_bis_taxo) %>%
        left_join(., metadata_i) %>%
        mutate(project = "Other")
    })
    
    # Bind 
    other_data_part2 <- bind_rows(list_other)
    
    # Bind all other projects together
    other_data <- rbind(other_data_part1, other_data_part2)
  } else if (project_i %ni% directories_multiples_projects & length(grep("Other(.*)", files_i, value=T)) == 0){
    ##### IF THERE IS A SINGLE PROJECT WITHIN THE DIRECTORY (CLASSIC CASE) & THERE IS NO OTHERS FILES 
    
    other_data <- data.frame()
    
  } else if (project_i %in% directories_multiples_projects & length(grep("Other(.*)", files_i, value=T)) == 0){
    ##### IF THERE IS A MULTIPLE PROJECTS WITHIN THE DIRECTORY (CLASSIC CASE) & THERE IS NO OTHERS FILES 
    
    # Other -- the other(s) project(s) also analyzed 
    other_projects <- unique(word(grep(paste0(project_i, "|Other"), files_i, ignore.case = T, value=T, invert=T), 1, sep="_"))
    
    list_other <- lapply(other_projects, function(x){
      
      # x = "Fakarava"
      # Same thing
      other_bis_table <- fread(paste0(dir_i, "/", grep(paste0(x, "(.*)table"), files_i, value=T)), sep="\t", stringsAsFactors = F, h=T)
      other_bis_taxo <- fread(paste0(dir_i, "/", grep(paste0(x, "(.*)ecotag"), files_i, value=T)), sep="\t", stringsAsFactors = F, h=T)
      other_bis_data <- assemble_data(table_otu = other_bis_table, taxo_otu = other_bis_taxo) %>%
        left_join(., metadata_i) %>%
        mutate(project = "Other")
    })
    
    # Bind 
    other_data <- bind_rows(list_other)
  }
  
  
         ###       if(project_i %in% directories_multiples_projects){
         ###         ##### IF THERE ARE MULTIPLE PROJECTS WITHIN THE DIRECTORY
         ###         message(paste0("There is multiple projects within the ", project_i, " directory, so other projects will be condensed together"))
         ###         
         ###         # Other - normal
         ###         other_table <- fread(paste0(dir_i, "/", grep("Other(.*)table", files_i, value=T)), sep="\t", stringsAsFactors = F, h=T)
         ###         other_taxo <- fread(paste0(dir_i, "/", grep("Other(.*)ecotag", files_i, value=T)), sep="\t", stringsAsFactors = F, h=T)
         ###         other_data_part1 <- assemble_data(table_otu = other_table, taxo_otu = other_taxo) %>%
         ###           left_join(., metadata_i)
         ###         
         ###         # Other -- the other(s) project(s) also analyzed 
         ###         other_projects <- unique(word(grep(paste0(project_i, "|Other|Blank"), files_i, ignore.case = T, value=T, invert=T), 1, sep="_"))
         ###         
         ###         list_other <- lapply(other_projects, function(x){
         ###           
         ###           # x = "Fakarava"
         ###           # Same thing
         ###           other_bis_table <- fread(paste0(dir_i, "/", grep(paste0(x, "(.*)table"), files_i, value=T)), sep="\t", stringsAsFactors = F, h=T)
         ###           other_bis_taxo <- fread(paste0(dir_i, "/", grep(paste0(x, "(.*)ecotag"), files_i, value=T)), sep="\t", stringsAsFactors = F, h=T)
         ###           other_bis_data <- assemble_data(table_otu = other_bis_table, taxo_otu = other_bis_taxo) %>%
         ###             left_join(., metadata_i) %>%
         ###             mutate(project = "Other")
         ###         })
         ###         
         ###         # Bind 
         ###         other_data_part2 <- bind_rows(list_other)
         ###         
         ###         # Bind all other projects together
         ###         other_data <- rbind(other_data_part1, other_data_part2)
         ###         
         ###       } else {
         ###         ##### IF THERE IS A SINGLE PROJECT WITHIN THE DIRECTORY (CLASSIC CASE)
         ###         
         ###         # Other 
         ###         other_table <- fread(paste0(dir_i, "/", grep("Other(.*)table", files_i, value=T)), sep="\t", stringsAsFactors = F, h=T)
         ###         other_taxo <- fread(paste0(dir_i, "/", grep("Other(.*)ecotag", files_i, value=T)), sep="\t", stringsAsFactors = F, h=T)
         ###         
         ###         # Assemble
         ###         other_data <- assemble_data(table_otu = other_table, taxo_otu = other_taxo) %>%
         ###           left_join(., metadata_i)
         ###         
         ###       }
         ###       
  # ----- # Assemble project data 
  project_data <- assemble_data(table_otu = project_table, taxo_otu = project_taxo) %>%
    left_join(., metadata_i)
  
  # ----- # Blanks 
  # Blanks - not always present 
  # Apply the code using the blanks only if those are present 
  
  if(length(grep("Blank(.*)", files_i, value=T)) == 0){message(paste0("There is no blank files for the ", project_i, " data"))}

  if(length(grep("Blank(.*)", files_i, value=T)) != 0){
    
    # Read 
    blank_table <- try(fread(paste0(dir_i, "/", grep("Blank(.*)table", files_i, value=T)), sep="\t", stringsAsFactors = F, h=T))
    blank_taxo <- try(fread(paste0(dir_i, "/", grep("Blank(.*)ecotag", files_i, value=T)), sep="\t", stringsAsFactors = F, h=T))
    
    # Assemble
    blank_data <- assemble_data(table_otu = blank_table, taxo_otu = blank_taxo) %>%
      left_join(., metadata_i)
    
    # CHeck NAs 
    lesna <- blank_data %>% filter(is.na(plaque))
    
    # Clean index-hoping (inter-library tag-jump)
    project_data <- clean_index_hoping(file_edna = rbind(project_data, other_data), 
                                       file_blank = blank_data)[[1]]
    project_data_blanks_discarded <- clean_index_hoping(file_edna = rbind(project_data, other_data), 
                                                        file_blank = blank_data)[[2]]
    
    seuil <- clean_index_hoping(file_edna = rbind(project_data, other_data), 
                                file_blank = blank_data)[[3]]
    
    message(paste0("The Blank threhsold for the ", project_i, " project is ", seuil$seuil_blank))

  } else { project_data_blanks_discarded <- data.frame()}
  
  # Verifs - no NA values on the run
  verif_metadata <- !is.na(project_data$run)
  if( length(verif_metadata[verif_metadata==FALSE]) > 0 ) stop(paste(" error: some samples do not have metadata fields"))

  # Verify NAs
  les_na <- project_data %>% filter(is.na(run))
  
  # Clean tag-jump (output is a list, using [[2]] will output the discarded reads)
  project_data_tag <- clean_tag_jump(file_edna = project_data, file_other = other_data)[[1]]
  
  # Add project column - normally not necessary as present in metadata but just in case 
  project_data_tag$project_i <- project_i
  
  # Store in list 
  list_read_step1[[i]] <- project_data_tag
  list_clean_lot_discarded[[i]] <- project_data_blanks_discarded
  
  # Print
  print(paste0(i, "_",  project_i))
}

save(list_read_step1, list_clean_lot_discarded, file = "Rdata/all_df_step1.Rdata")
# load("Rdata/all_df_step1.Rdata")

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- # 
# Step 2: Clean and complete taxonomy
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- # 

load("Rdata/archive_class_ncbi.Rdata")

# Remove null elements from list
list_read_step1 <- list_read_step1[!sapply(list_read_step1,is.null)]

# Apply workflow
list_read_step2 <- lapply(list_read_step1, function(file){
  
  # Filter out null object
  if(is.null(file)){next}
  
  # Clean column names 
  columns_to_remove <- c("amplicon", "family", "genus", "order", "species", "taxid", "OTU", "total", 
                         "cloud", "length", "abundance", "spread", "identity", "taxonomy", "references", 
                         "quality")
  
  file_short <- file %>%
    select(-one_of(columns_to_remove))
    
  # Clean taxonomy 
  file_taxo <- clean_taxonomy(file_short)
  
  # Add class name column 
  list_output <- add_class_name_archive(file_taxo, archive_class_ncbi)
  file_taxo_all <- list_output[[1]]
  archive_class_ncbi <- list_output[[2]]
  
  # output
  return(file_taxo_all)
})

# Save the new archive file 
save(archive_class_ncbi, file = "Rdata/archive_class_ncbi.Rdata")

# Save data
save(list_read_step1, list_clean_lot_discarded, list_read_step2, file = "Rdata/01_read_data.Rdata")



