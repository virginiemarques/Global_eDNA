# This script reads eDNA metabarcoding data from the SWARM clustering pipeline 
# It cleans and formats data correctly 


# Lib 
library(tidyverse)
library(data.table)
library(lulu)
'%ni%' <- Negate("%in%")

# Source functions
source("scripts/01_read_cleaning/00_functions.R")

# Load data
load("Rdata/01_read_data.Rdata")

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- # 
# Step 3 - data cleaning 
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- # 

# Faire un bout de manual cleaning pour enlever les trucs bizarre de freshwater 

# Combien d'sp sont perdues par projet 

# Detail technique: une seule PCR sur la totalité du jeu de données, ou par projet (comme c'est le cas en ce moment) ? 
# Clean with default values
list_read_step3 <- lapply(list_read_step2, function(x){
  clean_motus(x, min_PCR = 1)[[1]]
})
  
# Extract the data present in only one PCR and which is discarded 
list_read_step3_0PCR <- lapply(list_read_step2, function(x){
  clean_motus(x, min_PCR = 0)[[1]] %>%
    filter(n_PCR == 1)
})

# 
fish_species_0_PCR <- list_read_step3_0PCR %>%
  bind_rows() %>%
  filter(class_name %in% c("Actinopteri", "Chondrichthyes")) %>%
  # filter(best_identity_database == 1) %>%
  distinct(sequence,best_identity_database ,scientific_name_ncbi_corrected, project)

fish_species_0_PCR_explo <- list_read_step3_0PCR %>%
  bind_rows() %>%
  filter(class_name %in% c("Actinopteri", "Chondrichthyes")) %>%
  filter(best_identity_database == 1) %>%
  distinct(sequence, scientific_name_ncbi_corrected, project)

message(paste0("There is ", nrow(fish_species_0_PCR), " lost fish species during the PCR cleaning step (for each project separatly), including ", length(unique(fish_species_0_PCR$scientific_name_ncbi_corrected)), " unique taxa"))

# all_df_before <- bind_rows(list_read_step2)
# all_df_lost <- bind_rows(list_read_step3_0PCR)
# all_df_clean <- bind_rows(list_read_step3)


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- # 
# Step 4 - LULU post-clustering
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- # 

# Set the path 
path_lulu <- "data/LULU/"

# To apply this function, you need to have a UNIX OS system (The system function doesn't work on windows I think)
# And you need to install the blastn tools in your local machine 

list_read_step4 <- lapply(list_read_step3, function(x){
  # Apply LULU
  list_LULU <- apply_lulu(x, path_lulu)
  # MOTUs to keep
  lulu_motus_keep <- list_LULU$curated_otus
  # Filter out discarded MOTUs
  x_filtered <- x %>%
    filter(definition %in% lulu_motus_keep)
  # End 
  return(x_filtered)
})

# Cout the number of MOTUs before LULU
count_motu_before_lulu <- lapply(list_read_step3, function(x){
  length(unique(x$definition))
})

# Cout the number of MOTUs after LULU
count_motus_after_lulu <- lapply(list_read_step4, function(x){
  length(unique(x$definition))
})

# Lost MOTUs due to LULU 
Map('-', count_motu_before_lulu, count_motus_after_lulu)


## Join metadata
columns_delete_field_metadata <- c("turbidity", "gps_start", "gps_b", "lat_gps_b", "long_gps_b", "gps_c", "long_gps_c", "lat_gps_d", "gps_half_turn", "longitude_turn", "latitude_end", "longitude_end", 
                                   "gps_end", "long_gps_d", "gps_d", "lat_gps_c", "latitude_turn", "data_manager", "gps_owner", "project")

metadata_field <- read.csv("metadata/Metadata_eDNA_global_V6.csv", sep=";", stringsAsFactors = F)
metadata_field <- select(metadata_field, -c(columns_delete_field_metadata))

list_read_step4 <- lapply(list_read_step4, function(x){
  x$sample_name_all_pcr <- substr(x$sample_name,1, nchar(x$sample_name) -3)
  x <- left_join(x, metadata_field, by=c("sample_name_all_pcr" = "code_spygen"))
})

# Clean wrong family assigment

# BUG chez moi ici 
list_read_step4 <- lapply(list_read_step4, function(x){
  scaridae <- c("Scarus", "Bolbometopon", "Cetoscarus", "Chlorurus", "Hipposcarus", "Calotomus", "Cryptotomus", "Leptoscarus", "Nicholsina", "Sparisoma")
  x_scaridae <- x%>%
  filter(genus_name_corrected%in%scaridae)
  x <- x%>%
  filter(genus_name_corrected%ni%scaridae)
  x_scaridae$family_name_corrected <- "Scaridae"
  x <- rbind(x, x_scaridae)
  x_taeniura <- x%>%
    filter(genus_name_corrected=="Taeniura")
  x <- x%>%
    filter(genus_name_corrected!="Taeniura")
  x_taeniura$family_name_corrected <- "Dasyatidae"
  x <- rbind(x, x_taeniura)
})

# même chose codé différement pour que ça marche
scaridae <- c("Scarus", "Bolbometopon", "Cetoscarus", "Chlorurus", "Hipposcarus", "Calotomus", "Cryptotomus", "Leptoscarus", "Nicholsina", "Sparisoma")
list_read_step4 <- lapply(list_read_step4, function(x){
  
  x <- x %>%
    mutate(family_name_corrected = case_when(
      # Correct Scaridae
      genus_name_corrected %in% scaridae ~ "Scaridae", 
      genus_name_corrected %in% c("Taeniura") ~ "Dasyatidae", 
      TRUE ~ family_name_corrected
    )) 
    
  return(x)
})
 
df_all_filters <- bind_rows(list_read_step4)

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- # 
# save files
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- # 

save(list_read_step3, list_read_step3_0PCR, list_read_step4, df_all_filters, file = "Rdata/02-clean-data.Rdata")
