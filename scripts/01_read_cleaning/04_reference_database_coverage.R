# Extract all the species information from the reference database

# lib
library(seqRFLP)

# functions 
source("scripts/00_functions.R")
conflict_prefer("filter", "dplyr")

# --------------------------------------------------------------------------------------- # 
#                             Functions: 
# --------------------------------------------------------------------------------------- # 

# Extract species name from .fasta database
gb_to_txt <- function(entry){
  
  sp_gb <- sub('.*species_name=', '', entry)
  sp_gb2 <- sub('; family=.*', '', sp_gb)
  
  # Extract only species names (odd number lines)
  numbers <- seq(from=1, to=length(sp_gb2), by=2)
  
  # species from genbank
  sp_gb3 <- str_trim(sp_gb2[numbers])
  
  # reformate correctly
  sp_teleo <- gsub('\\ ', '_', sp_gb3)
  sp_teleo <- unique(sp_teleo)
  return(sp_teleo)
}

# filter fish from a vector, using the class from ncbi
class_vector <-function(x){ 
  
  require(taxize)
  
  # the classification
  liste_classification <- classification(x, db = "ncbi", rows=1)
  
  # Transform into data frame
  liste_classification <- map_df(liste_classification, ~as.data.frame(.x), .id="initial_name")
  
  # Extract useful data
  liste_classification_class <- liste_classification %>%
    filter(rank == "class") %>%
    rename(class_name = name) %>%
    rename(species_name = initial_name) %>%
    select(species_name, class_name) 
}

# --------------------------------------------------------------------------------------- # 
# End of functions 
# --------------------------------------------------------------------------------------- # 

# Open complete database - including all taxa (fish, non-fish) and those not resolutive at species level
all_database <- read.fasta("data/reference_database/teleo/v_embl_std_clean.fasta")
# Open resolutive database - including all taxa (fish, non-fish) and only the resolutive species level (the non-resolutive are present but no species name are visible so we can parse it)
resolutive_database <- read.fasta("data/reference_database/teleo/db_embl_std.fasta")

# Parse species_name
all_database_sp <- gb_to_txt(all_database)
resolutive_database_sp <- gb_to_txt(resolutive_database)

# Check who is a fish - long operation
# need to cut in two as too heavy in one shot
resolutive_database_sp_part1 <- resolutive_database_sp[seq(1, round(length(resolutive_database_sp) /2, 0))]
resolutive_database_sp_part2 <- resolutive_database_sp[seq(length(resolutive_database_sp_part1), length(resolutive_database_sp))]

class_resolutive_sp_part1 <- class_vector(resolutive_database_sp_part1) %>%
  filter(class_name %in% c("Actinopteri", "Chondrichthyes"))

class_resolutive_sp_part2 <- class_vector(resolutive_database_sp_part2) %>%
  filter(class_name %in% c("Actinopteri", "Chondrichthyes"))

# Rbind
class_resolutive_sp <- rbind(class_resolutive_sp_part1, class_resolutive_sp_part2)

# need to cut in two as too heavy in one shot
all_database_sp_part1 <- all_database_sp[seq(1, round(length(all_database_sp) /2, 0))]
all_database_sp_part2 <- all_database_sp[seq(length(all_database_sp_part1), length(all_database_sp))]

class_all_sp_part1 <- class_vector(all_database_sp_part1) %>%
  filter(class_name %in% c("Actinopteri", "Chondrichthyes"))

class_all_sp_part2 <- class_vector(all_database_sp_part2) %>%
  filter(class_name %in% c("Actinopteri", "Chondrichthyes"))

# Rbind
class_all_sp <- rbind(class_all_sp_part1, class_all_sp_part2)

# Save
save(class_resolutive_sp, class_all_sp, file = "data/reference_database/teleo/04_teleo_database.Rdata")

# Export csv
write.csv(class_resolutive_sp, "data/reference_database/teleo/csv/reference_ncbi_all_species_resolutive_teleo.csv", row.names = F)
write.csv(class_all_sp, "data/reference_database/teleo/csv/reference_ncbi_all_species_teleo.csv", row.names = F)
