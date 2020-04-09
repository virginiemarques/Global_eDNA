# Extract all the species information from the reference database - add the family level
# Use the previous script

# lib
library(seqRFLP)

# functions 
source("scripts/00_functions.R")
conflict_prefer("filter", "dplyr")

# data
load("data/reference_database/teleo/04_teleo_database.Rdata")

# --------------------------------------------------------------------------------------- # 
#                             Functions: 
# --------------------------------------------------------------------------------------- # 

# filter fish from a vector, using the class from ncbi
class_vector_get_taxo <- function(x, level = "family"){ 
  
  require(taxize)
  
  # Conditions
  if (!level %in% c("kingdom", "phylum", "class", "order", "family", "genus", "species")) {
    stop("this level does not exist, please choose an existing taxonomic rank")
  }
  
  # the classification
  liste_classification <- classification(x, db = "ncbi", rows=1)
  
  # Transform into data frame
  liste_classification <- map_df(liste_classification, ~as.data.frame(.x), .id="initial_name")
  
  # Extract useful data
  liste_classification_taxo <- liste_classification %>%
    filter(rank == level) %>%
    rename(!!level := name) %>%
    rename(species_name = initial_name) %>%
    select(species_name, !!level) 
  
  # output
  return(liste_classification_taxo)
}

# --------------------------------------------------------------------------------------- # 
# End of functions 
# --------------------------------------------------------------------------------------- # 

# First on all species
class_all_sp_family <- class_vector_get_taxo(class_all_sp$species_name, level = "family")

# Get the class also
class_all_sp_family_class <- class_all_sp_family %>%
  left_join(., class_all_sp)

# Then on resolutive species 
class_resolutive_sp_family_class <- class_resolutive_sp %>%
  left_join(., class_all_sp_family) 

# Save
save(class_all_sp_family_class, class_resolutive_sp_family_class, file = "data/reference_database/teleo/05_reference_database_family.Rdata")





