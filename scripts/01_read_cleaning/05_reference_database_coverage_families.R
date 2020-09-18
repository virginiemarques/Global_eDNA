# Extract all the species information from the reference database - add the family level
# Use the previous script

# lib
library(seqRFLP)
library(tidyverse)

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


# The same function, but cutting in several tasks because its annoying when it stops 
# Set parameters: maximum number of iterations in a single demand 
class_vector_get_taxo <- function(x, max_iterations = 1000, level = "family"){ 
  # Lib
  require(taxize)
  
  # Conditions
  if (!level %in% c("kingdom", "phylum", "class", "order", "family", "genus", "species")) {
    stop("this level does not exist, please choose an existing taxonomic rank")
  }
  
  # cut in n pieces 
  n <- round(length(x) / max_iterations, 0)
  
  # Store in a list 
  list_sp <- split(x, sort(rep_len(1:n, length(x))))
  
  # Apply on each element of the list
  classif_list_sp <- lapply(list_sp, function(z){
    # the classification
    liste_classification <- classification(z, db = "ncbi", rows=1)
    
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
  })
  
  # unlist
  liste_classification_taxo <- map_df(classif_list_sp, ~as.data.frame(.x))
  
  # output
  return(liste_classification_taxo)
}


# --------------------------------------------------------------------------------------- # 
# End of functions 
# --------------------------------------------------------------------------------------- # 

# First on all species
class_all_sp_family <- class_vector_get_taxo(class_all_sp$species_name, max_iterations = 500, level = "family")

# Get the class also
class_all_sp_family_class <- class_all_sp %>%
  left_join(., class_all_sp_family)

# Then on resolutive species 
class_resolutive_sp_family_class <- class_resolutive_sp %>%
  left_join(., class_all_sp_family) 

# Save
save(class_all_sp_family_class, class_resolutive_sp_family_class, class_all_sp_family, file = "data/reference_database/teleo/05_reference_database_family.Rdata")


# --------------------------------------------------------------------------------------- # 
# Indice of resolution per family 
# --------------------------------------------------------------------------------------- # 

library(rfishbase)

# Load all fishbase data
all_fishbase <- load_taxa() %>%
  mutate(Species_name = gsub("\\ ", "_", Species))

# correct species name - all species
all_species_sequenced_ncbi <- class_all_sp_family_class$species_name
all_species_sequenced_fishbase <- gsub("\\ ", "_", validate_names(gsub("\\_", " ", class_all_sp_family_class$species_name)))

# correct species name - resolutive species
resolutive_species_sequenced_ncbi <- class_resolutive_sp_family_class$species_name
resolutive_species_sequenced_fishbase <- gsub("\\ ", "_", validate_names(gsub("\\_", " ", class_resolutive_sp_family_class$species_name)))

# Family stats
# ATTENTION! WE NEED A BETTER WAY TO CALCULATE RESOLUTION 
# HERE WE CALCULATE JUST THE NUMBER OF RESOLUTIVE SPECIES, WHAT WE WOULD NEED IS THE NUMBER OF DIFFERENT SEQUENCES (?) PER FAMILY, BUT EXCLUDING INTRA-SPECIFIC VARIABILITY

family_stat <- all_fishbase %>%
  group_by(Family) %>%
  mutate(sequenced_teleo_all = ifelse(Species_name %in% c(all_species_sequenced_ncbi, all_species_sequenced_fishbase), 1, 0)) %>%
  mutate(sequenced_teleo_resolutive = ifelse(Species_name %in% c(resolutive_species_sequenced_ncbi, resolutive_species_sequenced_fishbase), 1, 0)) %>%
  summarise(n_sp_tot = n_distinct(Species_name), 
            n_sequenced = sum(sequenced_teleo_all), 
            n_sequenced_resolutive = sum(sequenced_teleo_resolutive), 
            percent_sequenced = round(n_sequenced / n_sp_tot, 2), 
            percent_sequenced_resolutive = round(n_sequenced_resolutive / n_sp_tot, 2), 
            percent_sequenced_resolutive_among_sequenced = case_when(
              n_sequenced != 0 & n_sequenced_resolutive == 0 ~ round((n_sequenced_resolutive+1) / n_sequenced, 2), 
              n_sequenced != 0 & n_sequenced_resolutive != 0 ~ round((n_sequenced_resolutive) / n_sequenced, 2), 
              n_sequenced == 0 ~ 0
            ))

# Prepare the resume table 
family_reference_coef <- family_stat %>%
  rename(coef_sequencing = percent_sequenced, 
         coef_resolution = percent_sequenced_resolutive_among_sequenced) %>%
  select(Family, coef_sequencing, coef_resolution)

# Among the families present in the eDNA dataset
write.csv(family_reference_coef, "outputs/01_read_data_stats/family_resolution_coefs.csv", row.names = F)




