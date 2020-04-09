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


## calculate indice of resolution per family

family_tot <- unique(class_all_sp_family_class$family)

n_total <- data.frame(family=character(451), n_total=numeric(451), stringsAsFactors = FALSE)

for (i in 1:length(family_tot)) {
  f <- family_tot[[i]]
  tot <- class_all_sp_family_class %>%
    filter(family==family_tot[i]) %>%
    nrow
  n_total[i,1] <- f
  n_total[i,2] <- tot
  
}


family_res <- unique(class_resolutive_sp_family_class$family)
n_res <- data.frame(family=character(441), n_res=numeric(441), stringsAsFactors = FALSE)

for (i in 1:length(family_res)) {
  f <- family_res[[i]]
  res <- class_resolutive_sp_family_class %>%
    filter(family==family_res[i]) %>%
    nrow
  n_res[i,1] <- f
  n_res[i,2] <- res
}

resolution <- left_join(n_total, n_res, by="family")
resolution$resolution <- resolution$n_res / resolution$n_total
resolution <- resolution[order(resolution$resolution, decreasing = TRUE),]
resolution[is.na(resolution)] <- 0
write.csv(resolution, "outputs/01_read_data_stats/family_resolution.csv")
