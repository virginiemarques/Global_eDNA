# This is code to control all quality thresholds applied on the ouput of 01 script + apply LULU + combine list into dataframe
# Code developed by Virginie Marques
# Last updated: 28 Jan. 2020

# Lib 
library(tidyverse)
library(conflicted)
library(lulu)

# Source functions
source('scripts/00_functions.R')

# Solve the conflicts with dplyr
conflict_prefer("filter", "dplyr")

# Load data
load('Rdata/01_liste_all_read_edna.Rdata')

# ------------------------------------------------------------------------------------- # 
# Cleaning thresholds 

# Apply all the default filters - notably all the marine filters - this excludes filters from resurgences, estuaries and rivers
liste_all_filters <- lapply(liste_read_edna, clean_data)

# Cout the number of MOTUs before LULU
lapply(liste_all_filters, function(x){
  length(unique(x$amplicon))
})

# ------------------------------------------------------------------------------------- # 
# Apply LULU to clean post-clustering

# Set the path 
path_lulu <- "data/LULU/"

# To apply this function, you need to have a UNIX OS system (The system function doesn't work on windows I think)
# And you need to install the blastn tools in your local machine 

liste_read_edna_LULU <- lapply(liste_all_filters, function(x){
  # Apply LULU
  list_LULU <- apply_lulu(x, path_lulu)
  # MOTUs to keep
  lulu_motus_keep <- list_LULU$curated_otus
  # Filter out discarded MOTUs
  x_filtered <- x %>%
    filter(amplicon %in% lulu_motus_keep)
  # End 
  return(x_filtered)
})

# Cout the number of MOTUs after LULU
lapply(liste_read_edna_LULU, function(x){
  length(unique(x$amplicon))
})


# ------------------------------------------------------------------------------------- # 
# Assemble the lists to obtain the big dataset
# simplify_sample_level simplifies the dataset at the sample level and not at the PCR level as before

# Combine list to one big dataset
df_all_filters <- do.call("rbind", liste_read_edna_LULU) 

# Simplify at the sample level, instead of PCR level?
df_all_filters_sample <- simplify_sample_level(df_all_filters)

# ------------------------------------------------------------------------------------- # 
# Save results
# List + df 

save(liste_all_filters, liste_read_edna_LULU, df_all_filters, file = 'Rdata/02_clean_all.Rdata')
