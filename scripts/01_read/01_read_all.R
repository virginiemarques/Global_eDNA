# This is code to read and combine all data from available eDNA metabarcoding outputs from SWARM clustering.
# Code developed by Virginie Marques
# Last updated: 24 Fev 2020

# Bugs to fix:
# 

# ----------------------------------------------------------------------------------------------------------------- # 
# Before running this script, if you add new metabarcoding data, make sure you have all the corresponding metadata 
# Otherwise, the script will fail
# ----------------------------------------------------------------------------------------------------------------- # 

# Lib 
library(tidyverse)
library(conflicted)
library(readxl)

# Source functions
source('scripts/00_functions.R')

# Solve the conflicts with dplyr
conflict_prefer("filter", "dplyr")

# Get the file names
liste_table <- list.files(path = "data/All/", pattern = ".table")
liste_taxo <- list.files(path = "data/All/", pattern = ".csv")

# Check 
if(length(liste_table) != length(liste_taxo)){
  message("Wrong number of input files, check in folder if all files are present without duplicates.")
} 

# Get the project name 
projects <- unique(word(liste_table, 1, sep="_"))

# Metadata - field. Import the csv; the excel files messes with dates and dates
# There is an issue with numerous empty columns being imported in the csv 
# metadata_sampling <- read_excel("metadata/Metadata_eDNA_Megafauna_EB_leng.xlsx", 1)
metadata_sampling <- read.csv("metadata/Metadata_eDNA_global_V3.csv", sep=";", stringsAsFactors = F)

# Clean the spaces before coordinates
metadata_sampling$longitude_start <- str_trim(metadata_sampling$longitude_start)
metadata_sampling$latitude_start <- str_trim(metadata_sampling$latitude_start)

# Metadata - sequencing
# For later: the generation of those files should be automated to always have the same columns?
# For later: add a column project to ease the check?
liste_run_metadata <- list.files(path = "data/metadata_sequencing/", pattern = ".table")

# Check metadata completness
if(length(projects) != length(liste_run_metadata)){
  message("Wrong number of metadata files, check in folder if all files are present.")
} 

# Open and clean
all.the.data <- lapply(paste("data/metadata_sequencing/", liste_run_metadata, sep=""),  read.table, header=TRUE, stringsAsFactors = F, sep=" ")

# Rbind
metadata_sequencing <- do.call("rbind", all.the.data) # All the data obitools

# Loop to open and store all data in a list 
# It is also possible to store each project in a distinct dataframe using 'assign' along the loop 
# Or transform for loop to apply 

liste_read_edna <- list()

# Loop to open all files 
for (i in 1:length(projects)){
  
  # For debug
  # i = 2
  
  # Get project name
  project_i <- projects[[i]]
  
  # Open files corresponding to said project
  table_otu <- read.csv(paste0('data/All/', liste_table[grepl(project_i, liste_table)]), h=T, sep="\t", stringsAsFactors = F)
  taxo_motu <- read.csv(paste0('data/All/', liste_taxo[grepl(project_i, liste_taxo)]), h=T, sep="\t", stringsAsFactors = F)
  
  # Assemble data - some projects have a long computation time (i.e. Lengguru)
  sw4 <- read_data(taxo_motu, table_otu, metadata_sampling, metadata_sequencing)
  
  # Add project name
  sw4 <- sw4 %>%
    mutate(project_name = project_i)
  
  # Add the class - can be slow as well - 15 min for 1200 queries (i.e. Lengguru)
  # You can choose no to do this step and just store sw4 in the list -- but that means you cannot use it for the fish filter on the next step
  sw4_with_class <- add_class_name(sw4)
  
  # Store in list 
  liste_read_edna[[i]] <- sw4_with_class
  
  # Follow progress
  print(i)
  
}

# Name the list elements
names(liste_read_edna) <- projects

# Save 
save(liste_read_edna, file="Rdata/01_liste_all_read_edna.Rdata")








