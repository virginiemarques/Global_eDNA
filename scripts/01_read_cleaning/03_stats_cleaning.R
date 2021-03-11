# R script
# Calculating the number of MOTUs at each step of the cleaning for each project / all the projects 

# ATTENTION! A RECODER
# POUR PRENDRE EN COMPTE LES FILTRES A L ECHELLE DES PROJETS ET NON PAS DES REGIONS
# POUR REFLETER CE QUI A ETE FAIT DANS LE MS PLUS FIDELEMENT (SAUF SI ON CHANGE LA FACON DE CALCULER POUR LE MS, C'EST AUSSI POSSIBLE; notamment pour les PCR)

# Lib
library(tidyverse)
library(conflicted)
library(purrr)

# Set preference
conflict_prefer("filter", "dplyr")

# Load data
#load('Rdata/01_liste_all_read_edna.Rdata')
#load("Rdata/02_clean_all.Rdata")

load("Rdata/01_read_data.Rdata")
load("Rdata/02-clean-data.Rdata")

## Join metadata
columns_delete_field_metadata <- c("turbidity", "gps_start", "gps_b", "lat_gps_b", "long_gps_b", "gps_c", "long_gps_c", "lat_gps_d", "gps_half_turn", "longitude_turn", "latitude_end", "longitude_end", 
                                   "gps_end", "long_gps_d", "gps_d", "lat_gps_c", "latitude_turn", "data_manager", "gps_owner", "project")

metadata_field <- read.csv("metadata/Metadata_eDNA_global_V6.csv", sep=";", stringsAsFactors = F)
metadata_field <- select(metadata_field, -c(columns_delete_field_metadata))

# Load functions
source("scripts/01_read_cleaning/00_functions.R")
"%ni%" <- Negate("%in%")

# Function
fct_filter_station <- function(x){
  x <- x %>%
    # Remove non-marine Colombia
    filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")) %>%
    filter(sample_method !="niskin" & province!="East_Pacific" & comment %ni% c("Distance decay 600m", "Distance decay 300m") & station!="glorieuse_distance_300m") %>%
    # Remove for NC
    filter(project != "SEAMOUNTS") %>% 
    filter(habitat_type %ni% c("BAIE", "Sommet")) %>%
    # Unused projects
    filter(project_i %ni% c("atlantique", "Med", "curacao", "malpelo"))
  
  return(x)
}

# ------------------------------- # 
# Reformate the lists by province 
# ------------------------------- # 

# ----- # After LULU 
df_all_filters_temp <- do.call("rbind", list_read_step4) %>%
  filter(province != "Tropical_East_Pacific") %>%
  # Filter marine 
  filter(habitat == "marine") %>%
  # Unused projects
  filter(project_i %ni% c("atlantique", "Med", "curacao", "malpelo"))

unique(df_all_filters_temp$project_i)
unique(df_all_filters_temp$habitat)
unique(df_all_filters_temp$province)

# Split by province
liste_read_edna_LULU <- split(df_all_filters_temp, df_all_filters_temp$province)

# -----# Before everything 
df_all_filters_temp <- do.call("rbind", list_read_step2) %>%
  mutate(sample_name_all_pcr = substr(sample_name,1, nchar(sample_name) -3)) %>%
  left_join(., metadata_field, by=c("sample_name_all_pcr" = "code_spygen")) %>%
  filter(province != "East_Pacific") %>%
  # Filter marine 
  filter(habitat == "marine") %>%
  # Unused projects
  filter(project_i %ni% c("atlantique", "Med", "curacao", "malpelo"))

unique(df_all_filters_temp$project_i)
unique(df_all_filters_temp$habitat)
unique(df_all_filters_temp$province)
# Split by province
liste_read_edna <- split(df_all_filters_temp, df_all_filters_temp$province)

# ----- # After cleaning but before LULU
df_all_filters_temp <- do.call("rbind", list_read_step3) %>%
  mutate(sample_name_all_pcr = substr(sample_name,1, nchar(sample_name) -3)) %>%
  left_join(., metadata_field, by=c("sample_name_all_pcr" = "code_spygen")) %>%
  # Filter marine 
  filter(habitat == "marine") %>%
  # Unused projects
  filter(project_i %ni% c("atlantique", "Med", "curacao", "malpelo"))

unique(df_all_filters_temp$project_i)
unique(df_all_filters_temp$habitat)
unique(df_all_filters_temp$province)


# Split by province
liste_read_edna_cleaned <- split(df_all_filters_temp, df_all_filters_temp$province)


# --------------------------------------------------------- # 
####             Count by province                     ----
# --------------------------------------------------------- #

# ---- # Before & unlist

before <- lapply(liste_read_edna, clean_motus, min_reads=0, remove_chimera=FALSE, remove_not_fish_taxize = FALSE,
                 min_size_seq = 0, max_size_seq = 500, min_PCR = 0, min_percentage_id = 0, remove_PCR_blanks = FALSE) %>%
  # Select first element of list 
  map(., 1) %>%
  # filter stations
  lapply(., fct_filter_station) %>%
  # Stat
  lapply(., infos_statistiques) %>%
  bind_rows(., .id = "province") %>%
  mutate(step = "before") %>% 
  select(step, colnames(.))


# ---- # 10 reads

tenreads <- lapply(liste_read_edna, clean_motus, min_reads=10, remove_chimera=FALSE, remove_not_fish_taxize = FALSE,
                   min_size_seq = 0, max_size_seq = 500, min_PCR = 0, min_percentage_id = 0, remove_PCR_blanks = FALSE) %>%
  # Select first element of list 
  map(., 1) %>%
  # filter stations
  lapply(., fct_filter_station) %>%
  # Stat
  lapply(., infos_statistiques) %>%
  bind_rows(., .id = "province") %>%
  mutate(step = "tenreads") %>% 
  select(step, colnames(.))


# ---- # PCR blanks & chimeras

blanks <- lapply(liste_read_edna, clean_motus, min_reads=10, remove_chimera=TRUE, remove_not_fish_taxize = FALSE,
                 min_size_seq = 0, max_size_seq = 500, min_PCR = 0, min_percentage_id = 0, remove_PCR_blanks = TRUE) %>%
  # Select first element of list 
  map(., 1) %>%
  # filter stations
  lapply(., fct_filter_station) %>%
  # Stat
  lapply(., infos_statistiques) %>%
  bind_rows(., .id = "province") %>%
  mutate(step = "PCR_blanks_chimeras") %>% 
  select(step, colnames(.))

# ---- # Fishes only 

fishonly <- lapply(liste_read_edna, clean_motus, min_reads=10, remove_chimera=TRUE, remove_not_fish_taxize = TRUE,
                   min_size_seq = 0, max_size_seq = 500, min_PCR = 0, min_percentage_id = 0, remove_PCR_blanks = TRUE) %>%
  # Select first element of list 
  map(., 1) %>%
  # filter stations
  lapply(., fct_filter_station) %>%
  # Stat
  lapply(., infos_statistiques) %>%
  bind_rows(., .id = "province") %>%
  mutate(step = "fishonly") %>% 
  select(step, colnames(.))


# ---- # Length

readlength <- lapply(liste_read_edna, clean_motus, min_reads=10, remove_chimera=TRUE, remove_not_fish_taxize = TRUE,
                     min_size_seq = 30, max_size_seq = 100, min_PCR = 0, min_percentage_id = 0, remove_PCR_blanks = TRUE) %>%
  # Select first element of list 
  map(., 1) %>%
  # filter stations
  lapply(., fct_filter_station) %>%
  # Stat
  lapply(., infos_statistiques) %>%
  bind_rows(., .id = "province") %>%
  mutate(step = "readlength") %>% 
  select(step, colnames(.))


# ---- # PCRfilter: on the whole dataset 

PCR_all <- liste_read_edna_cleaned %>%
  # filter stations
  lapply(., fct_filter_station) %>%
  # Stat
  lapply(., infos_statistiques) %>%
  bind_rows(., .id = "province") %>%
  mutate(step = "PCR_all") %>% 
  select(step, colnames(.))

# ---- # PCRfilter: on a pair of filters (1/24 PCR). This is the alternative of the previous filter, not to use at this scale (otherwise: deletion of too many true observations)

# I think this is useless

#    PCR_sample <- lapply(liste_read_edna, clean_data, remove_blanks = TRUE, min_reads=10, remove_chimera=TRUE, remove_not_fish_manual=FALSE, remove_not_fish_taxize = TRUE,
#                         min_size_seq = 30, max_size_seq = 90, tag_jump=TRUE, tag_jump_value = 0.001, min_PCR = 0,
#                         min_PCR_sample = 1, habitat_select = c("marine"), min_percentage_id = 0, delete_gps_col=TRUE) %>%
#      # filter stations
#      lapply(., fct_filter_station) %>%
#      # Stat
#      lapply(., infos_statistiques) %>%
#      bind_rows(., .id = "province") %>%
#      mutate(step = "PCR_sample") %>% 
#      select(step, colnames(.))

# ---- # LULU  

LULU <- liste_read_edna_LULU %>%
  lapply(., fct_filter_station) %>%
  lapply(., infos_statistiques) %>%
  bind_rows(., .id = "province") %>%
  mutate(step = "LULU") %>% 
  select(step, colnames(.))

# ---- # LULU family (no assignation above family rank)

LULU_family <- lapply(liste_read_edna_LULU, function(x){
  x <- x %>% 
    filter(rank_ncbi_corrected != "higher")
  return(x)
  }) %>%
  lapply(., fct_filter_station) %>%
  lapply(., infos_statistiques) %>%
  bind_rows(., .id = "province") %>%
  mutate(step = "LULU_family") %>% 
  select(step, colnames(.))

# ---------------------- # 
#  REGROUP
# ---------------------- # 

# Combine
stat_by_project <- rbind(before, tenreads, blanks, fishonly, readlength, PCR_all, LULU, LULU_family) %>%
  arrange(province) %>%
  select(province, step, colnames(.))

# Write
write.csv(stat_by_project, "outputs/01_read_data_stats/stats_by_project.csv", row.names = FALSE)
write.csv(stat_by_project, "outputs/00_Figures_for_paper/Extended_Data/ED_Table3.csv", row.names = FALSE)



# --------------------------------------------------------- # 
####        Count combining all projects              ----
# --------------------------------------------------------- # 

# Before
before <- lapply(liste_read_edna, clean_data, remove_blanks = FALSE, min_reads=0, remove_chimera=FALSE, remove_not_fish_manual=FALSE, remove_not_fish_taxize = FALSE,
                 min_size_seq = 0, max_size_seq = 500, tag_jump=FALSE, tag_jump_value = 0.001, min_PCR = 0,
                 min_PCR_sample = 0, habitat_select = c("marine"), min_percentage_id = 0, delete_gps_col=TRUE) %>%
  # filter stations
  lapply(., fct_filter_station) %>%
  # Stat
  bind_rows(.) %>%
  infos_statistiques() %>%
  mutate(step = "before") %>% 
  select(step, colnames(.))

# ---- # 10 reads

tenreads <- lapply(liste_read_edna, clean_data, remove_blanks = FALSE, min_reads=10, remove_chimera=TRUE, remove_not_fish_manual=FALSE, remove_not_fish_taxize = FALSE,
                   min_size_seq = 0, max_size_seq = 500, tag_jump=FALSE, tag_jump_value = 0.001, min_PCR = 0,
                   min_PCR_sample = 0, habitat_select = c("marine"), min_percentage_id = 0, delete_gps_col=TRUE) %>%
  # filter stations
  lapply(., fct_filter_station) %>%
  # Stat
  bind_rows(.) %>%
  infos_statistiques() %>%
  mutate(step = "tenreads") %>% 
  select(step, colnames(.))

# ---- # PCR blanks & tag-jump

blanks <- lapply(liste_read_edna, clean_data, remove_blanks = TRUE, min_reads=10, remove_chimera=TRUE, remove_not_fish_manual=FALSE, remove_not_fish_taxize = FALSE,
                 min_size_seq = 0, max_size_seq = 500, tag_jump=TRUE, tag_jump_value = 0.001, min_PCR = 0,
                 min_PCR_sample = 0, habitat_select = c("marine"), min_percentage_id = 0, delete_gps_col=TRUE) %>%
  # filter stations
  lapply(., fct_filter_station) %>%
  # Stat
  bind_rows(.) %>%
  infos_statistiques() %>%
  mutate(step = "blanks") %>% 
  select(step, colnames(.))

# ---- # Fishes only 

fishonly <- lapply(liste_read_edna, clean_data, remove_blanks = TRUE, min_reads=10, remove_chimera=TRUE, remove_not_fish_manual=FALSE, remove_not_fish_taxize = TRUE,
                   min_size_seq = 0, max_size_seq = 500, tag_jump=TRUE, tag_jump_value = 0.001, min_PCR = 0,
                   min_PCR_sample = 0, habitat_select = c("marine"), min_percentage_id = 0, delete_gps_col=TRUE) %>%
  # filter stations
  lapply(., fct_filter_station) %>%
  # Stat
  bind_rows(.) %>%
  infos_statistiques() %>%
  mutate(step = "fishonly") %>% 
  select(step, colnames(.))

# ---- # Length 

readlength <- lapply(liste_read_edna, clean_data, remove_blanks = TRUE, min_reads=10, remove_chimera=TRUE, remove_not_fish_manual=FALSE, remove_not_fish_taxize = TRUE,
                     min_size_seq = 30, max_size_seq = 90, tag_jump=TRUE, tag_jump_value = 0.001, min_PCR = 0,
                     min_PCR_sample = 0, habitat_select = c("marine"), min_percentage_id = 0, delete_gps_col=TRUE) %>%
  # filter stations
  lapply(., fct_filter_station) %>%
  # Stat
  bind_rows(.) %>%
  infos_statistiques() %>%
  mutate(step = "readlength") %>% 
  select(step, colnames(.))

# ---- # PCRfilter: on the whole dataset 

PCR_all <- lapply(liste_read_edna, clean_data, remove_blanks = TRUE, min_reads=10, remove_chimera=TRUE, remove_not_fish_manual=FALSE, remove_not_fish_taxize = TRUE,
                  min_size_seq = 30, max_size_seq = 90, tag_jump=TRUE, tag_jump_value = 0.001, min_PCR = 1,
                  min_PCR_sample = 0, habitat_select = c("marine"), min_percentage_id = 0, delete_gps_col=TRUE) %>%
  # filter stations
  lapply(., fct_filter_station) %>%
  # Stat
  bind_rows(.) %>%
  infos_statistiques() %>%
  mutate(step = "PCR_all") %>% 
  select(step, colnames(.))

# ---- # PCRfilter: on a pair of filters (1/24 PCR). This is the alternative of the previous filter, not to use at this scale (otherwise: deletion of too many true observations)

PCR_sample <- lapply(liste_read_edna, clean_data, remove_blanks = TRUE, min_reads=10, remove_chimera=TRUE, remove_not_fish_manual=FALSE, remove_not_fish_taxize = TRUE,
                     min_size_seq = 30, max_size_seq = 90, tag_jump=TRUE, tag_jump_value = 0.001, min_PCR = 0,
                     min_PCR_sample = 1, habitat_select = c("marine"), min_percentage_id = 0, delete_gps_col=TRUE) %>%
  # filter stations
  lapply(., fct_filter_station) %>%
  # Stat
  bind_rows(.) %>%
  infos_statistiques() %>%
  mutate(step = "PCR_sample") %>% 
  select(step, colnames(.))

# ---- # LULU  

LULU <- liste_read_edna_LULU %>%
  # filter stations
  lapply(., fct_filter_station) %>%
  # Stat
  bind_rows(.) %>%
  infos_statistiques() %>%
  mutate(step = "LULU") %>% 
  select(step, colnames(.))

# ---- # LULU family (no assignation above family rank)

LULU_family <- lapply(liste_read_edna_LULU, function(x){
  x %>% 
    filter(new_rank_ncbi != "higher")}) %>%
  # filter stations
  lapply(., fct_filter_station) %>%
  # Stat
  bind_rows(.) %>%
  infos_statistiques() %>%
  mutate(step = "LULU_family") %>% 
  select(step, colnames(.))

# ---------------------- # 
#  REGROUP
# ---------------------- # 

# Combine
stat_all_project <- rbind(before, tenreads, blanks, fishonly, readlength, PCR_all, LULU, LULU_family)

# Write
write.csv(stat_all_project, "outputs/01_read_data_stats/stats_all_project.csv", row.names = FALSE)
write.csv(stat_all_project, "outputs/00_Figures_for_paper/Extended_Data/ED_Table1.csv", row.names = FALSE)





