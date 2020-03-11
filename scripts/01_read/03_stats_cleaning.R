# R script
# Calculating the number of MOTUs at each step of the cleaning for each project / all the projects 

# Lib
library(tidyverse)
library(conflicted)

# Set preference
conflict_prefer("filter", "dplyr")

# Load data
load('Rdata/01_liste_all_read_edna.Rdata')
load("Rdata/02_clean_all.Rdata")

# Load functions
source("scripts/01_read/00_functions.R")

# --------------------------------------------------------- # 
####             Count by project                     ----
# --------------------------------------------------------- # 

# ---- # Before & unlist

before <- lapply(liste_read_edna, clean_data, remove_blanks = FALSE, min_reads=0, remove_chimera=FALSE, remove_not_fish_manual=FALSE, remove_not_fish_taxize = FALSE,
                   min_size_seq = 0, max_size_seq = 500, tag_jump=FALSE, tag_jump_value = 0.001, min_PCR = 0,
                   min_PCR_sample = 0, habitat_select = c("marine"), min_percentage_id = 0, delete_gps_col=TRUE) %>%
  lapply(., infos_statistiques) %>%
  bind_rows(., .id = "project_name") %>%
  mutate(step = "before") %>% 
  select(step, colnames(.))

# ---- # 10 reads

tenreads <- lapply(liste_read_edna, clean_data, remove_blanks = FALSE, min_reads=10, remove_chimera=FALSE, remove_not_fish_manual=FALSE, remove_not_fish_taxize = FALSE,
                   min_size_seq = 0, max_size_seq = 500, tag_jump=FALSE, tag_jump_value = 0.001, min_PCR = 0,
                   min_PCR_sample = 0, habitat_select = c("marine"), min_percentage_id = 0, delete_gps_col=TRUE) %>%
  lapply(., infos_statistiques) %>%
  bind_rows(., .id = "project_name") %>%
  mutate(step = "tenreads") %>% 
  select(step, colnames(.))

# ---- # PCR blanks & tag-jump

blanks <- lapply(liste_read_edna, clean_data, remove_blanks = TRUE, min_reads=10, remove_chimera=FALSE, remove_not_fish_manual=FALSE, remove_not_fish_taxize = FALSE,
                   min_size_seq = 0, max_size_seq = 500, tag_jump=TRUE, tag_jump_value = 0.001, min_PCR = 0,
                   min_PCR_sample = 0, habitat_select = c("marine"), min_percentage_id = 0, delete_gps_col=TRUE) %>%
  lapply(., infos_statistiques) %>%
  bind_rows(., .id = "project_name") %>%
  mutate(step = "blanks") %>% 
  select(step, colnames(.))

# ---- # Fishes only 

fishonly <- lapply(liste_read_edna, clean_data, remove_blanks = TRUE, min_reads=10, remove_chimera=FALSE, remove_not_fish_manual=FALSE, remove_not_fish_taxize = TRUE,
                   min_size_seq = 0, max_size_seq = 500, tag_jump=TRUE, tag_jump_value = 0.001, min_PCR = 0,
                   min_PCR_sample = 0, habitat_select = c("marine"), min_percentage_id = 0, delete_gps_col=TRUE) %>%
  lapply(., infos_statistiques) %>%
  bind_rows(., .id = "project_name") %>%
  mutate(step = "fishonly") %>% 
  select(step, colnames(.))

# ---- # Length: A ne plus faire ? Inutile ? 

readlength <- lapply(liste_read_edna, clean_data, remove_blanks = TRUE, min_reads=10, remove_chimera=FALSE, remove_not_fish_manual=FALSE, remove_not_fish_taxize = TRUE,
                   min_size_seq = 30, max_size_seq = 90, tag_jump=TRUE, tag_jump_value = 0.001, min_PCR = 0,
                   min_PCR_sample = 0, habitat_select = c("marine"), min_percentage_id = 0, delete_gps_col=TRUE) %>%
  lapply(., infos_statistiques) %>%
  bind_rows(., .id = "project_name") %>%
  mutate(step = "readlength") %>% 
  select(step, colnames(.))

# ---- # PCRfilter: sur tout le jeu de données 

PCR_all <- lapply(liste_read_edna, clean_data, remove_blanks = TRUE, min_reads=10, remove_chimera=FALSE, remove_not_fish_manual=FALSE, remove_not_fish_taxize = TRUE,
                     min_size_seq = 30, max_size_seq = 90, tag_jump=TRUE, tag_jump_value = 0.001, min_PCR = 1,
                     min_PCR_sample = 0, habitat_select = c("marine"), min_percentage_id = 0, delete_gps_col=TRUE) %>%
  lapply(., infos_statistiques) %>%
  bind_rows(., .id = "project_name") %>%
  mutate(step = "PCR_all") %>% 
  select(step, colnames(.))

# ---- # PCRfilter: sur une paire de filtre (1/24 PR). This is the alternative of the previous filter, not to use at this scale (otherwise: deletion of too many true observations)

PCR_sample <- lapply(liste_read_edna, clean_data, remove_blanks = TRUE, min_reads=10, remove_chimera=FALSE, remove_not_fish_manual=FALSE, remove_not_fish_taxize = TRUE,
                         min_size_seq = 30, max_size_seq = 90, tag_jump=TRUE, tag_jump_value = 0.001, min_PCR = 0,
                         min_PCR_sample = 1, habitat_select = c("marine"), min_percentage_id = 0, delete_gps_col=TRUE) %>%
  lapply(., infos_statistiques) %>%
  bind_rows(., .id = "project_name") %>%
  mutate(step = "PCR_sample") %>% 
  select(step, colnames(.))

# ---- # LULU  

LULU <- liste_read_edna_LULU %>%
  lapply(., infos_statistiques) %>%
  bind_rows(., .id = "project_name") %>%
  mutate(step = "LULU") %>% 
  select(step, colnames(.))

# ---- # LULU family (no assignation above family rank)

LULU_family <- lapply(liste_read_edna_LULU, function(x){
  x %>% 
    filter(new_rank_ncbi != "higher")}) %>%
  lapply(., infos_statistiques) %>%
  bind_rows(., .id = "project_name") %>%
  mutate(step = "LULU_family") %>% 
  select(step, colnames(.))


# ---------------------- # 
#  REGROUP
# ---------------------- # 

# Combine
stat_by_project <- rbind(before, tenreads, blanks, fishonly, readlength, PCR_all, LULU, LULU_family) %>%
  arrange(project_name) %>%
  select(project_name, step, colnames(.))

# Write
write.csv(stat_by_project, "outputs/01_read_data/stats_by_project.csv", row.names = FALSE)


# --------------------------------------------------------- # 
####        Count combining all projects              ----
# --------------------------------------------------------- # 

# Before
before <- lapply(liste_read_edna, clean_data, remove_blanks = FALSE, min_reads=0, remove_chimera=FALSE, remove_not_fish_manual=FALSE, remove_not_fish_taxize = FALSE,
                 min_size_seq = 0, max_size_seq = 500, tag_jump=FALSE, tag_jump_value = 0.001, min_PCR = 0,
                 min_PCR_sample = 0, habitat_select = c("marine"), min_percentage_id = 0, delete_gps_col=TRUE) %>%
  bind_rows(.) %>%
  infos_statistiques() %>%
  mutate(step = "before") %>% 
  select(step, colnames(.))

# ---- # 10 reads

tenreads <- lapply(liste_read_edna, clean_data, remove_blanks = FALSE, min_reads=10, remove_chimera=FALSE, remove_not_fish_manual=FALSE, remove_not_fish_taxize = FALSE,
                   min_size_seq = 0, max_size_seq = 500, tag_jump=FALSE, tag_jump_value = 0.001, min_PCR = 0,
                   min_PCR_sample = 0, habitat_select = c("marine"), min_percentage_id = 0, delete_gps_col=TRUE) %>%
  bind_rows(.) %>%
  infos_statistiques() %>%
  mutate(step = "tenreads") %>% 
  select(step, colnames(.))

# ---- # PCR blanks & tag-jump

blanks <- lapply(liste_read_edna, clean_data, remove_blanks = TRUE, min_reads=10, remove_chimera=FALSE, remove_not_fish_manual=FALSE, remove_not_fish_taxize = FALSE,
                 min_size_seq = 0, max_size_seq = 500, tag_jump=TRUE, tag_jump_value = 0.001, min_PCR = 0,
                 min_PCR_sample = 0, habitat_select = c("marine"), min_percentage_id = 0, delete_gps_col=TRUE) %>%
  bind_rows(.) %>%
  infos_statistiques() %>%
  mutate(step = "blanks") %>% 
  select(step, colnames(.))

# ---- # Fishes only 

fishonly <- lapply(liste_read_edna, clean_data, remove_blanks = TRUE, min_reads=10, remove_chimera=FALSE, remove_not_fish_manual=FALSE, remove_not_fish_taxize = TRUE,
                   min_size_seq = 0, max_size_seq = 500, tag_jump=TRUE, tag_jump_value = 0.001, min_PCR = 0,
                   min_PCR_sample = 0, habitat_select = c("marine"), min_percentage_id = 0, delete_gps_col=TRUE) %>%
  bind_rows(.) %>%
  infos_statistiques() %>%
  mutate(step = "fishonly") %>% 
  select(step, colnames(.))

# ---- # Length: A ne plus faire ? Inutile ? 

readlength <- lapply(liste_read_edna, clean_data, remove_blanks = TRUE, min_reads=10, remove_chimera=FALSE, remove_not_fish_manual=FALSE, remove_not_fish_taxize = TRUE,
                     min_size_seq = 30, max_size_seq = 90, tag_jump=TRUE, tag_jump_value = 0.001, min_PCR = 0,
                     min_PCR_sample = 0, habitat_select = c("marine"), min_percentage_id = 0, delete_gps_col=TRUE) %>%
  bind_rows(.) %>%
  infos_statistiques() %>%
  mutate(step = "readlength") %>% 
  select(step, colnames(.))

# ---- # PCRfilter: sur tout le jeu de données 

PCR_all <- lapply(liste_read_edna, clean_data, remove_blanks = TRUE, min_reads=10, remove_chimera=FALSE, remove_not_fish_manual=FALSE, remove_not_fish_taxize = TRUE,
                  min_size_seq = 30, max_size_seq = 90, tag_jump=TRUE, tag_jump_value = 0.001, min_PCR = 1,
                  min_PCR_sample = 0, habitat_select = c("marine"), min_percentage_id = 0, delete_gps_col=TRUE) %>%
  bind_rows(.) %>%
  infos_statistiques() %>%
  mutate(step = "PCR_all") %>% 
  select(step, colnames(.))

# ---- # PCRfilter: sur une paire de filtre (1/24 PR). This is the alternative of the previous filter, not to use at this scale (otherwise: deletion of too many true observations)

PCR_sample <- lapply(liste_read_edna, clean_data, remove_blanks = TRUE, min_reads=10, remove_chimera=FALSE, remove_not_fish_manual=FALSE, remove_not_fish_taxize = TRUE,
                     min_size_seq = 30, max_size_seq = 90, tag_jump=TRUE, tag_jump_value = 0.001, min_PCR = 0,
                     min_PCR_sample = 1, habitat_select = c("marine"), min_percentage_id = 0, delete_gps_col=TRUE) %>%
  bind_rows(.) %>%
  infos_statistiques() %>%
  mutate(step = "PCR_sample") %>% 
  select(step, colnames(.))

# ---- # LULU  

LULU <- liste_read_edna_LULU %>%
  bind_rows(.) %>%
  infos_statistiques() %>%
  mutate(step = "LULU") %>% 
  select(step, colnames(.))

# ---- # LULU family (no assignation above family rank)

LULU_family <- lapply(liste_read_edna_LULU, function(x){
  x %>% 
    filter(new_rank_ncbi != "higher")}) %>%
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
write.csv(stat_all_project, "outputs/01_read_data/stats_all_project.csv", row.names = FALSE)





