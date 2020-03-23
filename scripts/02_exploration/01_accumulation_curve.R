# Exploration script to vizualize accumulation curves from all projects
# Virginie Marques
# Last up: 27.02.2020

# To do:
# Remove the deep water filters. Question: even in Lengguru?

# Lib 
library(tidyverse)
library(ggplot2)

# data
load("Rdata/02_clean_all.Rdata")

# 
'%ni%' <- Negate("%in%")

# Functions
source('scripts/02_exploration/00_functions.R')

# Ajouter une colonne nom assigné unique lié aux MOTUs? Pour éviter d avoir des noms de colonnes avec des séquences?

# First, remove the assignations above family (we can remove it in the 02_read script if necessary)
liste_read_edna_LULU <- lapply(liste_read_edna_LULU, function(x){
  x %>%
    filter(new_rank_ncbi != "higher")
})


lapply(liste_read_edna_LULU, function(x){
  length(unique(x$amplicon))
})

lapply(liste_read_edna_LULU, function(x){
  length(unique(x$station))
})

leg <- df_all_filters %>%
  filter(project_name == "Lengguru") %>%
  distinct(amplicon, sequence, new_rank_ncbi)


table(leg$new_rank_ncbi)

length(unique(leg$sequence))

# ------------------------------------------------------------------------------- # 
# On individual filters
# ------------------------------------------------------------------------------- # 

# accumlation all plots
liste_accumulation <- lapply(liste_read_edna_LULU, accumulation_curve_df)

# Asymptote of all plots 
liste_asymptote <- lapply(liste_read_edna_LULU, asymptote_mm)

# Unlist
df_accumulation <- bind_rows(liste_accumulation, .id = "project_name")
df_asymptote <- bind_rows(liste_asymptote, .id = "project_name")

# Dataset for plot
df_accumulation_all <- left_join(df_accumulation, df_asymptote, by = "project_name") %>%
  group_by(project_name) %>%
  mutate(position_asymptote_y = 0.20 * max(asymptote), 
         position_asymptote_x = max(sites))

# Plots
plot_acc_all <- ggplot(df_accumulation_all) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd), fill = "lightblue", alpha = 0.8) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 1, col = "lightblue") +
  geom_text(aes(x = position_asymptote_x, y =position_asymptote_y, hjust = 1, label = paste("asymptote =", round(asymptote, 1), "MOTUs")), col = "black") +
  facet_wrap(~project_name, scales = "free") +
  ylab("Number of MOTUs") +
  xlab("Samples (filter)") +
  theme_bw()

plot_acc_all

# Save
ggsave("outputs/03_accumulation_curves/01_accumulation_curve_all_projects_filters.png", plot_acc_all, width = 12, height = 8)
  


# ------------------------------------------------------------------------------- # 
# On individual sites (couple of filters)
# ------------------------------------------------------------------------------- # 

# accumlation all plots
# We have a warning?
liste_accumulation <- lapply(liste_read_edna_LULU, accumulation_curve_df, column_station = "station")

# Asymptote of all plots 
liste_asymptote <- lapply(liste_read_edna_LULU, asymptote_mm)

# Unlist
df_accumulation <- bind_rows(liste_accumulation, .id = "project_name")
df_asymptote <- bind_rows(liste_asymptote, .id = "project_name")

# Dataset for plot
df_accumulation_all <- left_join(df_accumulation, df_asymptote, by = "project_name") %>%
  group_by(project_name) %>%
  mutate(position_asymptote_y = 0.20 * max(asymptote), 
         position_asymptote_x = max(sites))

# Plots
plot_acc_all <- ggplot(df_accumulation_all) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd), fill = "lightgreen", alpha = 0.8) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 1, col = "lightgreen") +
  geom_text(aes(x = position_asymptote_x, y =position_asymptote_y, hjust = 1, label = paste("asymptote =", round(asymptote, 1), "MOTUs")), col = "black") +
  facet_wrap(~project_name, scales = "free") +
  ylab("Number of MOTUs") +
  xlab("Station (pair of filters") +
  theme_bw()

plot_acc_all

# Save
ggsave("outputs/03_accumulation_curves/01_accumulation_curve_all_projects_station.png", plot_acc_all, width = 12, height = 8)


# ----------------------------------------------------------------------------- # 
# All samples together 

# Add a global saturation curve, i.e. all samples together?
all_accumulation <- accumulation_curve_df(df_all_filters) %>%
  mutate(project_name = "All") %>%
  select(project_name, richness, sd, sites)

# Asymptote of all plots 
all_asymptote <- asymptote_mm(df_all_filters) %>%
  mutate(project_name = "All") %>%
  select(project_name, asymptote)

# Bind together
df_all_accumulation <- rbind(df_accumulation, all_accumulation)
df_all_asymptote <- rbind(df_asymptote, all_asymptote)

# 
df_join_all <- df_all_accumulation %>%
  left_join(., df_all_asymptote, by = "project_name") 

# Plots
plot_acc_all <- ggplot(df_join_all 
                       %>% filter(project_name != "All"), 
                       aes(fill = project_name, col = project_name)) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd), alpha = 0.5) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote, col = project_name), linetype = "dashed", size = 0.5) +
  ylab("Number of MOTUs") +
  xlab("Samples (filter)") +
  theme_bw() # + scale_x_continuous(trans='log2') 

plot_acc_all

ggsave("outputs/03_accumulation_curves/01_accumulation_curve_all_projects_combination_no_facet.png", plot_acc_all, width = 12, height = 8)

# Plot with facet
plot_acc_all <- ggplot(df_join_all, aes(fill = project_name, col = project_name)) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.5) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 1) +
  facet_wrap(~project_name, scales = "free") +
  ylab("Number of MOTUs") +
  xlab("Samples (filter)") +
  theme_bw()

plot_acc_all

ggsave("outputs/03_accumulation_curves/01_accumulation_curve_all_projects_combination.png", plot_acc_all, width = 12, height = 8)

# ----------------------------------------------------------------------------- # 
# Counts

# Attention a la sur-estimations! 
# Des erreurs peuvent se glisser dans un jeu de données, et rester plus loin et participent a sur-evaluer la diversite gamma
# Approche par nombre de hill? q=1 ou q = 0.5?

otu_names <- df_all_filters %>%
  distinct(sequence, new_rank_ncbi, new_scientific_name_ncbi, best_identity_database)

proportions <- table(otu_names$new_rank_ncbi)  
pie(proportions)  

hist(otu_names$best_identity_database, breaks=100)

otu_names_family <- df_all_filters %>%
  distinct(sequence, new_rank_ncbi, new_family_name, new_scientific_name_ncbi, best_identity_database)


count_families <- data.frame(table(otu_names_family$new_family_name))

ggplot(count_families, aes(x=reorder(Var1, Freq), y = Freq, fill = Freq)) + 
  geom_bar(stat="identity") + 
  theme_bw() +
  theme(axis.text.x=element_text(angle = 0, hjust = 0)) + 
  coord_flip()

ggsave("outputs/03_accumulation_curves/01_number_motus_family.png", width=6, height=16)

# Leng
leng <- df_all_filters %>%
  filter(project_name == "Lengguru") %>%
  distinct(sequence, new_rank_ncbi, new_family_name, new_scientific_name_ncbi, best_identity_database)


# Most of the gobies are from Lengguru, i.e. 66 sp among the 95
