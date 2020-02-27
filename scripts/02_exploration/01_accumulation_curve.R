# Exploration script to vizualize accumulation curves from all projects
# Virginie Marques
# Last up: 27.02.2020

# Lib 
library(tidyverse)
library(ggplot2)

# data
load("Rdata/02_clean_all.Rdata")

# 
'%ni%' <- Negate("%in%")

# Ajouter une colonne nom assigné unique lié aux MOTUs? Pour éviter d avoir des noms de colonnes avec des séquences?

# First, remove the assignations above family (we can remove it in the 02_read script if necessary)
liste_read_edna_LULU <- lapply(liste_read_edna_LULU, function(x){
  x %>%
    filter(new_rank_ncbi != "higher")
})

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
ggsave("plots/01_exploration/01_accumulation_curve_all_projects_filters.png", plot_acc_all, width = 12, height = 8)
  


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
  xlab("Samples (filter)") +
  theme_bw()

plot_acc_all

# Save
ggsave("plots/01_exploration/01_accumulation_curve_all_projects_station.png", plot_acc_all, width = 12, height = 8)








