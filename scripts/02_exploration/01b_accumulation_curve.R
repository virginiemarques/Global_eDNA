# Accumulation curves on genus, family, order 
# Lib 
library(tidyverse)
library(ggplot2)
library(ggpubr)

# data
load("Rdata/02_clean_all.Rdata")

# 
'%ni%' <- Negate("%in%")

# Functions
source('scripts/02_exploration/00_functions.R')

# Ajouter une colonne nom assigné unique lié aux MOTUs? Pour éviter d avoir des noms de colonnes avec des séquences?

# First, remove the assignations above family (we can remove it in the 02_read script if necessary) --> keep it, put the families only in the supl. mat. 
liste_read_edna_LULU <- lapply(liste_read_edna_LULU, function(x){
  x %>%
    #filter(new_rank_ncbi != "higher") %>%
    filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")) %>%
    filter(sample_method !="niskin" & region!="East_Pacific" & comment %ni% c("Distance decay 600m", "Distance decay 300m") & station!="glorieuse_distance_300m")%>%
    filter(project != "SEAMOUNTS") %>% 
    filter(habitat_type %ni% c("BAIE", "Sommet"))
})

# On the df as well
df_all_filters <- df_all_filters %>%
  #filter(new_rank_ncbi != "higher") %>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")) %>%
  filter(sample_method !="niskin" & region!="East_Pacific" & comment %ni% c("Distance decay 600m", "Distance decay 300m") & station!="glorieuse_distance_300m")%>%
  filter(project != "SEAMOUNTS") %>% 
  filter(habitat_type %ni% c("BAIE", "Sommet"))

# Re-format at region scale
# -----# After LULU 
df_all_filters_temp <- do.call("rbind", liste_read_edna_LULU) %>%
  filter(region != "East_Pacific")

# Split by region
liste_read_edna_LULU <- split(df_all_filters_temp, df_all_filters_temp$region)

# Counts
lapply(liste_read_edna_LULU, function(x){
  length(unique(x$amplicon))
})

lapply(liste_read_edna_LULU, function(x){
  length(unique(x$station))
})

lapply(liste_read_edna_LULU, function(x){
  length(unique(x$new_family_name))
})



# ------------------------------------------------------------------------------- # 
#### On family ----
# ------------------------------------------------------------------------------- # 

# rank_specify
rank_choice = 'new_family_name'

# accumlation all plots
liste_accumulation <- lapply(liste_read_edna_LULU, accumulation_curve_df, species_unit = rank_choice)

# Asymptote of all plots 
liste_asymptote <- lapply(liste_read_edna_LULU, asymptote_mm, species_unit = rank_choice)

# Unlist
df_accumulation <- bind_rows(liste_accumulation, .id = "project_name")
df_asymptote <- bind_rows(liste_asymptote, .id = "project_name")

# Dataset for plot
df_accumulation_all <- left_join(df_accumulation, df_asymptote, by = "project_name") %>%
  group_by(project_name) %>%
  mutate(position_asymptote_y = 0.20 * max(asymptote), 
         position_asymptote_x = max(sites))

# Plots
plot_acc_family <- ggplot(df_accumulation_all) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd), fill = "lightblue", alpha = 0.8) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 1, col = "lightblue") +
  geom_text(aes(x = position_asymptote_x, y =position_asymptote_y, hjust = 1, label = paste("asymptote =", round(asymptote, 1), "MOTUs")), col = "black") +
  facet_wrap(~project_name, scales = "free") +
  ylab("Number of MOTUs") +
  xlab("Samples (filter)") +
  theme_bw()

plot_acc_family

# Save
# ggsave("outputs/03_accumulation_curves/01b_accumulation_curve_all_projects_family.png", plot_acc_all, width = 12, height = 8)

# Add All samples
# Add a global saturation curve, i.e. all samples together?
all_accumulation <- accumulation_curve_df(df_all_filters, species_unit = rank_choice) %>%
  mutate(project_name = "All") %>%
  select(project_name, richness, sd, sites)


# Asymptote of all plots 
all_asymptote <- asymptote_mm(df_all_filters, species_unit = rank_choice) %>%
  mutate(project_name = "All") %>%
  select(project_name, asymptote, slope)


# Bind together
df_all_accumulation <- rbind(df_accumulation, all_accumulation)
df_all_asymptote <- rbind(df_asymptote, all_asymptote)

# 
df_join_all <- df_all_accumulation %>%
  left_join(., df_all_asymptote, by = "project_name") %>%
  group_by(project_name) %>%
  mutate(position_asymptote_y = 1.05*asymptote, 
         position_asymptote_x = max(sites),
         position_slope_y = 0.50 * max(asymptote))

# Plots
plot_acc_all <- ggplot(df_join_all, 
                       aes(fill = project_name, col = project_name)) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd), alpha = 0.5) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote, col = project_name), linetype = "dashed", size = 0.5) +
  ylab("Number of MOTUs") +
  xlab("Samples (filter)") +
  theme_bw() # + scale_x_continuous(trans='log2') 

plot_acc_all

# ggsave("outputs/03_accumulation_curves/01_accumulation_curve_all_projects_combination_no_facet_family.png", plot_acc_all, width = 12, height = 8)

# Plot with facet
plot_acc_family <- ggplot(df_join_all, aes(fill = project_name)) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.7) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 1) +
  scale_fill_manual(values=c("#A41D1A", "#E5A729", "#8AAE8A", "#4F4D1D", "#863b34", "#C67052"))+ 
  facet_wrap(~project_name, scales = "free") +
  ylab("Number of families") +
  xlab("Samples (filter)") +
  theme_bw() + 
  ggtitle("Families") + 
  geom_text(aes(x = position_asymptote_x, y =position_asymptote_y, hjust = 1, label = paste(round(asymptote, 1), "Families")), col = "black", size=3)+
  geom_text(aes(x = position_asymptote_x, y =position_slope_y, hjust = 1, label = paste("Lomolino slope=",round(slope, 1))), col = "black", size=3)

plot_acc_family

ggsave("outputs/03_accumulation_curves/01b_accumulation_curve_all_projects_combination_no_facet_family.png", plot_acc_family, width = 12, height = 8)

# Table
stats_family <- df_join_all %>%
  group_by(project_name) %>%
  summarise(richness = max(richness), 
            asymptote = round(max(asymptote), 0)) %>%
  mutate(rank = 'Family')

# Save
write.csv(stats_family, "outputs/03_accumulation_curves/asymptotes_family.csv", row.names = F)

# Simple plot on all 
df_fam <- df_join_all %>% 
  filter(project_name == "All") %>%
  mutate(level = "family")

save(df_fam, file = "Rdata/accumulation_asymptote_families_all.rdata")

family <- ggplot(df_fam, aes(fill = project_name, col = project_name)) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.5) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 1) +
  ylab("Number of families") +
  xlab("Samples (filter)") +
  theme_bw() + 
  ggtitle("Families") + 
  geom_text(aes(x = position_asymptote_x, y =position_asymptote_y, hjust = 1, label = paste("asymptote =", round(asymptote, 0), "Families")), col = "black") 

family



# ------------------------------------------------------------------------------- # 
#### On MOTUs ----
# ------------------------------------------------------------------------------- # 

# rank_specify
rank_choice = 'sequence'

# accumlation all plots
liste_accumulation <- lapply(liste_read_edna_LULU, accumulation_curve_df, species_unit = rank_choice)

# Asymptote of all plots 
liste_asymptote <- lapply(liste_read_edna_LULU, asymptote_mm, species_unit = rank_choice)

# Unlist
df_accumulation <- bind_rows(liste_accumulation, .id = "project_name")
df_asymptote <- bind_rows(liste_asymptote, .id = "project_name")

# Dataset for plot
df_accumulation_all <- left_join(df_accumulation, df_asymptote, by = "project_name") %>%
  group_by(project_name) %>%
  mutate(position_asymptote_y = 0.20 * max(asymptote), 
         position_asymptote_x = max(sites))

# Add All samples
# Add a global saturation curve, i.e. all samples together?
all_accumulation_motu <- accumulation_curve_df(df_all_filters, species_unit = rank_choice) %>%
  mutate(project_name = "All") %>%
  select(project_name, richness, sd, sites)


# Asymptote of all plots 
all_asymptote_motu <- asymptote_mm(df_all_filters, species_unit = rank_choice) %>%
  mutate(project_name = "All") %>%
  select(project_name, asymptote, slope)


# Bind together
df_all_accumulation <- rbind(df_accumulation, all_accumulation_motu)
df_all_asymptote <- rbind(df_asymptote, all_asymptote_motu)

# 
df_join_all <- df_all_accumulation %>%
  left_join(., df_all_asymptote, by = "project_name") %>%
  group_by(project_name) %>%
  mutate(position_asymptote_y = 1.05*asymptote, 
         position_asymptote_x = max(sites),
         position_slope_y = 0.30 * max(asymptote))

# Plot with facet
plot_acc_motus <- ggplot(df_join_all, aes(fill = project_name)) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.7) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 1) +
  facet_wrap(~project_name, scales = "free") +
  scale_fill_manual(values=c("#A41D1A", "#E5A729", "#8AAE8A", "#4F4D1D", "#863b34", "#C67052"))+ 
  ylab("Number of genus") +
  xlab("Samples (filter)") +
  theme_bw() + 
  ggtitle("MOTUs") + 
  geom_text(aes(x = position_asymptote_x, y =position_asymptote_y, hjust = 1, label = paste(round(asymptote, 1), "MOTUs")), col = "black", size=3)+
  geom_text(aes(x = position_asymptote_x, y =position_slope_y, hjust = 1, label = paste("Lomolino slope=",round(slope, 1))), col = "black", size=3)

plot_acc_motus

ggsave("outputs/03_accumulation_curves/01b_accumulation_curve_all_projects_combination_no_facet_motus.png", width = 12, height = 8)

# Table
stats<- df_join_all %>%
  group_by(project_name) %>%
  summarise(richness = max(richness), 
            asymptote = round(max(asymptote), 0)) %>%
  mutate(rank = 'MOTUs')

# Save
write.csv(stats, "outputs/03_accumulation_curves/asymptotes_motus.csv", row.names = F)

# Simple plot on all 
df_motus <- df_join_all %>% 
  filter(project_name == "All") %>%
  mutate(level = "MOTUs")

save(df_motus, file = "Rdata/accumulation_asymptote_motus_all.rdata")

motus <- ggplot(df_motus, aes(fill = project_name, col = project_name)) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.5) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 1) +
  ylab("Number of MOTUs") +
  xlab("Samples (filter)") +
  theme_bw() + 
  ggtitle("MOTUs") + 
  geom_text(aes(x = position_asymptote_x, y =position_asymptote_y, hjust = 1, label = paste("asymptote =", round(asymptote, 0), "MOTUs")), col = "black") 

motus

# --------------------------------------------------------------------- # 
#### Final figure - combine all levels  ----
# --------------------------------------------------------------------- # 


ggarrange(plot_acc_motus, plot_acc_family, nrow=2, labels=c("a", "b"))
ggsave("outputs/Figures papier/ED_figure5.png", width=10, height = 12)















# ------------------------------------------------------------------------------- # 
#### On Order ----
# ------------------------------------------------------------------------------- # 

# rank_specify
rank_choice = 'order_name'

# accumlation all plots
liste_accumulation <- lapply(liste_read_edna_LULU, accumulation_curve_df, species_unit = rank_choice)

# Asymptote of all plots 
liste_asymptote <- lapply(liste_read_edna_LULU, asymptote_mm, species_unit = rank_choice)

# Unlist
df_accumulation <- bind_rows(liste_accumulation, .id = "project_name")
df_asymptote <- bind_rows(liste_asymptote, .id = "project_name")

# Dataset for plot
df_accumulation_all <- left_join(df_accumulation, df_asymptote, by = "project_name") %>%
  group_by(project_name) %>%
  mutate(position_asymptote_y = 0.20 * max(asymptote), 
         position_asymptote_x = max(sites))

# Add All samples
# Add a global saturation curve, i.e. all samples together?
all_accumulation <- accumulation_curve_df(df_all_filters, species_unit = rank_choice) %>%
  mutate(project_name = "All") %>%
  select(project_name, richness, sd, sites)

# Asymptote of all plots 
all_asymptote <- asymptote_mm(df_all_filters, species_unit = rank_choice) %>%
  mutate(project_name = "All") %>%
  select(project_name, asymptote)

# Bind together
df_all_accumulation <- rbind(df_accumulation, all_accumulation)
df_all_asymptote <- rbind(df_asymptote, all_asymptote)

# 
df_join_all <- df_all_accumulation %>%
  left_join(., df_all_asymptote, by = "project_name") %>%
  group_by(project_name) %>%
  mutate(position_asymptote_y = 0.20 * max(asymptote), 
         position_asymptote_x = max(sites))

# Plot with facet
plot_acc_order <- ggplot(df_join_all, aes(fill = project_name, col = project_name)) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.5) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 1) +
  facet_wrap(~project_name, scales = "free") +
  ylab("Number of orders") +
  xlab("Samples (filter)") +
  theme_bw() + 
  ggtitle("Orders") + 
  geom_text(aes(x = position_asymptote_x, y =position_asymptote_y, hjust = 1, label = paste("asymptote =", round(asymptote, 1), "Orders")), col = "black") 

plot_acc_order

# ggsave("outputs/03_accumulation_curves/01b_accumulation_curve_all_projects_combination_no_facet_order.png", width = 12, height = 8)

# Table
stats<- df_join_all %>%
  group_by(project_name) %>%
  summarise(richness = max(richness), 
            asymptote = round(max(asymptote), 0)) %>%
  mutate(rank = 'Order')

# Save
write.csv(stats, "outputs/03_accumulation_curves/asymptotes_order.csv", row.names = F)

# Simple plot on all 
df_order <- df_join_all %>% 
  filter(project_name == "All") %>%
  mutate(level = "order")

order <- ggplot(df_order, aes(fill = project_name, col = project_name)) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.5) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 1) +
  ylab("Number of orders") +
  xlab("Samples (filter)") +
  theme_bw() + 
  ggtitle("Orders") + 
  geom_text(aes(x = position_asymptote_x, y =position_asymptote_y, hjust = 1, label = paste("asymptote =", round(asymptote, 0), "Orders")), col = "black") 

order

# ------------------------------------------------------------------------------- # 
#### On genus ----
# ------------------------------------------------------------------------------- # 

# rank_specify
rank_choice = 'new_genus_name'

# accumlation all plots
liste_accumulation <- lapply(liste_read_edna_LULU, accumulation_curve_df, species_unit = rank_choice)

# Asymptote of all plots 
liste_asymptote <- lapply(liste_read_edna_LULU, asymptote_mm, species_unit = rank_choice)

# Unlist
df_accumulation <- bind_rows(liste_accumulation, .id = "project_name")
df_asymptote <- bind_rows(liste_asymptote, .id = "project_name")

# Dataset for plot
df_accumulation_all <- left_join(df_accumulation, df_asymptote, by = "project_name") %>%
  group_by(project_name) %>%
  mutate(position_asymptote_y = 0.20 * max(asymptote), 
         position_asymptote_x = max(sites))

# Add All samples
# Add a global saturation curve, i.e. all samples together?
all_accumulation <- accumulation_curve_df(df_all_filters, species_unit = rank_choice) %>%
  mutate(project_name = "All") %>%
  select(project_name, richness, sd, sites)

# Asymptote of all plots 
all_asymptote <- asymptote_mm(df_all_filters, species_unit = rank_choice) %>%
  mutate(project_name = "All") %>%
  select(project_name, asymptote)

# Bind together
df_all_accumulation <- rbind(df_accumulation, all_accumulation)
df_all_asymptote <- rbind(df_asymptote, all_asymptote)

# 
df_join_all <- df_all_accumulation %>%
  left_join(., df_all_asymptote, by = "project_name") %>%
  group_by(project_name) %>%
  mutate(position_asymptote_y = 0.20 * max(asymptote), 
         position_asymptote_x = max(sites))

# Plot with facet
plot_acc_genus <- ggplot(df_join_all, aes(fill = project_name)) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.7) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 1) +
  scale_fill_manual(values=c("#A41D1A", "#E5A729", "#8AAE8A", "#4F4D1D", "#863b34", "#C67052"))+ 
  facet_wrap(~project_name, scales = "free") +
  ylab("Number of genus") +
  xlab("Samples (filter)") +
  theme_bw() + 
  ggtitle("Genus") + 
  geom_text(aes(x = position_asymptote_x, y =position_asymptote_y, hjust = 1, label = paste("asymptote =", round(asymptote, 1), "Genus")), col = "black") 

plot_acc_genus

ggsave("outputs/03_accumulation_curves/01b_accumulation_curve_all_projects_combination_no_facet_genus.png", width = 12, height = 8)

# Table
stats<- df_join_all %>%
  group_by(project_name) %>%
  summarise(richness = max(richness), 
            asymptote = round(max(asymptote), 0)) %>%
  mutate(rank = 'Genus')

# Save
write.csv(stats, "outputs/03_accumulation_curves/asymptotes_genus.csv", row.names = F)

# Simple plot on all 
df_genus <- df_join_all %>% 
  filter(project_name == "All") %>%
  mutate(level = "genus")

genus <- ggplot(df_genus, aes(fill = project_name, col = project_name)) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.5) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 1) +
  ylab("Number of genus") +
  xlab("Samples (filter)") +
  theme_bw() + 
  ggtitle("Genus") + 
  geom_text(aes(x = position_asymptote_x, y =position_asymptote_y, hjust = 1, label = paste("asymptote =", round(asymptote, 0), "Genus")), col = "black") 

genus