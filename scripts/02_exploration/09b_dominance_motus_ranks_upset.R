# Lib
library(UpSetR)
library(devtools)
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)


# ---------------------------------------------------------------- # 

# Lib
library(tidyverse)
library(reshape2)
library(vegan)
library(betapart)
library(ggpubr)

# fct
"%ni%" <- Negate("%in%")

# data
load("Rdata/02_clean_all.Rdata")

# Remove estuary stations and deep niskin station
df_all_filters <- df_all_filters %>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")) %>%
  filter(sample_method !="niskin" & region!="East_Pacific" & comment %ni% c("Distance decay 600m", "Distance decay 300m") & station!="glorieuse_distance_300m")%>%
  filter(project != "SEAMOUNTS") %>% 
  filter(habitat_type %ni% c("BAIE", "Sommet"))

# on est d'accord qu'on enlève toujours les MOTUs non assignés au dessus de la famille dans ce MS? 
# -> faire les deux et mettre les familles seulement en SI. 
#df_all_filters <- df_all_filters %>%
#  filter(!is.na(new_family_name))

# Des que le style des figures est bien posé, faire une fonction pour simplifier le code 

unique(df_all_filters$region)

# the motus 
df_regions <- split(df_all_filters, df_all_filters$region)
df_site <- split(df_all_filters, df_all_filters$site)

# ------------------------------------------------------------------------------- # 
#### REGION - MOTUs ---- 

# Color panel 
#pal <- park_palette("Everglades", 3)
# order des regions: 
# 1) central indo pacific (lengguru)
# 2) Caraibes
# 3) West Indian 
# 4) Central pacific

pal <- c("#8AAE8A", "#E5A729", "#C67052", "#4F4D1D", "#863b34") 
pal

# Construct the inital matrix - add some important infos as a side
motus_fam_region <-  df_all_filters %>%
  group_by(region) %>%
  summarize(n_motus = n_distinct(sequence), 
            n_family = n_distinct(new_family_name)) %>%
  as.data.frame() %>%
  arrange(n_motus)

# metadata - needed to control colors?
metadata1 <- df_all_filters %>%
  distinct(region) %>% 
  mutate(sets = region) %>%
  select(sets, region) %>%
  mutate(color = case_when(
    region == "South_West_Pacific" ~ pal[5],
    region == "Central_Pacific" ~ pal[4],
    region == "West_Indian" ~ pal[3],
    region == "Caribbean" ~ pal[2], 
    region == "Central_Indo_Pacific" ~ pal[1]
  )) %>%
  as.data.frame() %>%
  left_join(., motus_fam_region)

# MOTUs
matrix_motus <- df_all_filters %>%
  distinct(sequence, new_scientific_name_ncbi) %>%
  mutate(`Central_Indo_Pacific` = ifelse(sequence %in% df_regions$`Central_Indo_Pacific`$sequence, 1, 0), 
         Central_Pacific = ifelse(sequence %in% df_regions$Central_Pacific$sequence, 1, 0),
         Caribbean = ifelse(sequence %in% df_regions$Caribbean$sequence, 1, 0),
         West_Indian = ifelse(sequence %in% df_regions$West_Indian$sequence, 1, 0),
         South_West_Pacific = ifelse(sequence %in% df_regions$South_West_Pacific$sequence, 1, 0)) %>%
  as.data.frame()

# Supp Settings
# Ideas: rajouter le nombre de familles/genres sur le cote, par region? en boxplot. Le nombre de samples peut être ? 

# Construct a nicer plot
# Base plot
upset(matrix_motus, 
      order.by = c("freq"),
      mainbar.y.label = "Number of MOTUs", 
      sets.x.label = "Number of MOTUs", 
      text.scale = c(1.2, 1.2, 1.2,1.2,1.2,1.2))

# Trials ameliorations: color in set bar, points and lines and alpha heatmap by region groups
# Work in progress
p1 <- upset(matrix_motus, 
            mb.ratio = c(0.7, 0.3),
      order.by = c("freq"),
      mainbar.y.label = "Number of MOTUs", 
      sets.x.label = "Number of MOTUs", 
      text.scale = c(1.2, 1.2, 1.2,1.2,1.2,1.2), 
      # Color bar 
      sets.bar.color=c("#8AAE8A", "#863b34", "#E5A729", "#C67052", "#4F4D1D"), 
      # Color matrix
      set.metadata = list(
        data = metadata1,
        plots = list(list(
          type = "matrix_rows",
          column = "region", 
          colors = c(`Central_Indo_Pacific` = pal[1], Caribbean =  pal[2], West_Indian = pal[3], Central_Pacific =  pal[4], South_West_Pacific = pal[5]),
          alpha = 0.3
        ))
      ))
p1

save(p1, file = "Rdata/upset_plot_motus_region.rdata")
png('outputs/09_dominance_motus_ranks/upset_plot_region_motus.png', width = 6, height=3, units = "in", res=300)
p1
grid.text("Regions - MOTUs",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# ------------------------------------------------------------------------------- # 
#### REGION - FAMILY ---- 

# rarity in samples - by families 
family_samples_rarity <- df_all_filters %>%
  filter(!is.na(new_family_name)) %>%
  group_by(new_family_name) %>%
  summarise(n_samples = n_distinct(sample_name_all_pcr))

# Family
matrix_family <- df_all_filters %>%
  filter(!is.na(new_family_name)) %>%
  distinct(new_family_name) %>%
  mutate(`Central_Indo_Pacific` = ifelse(new_family_name %in% df_regions$`Central_Indo_Pacific`$new_family_name, 1, 0), 
         Central_Pacific = ifelse(new_family_name %in% df_regions$Central_Pacific$new_family_name, 1, 0),
         Caribbean = ifelse(new_family_name %in% df_regions$Caribbean$new_family_name, 1, 0),
         West_Indian = ifelse(new_family_name %in% df_regions$West_Indian$new_family_name, 1, 0),
         South_West_Pacific = ifelse(new_family_name %in% df_regions$South_West_Pacific$new_family_name, 1, 0)) %>%
  left_join(., family_samples_rarity) %>%
  as.data.frame()

# Simple plot
upset(matrix_family, 
      order.by = c("freq"),
      mainbar.y.label = "Number of families", 
      sets.x.label = "Number of families", 
      text.scale = c(1.2, 1.2, 1.2,1.2,1.2,1.2))

# Alternative plot
p2 <- upset(matrix_family, 
            order.by = c("freq"),
            mainbar.y.label = "Number of families", 
            sets.x.label = "Number of families", 
            text.scale = c(1.2, 1.2, 1.2,1.2,1.2,1.2), 
            # Color bar 
            sets.bar.color=c("#8AAE8A", "#863b34", "#E5A729", "#C67052", "#4F4D1D"), 
            # Color matrix
            set.metadata = list(
              data = metadata1,
              plots = list(list(
                type = "matrix_rows",
                column = "region", 
                colors = c(`Central_Indo_Pacific` = pal[1], Caribbean =  pal[2], West_Indian = pal[3], Central_Pacific =  pal[4], South_West_Pacific = pal[5]),
                alpha = 0.3
              ))
            ))
p2

png('outputs/09_dominance_motus_ranks/upset_plot_region_family.png', width = 6, height=3, units = "in", res=300)
p2
grid.text("Regions - Family",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# ------------------------------------------------------------------------------- # 
#### SITE - MOTUs ---- 

all_sites <- names(df_site)

# Count MOTUs per site
motus_sites <- df_all_filters %>%
  group_by(site) %>%
  dplyr::summarize(n_motus = n_distinct(sequence)) %>%
  arrange(n_motus)

# metadata - needed to control colors?
metadata1 <- df_all_filters %>%
  distinct(site, region) %>% 
  mutate(sets = site) %>%
  select(sets, region) %>%
  mutate(color = case_when(
    region == "South_West_Pacific" ~ pal[5],
    region == "Central_Pacific" ~ pal[4], 
    region == "West_Indian" ~ pal[3],
    region == "Caribbean" ~ pal[2], 
    region == "Central_Indo_Pacific" ~ pal[1]
  )) %>%
  left_join(., motus_sites, by = c("sets" = "site")) %>%
  # Arrange by n_motu frequency to get the color right in the plot
  arrange(-n_motus) %>%
  as.data.frame()

# MOTUs
matrix_motus <- df_all_filters %>%
  distinct(sequence) %>%
  as.data.frame()

# Variables
rank <- "sequence"
dataset <- df_site

# Fill the sites - automatic
for (i in all_sites){
  matrix_motus[,i] <- ifelse(matrix_motus[,rank] %in% dataset[[i]][[rank]], 1, 0)
  print(i)
}

# Plotter les 40 premieres intersections pour la lisibilité
p3 <- upset(matrix_motus, 
            nsets = 25,
            order.by = c("freq"),
            mainbar.y.label = "Number of MOTUs", 
            sets.x.label = "Number of MOTUs", 
            nintersects = 40)

p3

# Save
#   png('outputs/09_dominance_motus_ranks/upset_plot_sites_motus.png', width = 12, height=8, units = "in", res=300)
#   p3
#   grid.text("Sites - MOTUs",x = 0.65, y=0.95, gp=gpar(fontsize=15))
#   dev.off()

# Alternative plot 
# Work in progress

# Trials ameliorations: color in set bar, points and lines and alpha heatmap by region groups
p3_color <- upset(matrix_motus, 
      nsets = 25,
      order.by = c("freq"),
      mainbar.y.label = "Number of MOTUs", 
      sets.x.label = "Number of MOTUs", 
      text.scale = c(1.2, 1.2, 1.2,1.2,1.2,1.7), 
      # Color bar 
      sets.bar.color=metadata1$color, 
      # Color matrix
      set.metadata = list(
        data = metadata1,
        plots = list(list(
          type = "matrix_rows",
          column = "region", 
          colors = c(`Central_Indo_Pacific` = pal[1], Caribbean =  pal[2], West_Indian = pal[3], Central_Pacific =  pal[4], South_West_Pacific = pal[5]),
          alpha = 0.3
        ))
      ))

p3_color

# Save - w/ colors
png('outputs/09_dominance_motus_ranks/upset_plot_sites_motus_colors.png', width = 12, height=12, units = "in", res=300)
p3_color
grid.text("Sites - MOTUs",x = 0.65, y=0.95, gp=gpar(fontsize=15))
dev.off()


# ------------------------------------------------------------------------------- # 
#### SITE - FAMILY ---- 

all_sites <- names(df_site)

# Count MOTUs per site
motus_family <- df_all_filters %>%
  group_by(site) %>%
  dplyr::summarize(n_family = n_distinct(new_family_name)) %>%
  arrange(n_family)

# metadata - needed to control colors?
metadata1 <- df_all_filters %>%
  distinct(site, region) %>% 
  mutate(sets = site) %>%
  select(sets, region) %>%
  mutate(color = case_when(
    region == "South_West_Pacific" ~ pal[5],
    region == "Central_Pacific" ~ pal[4],
    region == "West_Indian" ~ pal[3],
    region == "Caribbean" ~ pal[2], 
    region == "Central_Indo_Pacific" ~ pal[1]
  )) %>%
  left_join(., motus_family, by = c("sets" = "site")) %>%
  # Arrange by n_motu frequency to get the color right in the plot
  arrange(-n_family) %>%
  as.data.frame()

# MOTUs
matrix_family <- df_all_filters %>%
  filter(!is.na(new_family_name)) %>%
  distinct(new_family_name) %>%
  as.data.frame()

# Variables
rank <- "new_family_name"
dataset <- df_site

# Fill the sites - automatic
for (i in all_sites){
  matrix_family[,i] <- ifelse(matrix_family[,rank] %in% dataset[[i]][[rank]], 1, 0)
  print(i)
}

# Plot color
p4_color <- upset(matrix_family, 
                  nsets = 25,
                  nintersects = 13,
                  order.by = c("freq"),
                  mainbar.y.label = "Number of families", 
                  sets.x.label = "Number of families", 
                  text.scale = c(1.2, 1.2, 1.2,1.2,1.2,1.7), 
                  # Color bar 
                  sets.bar.color=metadata1$color, 
                  # Color matrix
                  set.metadata = list(
                    data = metadata1,
                    plots = list(list(
                      type = "matrix_rows",
                      column = "region", 
                      colors = c(`Central_Indo_Pacific` = pal[1], Caribbean =  pal[2], West_Indian = pal[3], Central_Pacific =  pal[4], South_West_Pacific = pal[5]),
                      alpha = 0.3
                    ))
                  ))

p4_color

# Save
png('outputs/09_dominance_motus_ranks/upset_plot_sites_family_color.png', width = 10, height=10, units = "in", res=300)
p4_color
grid.text("Sites - families",x = 0.65, y=0.95, gp=gpar(fontsize=15))
dev.off()


