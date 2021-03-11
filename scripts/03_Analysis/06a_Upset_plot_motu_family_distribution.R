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
library(conflicted)

# fct
"%ni%" <- Negate("%in%")
conflict_prefer("arrange", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("filter", "dplyr")
# data
load("Rdata/02-clean-data.Rdata")

# Remove estuary stations and deep niskin station
df_all_filters <- df_all_filters %>%
  filter(province %in% c("Western_Indian_Ocean", "Southeast_Polynesia", "Tropical_Northwestern_Atlantic", "Western_Coral_Triangle", "Tropical_Southwestern_Pacific"))%>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3", "glorieuse_distance_300m")) %>%
  filter(sample_method !="niskin" & comment %ni% c("Distance decay 600m", "Distance decay 300m"))%>%
  filter(project.y != "SEAMOUNTS") %>%
  filter(habitat=="marine")%>%
  filter(habitat_type %ni% c("BAIE"))

unique(df_all_filters$province)

# the motus 
df_provinces <- split(df_all_filters, df_all_filters$province)
df_site <- split(df_all_filters, df_all_filters$site35)

# ------------------------------------------------------------------------------- # 
#### REGION - MOTUs ---- 

# Color panel 
# order of regions: 
# 1) central indo pacific (lengguru)
# 2) Caribbean
# 3) West Indian 
# 4) Central pacific
# 5) South-west Pacific

pal <- c("#80cdc1", "#E5A729", "#015462", "#a6611a", "#b2182b") 


# Construct the inital matrix - add some important infos as a side
motus_fam_region <-  df_all_filters %>%
  group_by(province) %>%
  dplyr::summarize(n_motus = n_distinct(sequence), 
            n_family = n_distinct(family_name_corrected)) %>%
  as.data.frame() %>%
  dplyr::arrange(n_motus)

# metadata - 
metadata1 <- df_all_filters %>%
  distinct(province) %>% 
  dplyr::mutate(sets = province) %>%
  select(sets, province) %>%
  dplyr::mutate(color = case_when(
    province == "Tropical_Southwestern_Pacific" ~ pal[5],
    province == "Southeast_Polynesia" ~ pal[4],
    province == "Western_Indian_Ocean" ~ pal[3],
    province == "Tropical_Northwestern_Atlantic" ~ pal[2], 
    province == "Western_Coral_Triangle" ~ pal[1]
  )) %>%
  as.data.frame() %>%
  left_join(., motus_fam_region)

# MOTUs
matrix_motus <- df_all_filters %>%
  distinct(sequence, new_scientific_name_ncbi) %>%
  dplyr::mutate(`Western_Coral_Triangle` = ifelse(sequence %in% df_provinces$`Western_Coral_Triangle`$sequence, 1, 0), 
         Southeast_Polynesia = ifelse(sequence %in% df_provinces$Southeast_Polynesia$sequence, 1, 0),
         Tropical_Northwestern_Atlantic = ifelse(sequence %in% df_provinces$Tropical_Northwestern_Atlantic$sequence, 1, 0),
         Western_Indian_Ocean = ifelse(sequence %in% df_provinces$Western_Indian_Ocean$sequence, 1, 0),
         Tropical_Southwestern_Pacific = ifelse(sequence %in% df_provinces$Tropical_Southwestern_Pacific$sequence, 1, 0)) %>%
  as.data.frame()



# Plot MOTUs distribution in regions
p1 <- upset(matrix_motus, 
            mb.ratio = c(0.7, 0.3),
      order.by = c("freq"),
      mainbar.y.label = "Number of MOTUs", 
      sets.x.label = "Number of MOTUs", 
      text.scale = c(1.2, 1.2, 1.2,1.2,1.2,1.2), 
      # Color bar 
      sets.bar.color=c("#80cdc1", "#b2182b", "#E5A729", "#015462", "#a6611a"), 
      # Color matrix
      set.metadata = list(
        data = metadata1,
        plots = list(list(
          type = "matrix_rows",
          column = "province", 
          colors = c(`Western_Coral_Triangle` = pal[1], Tropical_Northwestern_Atlantic =  pal[2], Western_Indian_Ocean = pal[3], Southeast_Polynesia =  pal[4], Tropical_Southwestern_Pacific = pal[5]),
          alpha = 0.3
        ))
      ))
p1

save(p1, file = "Rdata/upset_plot_motus_region.rdata")
png('outputs/00_Figures_for_paper/Figure3a.png', width = 7, height=3.6, units = "in", res=300)
p1
grid.text("Provinces - MOTUs",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# ------------------------------------------------------------------------------- # 
#### REGION - FAMILY ---- 

# rarity in samples - by families 
family_samples_rarity <- df_all_filters %>%
  filter(!is.na(family_name_corrected)) %>%
  group_by(family_name_corrected) %>%
  summarise(n_samples = n_distinct(sample_name_all_pcr))

# Family
matrix_family <- df_all_filters %>%
  filter(!is.na(family_name_corrected)) %>%
  distinct(family_name_corrected) %>%
  mutate(`Western_Coral_Triangle` = ifelse(family_name_corrected %in% df_provinces$`Western_Coral_Triangle`$family_name_corrected, 1, 0), 
         Southeast_Polynesia = ifelse(family_name_corrected %in% df_provinces$Southeast_Polynesia$family_name_corrected, 1, 0),
         Tropical_Northwestern_Atlantic = ifelse(family_name_corrected %in% df_provinces$Tropical_Northwestern_Atlantic$family_name_corrected, 1, 0),
         Western_Indian_Ocean = ifelse(family_name_corrected %in% df_provinces$Western_Indian_Ocean$family_name_corrected, 1, 0),
         Tropical_Southwestern_Pacific = ifelse(family_name_corrected %in% df_provinces$Tropical_Southwestern_Pacific$family_name_corrected, 1, 0)) %>%
  left_join(., family_samples_rarity) %>%
  as.data.frame()


# Plot families distrbution across regions
p2 <- upset(matrix_family, 
            order.by = c("freq"),
            mainbar.y.label = "Number of families", 
            sets.x.label = "Number of families", 
            text.scale = c(1.2, 1.2, 1.2,1.2,1.2,1.2), 
            # Color bar 
            sets.bar.color=c("#80cdc1", "#b2182b", "#E5A729", "#015462", "#a6611a"), 
            # Color matrix
            set.metadata = list(
              data = metadata1,
              plots = list(list(
                type = "matrix_rows",
                column = "province", 
                colors = c(`Western_Coral_Triangle` = pal[1], Tropical_Northwestern_Atlantic =  pal[2], Western_Indian_Ocean = pal[3], Southeast_Polynesia =  pal[4], Tropical_Southwestern_Pacific = pal[5]),
                alpha = 0.3
              ))
            ))
p2

png('outputs/06_Upset_plots_Histograms_motus_family/upset_plot_region_family.png', width = 6, height=3, units = "in", res=300)
p2
grid.text("provinces - Family",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# ------------------------------------------------------------------------------- # 
#### SITE - MOTUs ---- 

all_sites <- names(df_site)

# Count MOTUs per site
motus_sites <- df_all_filters %>%
  group_by(site35) %>%
  dplyr::summarize(n_motus = n_distinct(sequence)) %>%
  arrange(n_motus)

# metadata - needed to control colors?
metadata1 <- df_all_filters %>%
  distinct(site35, province) %>% 
  mutate(sets = site35) %>%
  select(sets, province) %>%
  mutate(color = case_when(
    province == "Tropical_Southwestern_Pacific" ~ pal[5],
    province == "Southeast_Polynesia" ~ pal[4], 
    province == "Western_Indian_Ocean" ~ pal[3],
    province == "Tropical_Northwestern_Atlantic" ~ pal[2], 
    province == "Western_Coral_Triangle" ~ pal[1]
  )) %>%
  left_join(., motus_sites, by = c("sets" = "site35")) %>%
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

# Plot MOTUs distribution across sites
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
          column = "province", 
          colors = c(`Western_Coral_Triangle` = pal[1], Tropical_Northwestern_Atlantic =  pal[2], Western_Indian_Ocean = pal[3], Southeast_Polynesia =  pal[4], Tropical_Southwestern_Pacific = pal[5]),
          alpha = 0.3
        ))
      ))

p3_color

# Save - w/ colors
png('outputs/00_Figures_for_paper/Extended_Data/ED_Figure5.png', width = 12, height=12, units = "in", res=300)
p3_color
grid.text("Sites - MOTUs",x = 0.65, y=0.95, gp=gpar(fontsize=15))
dev.off()


# ------------------------------------------------------------------------------- # 
#### SITE - FAMILY ---- 

all_sites <- names(df_site)

# Count MOTUs per site
motus_family <- df_all_filters %>%
  group_by(site35) %>%
  dplyr::summarize(n_family = n_distinct(family_name_corrected)) %>%
  arrange(n_family)

# metadata
metadata1 <- df_all_filters %>%
  distinct(site35, province) %>% 
  mutate(sets = site35) %>%
  select(sets, province) %>%
  mutate(color = case_when(
    province == "Tropical_Southwestern_Pacific" ~ pal[5],
    province == "Southeast_Polynesia" ~ pal[4],
    province == "Western_Indian_Ocean" ~ pal[3],
    province == "Tropical_Northwestern_Atlantic" ~ pal[2], 
    province == "Western_Coral_Triangle" ~ pal[1]
  )) %>%
  left_join(., motus_family, by = c("sets" = "site35")) %>%
  # Arrange by n_motu frequency to get the color right in the plot
  arrange(-n_family) %>%
  as.data.frame()

# MOTUs
matrix_family <- df_all_filters %>%
  filter(!is.na(family_name_corrected)) %>%
  distinct(family_name_corrected) %>%
  as.data.frame()

# Variables
rank <- "family_name_corrected"
dataset <- df_site

# Fill the sites - automatic
for (i in all_sites){
  matrix_family[,i] <- ifelse(matrix_family[,rank] %in% dataset[[i]][[rank]], 1, 0)
  print(i)
}

# Plot families distribution across sites
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
                      column = "province", 
                      colors = c(`Western_Coral_Triangle` = pal[1], Tropical_Northwestern_Atlantic =  pal[2], Western_Indian_Ocean = pal[3], Southeast_Polynesia =  pal[4], Tropical_Southwestern_Pacific = pal[5]),
                      alpha = 0.3
                    ))
                  ))

p4_color

# Save
png('outputs/06_Upset_plots_Histograms_motus_family/upset_plot_sites_family_color.png', width = 10, height=10, units = "in", res=300)
p4_color
grid.text("Sites - families",x = 0.65, y=0.95, gp=gpar(fontsize=15))
dev.off()


