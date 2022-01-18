# Sensitivity analysis to rarefied stations
# Supl. Mat. R1 manuscript

# Produce the following numbers:
# - Recalculer les richesses MOTUs et familles globales.
#  Juste une phrase dans la réponse au reviewer avec les +- SD suffira
# - Une figure avec les richesse par regions de chaque simu et richesse moyenne (violin plot + boxplot)

# Lib 
library(tidyverse)
library(ggplot2)
library(patchwork)

# Functions
"%ni%" <- Negate("%in%")

# Data - load rarefied data 
load("Rdata/rarefied_regions.rdata")
load("Rdata/rarefied_sites.rdata")

# The full dataset
load("Rdata/02-clean-data.Rdata")

df_all_filters <- df_all_filters %>%
  filter(province %in% c("Western_Indian_Ocean", "Southeast_Polynesia", "Tropical_Northwestern_Atlantic", "Western_Coral_Triangle", "Tropical_Southwestern_Pacific"))%>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3", "glorieuse_distance_300m")) %>%
  filter(sample_method !="niskin" & comment %ni% c("Distance decay 600m", "Distance decay 300m"))%>%
  filter(project != "Curacao") %>%
  filter(habitat=="marine")%>%
  filter(habitat_type %ni% c("BAIE"))%>%
  filter(site35!="") %>%
  filter(depth<40) %>%
  filter(family_name_corrected %ni% "Salmonidae")

# Use lapply to compute the metrics for each element of the list 

# ------------------------------------------------------------------ # 
# Rarefaction to sites - Richness MOTUs and number of families 
# ------------------------------------------------------------------ # 

k <- df_all_filters %>%
  group_by(province, site35) %>%
  summarise(n_station = n_distinct(station))

dist_per_region <- df_all_filters %>%
  group_by(province) %>%
  summarise(mean_dist_to_CT = mean(dist_to_CT)) %>%
  mutate(mean_dist_to_CT = case_when(
    province == "Western_Coral_Triangle" ~ 0 ,
    TRUE ~ mean_dist_to_CT
  ))

# Output: dataframe, nb of MOTUs and nb of families
list_global_numbers <- lapply(rarefied_sites, function(x){
  
  # x<-rarefied_sites[[1]]
  
  nb.motus <- x %>% distinct(definition) %>% nrow()
  nb.family <- x %>% distinct(family_name_corrected) %>% nrow()
  
  y <- data.frame(nb_motus = nb.motus, 
             nb_family= nb.family)
  return(y)
})
df_global_numbers <- bind_rows(list_global_numbers)

list_region <- lapply(rarefied_sites, function(x){
  # x<-rarefied_sites[[1]]
  x %>%
    group_by(province) %>%
    summarise(motus_richness = n_distinct(definition), 
              family_richness = n_distinct(family_name_corrected))
})

df_region_numbers <- bind_rows(list_region, .id = "groups") %>%
 left_join(., dist_per_region)

# Compute numbers - MOTUs
round(mean(df_global_numbers$nb_motus), 0)
round(sd(df_global_numbers$nb_motus), 0)

df_all_filters %>% distinct(definition) %>% nrow()

# Families 
round(mean(df_global_numbers$nb_family), 0)
round(sd(df_global_numbers$nb_family), 0)

df_all_filters %>% distinct(family_name_corrected) %>% nrow()

count_all <- data.frame(type=c("nb_motus", "nb_family"), 
                        count=c(2023, 127))

# ----------------------------- # 
# Plot - nb of MOTUs and families global 

plot_counts <- df_global_numbers %>%
  mutate(j = rownames(.)) %>%
  pivot_longer(!j, names_to = "type", values_to = "count")

plot.motus <- ggplot(data = plot_counts %>% filter(type=="nb_motus"), aes(x=type, y = count)) + 
  geom_point(alpha=0.2) + 
  geom_violin(alpha=0.4) +
  geom_point(data=count_all%>% filter(type=="nb_motus"), aes(x=type, y = count), col="red") + 
  ylim(c(0,2100)) + theme_bw()
plot.motus

plot.family <- ggplot(data = plot_counts %>% filter(type=="nb_family"), aes(x=type, y = count)) + 
  geom_point(alpha=0.2) + 
  geom_violin(alpha=0.4) +
  geom_point(data=count_all%>% filter(type=="nb_family"), aes(x=type, y = count), col="red") + 
  ylim(c(0,130)) + theme_bw()

plot.motus+plot.family

ggsave("/Users/virginiemarques/Desktop/supl_rar_sites_global.png")

# ----------------------------- # 
# Plot - distance to coral triangle

for_plot_all <- df_all_filters %>%
  group_by(province) %>%
  summarise(motus_richness = n_distinct(definition), 
            family_richness = n_distinct(family_name_corrected),
            mean_dist_to_CT = mean(dist_to_CT)) %>%
  mutate(mean_dist_to_CT = case_when(
    province == "Western_Coral_Triangle" ~ 0 ,
    TRUE ~ mean_dist_to_CT
  ))

# Distance to CT 
CT_motus <- ggplot(data=df_region_numbers, aes(x=mean_dist_to_CT, y = motus_richness, col = province)) + 
  geom_point(shape=20, size=3, alpha=0.7) +
  geom_boxplot(alpha=0.5) +
  geom_point(data = for_plot_all, aes(x=mean_dist_to_CT, y = motus_richness), size=3,  col = "grey30", alpha=0.5, shape=17) +
  theme_bw() + 
  scale_color_manual(values = c("#a6611a", "#E5A729","#b2182b", "#80cdc1", "#015462"),
                    name = "province", labels = c("Southeast_Polynesia", "Tropical_Northwestern_Atlantic", "Tropical_Southwestern_Pacific", "Western_Coral_Triangle", "Western_Indian_Ocean")) +
  theme(legend.position = "none")
CT_motus

CT_families <- ggplot(data=df_region_numbers, aes(x=mean_dist_to_CT, y = family_richness, col = province)) + 
  geom_point(shape=20, size=3, alpha=0.7) +
  geom_boxplot(alpha=0.5) +
  geom_point(data = for_plot_all, aes(x=mean_dist_to_CT, y = family_richness), size=3,  col = "grey30", alpha=0.5, shape=17) +
  theme_bw() + 
  scale_color_manual(values = c("#a6611a", "#E5A729","#b2182b", "#80cdc1", "#015462"),
                     name = "province", labels = c("Southeast_Polynesia", "Tropical_Northwestern_Atlantic", "Tropical_Southwestern_Pacific", "Western_Coral_Triangle", "Western_Indian_Ocean")) +
  theme(legend.position = "none")

CT_motus + CT_families

ggsave("/Users/virginiemarques/Desktop/dist_coral_triangle_rar_sites.png")

# ------------------------------------------------------------------ # 
# Rarefaction to regions - Richness MOTUs and number of families 
# ------------------------------------------------------------------ # 

# Output: dataframe, nb of MOTUs and nb of families
list_global_numbers <- lapply(rarefied_regions, function(x){
  
  # x<-rarefied_sites[[1]]
  
  nb.motus <- x %>% distinct(definition) %>% nrow()
  nb.family <- x %>% distinct(family_name_corrected) %>% nrow()
  
  y <- data.frame(nb_motus = nb.motus, 
                  nb_family= nb.family)
  return(y)
})

df_global_numbers_regions <- bind_rows(list_global_numbers)

list_region <- lapply(rarefied_regions, function(x){
  # x<-rarefied_sites[[1]]
  x %>%
    group_by(province) %>%
    summarise(motus_richness = n_distinct(definition), 
              family_richness = n_distinct(family_name_corrected))
})

df_region_numbers_regions <- bind_rows(list_region, .id = "groups") %>%
  left_join(., dist_per_region)

# Compute numbers - MOTUs
round(mean(df_global_numbers_regions$nb_motus), 0)
round(sd(df_global_numbers_regions$nb_motus), 0)

df_all_filters %>% distinct(definition) %>% nrow()

# Families 
round(mean(df_global_numbers_regions$nb_family), 0)
round(sd(df_global_numbers_regions$nb_family), 0)

df_all_filters %>% distinct(family_name_corrected) %>% nrow()

count_all <- data.frame(type=c("nb_motus", "nb_family"), 
                        count=c(2023, 127))

# Plot 
plot_counts <- df_global_numbers_regions %>%
  mutate(j = rownames(.)) %>%
  pivot_longer(!j, names_to = "type", values_to = "count")

# Plot 

# ggplot(data = df_global_numbers) +
#   geom_point(aes(x="MOTUs", y = 2023), col="red", alpha=0.5) + 
#   geom_point(aes(x="MOTUs", y = nb_motus), alpha=0.5) + 
#   theme_bw() + 
#   ylim(c(0,2100))

plot.motus <- ggplot(data = plot_counts %>% filter(type=="nb_motus"), aes(x=type, y = count)) + 
  geom_point(alpha=0.2) + 
  geom_violin(alpha=0.4) +
  geom_point(data=count_all%>% filter(type=="nb_motus"), aes(x=type, y = count), col="red") + 
  ylim(c(0,2100)) + theme_bw()
plot.motus

plot.family <- ggplot(data = plot_counts %>% filter(type=="nb_family"), aes(x=type, y = count)) + 
  geom_point(alpha=0.2) + 
  geom_violin(alpha=0.4) +
  geom_point(data=count_all%>% filter(type=="nb_family"), aes(x=type, y = count), col="red") + 
  ylim(c(0,130)) + theme_bw()

plot.motus+plot.family

ggsave("/Users/virginiemarques/Desktop/supl_rar_regions_global.png")

# ----------------------------- # 
# Plot - distance to coral triangle

# Distance to CT 
CT_motus <- ggplot(data=df_region_numbers_regions, aes(x=mean_dist_to_CT, y = motus_richness, col = province)) + 
  geom_point(shape=20, size=3, alpha=0.7) +
  geom_boxplot(alpha=0.5) +
  geom_point(data = for_plot_all, aes(x=mean_dist_to_CT, y = motus_richness), size=3,  col = "grey30", alpha=0.5, shape=17) +
  theme_bw() + 
  scale_color_manual(values = c("#a6611a", "#E5A729","#b2182b", "#80cdc1", "#015462"),
                     name = "province", labels = c("Southeast_Polynesia", "Tropical_Northwestern_Atlantic", "Tropical_Southwestern_Pacific", "Western_Coral_Triangle", "Western_Indian_Ocean")) +
  theme(legend.position = "none")
CT_motus

CT_families <- ggplot(data=df_region_numbers_regions, aes(x=mean_dist_to_CT, y = family_richness, col = province)) + 
  geom_point(shape=20, size=3, alpha=0.7) +
  geom_boxplot(alpha=0.5) +
  geom_point(data = for_plot_all, aes(x=mean_dist_to_CT, y = family_richness), size=3,  col = "grey30", alpha=0.5, shape=17) +
  theme_bw() + 
  scale_color_manual(values = c("#a6611a", "#E5A729","#b2182b", "#80cdc1", "#015462"),
                     name = "province", labels = c("Southeast_Polynesia", "Tropical_Northwestern_Atlantic", "Tropical_Southwestern_Pacific", "Western_Coral_Triangle", "Western_Indian_Ocean")) +
  theme(legend.position = "none")

CT_motus + CT_families

ggsave("/Users/virginiemarques/Desktop/dist_coral_triangle_rar_regions.png")





