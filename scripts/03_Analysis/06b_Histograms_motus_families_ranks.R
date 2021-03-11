# This script will explore the dominance of MOTUs across its taxonomic ranks (genus, families, species-proxy)

# Lib
library(tidyverse)
library(reshape2)
library(vegan)
library(betapart)
library(ggpubr)
library(conflicted)

# 
conflict_prefer("filter", "dplyr")
conflict_prefer("summarise", "dplyr")

# fct
"%ni%" <- Negate("%in%")

# data
load("Rdata/02-clean-data.Rdata")

#Remove estuary stations and deep niskin station
df_all_filters <- df_all_filters %>%
  filter(province %in% c("Western_Indian_Ocean", "Southeast_Polynesia", "Tropical_Northwestern_Atlantic", "Western_Coral_Triangle", "Tropical_Southwestern_Pacific"))%>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3", "glorieuse_distance_300m")) %>%
  filter(sample_method !="niskin" & comment %ni% c("Distance decay 600m", "Distance decay 300m"))%>%
  filter(project.y != "SEAMOUNTS") %>%
  filter(habitat=="marine")%>%
  filter(habitat_type %ni% c("BAIE"))


# N global 
Nmotus <- length(unique(df_all_filters$sequence))
Nmotus
# 2160 MOTUs in total

# N families
Nfamily <- df_all_filters %>%
  filter(!is.na(family_name_corrected)) %>%
  distinct(family_name_corrected) %>% pull() %>% length()
Nfamily
# 145

# N genus
Ngenus <- df_all_filters %>%
  filter(!is.na(new_genus_name)) %>%
  distinct(new_genus_name) %>% pull() %>% length()
Ngenus
#373

# --------------------------------------------------------------------- # 
#### Functions ----

# rank_x : the spatial scale (region, site, station, sample)
# rank_y : the taxonomic unit (MOTUs, family, genus)
# n_tot : total number of taxonomic units in the dataset 

fct_barplot <- function(dataset, rank_x, rank_y = "MOTUs", color = "#56B4E9", n_tot = Nmotus){
  # Plot
  ggplot(dataset, aes(x = n, y=n_motus)) + 
    geom_bar(stat="identity", color = "black", fill = color) + 
    geom_text(aes(label=n_motus), vjust=-0.5, color="black", size=3.5) +
    theme_classic() + 
    labs(x= paste0("Number of ", rank_x), y = paste0("Number of ", rank_y)) + 
    ggtitle(paste0("By ", rank_x), subtitle = paste0(round(pull(dataset[1,2])/n_tot*100,0), "% of ", rank_y, " are present in only one ", rank_x))
}

# End of functions 
# --------------------------------------------------------------------- # 


# --------------------------------------------------------------------- # 
#### MOTUs unit  ----

# --------------------- # 
# Region scale 

length(unique(df_all_filters$province))

# Count
motu_province <- df_all_filters %>%
  group_by(sequence) %>%
  summarise(n = n_distinct(province)) %>%
  ungroup() %>%
  group_by(n) %>%
  summarise(n_motus = n_distinct(sequence)) %>%
  mutate(rank = "province")

# plot
p_motu_province <- fct_barplot(motu_province, rank_x="province", rank_y = "MOTUs")
p_motu_province

# 2 MOTUs are present in all 5 regions:


# --------------------- # 
# site scale 

length(unique(df_all_filters$site35))

# Count
motu_site <- df_all_filters %>%
  group_by(sequence) %>%
  summarise(n = n_distinct(site35)) %>%
  ungroup() %>%
  group_by(n) %>%
  summarise(n_motus = n_distinct(sequence))%>%
  mutate(rank = "site35")
save(motu_site, file="Rdata/rarete_motu_site.rdata")

# plot
p_motu_site <- fct_barplot(motu_site, rank_x="site35", rank_y = "MOTUs")
p_motu_site

ggplot(motu_site)+
  geom_bar(aes(x=n, y=n_motus), stat= 'identity', fill="#d2981a")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA))+
  labs(x="Number of sites",y="Number of MOTUs")

# --------------------- # 
# station scale 

length(unique(df_all_filters$station))

# Count
motu_station <- df_all_filters %>%
  group_by(sequence) %>%
  summarise(n = n_distinct(station)) %>%
  ungroup() %>%
  group_by(n) %>%
  summarise(n_motus = n_distinct(sequence))%>%
  mutate(rank = "station")
save(motu_station, file="Rdata/rarete_motu_station.rdata")

# plot
p_motu_station <- fct_barplot(motu_station, rank_x="station", rank_y = "MOTUs", color = "#d2981a")
p_motu_station

ggplot(motu_station)+
  geom_bar(aes(x=n, y=n_motus), stat= 'identity', fill="#d2981a")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA))+
  labs(x="Number of stations",y="Number of MOTUs")

# --------------------- # 
# sample scale 

length(unique(df_all_filters$sample_name_all_pcr))

# Count
motu_sample <- df_all_filters %>%
  group_by(sequence) %>%
  summarise(n = n_distinct(sample_name_all_pcr)) %>%
  ungroup() %>%
  group_by(n) %>%
  summarise(n_motus = n_distinct(sequence))%>%
  mutate(rank = "sample")

# plot
p_motu_sample <- fct_barplot(motu_sample, rank_x="sample", rank_y = "MOTUs", color = "grey")
p_motu_sample

# --------------------- # 
# all scales

ggarrange(p_motu_province, 
          p_motu_site + rremove("ylab"), 
          p_motu_station, 
          p_motu_sample+ rremove("ylab"),
          nrow=1, ncol=2)

ggsave("outputs/06_Upset_plots_Histograms_motus_family/repartition_MOTUs_regions_site_station.png", height = 10, width = 15)





# --------------------------------------------------------------------- # 
#### Family unit  ----

# --------------------- # 
# Region scale 

length(unique(df_all_filters$province))

# Count
family_province <- df_all_filters %>%
  filter(!is.na(family_name_corrected)) %>%
  group_by(family_name_corrected) %>%
  summarise(n = n_distinct(province)) %>%
  ungroup() %>%
  group_by(n) %>%
  summarise(n_motus = n_distinct(family_name_corrected)) %>%
  mutate(rank = "province")

# plot
p_family_province <- fct_barplot(family_province, rank_x="province", rank_y = "families", n_tot = Nfamily)
p_family_province


# --------------------- # 
# site scale 

length(unique(df_all_filters$site35))

# Count
family_site <- df_all_filters %>%
  filter(!is.na(family_name_corrected)) %>%
  group_by(family_name_corrected) %>%
  summarise(n = n_distinct(site35)) %>%
  ungroup() %>%
  group_by(n) %>%
  summarise(n_motus = n_distinct(family_name_corrected)) %>%
  mutate(rank = "site35")
save(family_site, file="Rdata/rarete_family_site.rdata")

# plot
p_family_site <- fct_barplot(dataset, rank_x="site35", rank_y = "families", n_tot = Nfamily)
p_family_site

ggplot(family_site)+
  geom_bar(aes(x=n, y=n_motus), stat= 'identity', fill="#457277")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA))+
  labs(x="Number of sites",y="Number of families")

# --------------------- # 
# station scale 

length(unique(df_all_filters$station))

# Count
family_station <- df_all_filters %>%
  filter(!is.na(family_name_corrected)) %>%
  group_by(family_name_corrected) %>%
  summarise(n = n_distinct(station)) %>%
  ungroup() %>%
  group_by(n) %>%
  summarise(n_motus = n_distinct(family_name_corrected)) %>%
  mutate(rank = "station")
save(family_station, file="Rdata/rarete_family_station.rdata")

# plot
p_family_station <- fct_barplot(dataset, rank_x="station", rank_y = "families", n_tot = Nfamily, color = "grey")
p_family_station

ggplot(family_station)+
  geom_bar(aes(x=n, y=n_motus), stat= 'identity', fill="#457277")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA))+
  labs(x="Number of stations",y="Number of families")

# --------------------- # 
# sample scale 

length(unique(df_all_filters$sample_name_all_pcr))

# Count
dataset <- df_all_filters %>%
  filter(!is.na(family_name_corrected)) %>%
  group_by(family_name_corrected) %>%
  summarise(n = n_distinct(sample_name_all_pcr)) %>%
  ungroup() %>%
  group_by(n) %>%
  summarise(n_motus = n_distinct(family_name_corrected)) %>%
  mutate(rank = "sample")

# plot
p_family_sample<- fct_barplot(dataset, rank_x="sample", rank_y = "families", n_tot = Nfamily, color = "grey")
p_family_sample

# --------------------- # 
# all scales


ggarrange(#p_motu_region, 
          #p_motu_site + rremove("ylab"), 
          p_family_station, 
          p_family_sample+ rremove("ylab"),
          nrow=1, ncol=2)

ggsave("outputs/06_Upset_plots_Histograms_motus_family/repartition_family_regions_site_station.png", height = 10, width = 15)


ggarrange(p_motu_station, 
          p_motu_sample+ rremove("ylab"),
          p_family_station, 
          p_family_sample+ rremove("ylab"),
          nrow=2, ncol=2)

ggsave("outputs/06_Upset_plots_Histograms_motus_family/histogrammes_rarete.png", height = 15, width = 15)


# --------------------------------------------------------------------- # 
#### Genus unit  ----

# --------------------- # 
# Region scale 

length(unique(df_all_filters$province))

# Count
dataset <- df_all_filters %>%
  filter(!is.na(new_genus_name))  %>%
  group_by(new_genus_name) %>%
  summarise(n = n_distinct(province)) %>%
  ungroup() %>%
  group_by(n) %>%
  summarise(n_motus = n_distinct(new_genus_name)) %>%
  mutate(rank = "province")

# plot
p_motu_province <- fct_barplot(dataset, rank_x="province", rank_y = "genus", n_tot = Ngenus)
p_motu_province

# --------------------- # 
# site scale 

length(unique(df_all_filters$site35))

# Count
dataset <- df_all_filters %>%
  filter(!is.na(new_genus_name))  %>%
  group_by(new_genus_name) %>%
  summarise(n = n_distinct(site35)) %>%
  ungroup() %>%
  group_by(n) %>%
  summarise(n_motus = n_distinct(new_genus_name)) %>%
  mutate(rank = "site35")

# plot
p_motu_site <- fct_barplot(dataset, rank_x="site35", rank_y = "genus", n_tot = Ngenus)
p_motu_site

# --------------------- # 
# station scale 

length(unique(df_all_filters$station))

# Count
dataset <- df_all_filters %>%
  filter(!is.na(new_genus_name))  %>%
  group_by(new_genus_name) %>%
  summarise(n = n_distinct(station)) %>%
  ungroup() %>%
  group_by(n) %>%
  summarise(n_motus = n_distinct(new_genus_name)) %>%
  mutate(rank = "station")

# plot
p_motu_station <- fct_barplot(dataset, rank_x="station", rank_y = "genus", n_tot = Ngenus)
p_motu_station


# --------------------- # 
# sample scale 

length(unique(df_all_filters$sample_name_all_pcr))

# Count
dataset <- df_all_filters %>%
  filter(!is.na(new_genus_name))  %>%
  group_by(new_genus_name) %>%
  summarise(n = n_distinct(sample_name_all_pcr)) %>%
  ungroup() %>%
  group_by(n) %>%
  summarise(n_motus = n_distinct(new_genus_name)) %>%
  mutate(rank = "sample_name")

# plot
p_motu_sample <- fct_barplot(dataset, rank_x="sample", rank_y = "genus", n_tot = Ngenus)
p_motu_sample

# --------------------- # 
# all scales

ggarrange(p_motu_province, 
          p_motu_site + rremove("ylab"), 
          p_motu_station, 
          p_motu_sample+ rremove("ylab"),
          nrow=2, ncol=2)

ggsave("outputs/06_Upset_plots_Histograms_motus_family/repartition_genus_regions_site_station.png", height = 12, width = 12)



