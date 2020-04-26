# This script will explore the dominance of MOTUs across its taxonomic ranks (genus, families, species-proxy)
# Go to read numbers as well ? Rarefaction by lowest number of reads per station? 

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
load("Rdata/02_clean_all.Rdata")

#Remove estuary stations and deep niskin station
df_all_filters <- subset(df_all_filters, !(station %in% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")))
df_all_filters <- subset(df_all_filters, sample_method!="niskin")
df_all_filters <- subset(df_all_filters, region!="East_Pacific")

# on est d'accord qu'on enlève toujours les MOTUs non assignés au dessus de la famille dans ce MS? Faire les deux
#df_all_filters <- df_all_filters %>%
#  filter(!is.na(new_family_name))

# N global 
Nmotus <- length(unique(df_all_filters$sequence))
Nmotus
# 1047 MOTUs in total

# N families
Nfamily <- df_all_filters %>%
  filter(!is.na(new_family_name)) %>%
  distinct(new_family_name) %>% pull() %>% length()
Nfamily
# 125

# N genus
Ngenus <- df_all_filters %>%
  filter(!is.na(new_genus_name)) %>%
  distinct(new_genus_name) %>% pull() %>% length()
Ngenus
#300

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

length(unique(df_all_filters$region))

# Count
motu_region <- df_all_filters %>%
  group_by(sequence) %>%
  summarise(n = n_distinct(region)) %>%
  ungroup() %>%
  group_by(n) %>%
  summarise(n_motus = n_distinct(sequence)) %>%
  mutate(rank = "region")

# plot
p_motu_region <- fct_barplot(motu_region, rank_x="region", rank_y = "MOTUs")
p_motu_region

# 3 MOTUs are present in all 3 regions:
# * One non-ID Myctophidae 
# * One Katsuwonus pelamis	
# * Decapterus macarellus	

# --------------------- # 
# site scale 

length(unique(df_all_filters$site))

# Count
motu_site <- df_all_filters %>%
  group_by(sequence) %>%
  summarise(n = n_distinct(site)) %>%
  ungroup() %>%
  group_by(n) %>%
  summarise(n_motus = n_distinct(sequence))%>%
  mutate(rank = "site")

# plot
p_motu_site <- fct_barplot(motu_site, rank_x="site", rank_y = "MOTUs")
p_motu_site

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

# plot
p_motu_station <- fct_barplot(motu_station, rank_x="station", rank_y = "MOTUs")
p_motu_station

# --------------------- # 
# station scale 

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
p_motu_sample <- fct_barplot(motu_sample, rank_x="sample", rank_y = "MOTUs")
p_motu_sample

# --------------------- # 
# all scales

# Trial 1

ggarrange(p_motu_region, 
          p_motu_site + rremove("ylab"), 
          p_motu_station, 
          p_motu_sample+ rremove("ylab"),
          nrow=2, ncol=2)

ggsave("outputs/09_dominance_motus_ranks/repartition_MOTUs_regions_site_station.png", height = 12, width = 12)

# Trial 2
# Combine data 
all_motus <- rbind(motu_region, motu_site, motu_station, motu_sample)

# Plot
p_motu_all <- ggplot(all_motus, aes(x = n, y=n_motus)) + 
  geom_bar(stat="identity", color = "black", fill = "#56B4E9") + 
  geom_text(aes(label=n_motus), vjust=-0.5, color="black", size=3.5) +
  theme_classic() + 
  labs(x= "Number of spatial units", y = "Number of MOTUs") + 
  facet_wrap(~rank, scales = "free")

p_motu_all

# --------------------------------------------------------------------- # 
#### Family unit  ----

# --------------------- # 
# Region scale 

length(unique(df_all_filters$region))

# Count
motu_region <- df_all_filters %>%
  filter(!is.na(new_family_name)) %>%
  group_by(new_family_name) %>%
  summarise(n = n_distinct(region)) %>%
  ungroup() %>%
  group_by(n) %>%
  summarise(n_motus = n_distinct(new_family_name)) %>%
  mutate(rank = "region")

# plot
p_motu_region <- fct_barplot(motu_region, rank_x="region", rank_y = "families", n_tot = Nfamily)
p_motu_region

# --------------------- # 
# site scale 

length(unique(df_all_filters$site))

# Count
dataset <- df_all_filters %>%
  filter(!is.na(new_family_name)) %>%
  group_by(new_family_name) %>%
  summarise(n = n_distinct(site)) %>%
  ungroup() %>%
  group_by(n) %>%
  summarise(n_motus = n_distinct(new_family_name)) %>%
  mutate(rank = "site")

# plot
p_motu_site <- fct_barplot(dataset, rank_x="site", rank_y = "families", n_tot = Nfamily)
p_motu_site

# --------------------- # 
# station scale 

length(unique(df_all_filters$station))

# Count
dataset <- df_all_filters %>%
  filter(!is.na(new_family_name)) %>%
  group_by(new_family_name) %>%
  summarise(n = n_distinct(station)) %>%
  ungroup() %>%
  group_by(n) %>%
  summarise(n_motus = n_distinct(new_family_name)) %>%
  mutate(rank = "station")

# plot
p_motu_station <- fct_barplot(dataset, rank_x="station", rank_y = "families", n_tot = Nfamily)
p_motu_station

# --------------------- # 
# sample scale 

length(unique(df_all_filters$sample_name_all_pcr))

# Count
dataset <- df_all_filters %>%
  filter(!is.na(new_family_name)) %>%
  group_by(new_family_name) %>%
  summarise(n = n_distinct(sample_name_all_pcr)) %>%
  ungroup() %>%
  group_by(n) %>%
  summarise(n_motus = n_distinct(new_family_name)) %>%
  mutate(rank = "sample")

# plot
p_motu_sample<- fct_barplot(dataset, rank_x="sample", rank_y = "families", n_tot = Nfamily)
p_motu_sample

# --------------------- # 
# all scales

# Trial 1

ggarrange(p_motu_region, 
          p_motu_site + rremove("ylab"), 
          p_motu_station, 
          p_motu_sample+ rremove("ylab"),
          nrow=2, ncol=2)

ggsave("outputs/09_dominance_motus_ranks/repartition_family_regions_site_station.png", height = 12, width = 12)

# --------------------------------------------------------------------- # 
#### Genus unit  ----

# --------------------- # 
# Region scale 

length(unique(df_all_filters$region))

# Count
dataset <- df_all_filters %>%
  filter(!is.na(new_genus_name))  %>%
  group_by(new_genus_name) %>%
  summarise(n = n_distinct(region)) %>%
  ungroup() %>%
  group_by(n) %>%
  summarise(n_motus = n_distinct(new_genus_name)) %>%
  mutate(rank = "region")

# plot
p_motu_region <- fct_barplot(dataset, rank_x="region", rank_y = "genus", n_tot = Ngenus)
p_motu_region

# --------------------- # 
# site scale 

length(unique(df_all_filters$site))

# Count
dataset <- df_all_filters %>%
  filter(!is.na(new_genus_name))  %>%
  group_by(new_genus_name) %>%
  summarise(n = n_distinct(site)) %>%
  ungroup() %>%
  group_by(n) %>%
  summarise(n_motus = n_distinct(new_genus_name)) %>%
  mutate(rank = "site")

# plot
p_motu_site <- fct_barplot(dataset, rank_x="site", rank_y = "genus", n_tot = Ngenus)
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

ggarrange(p_motu_region, 
          p_motu_site + rremove("ylab"), 
          p_motu_station, 
          p_motu_sample+ rremove("ylab"),
          nrow=2, ncol=2)

ggsave("outputs/09_dominance_motus_ranks/repartition_genus_regions_site_station.png", height = 12, width = 12)



