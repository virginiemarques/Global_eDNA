
## calculate number of unique motus and families in each site and station of Lengguru and Caribbean

library(tidyverse)
library(reshape2)
library(lisa)
library(gridExtra)
library(grid)
library(cowplot)
library(ggplot2)
library(ggpubr)


load("Rdata/02_clean_all.Rdata")

'%ni%' <- Negate("%in%")

df_all_filters <- df_all_filters %>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")) %>%
  filter(sample_method !="niskin" & province!="Tropical_East_Pacific" & comment %ni% c("Distance decay 600m", "Distance decay 300m") & station!="glorieuse_distance_300m")%>%
  filter(project != "SEAMOUNTS") %>% 
  filter(habitat_type %ni% c("BAIE", "Sommet"))

# calculate nb of reads per family at global scale

df_all_filters <- df_all_filters %>%
  filter(!is.na(new_family_name))

family <- unique(df_all_filters$new_family_name)
nb_reads_family <- data.frame(family=character(), reads=numeric(), stringsAsFactors = FALSE)

for (i in 1:length(family)) {
  df <- subset(df_all_filters, new_family_name==family[i])
  nb_reads_family[i,1] <- family[i]
  nb_reads_family[i,2] <- sum(df$count_reads_all_pcr)
}


nb_reads_family <- arrange(nb_reads_family, desc(reads))
save(nb_reads_family, file = "Rdata/nb_reads_per_family_global.Rdata")


## Lengguru data

lengguru <- df_all_filters %>%
  filter(province=="Western_Coral_Triangle")

  # total MOTUs and family richness in Lengguru region
rich_tot_lengguru <- data.frame(motu=numeric(1), genus=numeric(1), family=numeric(1), province="Central_IndoPacific")
rich_tot_lengguru$motu <- lengguru %>% 
  summarise(n = n_distinct(sequence))
rich_tot_lengguru$genus <- lengguru %>% 
  summarise(n = n_distinct(new_genus_name))
rich_tot_lengguru$family <- lengguru %>% 
  summarise(n = n_distinct(new_family_name))


  # calculate unique motus and families at each station   
station <- c(unique(lengguru$station))

rich_station_lengguru <- data.frame(province="Central_IndoPacific", site=character(length(station)), station=character(length(station)), motu=numeric(length(station)), genus=numeric(length(station)), family=numeric(length(station)), stringsAsFactors = FALSE)

for (i in 1:length(station)) {
  s <- unique(lengguru[lengguru$station == station[i],]$site35)
  st <- station[i]
  motu <- lengguru[lengguru$station == station[i],] %>%
    summarise(n = n_distinct(sequence))
  gen <- lengguru[lengguru$station == station[i],] %>%
    summarise(n = n_distinct(new_genus_name))
  fam <- lengguru[lengguru$station == station[i],] %>%
    summarise(n = n_distinct(new_family_name))
  rich_station_lengguru[i,2] <- s
  rich_station_lengguru[i,3] <- st
  rich_station_lengguru[i,4] <- motu
  rich_station_lengguru[i,5] <- gen
  rich_station_lengguru[i,6] <- fam
}

  # calculate unique motus and families in each sample  
sample <- c(unique(lengguru$sample_name_all_pcr))

rich_sample_lengguru <- data.frame(province="Central_IndoPacific", site=character(length(sample)), sample=character(length(sample)), motu=numeric(length(sample)), genus=numeric(length(sample)), family=numeric(length(sample)), stringsAsFactors = FALSE)

for (i in 1:length(sample)) {
  s <- unique(lengguru[lengguru$sample_name_all_pcr == sample[i],]$site35)
  sa <- sample[i]
  motu <- lengguru[lengguru$sample_name_all_pcr == sample[i],] %>%
    summarise(n = n_distinct(sequence))
  gen <- lengguru[lengguru$sample_name_all_pcr == sample[i],] %>%
    summarise(n = n_distinct(new_genus_name))
  fam <- lengguru[lengguru$sample_name_all_pcr == sample[i],] %>%
    summarise(n = n_distinct(new_family_name))
  rich_sample_lengguru[i,2] <- s
  rich_sample_lengguru[i,3] <- sa
  rich_sample_lengguru[i,4] <- motu
  rich_sample_lengguru[i,5] <- gen
  rich_sample_lengguru[i,6] <- fam
}

  # calculate unique motus and families at each site (calculate mean and sd per station or sample)
site <- c(unique(lengguru$site35))

rich_site_lengguru <- data.frame(province="Western_Coral_Triangle", site=character(length(site)), motu=numeric(length(site)), genus=numeric(length(site)), family=numeric(length(site)), mean_motu=numeric(length(site)), sd_motu=numeric(length(site)), mean_genus=numeric(length(site)), sd_genus=numeric(length(site)), mean_family=numeric(length(site)), sd_family=numeric(length(site)), stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  s <- site[i]
  motu <- lengguru[lengguru$site35 == site[i],] %>%
    summarise(n = n_distinct(sequence))
  mm <- mean(rich_station_lengguru[rich_station_lengguru$site == site[i],]$motu)
  sdm <- sd(rich_station_lengguru[rich_station_lengguru$site == site[i],]$motu)
  gen <- lengguru[lengguru$site35 == site[i],] %>%
    summarise(n = n_distinct(new_genus_name))
  mg <- mean(rich_station_lengguru[rich_station_lengguru$site == site[i],]$genus)
  sdg <- sd(rich_station_lengguru[rich_station_lengguru$site == site[i],]$genus)
  fam <- lengguru[lengguru$site35 == site[i],] %>%
    summarise(n = n_distinct(new_family_name))
  mf <- mean(rich_station_lengguru[rich_station_lengguru$site == site[i],]$family)
  sdf <- sd(rich_station_lengguru[rich_station_lengguru$site == site[i],]$family)
  rich_site_lengguru[i,2] <- s
  rich_site_lengguru[i,3] <- motu
  rich_site_lengguru[i,4] <- gen
  rich_site_lengguru[i,5] <- fam
  rich_site_lengguru[i,6] <- mm
  rich_site_lengguru[i,7] <- sdm
  rich_site_lengguru[i,8] <- mg
  rich_site_lengguru[i,9] <- sdg
  rich_site_lengguru[i,10] <- mf
  rich_site_lengguru[i,11] <- sdf
}

write.csv(rich_site_lengguru, "outputs/05_richness_motus_families/richness_lengguru.csv", row.names = FALSE)



## Caribbean data

caribbean <- df_all_filters %>%
  filter(province=="Tropical_Northwestern_Atlantic")

  # total MOTUs and family richness in Caribbean region
rich_tot_caribbean <- data.frame(province="Tropical_Northwestern_Atlantic", motu=numeric(1), genus=numeric(1), family=numeric(1))
rich_tot_caribbean$motu <- caribbean %>% 
  summarise(n = n_distinct(sequence))
rich_tot_caribbean$genus <- caribbean %>% 
  summarise(n = n_distinct(new_genus_name))
rich_tot_caribbean$family <- caribbean %>% 
  summarise(n = n_distinct(new_family_name))


  # calculate unique motus and families at each station   
station <- c(unique(caribbean$station))

rich_station_caribbean <- data.frame(province="Tropical_Northwestern_Atlantic", site=character(length(station)), station=character(length(station)), motu=numeric(length(station)), genus=numeric(length(station)), family=numeric(length(station)), stringsAsFactors = FALSE)

for (i in 1:length(station)) {
  s <- unique(caribbean[caribbean$station == station[i],]$site35)
  st <- station[i]
  motu <- caribbean[caribbean$station == station[i],] %>%
    summarise(n = n_distinct(sequence))
  gen <- caribbean[caribbean$station == station[i],] %>%
    summarise(n = n_distinct(new_genus_name))
  fam <- caribbean[caribbean$station == station[i],] %>%
    summarise(n = n_distinct(new_family_name))
  rich_station_caribbean[i,2] <- s
  rich_station_caribbean[i,3] <- st
  rich_station_caribbean[i,4] <- motu
  rich_station_caribbean[i,5] <- gen
  rich_station_caribbean[i,6] <- fam
}

  # calculate unique motus and families at each sample   
sample <- c(unique(caribbean$sample_name_all_pcr))

rich_sample_caribbean <- data.frame(province="Tropical_Northwestern_Atlantic", site=character(length(sample)), sample=character(length(sample)), motu=numeric(length(sample)), genus=numeric(length(sample)), family=numeric(length(sample)), stringsAsFactors = FALSE)

for (i in 1:length(sample)) {
  s <- unique(caribbean[caribbean$sample_name_all_pcr == sample[i],]$site35)
  sa <- sample[i]
  motu <- caribbean[caribbean$sample_name_all_pcr == sample[i],] %>%
    summarise(n = n_distinct(sequence))
  gen <- caribbean[caribbean$sample_name_all_pcr == sample[i],] %>%
    summarise(n = n_distinct(new_genus_name))
  fam <- caribbean[caribbean$sample_name_all_pcr == sample[i],] %>%
    summarise(n = n_distinct(new_family_name))
  rich_sample_caribbean[i,2] <- s
  rich_sample_caribbean[i,3] <- sa
  rich_sample_caribbean[i,4] <- motu
  rich_sample_caribbean[i,5] <- gen
  rich_sample_caribbean[i,6] <- fam
}
  
  # calculate unique motus and families at each site (calculate mean and sd per station or sample) 
site <- c(unique(caribbean$site35))

rich_site_caribbean <- data.frame(province="Tropical_Northwestern_Atlantic", site=character(length(site)), motu=numeric(length(site)), genus=numeric(length(site)), family=numeric(length(site)), mean_motu=numeric(length(site)), sd_motu=numeric(length(site)), mean_genus=numeric(length(site)), sd_genus=numeric(length(site)), mean_family=numeric(length(site)), sd_family=numeric(length(site)), stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  s <- site[i]
  motu <- caribbean[caribbean$site35 == site[i],] %>%
    summarise(n = n_distinct(sequence))
  mm <- mean(rich_station_caribbean[rich_station_caribbean$site == site[i],]$motu)
  sdm <- sd(rich_station_caribbean[rich_station_caribbean$site == site[i],]$motu)
  gen <- caribbean[caribbean$site35 == site[i],] %>%
    summarise(n = n_distinct(new_genus_name))
  mg <- mean(rich_station_caribbean[rich_station_caribbean$site == site[i],]$genus)
  sdg <- sd(rich_station_caribbean[rich_station_caribbean$site == site[i],]$genus)
  fam <- caribbean[caribbean$site35 == site[i],] %>%
    summarise(n = n_distinct(new_family_name))
  mf <- mean(rich_station_caribbean[rich_station_caribbean$site == site[i],]$family)
  sdf <- sd(rich_station_caribbean[rich_station_caribbean$site == site[i],]$family)
  rich_site_caribbean[i,2] <- s
  rich_site_caribbean[i,3] <- motu
  rich_site_caribbean[i,4] <- gen
  rich_site_caribbean[i,5] <- fam
  rich_site_caribbean[i,6] <- mm
  rich_site_caribbean[i,7] <- sdm
  rich_site_caribbean[i,8] <- mg
  rich_site_caribbean[i,9] <- sdg
  rich_site_caribbean[i,10] <- mf
  rich_site_caribbean[i,11] <- sdf
}

write.csv(rich_site_caribbean, "outputs/05_richness_motus_families/richness_caribbean.csv", row.names = FALSE)



## Fakarava data

fakarava <- df_all_filters %>%
  filter(province=="Southeast_Polynesia")

  # calculate unique motus and families at each station   
station <- c(unique(fakarava$station))

rich_station_fakarava <- data.frame(province="Southeast_Polynesia", site="fakarava", station=character(length(station)), motu=numeric(length(station)), genus=numeric(length(station)), family=numeric(length(station)), stringsAsFactors = FALSE)

for (i in 1:length(station)) {
  st <- station[i]
  motu <- fakarava[fakarava$station == station[i],] %>%
    summarise(n = n_distinct(sequence))
  gen <- fakarava[fakarava$station == station[i],] %>%
    summarise(n = n_distinct(new_genus_name))
  fam <- fakarava[fakarava$station == station[i],] %>%
    summarise(n = n_distinct(new_family_name))
  rich_station_fakarava[i,3] <- st
  rich_station_fakarava[i,4] <- motu
  rich_station_fakarava[i,5] <- gen
  rich_station_fakarava[i,6] <- fam
}
colnames(rich_station_fakarava) <- c("province", "site", "station", "motu", "genus", "family")



# calculate unique motus and families in each sample  
sample <- c(unique(fakarava$sample_name_all_pcr))

rich_sample_fakarava <- data.frame(province="Southeast_Polynesia", site="fakarava", sample=character(length(sample)), motu=numeric(length(sample)), genus=numeric(length(sample)), family=numeric(length(sample)), stringsAsFactors = FALSE)

for (i in 1:length(sample)) {
  sa <- sample[i]
  motu <- fakarava[fakarava$sample_name_all_pcr == sample[i],] %>%
    summarise(n = n_distinct(sequence))
  gen <- fakarava[fakarava$sample_name_all_pcr == sample[i],] %>%
    summarise(n = n_distinct(new_genus_name))
  fam <- fakarava[fakarava$sample_name_all_pcr == sample[i],] %>%
    summarise(n = n_distinct(new_family_name))
  rich_sample_fakarava[i,3] <- sa
  rich_sample_fakarava[i,4] <- motu
  rich_sample_fakarava[i,5] <- gen
  rich_sample_fakarava[i,6] <- fam
}
colnames(rich_sample_fakarava) <- c("province", "site", "sample", "motu", "genus", "family")



  # total MOTUs and family richness in fakarava region (=site) (calculate mean and sd per station or sample) 
rich_site_fakarava <- data.frame(province="Southeast_Polynesia",
                                 site="fakarava",
                                 motu=fakarava %>% 
                                   summarise(n = n_distinct(sequence)),
                                 genus=fakarava %>% 
                                   summarise(n = n_distinct(new_genus_name)),
                                 family=fakarava %>% 
                                   summarise(n = n_distinct(new_family_name)),
                                 mean_motu=mean(rich_station_fakarava$motu),
                                 sd_motu=sd(rich_station_fakarava$motu),
                                 mean_genus=mean(rich_station_fakarava$genus),
                                 sd_genus=sd(rich_station_fakarava$genus),
                                 mean_family=mean(rich_station_fakarava$family),
                                 sd_family=sd(rich_station_fakarava$family),
                                 stringsAsFactors = FALSE)

colnames(rich_site_fakarava) <- c("province", "site", "motu", "genus", "family", "mean_motu", "sd_motu", "mean_genus", "sd_genus", "mean_family", "sd_family")

write.csv(rich_site_fakarava, "outputs/05_richness_motus_families/richness_fakarava.csv", row.names = FALSE)


## Eparses data

eparse <- df_all_filters %>%
  filter(province=="Western_Indian_Ocean")

  # total MOTUs and family richness in eparse region
rich_tot_eparse <- data.frame(motu=numeric(1), genus=numeric(1), family=numeric(1), province="Western_Indian_Ocean")
rich_tot_eparse$motu <- eparse %>% 
  summarise(n = n_distinct(sequence))
rich_tot_eparse$genus <- eparse %>% 
  summarise(n = n_distinct(new_genus_name))
rich_tot_eparse$family <- eparse %>% 
  summarise(n = n_distinct(new_family_name))


  # calculate unique motus and families at each station   
station <- c(unique(eparse$station))

rich_station_eparse <- data.frame(province="Western_Indian_Ocean", site=character(length(station)), station=character(length(station)), motu=numeric(length(station)), genus=numeric(length(station)), family=numeric(length(station)), stringsAsFactors = FALSE)

for (i in 1:length(station)) {
  s <- unique(eparse[eparse$station == station[i],]$site35)
  st <- station[i]
  motu <- eparse[eparse$station == station[i],] %>%
    summarise(n = n_distinct(sequence))
  gen <- eparse[eparse$station == station[i],] %>%
    summarise(n = n_distinct(new_genus_name))
  fam <- eparse[eparse$station == station[i],] %>%
    summarise(n = n_distinct(new_family_name))
  rich_station_eparse[i,2] <- s
  rich_station_eparse[i,3] <- st
  rich_station_eparse[i,4] <- motu
  rich_station_eparse[i,5] <- gen
  rich_station_eparse[i,6] <- fam
}


  # calculate unique motus and families in each sample   
sample <- c(unique(eparse$sample_name_all_pcr))

rich_sample_eparse <- data.frame(province="Western_Indian_Ocean", site=character(length(sample)), sample=character(length(sample)), motu=numeric(length(sample)), genus=numeric(length(sample)), family=numeric(length(sample)), stringsAsFactors = FALSE)

for (i in 1:length(sample)) {
  s <- unique(eparse[eparse$sample_name_all_pcr == sample[i],]$site35)
  sa <- sample[i]
  motu <- eparse[eparse$sample_name_all_pcr == sample[i],] %>%
    summarise(n = n_distinct(sequence))
  gen <- eparse[eparse$sample_name_all_pcr == sample[i],] %>%
    summarise(n = n_distinct(new_genus_name))
  fam <- eparse[eparse$sample_name_all_pcr == sample[i],] %>%
    summarise(n = n_distinct(new_family_name))
  rich_sample_eparse[i,2] <- s
  rich_sample_eparse[i,3] <- sa
  rich_sample_eparse[i,4] <- motu
  rich_sample_eparse[i,5] <- gen
  rich_sample_eparse[i,6] <- fam
}


  # calculate unique motus and families at each site
site <- c(unique(eparse$site35))

rich_site_eparse <- data.frame(province="Western_Indian_Ocean", site=character(length(site)), motu=numeric(length(site)), genus=numeric(length(site)), family=numeric(length(site)), mean_motu=numeric(length(site)), sd_motu=numeric(length(site)), mean_genus=numeric(length(site)), sd_genus=numeric(length(site)), mean_family=numeric(length(site)), sd_family=numeric(length(site)), stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  s <- site[i]
  motu <- eparse[eparse$site35 == site[i],] %>%
    summarise(n = n_distinct(sequence))
  mm <- mean(rich_station_eparse[rich_station_eparse$site == site[i],]$motu)
  sdm <- sd(rich_station_eparse[rich_station_eparse$site == site[i],]$motu)
  gen <- eparse[eparse$site35 == site[i],] %>%
    summarise(n = n_distinct(new_genus_name))
  mg <- mean(rich_station_eparse[rich_station_eparse$site == site[i],]$genus)
  sdg <- sd(rich_station_eparse[rich_station_eparse$site == site[i],]$genus)
  fam <- eparse[eparse$site35 == site[i],] %>%
    summarise(n = n_distinct(new_family_name))
  mf <- mean(rich_station_eparse[rich_station_eparse$site == site[i],]$family)
  sdf <- sd(rich_station_eparse[rich_station_eparse$site == site[i],]$family)
  rich_site_eparse[i,2] <- s
  rich_site_eparse[i,3] <- motu
  rich_site_eparse[i,4] <- gen
  rich_site_eparse[i,5] <- fam
  rich_site_eparse[i,6] <- mm
  rich_site_eparse[i,7] <- sdm
  rich_site_eparse[i,8] <- mg
  rich_site_eparse[i,9] <- sdg
  rich_site_eparse[i,10] <- mf
  rich_site_eparse[i,11] <- sdf
}

write.csv(rich_site_eparse, "outputs/05_richness_motus_families/richness_eparse.csv", row.names = FALSE)

## caledonia data

caledonia <- df_all_filters %>%
  filter(province=="Tropical_Southwestern_Pacific")

# total MOTUs and family richness in caledonia region
rich_tot_caledonia <- data.frame(motu=numeric(1), genus=numeric(1), family=numeric(1), province="Tropical_Southwestern_Pacific")
rich_tot_caledonia$motu <- caledonia %>% 
  summarise(n = n_distinct(sequence))
rich_tot_caledonia$genus <- caledonia %>% 
  summarise(n = n_distinct(new_genus_name))
rich_tot_caledonia$family <- caledonia %>% 
  summarise(n = n_distinct(new_family_name))


# calculate unique motus and families at each station   
station <- c(unique(caledonia$station))

rich_station_caledonia <- data.frame(province="Tropical_Southwestern_Pacific", site=character(length(station)), station=character(length(station)), motu=numeric(length(station)), genus=numeric(length(station)), family=numeric(length(station)), stringsAsFactors = FALSE)

for (i in 1:length(station)) {
  s <- unique(caledonia[caledonia$station == station[i],]$site35)
  st <- station[i]
  motu <- caledonia[caledonia$station == station[i],] %>%
    summarise(n = n_distinct(sequence))
  gen <- caledonia[caledonia$station == station[i],] %>%
    summarise(n = n_distinct(new_genus_name))
  fam <- caledonia[caledonia$station == station[i],] %>%
    summarise(n = n_distinct(new_family_name))
  rich_station_caledonia[i,2] <- s
  rich_station_caledonia[i,3] <- st
  rich_station_caledonia[i,4] <- motu
  rich_station_caledonia[i,5] <- gen
  rich_station_caledonia[i,6] <- fam
}

# calculate unique motus and families in each sample  
sample <- c(unique(caledonia$sample_name_all_pcr))

rich_sample_caledonia <- data.frame(province="Tropical_Southwestern_Pacific", site=character(length(sample)), sample=character(length(sample)), motu=numeric(length(sample)), genus=numeric(length(sample)), family=numeric(length(sample)), stringsAsFactors = FALSE)

for (i in 1:length(sample)) {
  s <- unique(caledonia[caledonia$sample_name_all_pcr == sample[i],]$site35)
  sa <- sample[i]
  motu <- caledonia[caledonia$sample_name_all_pcr == sample[i],] %>%
    summarise(n = n_distinct(sequence))
  gen <- caledonia[caledonia$sample_name_all_pcr == sample[i],] %>%
    summarise(n = n_distinct(new_genus_name))
  fam <- caledonia[caledonia$sample_name_all_pcr == sample[i],] %>%
    summarise(n = n_distinct(new_family_name))
  rich_sample_caledonia[i,2] <- s
  rich_sample_caledonia[i,3] <- sa
  rich_sample_caledonia[i,4] <- motu
  rich_sample_caledonia[i,5] <- gen
  rich_sample_caledonia[i,6] <- fam
}

# calculate unique motus and families at each site (calculate mean and sd per station or sample)
site <- c(unique(caledonia$site35))

rich_site_caledonia <- data.frame(province="Tropical_Southwestern_Pacific", site=character(length(site)), motu=numeric(length(site)), genus=numeric(length(site)), family=numeric(length(site)), mean_motu=numeric(length(site)), sd_motu=numeric(length(site)), mean_genus=numeric(length(site)), sd_genus=numeric(length(site)), mean_family=numeric(length(site)), sd_family=numeric(length(site)), stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  s <- site[i]
  motu <- caledonia[caledonia$site35 == site[i],] %>%
    summarise(n = n_distinct(sequence))
  mm <- mean(rich_station_caledonia[rich_station_caledonia$site == site[i],]$motu)
  sdm <- sd(rich_station_caledonia[rich_station_caledonia$site == site[i],]$motu)
  gen <- caledonia[caledonia$site35 == site[i],] %>%
    summarise(n = n_distinct(new_genus_name))
  mg <- mean(rich_station_caledonia[rich_station_caledonia$site == site[i],]$genus)
  sdg <- sd(rich_station_caledonia[rich_station_caledonia$site == site[i],]$genus)
  fam <- caledonia[caledonia$site35 == site[i],] %>%
    summarise(n = n_distinct(new_family_name))
  mf <- mean(rich_station_caledonia[rich_station_caledonia$site == site[i],]$family)
  sdf <- sd(rich_station_caledonia[rich_station_caledonia$site == site[i],]$family)
  rich_site_caledonia[i,2] <- s
  rich_site_caledonia[i,3] <- motu
  rich_site_caledonia[i,4] <- gen
  rich_site_caledonia[i,5] <- fam
  rich_site_caledonia[i,6] <- mm
  rich_site_caledonia[i,7] <- sdm
  rich_site_caledonia[i,8] <- mg
  rich_site_caledonia[i,9] <- sdg
  rich_site_caledonia[i,10] <- mf
  rich_site_caledonia[i,11] <- sdf
}

write.csv(rich_site_caledonia, "outputs/05_richness_motus_families/richness_caledonia.csv", row.names = FALSE)


# plot motu, genus and family richness per region


rich_station_total <- rbind(rich_station_caribbean, rich_station_eparse, rich_station_fakarava, rich_station_lengguru, rich_station_caledonia)
rich_station_total$province <- gsub("Southeast_Polynesia", "SP", rich_station_total$province)

ggplot(rich_station_total)+
  geom_violin(aes(site, motu), color="#D2981A", fill="#D2981A", position = position_nudge(-0.3, 0))+
  geom_boxplot(aes(site, motu), width=0.1, position = position_nudge(-0.3, 0))+
  geom_violin(aes(site, genus), color="#A53E1F", fill="#A53E1F")+
  geom_boxplot(aes(site, genus), width=0.1)+
  geom_violin(aes(site, family), color="#457277", fill="#457277", position = position_nudge(0.3, 0))+
  geom_boxplot(aes(site, family), width=0.1, position = position_nudge(0.3, 0))+
  facet_grid(~province, scales = "free_x", space = "free_x")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5))+
  labs(x="Sites", y = "Richness")

ggsave("outputs/05_richness_motus_families/violin_plot_all.png", width = 20, height = 8)


motu <- ggplot(rich_station_total)+
  geom_violin(aes(site, motu), color="#D2981A", fill="#D2981A")+
  geom_boxplot(aes(site, motu), width=0.1)+
  facet_grid(~province, scales = "free_x", space = "free_x")+
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5))+
  labs(x="", y = "MOTUs richness")+
  theme_bw()+
  theme(axis.text.x = element_blank())


genus <- ggplot(rich_station_total)+
  geom_violin(aes(site, genus), color="#A53E1F", fill="#A53E1F")+
  geom_boxplot(aes(site, genus), width=0.1)+
  facet_grid(~province, scales = "free_x", space = "free_x")+
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5))+
  labs(x="", y = "Genus richness")+
  theme_bw()+
  theme(axis.text.x = element_blank())

family <- ggplot(rich_station_total)+
  geom_violin(aes(site, family), color="#457277", fill="#457277")+
  geom_boxplot(aes(site, family), width=0.1)+
  facet_grid(~province, scales = "free_x", space = "free_x")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5))+
  labs(x="", y = "Family richness")

ggarrange(motu, family, ncol=1, nrow=2, labels=c("A", "B"), heights = c(1,1.3))

ggsave("outputs/05_richness_motus_families/violin_plot_per_rank.png", width = 12, height = 8)
ggsave("outputs/00_Figures_for_paper/Extended_Data/ED_Figure6.png", width = 7, height = 6)





## load useful metadata 
metadata <- read.csv("metadata/Metadata_eDNA_global_V6.csv", sep=",", stringsAsFactors = FALSE)
metadata <- metadata %>%
  filter(province%in%c("Western_Coral_Triangle", "Southeast_Polynesia", "Tropical_Northwestern_Atlantic", "Western_Indian_Ocean", "Tropical_Southwestern_Pacific"))%>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")) %>%
  filter(sample_method !="niskin" & province!="Tropical_East_Pacific" & comment %ni% c("Distance decay 600m", "Distance decay 300m") & station!="glorieuse_distance_300m")%>%
  filter(project != "SEAMOUNTS") %>% 
  filter(habitat_type %ni% c("BAIE", "Sommet"))%>%
  subset(habitat=="marine")



metadata <- metadata[,c("site35", "station", "latitude_start_clean", "longitude_start_clean", "dist_to_coast..m.", "dist_to_CT")]



rich_site <- rbind(rich_site_caribbean, rich_site_fakarava, rich_site_lengguru, rich_site_eparse, rich_site_caledonia)
rich_site <- left_join(rich_site, metadata[,c("site35", "dist_to_CT")], by=c("site"="site35"))

# Compile data - calculate mean distance to CT for each region 
to_merge <- rich_site %>%
  dplyr::select(province, dist_to_CT) %>%
  dplyr::group_by(province) %>%
  dplyr::summarise(dist_to_CT = mean(dist_to_CT))
save(rich_site, file="Rdata/richness_station_site_region.rdata")

### MOTUs 
# DF for plot - count number of motus per region and join precedent df to get distance to CT
all_province <- df_all_filters %>%
  group_by(province) %>%
  summarise(n_motus = n_distinct(sequence)) %>%
  left_join(., to_merge)

save(all_province, file = "Rdata/richness_motu_region.rdata")
# Plot
all_motus <- ggplot(rich_site, aes(col=province))+
  geom_jitter(aes(x=dist_to_CT, y=motu), shape=17, size=2, alpha=0.7, show.legend = FALSE) +
  geom_jitter(aes(x=dist_to_CT, y=mean_motu), shape=20, size=4, alpha=0.7, show.legend = FALSE) +
  geom_errorbar(aes(x=dist_to_CT, ymin=mean_motu-sd_motu, ymax=mean_motu+sd_motu), show.legend = FALSE, alpha=0.7)+
  geom_bar(data=all_province, aes(x=dist_to_CT, y=n_motus), stat= 'identity', orientation = "x", alpha=0.05, show.legend = FALSE) + 
  scale_color_manual(values=c("#E5A729", "#8AAE8A", "#4F4D1D", "#863b34", "#C67052"))+ 
  theme(legend.position = "none")+
  theme_bw()+
  theme(axis.title.y = element_text(size = 10, face = "bold"), plot.margin=unit(c(0.2,0.4,0,0.1), "cm"), 
        text = element_text(size=12))+
  labs(x="",y="MOTU richness")



### FAMILIES
# DF for plot - count number of families per region and join precedent df to get distance to CT
all_province <- df_all_filters %>%
  group_by(province) %>%
  summarise(n_family = n_distinct(new_family_name)) %>%
  left_join(., to_merge)
save(all_province, file = "Rdata/richness_family_region.rdata")

# Plot
all_family <- ggplot(rich_site, aes(col=province))+
  geom_jitter(aes(x=dist_to_CT, y=family), shape=17, size=2, alpha=0.7, show.legend = FALSE) +
  geom_jitter(aes(x=dist_to_CT, y=mean_family), shape=20, size=4, alpha=0.7, show.legend = FALSE) +
  geom_errorbar(aes(x=dist_to_CT, ymin=mean_family-sd_family, ymax=mean_family+sd_family), show.legend = FALSE, alpha=0.7)+
  geom_bar(data=all_province, aes(x=dist_to_CT, y=n_family), stat= 'identity', orientation = "x", alpha=0.05, show.legend = FALSE) + 
  scale_color_manual(values=c("#E5A729", "#8AAE8A", "#4F4D1D", "#863b34", "#C67052"))+ 
  theme(legend.position = "none")+
  theme_bw()+
  theme(axis.title.y = element_text(size = 10, face = "bold"), plot.margin=unit(c(0.2,0.4,0,0.1), "cm"), 
        text = element_text(size=12))+
  labs(x="",y="Family richness")

 

plot <- ggarrange(all_motus, all_family, ncol = 2, nrow=1)
x.grob <- textGrob("Distance to Coral Triangle (km, W-E)", 
                   gp=gpar(fontface="bold", col="black", fontsize=12), vjust = -0.5)

plot_all_rich_site <- grid.arrange(plot, bottom=x.grob)
save(plot_all_rich_site, file = "Rdata/plot_richness_site~dist_CT.rdata")


library(pgirmess)
library(FSA)

kruskal.test(motu~province, data = rich_site)
kruskalmc(rich_site$motu~rich_site$province)
dt <- dunnTest(motu~province, data = rich_site)
kruskal(rich_site$mean_motu, rich_site$site, group=TRUE)$groups
kruskal(rich_site$motu, rich_site$province, group=TRUE)$groups


