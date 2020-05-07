## calculate number of unique motus and families in each site and station of Lengguru and Caribbean

library(tidyverse)
library(reshape2)
library(lisa)
library(gridExtra)
library(grid)
library(cowplot)
library(ggplot2)
library(ggpubr)

setwd("c:/Users/mathon/Desktop/linux/Global_eDNA/")
load("Rdata/02_clean_all.Rdata")

df_all_filters <- subset(df_all_filters, !(station %in% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")))
df_all_filters <- subset(df_all_filters, sample_method!="niskin")
df_all_filters <- subset(df_all_filters, region!="East_Pacific")
df_all_filters <- subset(df_all_filters, !(comment %in% c("Distance decay 600m", "Distance decay 300m")))
df_all_filters <- subset(df_all_filters, station!="glorieuse_distance_300m")


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
  filter(region=="Central_IndoPacific")

  # total MOTUs and family richness in Lengguru region
rich_tot_lengguru <- data.frame(motu=numeric(1), genus=numeric(1), family=numeric(1), region="Central_IndoPacific")
rich_tot_lengguru$motu <- lengguru %>% 
  summarise(n = n_distinct(sequence))
rich_tot_lengguru$genus <- lengguru %>% 
  summarise(n = n_distinct(new_genus_name))
rich_tot_lengguru$family <- lengguru %>% 
  summarise(n = n_distinct(new_family_name))


  # calculate unique motus and families at each station   
station <- c(unique(lengguru$station))

rich_station_lengguru <- data.frame(region="Central_IndoPacific", site=character(46), station=character(46), motu=numeric(46), genus=numeric(46), family=numeric(46), stringsAsFactors = FALSE)

for (i in 1:length(station)) {
  s <- unique(lengguru[lengguru$station == station[i],]$site)
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


  # calculate unique motus and families at each site
site <- c(unique(lengguru$site))

rich_site_lengguru <- data.frame(region="Central_Indo_Pacific", site=character(11), motu=numeric(11), genus=numeric(11), family=numeric(11), mean_motu=numeric(11), sd_motu=numeric(11), mean_genus=numeric(11), sd_genus=numeric(11), mean_family=numeric(11), sd_family=numeric(11), stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  s <- site[i]
  motu <- lengguru[lengguru$site == site[i],] %>%
    summarise(n = n_distinct(sequence))
  mm <- mean(rich_station_lengguru[rich_station_lengguru$site == site[i],]$motu)
  sdm <- sd(rich_station_lengguru[rich_station_lengguru$site == site[i],]$motu)
  gen <- lengguru[lengguru$site == site[i],] %>%
    summarise(n = n_distinct(new_genus_name))
  mg <- mean(rich_station_lengguru[rich_station_lengguru$site == site[i],]$genus)
  sdg <- sd(rich_station_lengguru[rich_station_lengguru$site == site[i],]$genus)
  fam <- lengguru[lengguru$site == site[i],] %>%
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

write.csv(rich_site_lengguru, "outputs/04_exploration_richness/richness_lengguru.csv", row.names = FALSE)


rich_site_lengguru_melt <- melt(rich_site_lengguru[,1:5])
rich_site_lengguru_melt[1:11,"mean"] <- rich_site_lengguru$mean_motu
rich_site_lengguru_melt[12:22,"mean"] <- rich_site_lengguru$mean_genus
rich_site_lengguru_melt[23:33,"mean"] <- rich_site_lengguru$mean_family

rich_site_lengguru_melt[1:11,"sd"] <- rich_site_lengguru$sd_motu
rich_site_lengguru_melt[12:22,"sd"] <- rich_site_lengguru$sd_family
rich_site_lengguru_melt[23:33,"sd"] <- rich_site_lengguru$sd_family






## Caribbean data

caribbean <- df_all_filters %>%
  filter(region=="Caribbean")

  # total MOTUs and family richness in Caribbean region
rich_tot_caribbean <- data.frame(region="Caribbean", motu=numeric(1), genus=numeric(1), family=numeric(1))
rich_tot_caribbean$motu <- caribbean %>% 
  summarise(n = n_distinct(sequence))
rich_tot_caribbean$genus <- caribbean %>% 
  summarise(n = n_distinct(new_genus_name))
rich_tot_caribbean$family <- caribbean %>% 
  summarise(n = n_distinct(new_family_name))


  # calculate unique motus and families at each station   
station <- c(unique(caribbean$station))

rich_station_caribbean <- data.frame(region="Caribbean", site=character(31), station=character(31), motu=numeric(31), genus=numeric(31), family=numeric(31), stringsAsFactors = FALSE)

for (i in 1:length(station)) {
  s <- unique(caribbean[caribbean$station == station[i],]$site)
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

  # calculate unique motus and families at each site
site <- c(unique(caribbean$site))

rich_site_caribbean <- data.frame(region="Caribbean", site=character(3), motu=numeric(3), genus=numeric(3), family=numeric(3), mean_motu=numeric(3), sd_motu=numeric(3), mean_genus=numeric(3), sd_genus=numeric(3), mean_family=numeric(3), sd_family=numeric(3), stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  s <- site[i]
  motu <- caribbean[caribbean$site == site[i],] %>%
    summarise(n = n_distinct(sequence))
  mm <- mean(rich_station_caribbean[rich_station_caribbean$site == site[i],]$motu)
  sdm <- sd(rich_station_caribbean[rich_station_caribbean$site == site[i],]$motu)
  gen <- caribbean[caribbean$site == site[i],] %>%
    summarise(n = n_distinct(new_genus_name))
  mg <- mean(rich_station_caribbean[rich_station_caribbean$site == site[i],]$genus)
  sdg <- sd(rich_station_caribbean[rich_station_caribbean$site == site[i],]$genus)
  fam <- caribbean[caribbean$site == site[i],] %>%
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

write.csv(rich_site_caribbean, "outputs/04_exploration_richness/richness_caribbean.csv", row.names = FALSE)


rich_site_caribbean_melt <- melt(rich_site_caribbean[,1:5])
rich_site_caribbean_melt[1:3,"mean"] <- rich_site_caribbean$mean_motu
rich_site_caribbean_melt[4:6,"mean"] <- rich_site_caribbean$mean_genus
rich_site_caribbean_melt[7:9,"mean"] <- rich_site_caribbean$mean_family

rich_site_caribbean_melt[1:3,"sd"] <- rich_site_caribbean$sd_motu
rich_site_caribbean_melt[4:6,"sd"] <- rich_site_caribbean$sd_genus
rich_site_caribbean_melt[7:9,"sd"] <- rich_site_caribbean$sd_family



## Fakarava data

fakarava <- df_all_filters %>%
  filter(region=="Central_Pacific")

  # calculate unique motus and families at each station   
station <- c(unique(fakarava$station))

rich_station_fakarava <- data.frame(region="Central_Pacific", site="fakarava", station=character(4), motu=numeric(4), genus=numeric(4), family=numeric(4), stringsAsFactors = FALSE)

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
colnames(rich_station_fakarava) <- c("region", "site", "station", "motu", "genus", "family")

  # total MOTUs and family richness in fakarava region (=site)
rich_site_fakarava <- data.frame(region="Central_Pacific",
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

colnames(rich_site_fakarava) <- c("region", "site", "motu", "genus", "family", "mean_motu", "sd_motu", "mean_genus", "sd_genus", "mean_family", "sd_family")

write.csv(rich_site_fakarava, "outputs/04_exploration_richness/richness_fakarava.csv", row.names = FALSE)

rich_site_fakarava_melt <- melt(rich_site_fakarava[,1:5])
rich_site_fakarava_melt[1,"mean"] <- rich_site_fakarava$mean_motu
rich_site_fakarava_melt[2,"mean"] <- rich_site_fakarava$mean_genus
rich_site_fakarava_melt[3,"mean"] <- rich_site_fakarava$mean_family

rich_site_fakarava_melt[1,"sd"] <- rich_site_fakarava$sd_motu
rich_site_fakarava_melt[2,"sd"] <- rich_site_fakarava$sd_genus
rich_site_fakarava_melt[3,"sd"] <- rich_site_fakarava$sd_family



## Eparses data

eparse <- df_all_filters %>%
  filter(region=="West_Indian")

  # total MOTUs and family richness in eparse region
rich_tot_eparse <- data.frame(motu=numeric(1), genus=numeric(1), family=numeric(1), region="West_Indian")
rich_tot_eparse$motu <- eparse %>% 
  summarise(n = n_distinct(sequence))
rich_tot_eparse$genus <- eparse %>% 
  summarise(n = n_distinct(new_genus_name))
rich_tot_eparse$family <- eparse %>% 
  summarise(n = n_distinct(new_family_name))


  # calculate unique motus and families at each station   
station <- c(unique(eparse$station))

rich_station_eparse <- data.frame(region="West_Indian", site=character(16), station=character(16), motu=numeric(16), genus=numeric(16), family=numeric(16), stringsAsFactors = FALSE)

for (i in 1:length(station)) {
  s <- unique(eparse[eparse$station == station[i],]$site)
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


  # calculate unique motus and families at each site
site <- c(unique(eparse$site))

rich_site_eparse <- data.frame(region="West_Indian", site=character(4), motu=numeric(4), genus=numeric(4), family=numeric(4), mean_motu=numeric(4), sd_motu=numeric(4), mean_genus=numeric(4), sd_genus=numeric(4), mean_family=numeric(4), sd_family=numeric(4), stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  s <- site[i]
  motu <- eparse[eparse$site == site[i],] %>%
    summarise(n = n_distinct(sequence))
  mm <- mean(rich_station_eparse[rich_station_eparse$site == site[i],]$motu)
  sdm <- sd(rich_station_eparse[rich_station_eparse$site == site[i],]$motu)
  gen <- eparse[eparse$site == site[i],] %>%
    summarise(n = n_distinct(new_genus_name))
  mg <- mean(rich_station_eparse[rich_station_eparse$site == site[i],]$genus)
  sdg <- sd(rich_station_eparse[rich_station_eparse$site == site[i],]$genus)
  fam <- eparse[eparse$site == site[i],] %>%
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

write.csv(rich_site_eparse, "outputs/04_exploration_richness/richness_eparse.csv", row.names = FALSE)


rich_site_eparse_melt <- melt(rich_site_eparse[,1:5])
rich_site_eparse_melt[1:4,"mean"] <- rich_site_eparse$mean_motu
rich_site_eparse_melt[5:8,"mean"] <- rich_site_eparse$mean_genus
rich_site_eparse_melt[9:12,"mean"] <- rich_site_eparse$mean_family

rich_site_eparse_melt[1:4,"sd"] <- rich_site_eparse$sd_motu
rich_site_eparse_melt[5:8,"sd"] <- rich_site_eparse$sd_genus
rich_site_eparse_melt[9:12,"sd"] <- rich_site_eparse$sd_family



# plot motu and family richness per region

rich_total <- rbind(rich_site_eparse_melt, rich_site_lengguru_melt, rich_site_caribbean_melt, rich_site_fakarava_melt)

gg <- ggplot(rich_total, aes(fill=variable)) +
  geom_col(aes(site, value), position = position_dodge(), show.legend = TRUE)+
  geom_point(data=rich_total[rich_total$variable=="motu",], aes(site, mean), position = position_nudge(-0.3, 0))+
  geom_errorbar(data=rich_total[rich_total$variable=="motu",], aes(site, ymin=mean-sd, ymax=mean+sd), position = position_nudge(-0.3, 0), width=0.32)+
  geom_point(data=rich_total[rich_total$variable=="genus",], aes(site, mean))+
  geom_errorbar(data=rich_total[rich_total$variable=="genus",], aes(site, ymin=mean-sd, ymax=mean+sd), width=0.32)+
  geom_point(data=rich_total[rich_total$variable=="family",], aes(site, mean), position = position_nudge(0.3, 0))+
  geom_errorbar(data=rich_total[rich_total$variable=="family",], aes(site, ymin=mean-sd, ymax=mean+sd), position = position_nudge(0.3, 0), width=0.32)+
  facet_grid(~region, scales = "free_x", space = "free_x")+
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5))+
  theme(strip.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 11))+
  theme(axis.title.y = element_text(size = 15))+
  theme(axis.title.x = element_text(size = 15))+
  scale_fill_manual(values=c("#d2981a", "#a53e1f", "#457277"))+
  labs(x="Sites",
       y = "Richness",
       fill=NULL)
gg

ggsave("outputs/04_exploration_richness/richness_region.png", width = 18, height = 8)






## load useful metadata 
metadata <- read.csv("metadata/Metadata_eDNA_global_V6.csv", stringsAsFactors = TRUE)
metadata <- metadata %>%
  filter(region%in%c("Central_IndoPacific", "Central_Pacific", "Caribbean", "West_Indian"))
metadata <- subset(metadata, !(station %in% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")))
metadata <- subset(metadata, sample_method!="niskin")
metadata <- subset(metadata, habitat=="marine")
metadata <- metadata[,c("site", "station", "latitude_start_clean", "longitude_start_clean", "dist_to_coast..m.", "dist_to_CT")]




## plot station richness ~ distance au CT
metadata <- metadata %>%
  distinct(station, .keep_all = TRUE)


rich_station <- rbind(rich_station_caribbean, rich_station_fakarava, rich_station_lengguru, rich_station_eparse)

rich_station <- left_join(rich_station, metadata[,c("station", "dist_to_CT")], by="station")

rich_motu <- ggplot(rich_station, aes(dist_to_CT, motu, col=region))+
  geom_point(size=2)+ 
  xlim(-11000,18000)+
  ylim(0,310)+
  scale_color_manual(values=c("#E5A729", "#8AAE8A", "#4F4D1D", "#C67052"))+ 
  theme(legend.position = "none")+
  theme(axis.title.y = element_text(size = 10, face = "bold"), plot.margin=unit(c(0.2,0.1,0,0.1), "cm"))+
  labs(x="",
       y="MOTU richness")
rich_motu
ggsave("outputs/04_exploration_richness/richness_motu_distCT.png")

rich_genus <- ggplot(rich_station, aes(dist_to_CT, genus, col=region))+
  geom_point(size=2)+ 
  xlim(-11000,18000)+
  scale_color_manual(values=c("#E5A729", "#8AAE8A", "#4F4D1D", "#C67052"))+ 
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100))+
  theme(legend.position = "none")+
  theme(axis.title.y = element_text(size = 10, face = "bold"), plot.margin=unit(c(0.2,0.1,0,0.1), "cm"))+
  labs(x="",
       y="Genus richness")
rich_genus
ggsave("outputs/04_exploration_richness/richness_genus_distCT.png")


rich_family <- ggplot(rich_station, aes(dist_to_CT, family, col=region))+
  geom_point(size=2)+ 
  xlim(-11000,18000)+
  ylim(0,60)+
  scale_color_manual(values=c("#E5A729", "#8AAE8A", "#4F4D1D", "#C67052"))+ 
  theme(legend.position = "none")+
  theme(axis.title.y = element_text(size = 10, face = "bold"), plot.margin=unit(c(0.2,0.1,0,0.1), "cm"))+
  labs(x="",
       y="Family richness")
rich_family
ggsave("outputs/04_exploration_richness/richness_family_distCT.png")

plot <- ggarrange(rich_motu, rich_genus, rich_family, ncol = 3, nrow=1)
x.grob <- textGrob("Distance to Coral Triangle (km, W-E)", 
                   gp=gpar(fontface="bold", col="black", fontsize=12), vjust = -0.5)

plot_all_rich_station <- grid.arrange(plot, bottom=x.grob)
save(plot_all_rich_station, file = "Rdata/plot_richness~dist_CT.rdata")


## plot site richness ~ distance au CT
metadata <- metadata %>%
  distinct(site, .keep_all = TRUE)


rich_site <- rbind(rich_site_caribbean, rich_site_fakarava, rich_site_lengguru, rich_site_eparse)

rich_site <- left_join(rich_site, metadata[,c("site", "dist_to_CT")], by="site")

rich_motu_site <- ggplot(rich_site, aes(dist_to_CT, motu, col=region))+
  geom_point(size=2)+ 
  xlim(-11000,18000)+
  ylim(0,730)+
  scale_color_manual(values=c("#E5A729", "#8AAE8A", "#4F4D1D", "#C67052"))+ 
  theme(legend.position = "none")+
  theme(axis.title.y = element_text(size = 10, face = "bold"), plot.margin=unit(c(0.2,0.1,0,0.1), "cm"))+
  labs(x="",
       y="MOTU richness")
rich_motu_site
ggsave("outputs/04_exploration_richness/richness_motu_site_distCT.png")

rich_genus_site <- ggplot(rich_site, aes(dist_to_CT, genus, col=region))+
  geom_point(size=2)+ 
  xlim(-11000,18000)+
  scale_color_manual(values=c("#E5A729", "#8AAE8A", "#4F4D1D", "#C67052"))+ 
  #scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100))+
  theme(legend.position = "none")+
  theme(axis.title.y = element_text(size = 10, face = "bold"), plot.margin=unit(c(0.2,0.1,0,0.1), "cm"))+
  labs(x="",
       y="Genus richness")
rich_genus_site
ggsave("outputs/04_exploration_richness/richness_genus_site_distCT.png")


rich_family_site <- ggplot(rich_site, aes(dist_to_CT, family, col=region))+
  geom_point(size=2)+ 
  xlim(-11000,18000)+
  ylim(0,90)+
  scale_color_manual(values=c("#E5A729", "#8AAE8A", "#4F4D1D", "#C67052"))+ 
  theme(legend.position = "none")+
  theme(axis.title.y = element_text(size = 10, face = "bold"), plot.margin=unit(c(0.2,0.1,0,0.1), "cm"))+
  labs(x="",
       y="Family richness")
rich_family_site
ggsave("outputs/04_exploration_richness/richness_family_site_distCT.png")

plot <- ggarrange(rich_motu_site, rich_genus_site, rich_family_site, ncol = 3, nrow=1)
x.grob <- textGrob("Distance to Coral Triangle (km, W-E)", 
                   gp=gpar(fontface="bold", col="black", fontsize=12), vjust = -0.5)

plot_all_rich_site <- grid.arrange(plot, bottom=x.grob)
save(plot_all_rich_site, file = "Rdata/plot_richness_site~dist_CT.rdata")






## plot station richness ~ latitude
rich_station <- rbind(rich_station_caribbean, rich_station_eparse, rich_station_fakarava, rich_station_lengguru)

rich_station <- left_join(rich_station, metadata[,c("station", "latitude_start_clean")], by="station")

ggplot(rich_station, aes(latitude_start_clean, motu))+
  geom_point(color="blue")+
  labs(x="latitude",
       y="MOTU richness")

ggsave("outputs/04_exploration_richness/richness_motu_latitude.png")

ggplot(rich_station, aes(latitude_start_clean, family))+
  geom_point(color="red")+
  labs(x="latitude",
       y="Family richness")

ggsave("outputs/04_exploration_richness/richness_family_latitude.png")


## plot station richness ~ distance_to_coast

rich_station <- rbind(rich_station_caribbean, rich_station_eparse, rich_station_fakarava, rich_station_lengguru)

rich_station <- left_join(rich_station, metadata[,c("station", "dist_to_coast..m.")], by="station")

ggplot(rich_station, aes(dist_to_coast..m., motu))+
  geom_point(color="blue")+
  labs(x="distance to coast (m)",
       y="MOTU richness")

ggsave("outputs/04_exploration_richness/richness_motu_distance_to_coast.png")

ggplot(rich_station, aes(dist_to_coast..m., family))+
  geom_point(color="red")+
  labs(x="Distance to coast (m)",
       y="Family richness")

ggsave("outputs/04_exploration_richness/richness_family_distance_to_coast.png")

