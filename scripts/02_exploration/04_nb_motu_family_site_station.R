## calculate number of unique motus and families in each site and station of Lengguru and Caribbean

library(tidyverse)
library(reshape2)

setwd("c:/Users/mathon/Desktop/linux/Global_eDNA/")
load("Rdata/02_clean_all.Rdata")

df_all_filters <- subset(df_all_filters, !(station %in% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")))
df_all_filters <- subset(df_all_filters, sample_method!="niskin")
df_all_filters <- subset(df_all_filters, region!="East_Pacific")



## Lengguru data

lengguru <- df_all_filters %>%
  filter(region=="West_Papua")

  # total MOTUs and family richness in Lengguru region
rich_tot_lengguru <- data.frame(motu=numeric(1), family=numeric(1), region="West_Papua")
rich_tot_lengguru$motu <- lengguru %>% 
  summarise(n = n_distinct(sequence))
rich_tot_lengguru$family <- lengguru %>% 
  summarise(n = n_distinct(new_family_name))


  # calculate unique motus and families at each station   
station <- c(unique(lengguru$station))

rich_station_lengguru <- data.frame(region="West_Papua", site=character(46), station=character(46), motu=numeric(46), family=numeric(46), stringsAsFactors = FALSE)

for (i in 1:length(station)) {
  s <- unique(lengguru[lengguru$station == station[i],]$site)
  st <- station[i]
  motu <- lengguru[lengguru$station == station[i],] %>%
    summarise(n = n_distinct(sequence))
  fam <- lengguru[lengguru$station == station[i],] %>%
    summarise(n = n_distinct(new_family_name))
  rich_station_lengguru[i,2] <- s
  rich_station_lengguru[i,3] <- st
  rich_station_lengguru[i,4] <- motu
  rich_station_lengguru[i,5] <- fam
}


  # calculate unique motus and families at each site
site <- c(unique(lengguru$site))

rich_site_lengguru <- data.frame(region="West_Papua", site=character(11), motu=numeric(11), family=numeric(11), mean_motu=numeric(11), sd_motu=numeric(11), mean_family=numeric(11), sd_family=numeric(11), stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  s <- site[i]
  motu <- lengguru[lengguru$site == site[i],] %>%
    summarise(n = n_distinct(sequence))
  mm <- mean(rich_station_lengguru[rich_station_lengguru$site == site[i],]$motu)
  sdm <- sd(rich_station_lengguru[rich_station_lengguru$site == site[i],]$motu)
  fam <- lengguru[lengguru$site == site[i],] %>%
    summarise(n = n_distinct(new_family_name))
  mf <- mean(rich_station_lengguru[rich_station_lengguru$site == site[i],]$family)
  sdf <- sd(rich_station_lengguru[rich_station_lengguru$site == site[i],]$family)
  rich_site_lengguru[i,2] <- s
  rich_site_lengguru[i,3] <- motu
  rich_site_lengguru[i,4] <- fam
  rich_site_lengguru[i,5] <- mm
  rich_site_lengguru[i,6] <- sdm
  rich_site_lengguru[i,7] <- mf
  rich_site_lengguru[i,8] <- sdf
}

write.csv(rich_site_lengguru, "outputs/04_exploration_richness/richness_lengguru.csv", row.names = FALSE)


rich_site_lengguru_melt <- melt(rich_site_lengguru[,1:4])
rich_site_lengguru_melt[1:11,"mean"] <- rich_site_lengguru$mean_motu
rich_site_lengguru_melt[12:22,"mean"] <- rich_site_lengguru$mean_family

rich_site_lengguru_melt[1:11,"sd"] <- rich_site_lengguru$sd_motu
rich_site_lengguru_melt[12:22,"sd"] <- rich_site_lengguru$sd_family






## Caribbean data

caribbean <- df_all_filters %>%
  filter(region=="Caribbean")

  # total MOTUs and family richness in Caribbean region
rich_tot_caribbean <- data.frame(region="Caribbean", motu=numeric(1), family=numeric(1))
rich_tot_caribbean$motu <- caribbean %>% 
  summarise(n = n_distinct(sequence))
rich_tot_caribbean$family <- caribbean %>% 
  summarise(n = n_distinct(new_family_name))


  # calculate unique motus and families at each station   
station <- c(unique(caribbean$station))

rich_station_caribbean <- data.frame(region="Caribbean", site=character(31), station=character(31), motu=numeric(31), family=numeric(31), stringsAsFactors = FALSE)

for (i in 1:length(station)) {
  s <- unique(caribbean[caribbean$station == station[i],]$site)
  st <- station[i]
  motu <- caribbean[caribbean$station == station[i],] %>%
    summarise(n = n_distinct(sequence))
  fam <- caribbean[caribbean$station == station[i],] %>%
    summarise(n = n_distinct(new_family_name))
  rich_station_caribbean[i,2] <- s
  rich_station_caribbean[i,3] <- st
  rich_station_caribbean[i,4] <- motu
  rich_station_caribbean[i,5] <- fam
}

  # calculate unique motus and families at each site
site <- c(unique(caribbean$site))

rich_site_caribbean <- data.frame(region="Caribbean", site=character(3), motu=numeric(3), family=numeric(3), mean_motu=numeric(3), sd_motu=numeric(3), mean_family=numeric(3), sd_family=numeric(3), stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  s <- site[i]
  motu <- caribbean[caribbean$site == site[i],] %>%
    summarise(n = n_distinct(sequence))
  mm <- mean(rich_station_caribbean[rich_station_caribbean$site == site[i],]$motu)
  sdm <- sd(rich_station_caribbean[rich_station_caribbean$site == site[i],]$motu)
  fam <- caribbean[caribbean$site == site[i],] %>%
    summarise(n = n_distinct(new_family_name))
  mf <- mean(rich_station_caribbean[rich_station_caribbean$site == site[i],]$family)
  sdf <- sd(rich_station_caribbean[rich_station_caribbean$site == site[i],]$family)
  rich_site_caribbean[i,2] <- s
  rich_site_caribbean[i,3] <- motu
  rich_site_caribbean[i,4] <- fam
  rich_site_caribbean[i,5] <- mm
  rich_site_caribbean[i,6] <- sdm
  rich_site_caribbean[i,7] <- mf
  rich_site_caribbean[i,8] <- sdf
}

write.csv(rich_site_caribbean, "outputs/04_exploration_richness/richness_caribbean.csv", row.names = FALSE)


rich_site_caribbean_melt <- melt(rich_site_caribbean[,1:4])
rich_site_caribbean_melt[1:3,"mean"] <- rich_site_caribbean$mean_motu
rich_site_caribbean_melt[4:6,"mean"] <- rich_site_caribbean$mean_family

rich_site_caribbean_melt[1:3,"sd"] <- rich_site_caribbean$sd_motu
rich_site_caribbean_melt[4:6,"sd"] <- rich_site_caribbean$sd_family



## Fakarava data

fakarava <- df_all_filters %>%
  filter(region=="French_Polynesia")

  # calculate unique motus and families at each station   
station <- c(unique(fakarava$station))

rich_station_fakarava <- data.frame(region="French_Polynesia", site="fakarava", station=character(4), motu.n=numeric(4), family.n=numeric(4), stringsAsFactors = FALSE)

for (i in 1:length(station)) {
  st <- station[i]
  motu <- fakarava[fakarava$station == station[i],] %>%
    summarise(n = n_distinct(sequence))
  fam <- fakarava[fakarava$station == station[i],] %>%
    summarise(n = n_distinct(new_family_name))
  rich_station_fakarava[i,3] <- st
  rich_station_fakarava[i,4] <- motu
  rich_station_fakarava[i,5] <- fam
}

  # total MOTUs and family richness in fakarava region (=site)
rich_site_fakarava <- data.frame(region="French_Polynesia",
                                 site="fakarava",
                                 motu=fakarava %>% 
                                   summarise(n = n_distinct(sequence)),
                                 family=fakarava %>% 
                                   summarise(n = n_distinct(new_family_name)),
                                 mean_motu=mean(rich_station_fakarava$motu.n),
                                 sd_motu=sd(rich_station_fakarava$motu.n),
                                 mean_family=mean(rich_station_fakarava$family.n),
                                 sd_family=sd(rich_station_fakarava$family.n),
                                 stringsAsFactors = FALSE)

colnames(rich_site_fakarava) <- c("region", "site", "motu", "family", "mean_motu", "sd_motu", "mean_family", "sd_family")

write.csv(rich_site_fakarava, "outputs/04_exploration_richness/richness_fakarava.csv", row.names = FALSE)

rich_site_fakarava_melt <- melt(rich_site_fakarava[,1:4])
rich_site_fakarava_melt[1,"mean"] <- rich_site_fakarava$mean_motu
rich_site_fakarava_melt[2,"mean"] <- rich_site_fakarava$mean_family

rich_site_fakarava_melt[1,"sd"] <- rich_site_fakarava$sd_motu
rich_site_fakarava_melt[2,"sd"] <- rich_site_fakarava$sd_family



## Eparses data


# plot motu and family richness per region

rich_total <- rbind(rich_site_lengguru_melt, rich_site_caribbean_melt, rich_site_fakarava_melt)

gg <- ggplot(rich_total, aes(fill=variable)) +
  geom_col(aes(site, value), position = position_dodge(), show.legend = TRUE)+
  geom_point(data=rich_total[rich_total$variable=="motu",], aes(site, mean), position = position_nudge(-0.22, 0))+
  geom_errorbar(data=rich_total[rich_total$variable=="motu",], aes(site, ymin=mean-sd, ymax=mean+sd), position = position_nudge(-0.22, 0), width=0.45)+
  geom_point(data=rich_total[rich_total$variable=="family",], aes(site, mean), position = position_nudge(0.22, 0))+
  geom_errorbar(data=rich_total[rich_total$variable=="family",], aes(site, ymin=mean-sd, ymax=mean+sd), position = position_nudge(0.22, 0), width=0.45)+
  facet_grid(~region, scales = "free_x", space = "free_x")+
  theme(axis.text.x = element_text(size = 9, angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.title.y = element_text(size = 15))+
  labs(x=NULL,
       y = "richness",
       fill=NULL)


ggsave("outputs/04_exploration_richness/richness_region.png", width = 15, height = 8)


## plot station richness ~ latitude
colnames(rich_station_fakarava) <- c("region", "site", "station", "motu", "family")

rich_station <- rbind(rich_station_caribbean, rich_station_fakarava, rich_station_lengguru)

metadata <- read.csv("metadata/Metadata_eDNA_global_V4.csv", stringsAsFactors = TRUE)
metadata <- metadata %>%
  filter(region%in%c("West_Papua", "French_Polynesia", "Caribbean"))
metadata <- subset(metadata, !(station %in% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")))
metadata <- subset(metadata, sample_method!="niskin")
metadata <- subset(metadata, habitat=="marine")
metadata <- metadata[,c("station", "latitude_start_clean")]
metadata <- metadata %>%
  distinct(station, .keep_all = TRUE)

rich_station <- left_join(rich_station, metadata, by="station")

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



## plot station richness ~ longitude
colnames(rich_station_fakarava) <- c("region", "site", "station", "motu", "family")

rich_station <- rbind(rich_station_caribbean, rich_station_fakarava, rich_station_lengguru)

metadata <- read.csv("metadata/Metadata_eDNA_global_V4.csv", stringsAsFactors = TRUE)
metadata <- metadata %>%
  filter(region%in%c("West_Papua", "French_Polynesia", "Caribbean"))
metadata <- subset(metadata, !(station %in% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")))
metadata <- subset(metadata, sample_method!="niskin")
metadata <- subset(metadata, habitat=="marine")
metadata <- metadata[,c("station", "longitude_start_clean")]
metadata <- metadata %>%
  distinct(station, .keep_all = TRUE)

rich_station <- left_join(rich_station, metadata, by="station")

ggplot(rich_station, aes(longitude_start_clean, motu))+
  geom_point(color="blue")+
  xlim(-180, 180)+
  labs(x="longitude",
       y="MOTU richness")

ggsave("outputs/04_exploration_richness/richness_motu_longitude.png")

ggplot(rich_station, aes(longitude_start_clean, family))+
  geom_point(color="red")+
  xlim(-180, 180)+
  labs(x="longitude",
       y="Family richness")

ggsave("outputs/04_exploration_richness/richness_family_longitude.png")
