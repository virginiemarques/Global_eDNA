## calculate number of unique motus and families in each site and station of Lengguru and Caribbean

library(tidyverse)
library(reshape2)

setwd("c:/Users/mathon/Desktop/linux/Global_eDNA/")
load("Rdata/02_clean_all.Rdata")

df_all_filters <- subset(df_all_filters, !(station %in% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")))
df_all_filters <- subset(df_all_filters, sample_method!="niskin")


## Lengguru data

lengguru <- df_all_filters %>%
  filter(region=="West_Papua")

  # total MOTUs and family richness in Lengguru region
rich_tot_lengguru <- data.frame(motu=numeric(1), family=numeric(1))
rich_tot_lengguru$motu <- lengguru %>% 
  summarise(n = n_distinct(sequence))
rich_tot_lengguru$family <- lengguru %>% 
  summarise(n = n_distinct(new_family_name))


  # calculate unique motus and families at each site
site <- c(unique(lengguru$site))

rich_site_lengguru <- data.frame(site=character(11), station="total", motu=numeric(11), family=numeric(11), stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  s <- site[i]
  motu <- lengguru[lengguru$site == site[i],] %>%
    summarise(n = n_distinct(sequence))
  fam <- lengguru[lengguru$site == site[i],] %>%
    summarise(n = n_distinct(new_family_name))
  rich_site_lengguru[i,1] <- s
  rich_site_lengguru[i,3] <- motu
  rich_site_lengguru[i,4] <- fam
}

  
  # calculate unique motus and families at each station   
station <- c(unique(lengguru$station))

rich_station_lengguru <- data.frame(site=character(46), station=character(46), motu=numeric(46), family=numeric(46), stringsAsFactors = FALSE)

for (i in 1:length(station)) {
  s <- unique(lengguru[lengguru$station == station[i],]$site)
  st <- station[i]
  motu <- lengguru[lengguru$station == station[i],] %>%
    summarise(n = n_distinct(sequence))
  fam <- lengguru[lengguru$station == station[i],] %>%
    summarise(n = n_distinct(new_family_name))
  rich_station_lengguru[i,1] <- s
  rich_station_lengguru[i,2] <- st
  rich_station_lengguru[i,3] <- motu
  rich_station_lengguru[i,4] <- fam
}

  # combine tables for sites and stations
rich_lengguru <- rbind(rich_site_lengguru, rich_station_lengguru)


write.csv(rich_lengguru, "outputs/04_exploration_richness/richness_lengguru.csv", row.names = FALSE)


  # plot motu and family richness per site
rich_lengguru_melt <- melt(rich_lengguru)

plot_rich_lengguru <- ggplot(rich_lengguru_melt, aes(station, value, fill=variable)) +
    geom_col(position = position_dodge(), show.legend = TRUE)+
    facet_wrap(~site, scales = "free_x")+
    theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5))+
    theme(axis.title.y = element_text(size = 15))+
    labs(x=NULL,
        y = "richness",
        fill=NULL)

ggsave("outputs/04_exploration_richness/richness_lengguru.png")



## Caribbean data

caribbean <- df_all_filters %>%
  filter(region=="Caribbean")

  # total MOTUs and family richness in Caribbean region
rich_tot_caribbean <- data.frame(motu=numeric(1), family=numeric(1))
rich_tot_caribbean$motu <- caribbean %>% 
  summarise(n = n_distinct(sequence))
rich_tot_caribbean$family <- caribbean %>% 
  summarise(n = n_distinct(new_family_name))


  # calculate unique motus and families at each site
site <- c(unique(caribbean$site))

rich_site_caribbean <- data.frame(site=character(3), station="total", motu=numeric(3), family=numeric(3), stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  s <- site[i]
  motu <- caribbean[caribbean$site == site[i],] %>%
    summarise(n = n_distinct(sequence))
  fam <- caribbean[caribbean$site == site[i],] %>%
    summarise(n = n_distinct(new_family_name))
  rich_site_caribbean[i,1] <- s
  rich_site_caribbean[i,3] <- motu
  rich_site_caribbean[i,4] <- fam
}


  # calculate unique motus and families at each station   
station <- c(unique(caribbean$station))

rich_station_caribbean <- data.frame(site=character(31), station=character(31), motu=numeric(31), family=numeric(31), stringsAsFactors = FALSE)

for (i in 1:length(station)) {
  s <- unique(caribbean[caribbean$station == station[i],]$site)
  st <- station[i]
  motu <- caribbean[caribbean$station == station[i],] %>%
    summarise(n = n_distinct(sequence))
  fam <- caribbean[caribbean$station == station[i],] %>%
    summarise(n = n_distinct(new_family_name))
  rich_station_caribbean[i,1] <- s
  rich_station_caribbean[i,2] <- st
  rich_station_caribbean[i,3] <- motu
  rich_station_caribbean[i,4] <- fam
}

  # combine tables for sites and stations
rich_caribbean <- rbind(rich_site_caribbean, rich_station_caribbean)


write.csv(rich_caribbean, "outputs/04_exploration_richness/richness_caribbean.csv", row.names = FALSE)


  # plot motu and family richness per site
rich_caribbean_melt <- melt(rich_caribbean)

plot_rich_caribbean <- ggplot(rich_caribbean_melt, aes(station, value, fill=variable)) +
  geom_col(position = position_dodge(), show.legend = TRUE)+
  facet_wrap(~site, scales = "free_x")+
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.title.y = element_text(size = 15))+
  labs(x=NULL,
       y = "richness",
       fill=NULL)
ggsave("outputs/04_exploration_richness/richness_caribbean.png")


## Fakarava data

fakarava <- df_all_filters %>%
  filter(region=="French_Polynesia")

  # calculate unique motus and families at each station   
station <- c(unique(fakarava$station))

rich_station_fakarava <- data.frame(site="fakarava", station=character(4), motu.n=numeric(4), family.n=numeric(4), stringsAsFactors = FALSE)

for (i in 1:length(station)) {
  st <- station[i]
  motu <- fakarava[fakarava$station == station[i],] %>%
    summarise(n = n_distinct(sequence))
  fam <- fakarava[fakarava$station == station[i],] %>%
    summarise(n = n_distinct(new_family_name))
  rich_station_fakarava[i,2] <- st
  rich_station_fakarava[i,3] <- motu
  rich_station_fakarava[i,4] <- fam
}

  # total MOTUs and family richness in fakarava region (=site)
rich_fakarava <- rich_station_fakarava
rich_fakarava[5,1] <- "fakarava"
rich_fakarava[5,2] <- "total"
rich_fakarava[5,3] <- fakarava %>% 
  summarise(n = n_distinct(sequence))
rich_fakarava[5,4] <- fakarava %>% 
  summarise(n = n_distinct(new_family_name))


write.csv(rich_fakarava, "outputs/04_exploration_richness/richness_fakarava.csv", row.names = FALSE)


  # plot motu and family richness per site
rich_fakarava_melt <- melt(rich_fakarava)

plot_rich_fakarava <- ggplot(rich_fakarava_melt, aes(station, value, fill=variable)) +
  geom_col(position = position_dodge(), show.legend = TRUE)+
  facet_wrap(~site, scales = "free_x")+
  theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.title.y = element_text(size = 15))+
  labs(x=NULL,
       y = "richness",
       fill=NULL)

ggsave("outputs/04_exploration_richness/richness_fakarava.png")

## malpelo data

malpelo <- df_all_filters %>%
  filter(region=="East_Pacific")

  # calculate unique motus and families at each station   
station <- c(unique(malpelo$station))

rich_station_malpelo <- data.frame(site="malpelo", station=character(13), motu.n=numeric(13), family.n=numeric(13), stringsAsFactors = FALSE)

for (i in 1:length(station)) {
  st <- station[i]
  motu <- malpelo[malpelo$station == station[i],] %>%
    summarise(n = n_distinct(sequence))
  fam <- malpelo[malpelo$station == station[i],] %>%
    summarise(n = n_distinct(new_family_name))
  rich_station_malpelo[i,2] <- st
  rich_station_malpelo[i,3] <- motu
  rich_station_malpelo[i,4] <- fam
}

  # total MOTUs and family richness in malpelo region (=site)
rich_malpelo <- rich_station_malpelo
rich_malpelo[14,1] <- "malpelo"
rich_malpelo[14,2] <- "total"
rich_malpelo[14,3] <- malpelo %>% 
  summarise(n = n_distinct(sequence))
rich_malpelo[14,4] <- malpelo %>% 
  summarise(n = n_distinct(new_family_name))


write.csv(rich_malpelo, "outputs/04_exploration_richness/richness_malpelo.csv", row.names = FALSE)


# plot motu and family richness per site
rich_malpelo_melt <- melt(rich_malpelo)

plot_rich_malpelo <- ggplot(rich_malpelo_melt, aes(station, value, fill=variable)) +
  geom_col(position = position_dodge(), show.legend = TRUE)+
  facet_wrap(~site, scales = "free_x")+
  theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.title.y = element_text(size = 15))+
  labs(x=NULL,
       y = "richness",
       fill=NULL)

ggsave("outputs/04_exploration_richness/richness_malpelo.png")



## Mediterranean sea data
## Eparses data



## plot station richness ~ latitude
colnames(rich_station_fakarava) <- c("site", "station", "motu", "family")
colnames(rich_station_malpelo) <- c("site", "station", "motu", "family")

rich_station <- rbind(rich_station_caribbean, rich_station_fakarava, rich_station_lengguru, rich_station_malpelo)

metadata <- read.csv("metadata/Metadata_eDNA_global_V4.csv", stringsAsFactors = TRUE)
metadata <- metadata %>%
  filter(region%in%c("West_Papua", "French_Polynesia", "Caribbean", "East_Pacific"))
metadata <- subset(metadata, !(station %in% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")))
metadata <- subset(metadata, sample_method!="niskin")
metadata <- subset(metadata, habitat=="marine")
metadata <- metadata[,c("station", "latitude_start_clean")]
metadata <- metadata %>%
  distinct(station, .keep_all = TRUE)

rich_station <- left_join(rich_station, metadata, by="station")
rich_station <- rich_station[,c(-1)]

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
