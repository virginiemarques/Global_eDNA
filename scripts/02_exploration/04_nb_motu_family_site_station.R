## calculate number of unique motus and families in each site and station of Lengguru and Caribbean

library(tidyverse)
library(conflicted)
library(reshape2)

setwd("c:/Users/mathon/Desktop/linux/Global_eDNA/")
load("Rdata/02_clean_all.Rdata")

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

  
  # caluculate unique motus and families at each station   
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



## Caribbean data

caribbean <- df_all_filters %>%
  filter(region=="Caribbean")

  # total MOTUs and family richness in Lengguru region
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


  # caluculate unique motus and families at each station   
station <- c(unique(caribbean$station))

rich_station_caribbean <- data.frame(site=character(40), station=character(40), motu=numeric(40), family=numeric(40), stringsAsFactors = FALSE)

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

## Mediterranean sea data
## Eparses data