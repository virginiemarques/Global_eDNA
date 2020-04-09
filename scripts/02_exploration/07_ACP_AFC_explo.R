library(tidyverse)
library(reshape2)
library(vegan)
library(betapart)

setwd("c:/Users/mathon/Desktop/linux/Global_eDNA/")
load("Rdata/02_clean_all.Rdata")

#Remove estuary stations and deep niskin station
df_all_filters <- subset(df_all_filters, !(station %in% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")))
df_all_filters <- subset(df_all_filters, sample_method!="niskin")


## table presence-absence of each unique MOTU per station

station <- unique(df_all_filters$station)
pa_station <- data.frame(motu=character(), stringsAsFactors = FALSE)
for (j in 1:length(station)) {
  df2 <- df_all_filters[df_all_filters$station==station[j],] %>%
    distinct(sequence, station)
  colnames(df2) <- c("motu", station[j])
  pa_station <- full_join(pa_station, df2, by="motu")
}
rownames(pa_station) <- pa_station[,1]
pa_station <- decostand(pa_station[,c(-1)], "pa",na.rm = TRUE)
pa_station[is.na(pa_station)] <- 0
pa_station <- as.data.frame(t(pa_station))


## tableau metadata par station (region, site, prof, distance cote, latitude, methode, nb litre)
## faire ACP / AFC / dbRDA

## plot richesse spe par station en foncton de latitude
