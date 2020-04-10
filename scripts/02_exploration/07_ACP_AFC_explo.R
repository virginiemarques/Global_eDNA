library(tidyverse)
library(reshape2)
library(vegan)
library(betapart)
library(ade4)

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
pa_station$station <- rownames(pa_station)

## tableau metadata par station (region, site, prof, distance cote, latitude, methode, nb litre)
metadata <- read.csv("metadata/Metadata_eDNA_global_V4.csv", stringsAsFactors = TRUE)
metadata <- metadata %>%
  filter(region%in%c("West_Papua", "French_Polynesia", "Caribbean", "East_Pacific"))
metadata <- subset(metadata, !(station %in% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")))
metadata <- subset(metadata, !(sample_method %in% c("niskin", "control")))
metadata <- subset(metadata, habitat=="marine")
metadata <- metadata[,c("station", "region", "site", "latitude_start_clean", "volume", "sample_method", "depth", "dist_to_coast..m.")]
metadata <- metadata %>%
  distinct(station, .keep_all = TRUE)


pa_station <- left_join(pa_station, metadata, by="station")

## faire ACP / AFC / dbRDA

afc_pa_station <-dudi.coa(pa_station[,1:1603])
pourc=round((afc_pa_station$eig/sum(afc_pa_station$eig))*100,2)
pourc
cumsum(pourc)

s.label(afc_pa_station$co, clab=0.5, boxes=FALSE, sub="MOTUs AFC1") #representation des especes
s.label(afc_pa_station$li, clab=0.5, sub="Stations") #representation des transects
region<-as.factor(pa_station$region)
s.class(afc_pa_station$li,region, col=c(1:4))

