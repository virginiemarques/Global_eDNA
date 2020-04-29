library(tidyverse)
library(sp)
library(rgeos)
library(rgdal)
library(maptools)
library(gissr)
library(raster)


## open metadata uptodate
metadata_sampling <- read.csv("metadata/Metadata_eDNA_global_V5.csv", header = T, sep = ",", stringsAsFactors = F, na.strings=c("","NA"))

## select samples with GPS data
metadata_dist  <- metadata_sampling %>%
  filter(!is.na(longitude_start_clean))

# coordinates of the center of the Coral Triangle
center_CT <- data.frame(longitude=133.679826, latitude=-1.307436)


pts = metadata_dist[,43:44]
pts$latitude_start_clean <- as.numeric(pts$latitude_start_clean)
pts$longitude_start_clean <- as.numeric(pts$longitude_start_clean)


# calculate distance
dist_to_CT <- pointDistance(pts, center_CT, lonlat=TRUE)
dist_to_CT <- as.data.frame(dist_to_CT)

# transform in KM with 12 decimals
dist_to_CT <- dist_to_CT/1000

dist_to_CT <- dist_to_CT %>% 
  mutate_if(is.numeric, round, digits = 2)


dist_to_CT$code_spygen <- metadata_dist$code_spygen
colnames(dist_to_CT) <- c("dist_to_CT", "code_spygen")


# assemble with metadata
metadata <- left_join(metadata_sampling, dist_to_CT, by="code_spygen")

write.csv(metadata, "metadata/Metadata_eDNA_global_V6.csv")
