library(tidyverse)
library(sp)
library(rgeos)
library(rgdal)
library(maptools)
library(gissr)
library(raster)


## open metadata uptodate
metadata_sampling <- read.csv("metadata/Metadata_eDNA_global_V5.csv", header = T, sep = ";", stringsAsFactors = F, na.strings=c("","NA"))

## select samples with GPS data
metadata_dist  <- metadata_sampling %>%
  filter(!is.na(longitude_start_clean))
metadata_na  <- metadata_sampling %>%
  filter(is.na(longitude_start_clean))


# coordinates of the center of the Coral Triangle
center_CT <- data.frame(longitude=133.679826, latitude=-1.307436)


# calculate distance
metadata_dist$dist_to_CT <- pointDistance(metadata_dist[,49:50], center_CT, lonlat=TRUE)


# transform in KM with 2 decimals
metadata_dist$dist_to_CT <- metadata_dist$dist_to_CT/1000

metadata_dist$dist_to_CT <- round(metadata_dist$dist_to_CT, 2) 


## put Coral triangle in the middle
metadata_dist <- metadata_dist %>%
  mutate(dist_to_CT = case_when(
    province == "Southeast_Polynesia" ~ paste0(dist_to_CT),
    province == "Tropical_Northwestern_Atlantic" ~ paste0(dist_to_CT),
    province == "Tropical_Southwestern_Pacific" ~ paste0(dist_to_CT),
    province == "Tropical_East_Pacific" ~ paste0(dist_to_CT),
    province == "Western_Indian_Ocean" ~ paste0("-", dist_to_CT),
    province == "Western_Coral_Triangle" ~ paste0("-", dist_to_CT),
    province == "Lusitanian" ~ paste0("-", dist_to_CT),
    province == "Mediterranean_Sea" ~ paste0("-", dist_to_CT)))

metadata_na$dist_to_CT <- NA
metadata <- rbind(metadata_dist, metadata_na)

write.csv(metadata, "metadata/Metadata_eDNA_global_V6.csv", row.names = F)
