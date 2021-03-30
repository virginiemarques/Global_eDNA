library(tidyverse)
library(sp)
library(rgeos)
library(rgdal)
library(maptools)
library(gissr)


## open metadata uptodate
metadata_sampling <- read.csv("metadata/Metadata_eDNA_global_V4.csv", header = T, sep = ";", stringsAsFactors = F, na.strings=c("","NA"))
metadata_sampling$longitude_start_clean <- gsub('\\?', '', metadata_sampling$longitude_start)
metadata_sampling$latitude_start_clean <- gsub('\\?', '', metadata_sampling$latitude_start)

## select samples with GPS data
metadata_dist  <- metadata_sampling %>%
  filter(!is.na(longitude_start_clean))

metadata_na  <- metadata_sampling %>%
  filter(is.na(longitude_start_clean))
## distance for whole world (EPSG:3857)

  # Define projections
crswgs84=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

crsmerc=CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")

  # import coastlines and change projection
coastline <- readOGR("C:/Users/mathon/Desktop/PhD/Projets/Megafauna/Carto_megafauna/GSHHS_f_L1.shp", verbose=TRUE)
coastline <- spTransform(coastline, crsmerc)
  
  # Formate GPS points and project
pts = metadata_dist[,43:44]
pts$latitude_start_clean <- as.numeric(pts$latitude_start_clean)
pts$longitude_start_clean <- as.numeric(pts$longitude_start_clean)
pts_sp <- SpatialPoints(pts,proj4string = crswgs84)
pts_sp <- spTransform(pts_sp, crsmerc)

  # calculate and formate distance
dist <- gDistance(pts_sp, coastline, byid = T)
dist_min <- apply(dist,2,min)
dist_min <- replace(dist_min, dist_min< 10, 10)
dist_min <- round(dist_min, 0) 
metadata_dist$dist_to_coast <- dist_min


  # assemble with metadata

metadata_na$dist_to_coast <- NA
metadata <- rbind(metadata_dist, metadata_na)

write.csv(metadata, "metadata/Metadata_eDNA_global_V4bis.csv", row.names = F)




