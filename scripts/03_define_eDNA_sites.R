library(ade4)
library(ape)
library(cluster)
library(combinat)
library(maptree)
library(e1071)
library(maptree)
library(apTreeshape)
library(picante)
library(ecodist)
library(vegan)
library(dplyr)


# distances geographiques

geogdist <- function(lat1, lon1, lat2, lon2)
{
  # takes latitude and longitue for points 1 and 2
  # returns great circle distances in km
  # south and west are negative
  
  rad <- function(x) { (x/360) * 2 * pi }
  deg <- function(x) { (x/(2 * pi)) * 360 }
  
  gd <- (sin(rad(lat1)) * sin(rad(lat2))) + (cos(rad(lat1)) *
                                               cos(rad(lat2)) * cos(rad(abs(lon2-lon1))))
  
  gd <- deg(acos(gd))
  
  111.23 * gd
  
}


# importation data

data <- read.csv("metadata/Metadata_eDNA_global_V4bis.csv", sep=",")
coord <- data.frame(station=data$station, lat=data$latitude_start_clean, long=data$longitude_start_clean)
coord_uni <- coord %>%
  distinct(station, lat, long)%>%
  filter(!is.na(lat))

dgeo=matrix(0,nrow=509,ncol=509)

for (i in 1:508)
{
  
  for (j in (i+1):509)
    dgeo[i,j]=dgeo[j,i]=geogdist(coord_uni[i,2],coord_uni[i,3],coord_uni[j,2],coord_uni[j,3])
  
}

dgeo[is.na(dgeo)] <- 0
dgeo=as.dist(dgeo)
clust.transect=hclust(dgeo,method = "complete")
plot(clust.transect)

site35=cutree(clust.transect, h = 35)
site30=cutree(clust.transect, h = 30)
site25=cutree(clust.transect, h = 25)
site20=cutree(clust.transect, h = 20)
site10=cutree(clust.transect, h = 10)
site5=cutree(clust.transect, h = 5)

coord_uni$site35 <- paste("site",site35, sep="")
coord_uni$site30 <- paste("site",site30, sep="")
coord_uni$site25 <- paste("site",site25, sep="")
coord_uni$site20 <- paste("site",site20, sep="")
coord_uni$site10 <- paste("site",site10, sep="")
coord_uni$site5 <- paste("site",site5, sep="")



write.csv(coord_uni, "metadata/site_eDNA_distance.csv")

#metadata <- read.csv("metadata/Metadata_eDNA_global_V4bis.csv", sep=",")
metadata <- full_join(data, coord_uni[,-c(2,3)], by="station")
metadata <- metadata %>%
  distinct(station, code_spygen, .keep_all = T)


write.csv(metadata, "metadata/Metadata_eDNA_global_V5.csv", row.names = F)
