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
library(tidyverse)

'%ni%' <- Negate("%in%")

# remove temperate transects

RLS <- read.csv("data/RLS/RLS_species_NEW_all.csv", sep=";")
RLS<- RLS %>%
  filter(Realm%ni%c("Temperate Australasia", "Temperate Northern Atlantic", "Temperate Southern Africa", "Temperate Northern Pacific"))

#--------------------------------------------------------------------------------------------------------------------------------
# Select 1 year per survey RLS
RLS$Year <- substr(RLS$SurveyDate, 1, 4)

station <- unique(RLS$Station)
RLS_sub <- data.frame()
RLS_st <- vector("list")
for (i in 1:length(station)) {
  RLS_st[[i]] <- RLS %>%
    filter(Station==station[i]) %>%
    filter(Year == max(Year))
  RLS_sub <- rbind(RLS_sub, RLS_st[[i]])
}

RLS_sub <- RLS_sub[,-2167]  

#-------------------------------------------------------------------------------------------------------------------------------
# define RLS sites

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

coord <- data.frame(Station=RLS_sub$Station, lat=RLS_sub$SiteLat, long=RLS_sub$SiteLong)
coord_uni <- coord %>%
  distinct(Station, lat, long)%>%
  filter(!is.na(lat))

dgeo=matrix(0,nrow=1311,ncol=1311)

for (i in 1:1310)
{
  
  for (j in (i+1):1311)
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

RLS_sub_site <- full_join(RLS_sub, coord_uni[,-c(2,3)], by="Station")



#-------------------------------------------------------------------------------------------------------------------------------
# filter coral cover
RLS_CC <- read.csv("data/RLS/RLS PQ.csv")

survey <- unique(RLS_CC$SurveyID)
RLS_coral <- data.frame()
Survey_CC <- data.frame(SurveyID=numeric(), coral_cover=numeric())
for (i in 1:length(survey)) {
  RLS_coral <- RLS_CC %>%
    filter(SurveyID==survey[i])%>%
    filter(MajorCategory=="Coral")
  Survey_CC[i,1] <- survey[i]
  Survey_CC[i,2] <- sum(RLS_coral$Percent.Coverage)
}

RLS_sub_site_CC <- left_join(RLS_sub_site, Survey_CC, by="SurveyID")

RLS_fin <- RLS_sub_site_CC %>%
  filter(coral_cover >= 5)

write.csv(RLS_fin, "data/RLS/RLS_species_NEW.csv", row.names = F)

#------------------------------------------------------------------------------------------------------------------------------
# same filters on families file

RLS_fam <- read.csv("data/RLS/RLS_families_NEW_all.csv")
RLS_fam_sub <- RLS_fam %>%
  subset(SurveyID %in% RLS_fin$SurveyID)

RLS_fam_fin <- left_join(RLS_fam_sub, RLS_fin[,c("SurveyID","Province", "site35", "site30", "site25", "site20", "site10", "site5")], by="SurveyID")
write.csv(RLS_fam_fin, "data/RLS/RLS_families_NEW.csv", row.names = F)
