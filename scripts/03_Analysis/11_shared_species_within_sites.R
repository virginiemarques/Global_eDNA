library(tidyverse)
library(reshape2)
library(gridExtra)
library(grid)
library(cowplot)
library(ggplot2)
library(ggpubr)

load("Rdata/02-clean-data.Rdata")
'%ni%' <- Negate("%in%")

#-------------------------------------------------------------------------------------------------------------------------
# edna
#-------------------------------------------------------------------------------------------------------------------------

#Remove estuary stations and deep niskin station
df_all_filters <- df_all_filters %>%
  filter(province %in% c("Western_Indian_Ocean", "Southeast_Polynesia", "Tropical_Northwestern_Atlantic", "Western_Coral_Triangle", "Tropical_Southwestern_Pacific"))%>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3", "glorieuse_distance_300m")) %>%
  filter(sample_method !="niskin" & comment %ni% c("Distance decay 600m", "Distance decay 300m"))%>%
  filter(project != "Curacao") %>%
  filter(habitat=="marine")%>%
  filter(habitat_type %ni% c("BAIE"))%>%
  filter(depth<40) %>%
  filter(family_name_corrected %ni% "Salmonidae")


site <- unique(df_all_filters$site35)
df2 <- vector("list")
shared_site <- data.frame(site=character(), shared=numeric(), stringsAsFactors = F)

for (i in 1:length(site)) {
  df <- df_all_filters[df_all_filters$site35==site[i],]
  station <- unique(df$station)
  df_site <- data.frame(species=character())
  for (j in 1:length(station)) {
    df2[[j]] <- df[df$station==station[j],] %>%
      distinct(sequence)%>%
      as.data.frame()
      colnames(df2[[j]]) <- "species"
      df2[[j]]$station <- 1
    df_site <- full_join(df_site, df2[[j]], by="species")
  }
  shared_site[i,"site"] <- site[i]
  shared_site[i,"shared"] <- nrow(na.omit(df_site)) / nrow(df_site)
}

shared_site <- shared_site %>%
  filter(shared != 1)

mean(shared_site$shared)
sd(shared_site$shared)

#---------------------------------------------------------------------------------------------------------------------------
# RLS
#---------------------------------------------------------------------------------------------------------------------------

RLS_species <- read.csv("data/RLS/RLS_species_NEW.csv", sep = ";", stringsAsFactors = F, check.names = F)
RLS_species <- RLS_species[, c(1,2,11,7,18:2173)]
RLS_species <- reshape2::melt(RLS_species, id=c("SurveyID", "Station", "site35", "Province"))
RLS_species <- RLS_species%>%
  filter(value!=0)
RLS_species <- RLS_species[,-6]
colnames(RLS_species) <- c("SurveyID", "Station","site35", "Province", "Species")

site <- unique(RLS_species$site35)
df2 <- vector("list")
shared_site <- data.frame(site=character(), shared=numeric(), stringsAsFactors = F)

for (i in 1:length(site)) {
  df <- RLS_species[RLS_species$site35==site[i],]
  station <- unique(df$Station)
  df_site <- data.frame(species=character())
  for (j in 1:length(station)) {
    df2[[j]] <- df[df$Station==station[j],] %>%
      distinct(Species)%>%
      as.data.frame()
    colnames(df2[[j]]) <- "species"
    df2[[j]]$Station <- 1
    df_site <- full_join(df_site, df2[[j]], by="species")
  }
  shared_site[i,"site"] <- site[i]
  shared_site[i,"shared"] <- nrow(na.omit(df_site)) / nrow(df_site)
}

shared_site <- shared_site %>%
  filter(shared != 1)

mean(shared_site$shared)
sd(shared_site$shared)
