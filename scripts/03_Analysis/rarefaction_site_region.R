library(tidyverse)
'%ni%' <- Negate("%in%")

load("Rdata/02-clean-data.Rdata")

df_all_filters <- df_all_filters %>%
  filter(province %in% c("Western_Indian_Ocean", "Southeast_Polynesia", "Tropical_Northwestern_Atlantic", "Western_Coral_Triangle", "Tropical_Southwestern_Pacific"))%>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3", "glorieuse_distance_300m")) %>%
  filter(sample_method !="niskin" & comment %ni% c("Distance decay 600m", "Distance decay 300m"))%>%
  filter(project != "Curacao") %>%
  filter(habitat=="marine")%>%
  filter(habitat_type %ni% c("BAIE"))%>%
  filter(site35!="") %>%
  filter(depth<40) %>%
  filter(family_name_corrected %ni% "Salmonidae")

# ------------------------------------------------------------------------------------------------
# Rarefaction of 4 stations per site
# ------------------------------------------------------------------------------------------------

# je ne fais que 15 simulations car le plus grand nombre de station par site est de 10 
# donc je pense que ça suffit pour avoir toutes les combinaisons

rarefied_sites <- vector("list", 15)

for (j in 1:15) {
  
  sites <- unique(df_all_filters$site35)
  
  stations <- vector("list", length(sites))
  
  for (i in 1:length(sites)) {
    stations[[i]] <- df_all_filters %>%
      filter(site35==sites[i])%>%
      distinct(station)
    stations[[i]] <- as.character(stations[[i]]$station)
    if (length(stations[[i]])>4){
      stations[[i]] <- sample(stations[[i]], 4, replace = F)
    }
  }
  
  selected_stations <- unlist(stations)
  
  rarefied_sites[[j]] <- df_all_filters %>%
    filter(station %in% selected_stations)
}

# On verifie qu'on a toujours le meme nombre de stations

lapply(rarefied_sites, function(x){
  length(unique(x$station))
})

# on sauve
save(rarefied_sites, file="Rdata/rarefied_sites.rdata")




# ---------------------------------------------------------------------------------------------
# Rarefaction of 4 stations per region
# ---------------------------------------------------------------------------------------------

rarefied_regions <- vector("list", 50)

for (j in 1:50) {
  
  regions <- unique(df_all_filters$province)
  
  stations <- vector("list", length(regions))
  
  for (i in 1:length(regions)) {
    stations[[i]] <- df_all_filters %>%
      filter(province==regions[i])%>%
      distinct(station)
    stations[[i]] <- as.character(stations[[i]]$station)
    stations[[i]] <- sample(stations[[i]], 4, replace = F)
  }
  
  selected_stations <- unlist(stations)
  
  rarefied_regions[[j]] <- df_all_filters %>%
    filter(station %in% selected_stations)
}

# On verifie qu'on a toujours le meme nombre de stations

lapply(rarefied_regions, function(x){
  length(unique(x$station))
})

# on sauve
save(rarefied_regions, file="Rdata/rarefied_regions.rdata")


