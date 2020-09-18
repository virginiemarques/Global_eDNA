# Exploration script to map all sites to visualize the spatial extent of the study
# Virginie Marques
# Last up: 02.03.2020


# Lib 
library(tidyverse)
library(sf)
library(leaflet)
library(rnaturalearth)
library(ggplot2)
library(htmlwidgets)

setwd("c:/Users/mathon/Desktop/linux/Global_eDNA/")
'%ni%' <- Negate("%in%")

# data 
# WARNING: there is some empty lines in the .csv
RLS_sampling <- read.csv("data/RLS/Coral_reef_fish_composition.csv", sep=";", stringsAsFactors = F, na.strings=c("","NA"))
RLS_sampling <- RLS_sampling %>%
  filter(realm%ni%c("Temperate Australasia", "Temperate Southern Africa", "Temperate Northern Atlantic", "Temperate Northern Pacific"))

# ---------------------------------------------------------------------------------------------- #
#                            Some mapping 


# Load world data
world <- ne_countries(returnclass = 'sf')

# Convert point in dataset to sf object
RLS_map_sf = st_as_sf(RLS_sampling, coords = c("SiteLong", "SiteLat"), 
                             crs = 4326)

# Map 1
ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80") + 
  geom_sf(data= RLS_map_sf, shape=20, alpha=0.7, size=2) + 
  coord_sf(crs = "+proj=robin") + 
  theme_minimal() +
  theme(panel.grid.minor = element_line(linetype = "blank"),
        plot.background = element_rect(colour = NA), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        text = element_text(size=25),
        legend.key.size = unit(2, "mm"),
        legend.text = element_text(size=3),
        legend.title = element_blank(),
        panel.grid.major = element_line(colour = "gray70"),
        plot.title = element_text(lineheight=.8, face="bold")) 

# Save
ggsave("outputs/02_spatial_mapping/02_RLS_500km.png", width=8, height=6)




# ---------------------------------------------------------------------------------------------- #
#                            Interactive map


m <- leaflet(RLS_map_sf) %>%
  setView(lat=10, lng=0 , zoom=2) %>%
  addTiles() %>%
  addCircles(weight = 1,
             popup = 1, 
             radius = 900, 
             fillColor = "blue", 
             label = ~paste(SiteCode), 
             color = 'black', 
             opacity = 1) %>%
  # add ocean basemap
  # addProviderTiles(providers$Esri.OceanBasemap) %>%
  addProviderTiles(providers$GeoportailFrance.orthos) 

m

# I cant save the file where I want, for some reason. So I save it in root, then move it to wanted destination. 
saveWidget(m, file="02_interactive_RLS.html")
system("mv 02_interactive_map.html outputs/02_spatial_mapping/.")

# ---------------------------------------------------------------------------------------------- #
#                            bbox

# Longitude min et max
st_bbox(metadata_map_sf)
# -22.33944
# 46.89614 




