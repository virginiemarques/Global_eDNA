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


# data 
# WARNING: there is some empty lines in the .csv
metadata_sampling <- read.csv("metadata/Metadata_eDNA_global_V4.csv", sep=";", stringsAsFactors = F, na.strings=c("","NA"))

# Clean the spaces before coordinates
metadata_sampling$longitude_start_clean <- gsub('\\?', '', metadata_sampling$longitude_start)
metadata_sampling$latitude_start_clean <- gsub('\\?', '', metadata_sampling$latitude_start)


metadata_sampling <- subset(metadata_sampling, !(station %in% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3", "glorieuse_distance_300m")))
metadata_sampling <- subset(metadata_sampling, sample_method!="niskin")
metadata_sampling <- subset(metadata_sampling, !region %in% c("Mediterranean"))
metadata_sampling <- subset(metadata_sampling, !(comment %in% c("Distance decay 600m", "Distance decay 300m")))
metadata_sampling <- subset(metadata_sampling, habitat=="marine")

# ---------------------------------------------------------------------------------------------- #
#                            Some mapping 

# Mapping object
metadata_map <- metadata_sampling %>%
  filter(!is.na(longitude_start_clean))

# Load world data
world <- ne_countries(returnclass = 'sf')

# Convert point in dataset to sf object
metadata_map_sf = st_as_sf(metadata_map, coords = c("longitude_start_clean", "latitude_start_clean"), 
                             crs = 4326)

# Map 1
ggplot() + 
  geom_sf(aes(), data = world, fill = "grey80") + 
  geom_sf(fill = "blue", data= metadata_map_sf, shape=21, alpha=0.7) + 
  coord_sf(crs = "+proj=robin") + 
  theme_minimal() +
  theme(panel.grid.minor = element_line(linetype = "blank"),
        plot.background = element_rect(colour = NA), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        text = element_text(size=25),
        legend.key.size = unit(8, "mm"),
        panel.grid.major = element_line(colour = "gray70"),
        plot.title = element_text(lineheight=.8, face="bold")) 

# Save
ggsave("outputs/02_spatial_mapping/02_spatial_map_global_points_kept.png", width=8, height=6)




# ---------------------------------------------------------------------------------------------- #
#                            Interactive map


m <- leaflet(metadata_map_sf) %>%
  setView(lat=10, lng=0 , zoom=2) %>%
  addTiles() %>%
  addCircles(weight = 1,
             popup = 1, 
             radius = 900, 
             fillColor = "blue", 
             label = ~paste(station), 
             color = 'black', 
             opacity = 1) %>%
  # add ocean basemap
  # addProviderTiles(providers$Esri.OceanBasemap) %>%
  addProviderTiles(providers$GeoportailFrance.orthos) 

m

# I cant save the file where I want, for some reason. So I save it in root, then move it to wanted destination. 
saveWidget(m, file="02_interactive_map_coral_samples.html")
system("mv 02_interactive_map.html outputs/02_spatial_mapping/.")

# ---------------------------------------------------------------------------------------------- #
#                            bbox

# Longitude min et max
st_bbox(metadata_map_sf)
# -22.33944
# 46.89614 




