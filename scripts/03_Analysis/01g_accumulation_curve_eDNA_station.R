
# Accumulation curves by stations


# Lib 
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(conflicted)

# data
load("Rdata/02-clean-data.Rdata")

# 
'%ni%' <- Negate("%in%")
conflict_prefer("select", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("filter", "dplyr")

# Functions
source('scripts/03_Analysis/00_functions.R')



# On the df as well
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

#### ----------------------- Site scale ---------------------------------------------

# Split by region
list_read_step4 <- split(df_all_filters, df_all_filters$station)


lapply(list_read_step4, function(x){
  n_distinct(x$sample_name_all_pcr)
})

cond <- lapply(list_read_step4, function(x) n_distinct(x$sample_name_all_pcr) > 1)

list_read_step4 <- list_read_step4[unlist(cond)]


# rank_specify
rank_choice = 'sequence'

# accumulation all plots
liste_accumulation <- lapply(list_read_step4, accumulation_curve_df, species_unit = rank_choice)

for (i in 1:length(liste_accumulation)) {
  liste_accumulation[[i]]$station_id  <- paste("station_", i, sep="")
}

# Unlist
df_accumulation <- bind_rows(liste_accumulation, .id = "Station")


# Plot with facet
colnames(df_accumulation) <- c("Station", "richness", "sd", "samples", "station_id")
df_accumulation$id <- c(1:nrow(df_accumulation))
df_accumulation <- left_join(df_accumulation, df_all_filters[,c("site35", "station")], by=c("Station"="station"))
df_accumulation <- df_accumulation %>% 
  distinct(id, .keep_all=T)


station <- unique(df_accumulation$station_id)

df_plot1 <- df_accumulation %>%
  filter(station_id %in% c(station[1:36]))

plot1 <- ggplot(df_plot1) + 
  geom_ribbon(aes(x = samples, ymin = richness-sd, ymax = richness+sd),  alpha = 0.7, show.legend = F) +
  geom_line(aes(x = samples, y = richness), show.legend = F) +
  facet_wrap(~station_id, scales = "free") +
  scale_x_continuous(breaks=c(1,2,3,4,5,6))+
  ylab("Number of MOTUs") +
  xlab("Samples (filter)") +
  theme_bw() + 
  theme(strip.text.x = element_text(size=7))
  


df_plot2 <- df_accumulation %>%
  filter(station_id %in% c(station[37:72]))

plot2 <- ggplot(df_plot2) + 
  geom_ribbon(aes(x = samples, ymin = richness-sd, ymax = richness+sd),  alpha = 0.7, show.legend = F) +
  geom_line(aes(x = samples, y = richness), show.legend = F) +
  facet_wrap(~station_id, scales = "free") +
  scale_x_continuous(breaks=c(1,2,3,4,5,6))+
  ylab("Number of MOTUs") +
  xlab("Samples (filter)") +
  theme_bw()

df_plot3 <- df_accumulation %>%
  filter(station_id %in% c(station[73:93]))

plot3 <- ggplot(df_plot3) + 
  geom_ribbon(aes(x = samples, ymin = richness-sd, ymax = richness+sd),  alpha = 0.7, show.legend = F) +
  geom_line(aes(x = samples, y = richness), show.legend = F) +
  facet_wrap(~station_id, scales = "free") +
  scale_x_continuous(breaks=c(1,2,3,4,5,6))+
  ylab("Number of MOTUs") +
  xlab("Samples (filter)") +
  theme_bw()
