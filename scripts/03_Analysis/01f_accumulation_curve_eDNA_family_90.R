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

# Select only assignments of families > 90%
df_90 <- df_all_filters %>%
  filter(!is.na(family_name_corrected)) %>%
  filter(best_identity_database >= 0.90)

# rank_specify
rank_choice = 'family_name_corrected'


all_accumulation <- accumulation_curve_df(df_90, species_unit = rank_choice) %>%
  mutate(project_name = "All") %>%
  select(project_name, richness, sd, sites)


# Asymptote of all plots 
all_asymptote <- asymptote_mm(df_90, species_unit = rank_choice) %>%
  mutate(project_name = "All") %>%
  select(project_name, asymptote, slope)


# 
df_join_all <- all_accumulation %>%
  left_join(., all_asymptote, by = "project_name") %>%
  group_by(project_name) %>%
  mutate(position_asymptote_y = 1.05*asymptote, 
         position_asymptote_x = max(sites),
         position_slope_y = 0.50 * max(asymptote)) 

# Plots
colnames(df_join_all) <- c("Province", "richness", "sd", "sites", "asymptote", "slope", "position_asymptote_y", "position_asymptote_x", "position_slope_y")


family <- ggplot(df_join_all) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.8, fill="#457277") +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "solid", size = 1, col="#457277") +
  annotate(geom="text", x=250, y=141+10, label="Asymptote : 141",hjust=1,color="#457277", size=4) +
  annotate(geom="text", x=250, y=115+10, label="eDNA Families : 115",hjust=1,color="#457277", size=4) +
  annotate(geom="text", x=225, y=30, label="Slope = 1.8",hjust=1, alpha=0.7, size=4)+
  ylim(0,190)+
  labs(y="", x="Number of samples")+
  ylab("Number of families") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(size = 12, face="bold"),
        plot.margin=unit(c(0,0.1,0,0.1), "cm"))

