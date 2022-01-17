
# Accumulation curves by sites


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
list_read_step4 <- split(df_all_filters, df_all_filters$site35)


lapply(list_read_step4, function(x){
  length(unique(x$sample_name_all_pcr))
})


# rank_specify
rank_choice = 'sequence'

# accumlation all plots
liste_accumulation <- lapply(list_read_step4, accumulation_curve_df, species_unit = rank_choice)

# Asymptote of all plots 
liste_asymptote <- lapply(list_read_step4, asymptote_mm, species_unit = rank_choice)

# Unlist
df_accumulation <- bind_rows(liste_accumulation, .id = "project_name")
df_asymptote <- bind_rows(liste_asymptote, .id = "project_name")

# Dataset for plot
df_accumulation_all <- left_join(df_accumulation, df_asymptote, by = "project_name") %>%
  group_by(project_name) %>%
  mutate(position_asymptote_y = 0.20 * max(asymptote), 
         position_asymptote_x = max(sites))


# Plot with facet
colnames(df_accumulation) <- c("Sites", "richness", "sd", "samples")
df_accumulation$id <- c(1:nrow(df_accumulation))
df_accumulation <- left_join(df_accumulation, df_all_filters[,c("site35", "province")], by=c("Sites"="site35"))
df_accumulation <- df_accumulation %>% 
  distinct(id, .keep_all=T)

plot_acc_motus <- ggplot(df_accumulation, aes( fill=province)) + 
  geom_ribbon(aes(x = samples, ymin = richness-sd, ymax = richness+sd),  alpha = 0.7) +
  geom_line(aes(x = samples, y = richness)) +
  #geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 1) +
  facet_wrap(~province, scales = "free") +
  scale_fill_manual(values=c("#a6611a", "#E5A729","#b2182b", "#80cdc1", "#015462"))+ 
  ylab("Number of MOTUs") +
  xlab("Samples (filter)") +
  theme_bw() + 
  ggtitle("MOTUs") 
  #geom_text(aes(x = position_asymptote_x, y =position_asymptote_y, hjust = 1, label = paste(round(asymptote, 1), "MOTUs")), col = "black", size=3)
  #geom_text(aes(x = position_asymptote_x, y =position_slope_y, hjust = 1, label = paste("Lomolino slope=",round(slope, 1))), col = "black", size=3)

plot_acc_motus

ggsave("outputs/03_accumulation_curves/01b_accumulation_curve_all_projects_combination_no_facet_motus.png", width = 12, height = 8)

