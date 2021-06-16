library(tidyverse)
library(reshape2)
library(lisa)
library(gridExtra)
library(grid)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(scales)


load("Rdata/02-clean-data.Rdata")

'%ni%' <- Negate("%in%")

df_all_filters <- df_all_filters %>%
  filter(province %in% c("Western_Indian_Ocean", "Southeast_Polynesia", "Tropical_Northwestern_Atlantic", "Western_Coral_Triangle", "Tropical_Southwestern_Pacific"))%>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3", "glorieuse_distance_300m")) %>%
  filter(sample_method !="niskin" & comment %ni% c("Distance decay 600m", "Distance decay 300m"))%>%
  filter(project != "Curacao") %>%
  filter(habitat=="marine")%>%
  filter(habitat_type %ni% c("BAIE"))%>%
  filter(site35!="")%>%
  filter(depth<40) %>%
  filter(family_name_corrected %ni% "Salmonidae")

df_all_filters <- df_all_filters %>%
  filter(!is.na(family_name_corrected))

perc_id <- df_all_filters %>% distinct(sequence, .keep_all=T)
perc_id <- perc_id[,"best_identity_database"]

perc_id2 <- as.data.frame(table(perc_id$best_identity_database))
perc_id2$perc <- (perc_id2$Freq*100)/1500

ggplot(perc_id2)+
  geom_col(aes(x=Var1, y=perc))+
  scale_x_discrete(breaks=c(0.85, 0.90, 0.95, 1))+
  xlab("% of similarity with reference sequence")+
  ylab("% of MOTUs")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank())
  

ggsave(filename = "outputs/00_Figures_for_paper/Extended_Data/ED_Figure10.png")
