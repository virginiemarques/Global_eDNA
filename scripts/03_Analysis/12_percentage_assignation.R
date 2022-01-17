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


perc_id_class <- data.frame(percentage_identity=numeric(), number_of_MOTUs=numeric(), percentage_of_MOTUs=numeric())
perc_id_class[1, "percentage_identity"] <- "85-87%"
perc_id_class[1, "number_of_MOTUs"] <- perc_id %>%
  filter(best_identity_database >= 0.85 & best_identity_database < 0.87) %>%
  nrow()
perc_id_class[1, "percentage_of_MOTUs"] <- (perc_id_class[1, "number_of_MOTUs"]*100)/nrow(perc_id)
perc_id_class[2, "percentage_identity"] <- "87-89%"
perc_id_class[2, "number_of_MOTUs"] <- perc_id %>%
  filter(best_identity_database >= 0.87 & best_identity_database < 0.89) %>%
  nrow()
perc_id_class[2, "percentage_of_MOTUs"] <- (perc_id_class[2, "number_of_MOTUs"]*100)/nrow(perc_id)
perc_id_class[3, "percentage_identity"] <- "89-91%"
perc_id_class[3, "number_of_MOTUs"] <- perc_id %>%
  filter(best_identity_database >= 0.89 & best_identity_database < 0.91) %>%
  nrow()
perc_id_class[3, "percentage_of_MOTUs"] <- (perc_id_class[3, "number_of_MOTUs"]*100)/nrow(perc_id)
perc_id_class[4, "percentage_identity"] <- "91-93%"
perc_id_class[4, "number_of_MOTUs"] <- perc_id %>%
  filter(best_identity_database >= 0.91 & best_identity_database < 0.93) %>%
  nrow()
perc_id_class[4, "percentage_of_MOTUs"] <- (perc_id_class[4, "number_of_MOTUs"]*100)/nrow(perc_id)
perc_id_class[5, "percentage_identity"] <- "93-95%"
perc_id_class[5, "number_of_MOTUs"] <- perc_id %>%
  filter(best_identity_database >= 0.93 & best_identity_database < 0.95) %>%
  nrow()
perc_id_class[5, "percentage_of_MOTUs"] <- (perc_id_class[5, "number_of_MOTUs"]*100)/nrow(perc_id)
perc_id_class[6, "percentage_identity"] <- "95-97%"
perc_id_class[6, "number_of_MOTUs"] <- perc_id %>%
  filter(best_identity_database >= 0.95 & best_identity_database < 0.97) %>%
  nrow()
perc_id_class[6, "percentage_of_MOTUs"] <- (perc_id_class[6, "number_of_MOTUs"]*100)/nrow(perc_id)
perc_id_class[7, "percentage_identity"] <- "97-99%"
perc_id_class[7, "number_of_MOTUs"] <- perc_id %>%
  filter(best_identity_database >= 0.97 & best_identity_database < 0.99) %>%
  nrow()
perc_id_class[7, "percentage_of_MOTUs"] <- (perc_id_class[7, "number_of_MOTUs"]*100)/nrow(perc_id)
perc_id_class[8, "percentage_identity"] <- ">99%"
perc_id_class[8, "number_of_MOTUs"] <- perc_id %>%
  filter(best_identity_database >= 0.99) %>%
  nrow()
perc_id_class[8, "percentage_of_MOTUs"] <- (perc_id_class[8, "number_of_MOTUs"]*100)/nrow(perc_id)


write.csv(perc_id_class, "outputs/00_Figures_for_paper/Extended_Data/Table_S8.csv", row.names = F)

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
