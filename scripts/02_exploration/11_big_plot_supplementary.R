library(tidyverse)
library(dplyr)
library(reshape2)
library(vegan)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggpubr)
library(RColorBrewer)


# code for the figure like de Vargas 2015
load("Rdata/02_clean_all.Rdata")
'%ni%' <- Negate("%in%")
df_all_filters <- df_all_filters %>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")) %>%
  filter(sample_method !="niskin" & region!="East_Pacific" & comment %ni% c("Distance decay 600m", "Distance decay 300m") & station!="glorieuse_distance_300m")%>%
  filter(project != "SEAMOUNTS") %>% 
  filter(habitat_type %ni% c("BAIE", "Sommet"))



# load data to plot together

load("Rdata/nb_motus_per_family_global.Rdata")
load("Rdata/family_proportion_global.Rdata")
load("Rdata/proportion_motus_similarity_threshold.Rdata")
family_coverage <- read.csv("outputs/01_read_data_stats/family_resolution_coefs.csv")
family <- unique(df_all_filters$new_family_name)
family_coverage <- subset(family_coverage, Family%in%family)


# select families with > 5 motus for main and < 5 motus for ED

count_families_global <- arrange(count_families_global, n_motus)
count_families_global_main <- count_families_global%>%
  subset(n_motus >= 5)
count_families_global_ED <- count_families_global %>%
  subset(n_motus < 5)
families_main <- as.character(unique(count_families_global_main$family))
families_ED <- as.character(unique(count_families_global_ED$family))

# mise en forme data pour ordonner correctement
families_prop_global_main <- families_prop_global %>%
  subset(family%in%families_main)
families_prop_global_main <- reshape2::melt(families_prop_global_main)
colnames(families_prop_global_main) <- c("family", "Region", "prop")

families_prop_global_ED <- families_prop_global %>%
  subset(family%in%family_coverage_ED$Family)
families_prop_global_ED <- reshape2::melt(families_prop_global_ED)
colnames(families_prop_global_ED) <- c("family", "Region", "prop")

prop_similarity <- left_join(prop_similarity, families_prop_global, by="family")
prop_similarity <- melt(prop_similarity, id=c("family", "85-90%", "90-95%", "95-<100%", "100%"))
colnames(prop_similarity) <- c("family", "85-90%", "90-95%", "95-<100%", "100%", "region", "prop")
prop_similarity_main <- prop_similarity %>%
  subset(family%in%families_main)
prop_similarity_main <- reshape2::melt(prop_similarity_main, id=c("family", "region", "prop"))
colnames(prop_similarity_main) <- c("family", "region", "prop", "class", "percentage")

prop_similarity_ED <- prop_similarity %>%
  subset(family%in%family_coverage_ED$Family)
prop_similarity_ED <- reshape2::melt(prop_similarity_ED, id=c("family", "region", "prop"))
colnames(prop_similarity_ED) <- c("family", "region", "prop", "class", "percentage")

family_coverage <- left_join(family_coverage, families_prop_global, by=c("Family"="family"))
family_coverage <- melt(family_coverage, id=c("Family", "coef_sequencing", "coef_resolution"))
colnames(family_coverage) <- c("Family", "coef_sequencing", "coef_resolution", "region", "prop")
family_coverage_main <- family_coverage %>%
  subset(Family%in%families_main)
family_coverage_ED <- family_coverage %>%
  subset(Family%in%families_ED)




# plot chaque plot main
prop <- ggplot(families_prop_global_main, aes(x=reorder(family, prop), y = prop, fill = Region)) + 
  geom_bar(stat="identity", show.legend = TRUE) + 
  theme_bw() +
  scale_fill_manual(values =c("#8AAE8A", "#E5A729", "#4F4D1D", "#C67052", "#863b34"))+ 
  labs(title="Proportion of MOTUs at global scale, \nand their distribution in regions", x="", y="")+ 
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 6, face="bold"), plot.margin=unit(c(0.1,0.2,0.6,0), "cm"))+
  coord_flip()

coverage <- ggplot(family_coverage_main, aes(x=reorder(Family, prop), y = coef_sequencing)) + 
  geom_bar(stat="identity", show.legend = TRUE) + 
  theme_bw() +
  labs(title="Percentage of sequences \nknown in databases", x="", y="")+ 
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 6, face = "bold"), plot.margin=unit(c(0.1,0.1,0.6,0), "cm"))+
  theme(axis.text.y = element_blank())+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  coord_flip()

resolution <- ggplot(family_coverage_main, aes(x=reorder(Family, prop), y = coef_resolution)) + 
  geom_bar(stat="identity", show.legend = TRUE) + 
  theme_bw() +
  labs(title="Percentage of resolutive \nsequences in databases", x="", y="")+ 
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 6, face = "bold"), plot.margin=unit(c(0.1,0.1,0.6,0), "cm"))+
  theme(axis.text.y = element_blank())+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  coord_flip()


pal <- brewer.pal(9, "Greys")
similarity <- ggplot(prop_similarity_main, aes(x=class, y = reorder(family, prop))) + 
  geom_tile(aes(fill=percentage))+
  scale_fill_gradientn(colours = pal)+
  theme_bw() +
  labs(title="Percentage of reads \nby similarity class", x="", y="")+ 
  theme(legend.position = c(0.8,1.025), legend.direction = "horizontal", legend.text = element_text(size=4), legend.key.height =unit(0.25,"cm"), legend.key.width=unit(0.2,"cm"), legend.title = element_blank())+
  theme(plot.title = element_text(size = 6, face = "bold"), plot.margin=unit(c(0.1,0.1,-0.1,0), "cm"))+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.8))


plot_all_main <- ggarrange(prop, coverage, resolution, similarity, ncol=4, nrow=1, widths = c(2,1,1,1), labels = c("A", "B", "C", "D"))
plot_all_main
ggsave("outputs/Figures papier/Figure2.png", width = 7.8, height = 8)


## plot all for extended Data
prop <- ggplot(families_prop_global_ED, aes(x=reorder(family, prop), y = prop, fill = Region)) + 
  geom_bar(stat="identity", show.legend = TRUE) + 
  theme_bw() +
  scale_fill_manual(values =c("#8AAE8A", "#E5A729", "#4F4D1D", "#C67052", "#863b34"))+ 
  labs(title="Proportion of MOTUs at global scale, \nand their distribution in regions", x="", y="")+ 
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 6, face="bold"), plot.margin=unit(c(0.1,0.2,0.6,0), "cm"))+
  coord_flip()

coverage <- ggplot(family_coverage_ED, aes(x=reorder(Family, prop), y = coef_sequencing)) + 
  geom_bar(stat="identity", show.legend = TRUE) + 
  theme_bw() +
  labs(title="Percentage of sequences \nknown in databases", x="", y="")+ 
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 6, face = "bold"), plot.margin=unit(c(0.1,0.1,0.6,0), "cm"))+
  theme(axis.text.y = element_blank())+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  coord_flip()

resolution <- ggplot(family_coverage_ED, aes(x=reorder(Family, prop), y = coef_resolution)) + 
  geom_bar(stat="identity", show.legend = TRUE) + 
  theme_bw() +
  labs(title="Percentage of resolutive \nsequences in databases", x="", y="")+ 
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 6, face = "bold"), plot.margin=unit(c(0.1,0.1,0.6,0), "cm"))+
  theme(axis.text.y = element_blank())+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  coord_flip()


pal <- brewer.pal(9, "Greys")
similarity <- ggplot(prop_similarity_ED, aes(x=class, y = reorder(family, prop))) + 
  geom_tile(aes(fill=percentage))+
  scale_fill_gradientn(colours = pal)+
  theme_bw() +
  labs(title="Percentage of reads \nby similarity class", x="", y="")+ 
  theme(legend.position = c(0.8,1.025), legend.direction = "horizontal", legend.text = element_text(size=4), legend.key.height =unit(0.25,"cm"), legend.key.width=unit(0.2,"cm"), legend.title = element_blank())+
  theme(plot.title = element_text(size = 6, face = "bold"), plot.margin=unit(c(0.1,0.1,-0.1,0), "cm"))+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.8))


plot_all_ED <- ggarrange(prop, coverage, resolution, similarity, ncol=4, nrow=1, widths = c(2,1,1,1), labels = c("A", "B", "C", "D"))
plot_all_ED
ggsave("outputs/Figures papier/ED_Figure3.png", width = 7.8, height = 10)

