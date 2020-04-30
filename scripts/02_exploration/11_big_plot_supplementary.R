library(tidyverse)
library(reshape2)
library(vegan)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)


# code for the figure like de Vargas 2015
load("Rdata/02_clean_all.Rdata")
df_all_filters <- subset(df_all_filters, !(station %in% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")))
df_all_filters <- subset(df_all_filters, sample_method!="niskin")
df_all_filters <- subset(df_all_filters, region!="East_Pacific")


# load data to plot together

load("Rdata/nb_motus_per_family_global.Rdata")
load("Rdata/family_proportion_global.Rdata")
load("Rdata/proportion_motus_similarity_threshold.Rdata")
family_coverage <- read.csv("outputs/01_read_data_stats/family_resolution_coefs.csv")
family <- unique(df_all_filters$new_family_name)
family_coverage <- subset(family_coverage, Family%in%family)


# mise en forme data
colnames(family_coverage) <- c("family", "coef_sequencing", "coef_resolution")

all <- left_join(count_families_global, family_coverage, by="family")

prop_similarity <- melt(prop_similarity)


# plot chaque plot
prop <- ggplot(families_prop_global2, aes(x=reorder(family, prop), y = prop, fill = Region)) + 
  geom_bar(stat="identity", show.legend = TRUE) + 
  theme_bw() +
  scale_fill_manual(values =c("#d7191c", "#2c7bb6", "#fdae61"))+#, "#abd9e9"
  labs(title="Proportion of MOTUs at global scale, \nand their distribution in regions", x="", y="")+ 
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 6, face="bold"), plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"))+
  coord_flip()

coverage <- ggplot(all, aes(x=reorder(family, n_motus), y = coef_sequencing)) + 
  geom_bar(stat="identity", show.legend = TRUE) + 
  theme_bw() +
  labs(title="Percentage of sequences \nknown in databases", x="", y="")+ 
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 6, face = "bold"), plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"))+
  theme(axis.text.y = element_blank())+
  coord_flip()

resolution <- ggplot(all, aes(x=reorder(family, n_motus), y = coef_resolution)) + 
  geom_bar(stat="identity", show.legend = TRUE) + 
  theme_bw() +
  labs(title="Percentage of resolutive \nsequences in databases", x="", y="")+ 
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 6, face = "bold"), plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"))+
  theme(axis.text.y = element_blank())+
  coord_flip()


ggarrange(prop, coverage, resolution, ncol=3, nrow=1, widths = c(2,1,1))


## add proportion similarité mais complexe
similarity <- ggplot(prop_similarity, aes(x=family, y = value, fill=variable)) + 
  geom_tile() + 
  theme_bw() +
  
  labs(title="Percentage of resolutive sequences in databases", x="", y="")+ 
  theme(legend.position = "none")+
  coord_flip()
