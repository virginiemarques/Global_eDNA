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

count_families_global <- arrange(count_families_global, n_motus)

prop_similarity <- melt(prop_similarity)
prop_similarity$family <- factor(prop_similarity$family, levels = as.character(factor(count_families_global$family)))
colnames(prop_similarity) <- c("family", "class", "percentage")


family_coverage$Family <- factor(family_coverage$Family, levels = as.character(factor(count_families_global$family)))


# plot chaque plot
prop <- ggplot(families_prop_global2, aes(x=reorder(family, prop), y = prop, fill = Region)) + 
  geom_bar(stat="identity", show.legend = TRUE) + 
  theme_bw() +
  scale_fill_manual(values =c("#8AAE8A", "#E5A729", "#4F4D1D"))+#, "#C67052"
  labs(title="Proportion of MOTUs at global scale, \nand their distribution in regions", x="", y="")+ 
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 6, face="bold"), plot.margin=unit(c(0.1,0.2,0.6,0), "cm"))+
  coord_flip()

coverage <- ggplot(family_coverage, aes(x=Family, y = coef_sequencing)) + 
  geom_bar(stat="identity", show.legend = TRUE) + 
  theme_bw() +
  labs(title="Percentage of sequences \nknown in databases", x="", y="")+ 
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 6, face = "bold"), plot.margin=unit(c(0.1,0.1,0.6,0), "cm"))+
  theme(axis.text.y = element_blank())+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  coord_flip()

resolution <- ggplot(family_coverage, aes(x=Family, y = coef_resolution)) + 
  geom_bar(stat="identity", show.legend = TRUE) + 
  theme_bw() +
  labs(title="Percentage of resolutive \nsequences in databases", x="", y="")+ 
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 6, face = "bold"), plot.margin=unit(c(0.1,0.1,0.6,0), "cm"))+
  theme(axis.text.y = element_blank())+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  coord_flip()


pal <- brewer.pal(9, "Greys")
similarity <- ggplot(prop_similarity, aes(x=class, y = family)) + 
  geom_tile(aes(fill=percentage))+
  scale_fill_gradientn(colours = pal)+
  theme_bw() +
  labs(title="Percentage of reads \nby similarity class", x="", y="")+ 
  theme(legend.position = c(0.8,1.015), legend.direction = "horizontal", legend.text = element_text(size=4), legend.key.height =unit(0.25,"cm"), legend.key.width=unit(0.2,"cm"), legend.title = element_blank())+
  theme(plot.title = element_text(size = 6, face = "bold"), plot.margin=unit(c(0.1,0.1,-0.1,0), "cm"))+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.8))


plot_all <- ggarrange(prop, coverage, resolution, similarity, ncol=4, nrow=1, widths = c(2,1,1,1))

# save wavec heigth=1300 et width=800

