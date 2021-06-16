library(png)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggpubr)
library(scales)
library(dplyr)
library(conflicted)
library(gambin)
library(sads)
library(vegan)
library(plyr)

# 
conflict_prefer("filter", "dplyr")
conflict_prefer("summarise", "dplyr")
'%ni%' <- Negate("%in%")


## For panel a, b, c and d : Run scripts "scripts/03_Analysis/01a_accumulation_curve.R"  
#                                         "scripts/03_Analysis/01b_accumulation_curve_RLS.R"
## Or load Rdata :
load("Rdata/accumulation_asymptote_motus_all.rdata")
load("Rdata/accumulation_species_RLS.rdata")
load("Rdata/accumulation_asymptote_families_all.rdata")
load("Rdata/accumulation_families_RLS.rdata")


# plot panel a : accumulation curve MOTUs eDNA

a <- ggplot(df_motus) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.8, fill="#d2981a") +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "solid", size = 1, col="#d2981a") +
  annotate(geom="text", x=250, y=2650+150, label="Asymptote : 2650", hjust=1, color="#d2981a", size=3.2) +
  annotate(geom="text", x=250, y=2023+150, label="eDNA MOTUs : 2023", hjust=1, color="#d2981a", size=3.2) +
  annotate(geom="text", x=250, y=600, label="Slope = 2.27", hjust=1, alpha=0.7, size=3.2) +
  ylim(0,3000)+
  ylab("") +
  xlab("Number of samples")+
  ylab("Nb of species / MOTUs") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(size = 12, face="bold"),
        plot.margin=unit(c(0,0.1,0,0.1), "cm")) + 
  ggtitle("a")

# plot panel b : accumulation curve species visual survey

b <- ggplot(all_accumulation_species_RLS_df) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.5) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 1, alpha=0.7) +
  annotate(geom="text", x=2190, y=2268+150, label="Asymptote : 2268",hjust=1, alpha=0.7, size=3.2) +
  annotate(geom="text", x=2190, y=1818+150, label="Visual census Species : 1818",hjust=1, alpha=0.7, size=3.2) +
  annotate(geom="text", x=2190, y=600, label="Slope = 1.76",hjust=1, alpha=0.7, size=3.2) +
  ylim(0,3000)+
  xlab("Number of transects") +
  ylab("")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(size = 12, face="bold"),
        plot.margin=unit(c(0,0.15,0,0), "cm")) + 
  ggtitle("b")

# plot panel c : accumulation curve families eDNA

c <- ggplot(df_fam) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.8, fill="#457277") +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "solid", size = 1, col="#457277") +
  annotate(geom="text", x=250, y=147+10, label="Asymptote : 147",hjust=1,color="#457277", size=3.2) +
  annotate(geom="text", x=250, y=126+10, label="eDNA Families : 126",hjust=1,color="#457277", size=3.2) +
  annotate(geom="text", x=250, y=30, label="Slope = 1.9",hjust=1, alpha=0.7, size=3.2)+
  ylim(0,190)+
  labs(y="", x="Number of samples")+
  ylab("Number of families") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(size = 12, face="bold"),
        plot.margin=unit(c(0,0.1,0,0.1), "cm")) + 
  ggtitle("c")

# plot panel d : accumulation curve families visual census

d <- ggplot(all_accumulation_RLS_df) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.5) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 1, alpha=0.7) +
  annotate(geom="text", x=2190, y=118+10, label="Asymptote : 118",hjust=1, alpha=0.7, size=3.2) +
  annotate(geom="text", x=2190, y=96+10, label="Visual census Families : 96",hjust=1, alpha=0.7, size=3.2)+
  annotate(geom="text", x=2190, y=30, label="Slope = 1.44",hjust=1, alpha=0.7, size=3.2)+
  ylim(0,190)+
  xlab("Number of transects") +
  ylab("")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(size = 12, face="bold"),
        plot.margin=unit(c(0,0.15,0,0), "cm")) + 
  ggtitle("d")


## For panel e : Run script "scripts/03_Analysis/02_prediction_family_species_richness.R"
## Or load Rdata :
load("Rdata/all_predictions_motus_species.rdata")

# plot panel e : log nb of species in visual survey families ~ in eDNA families + prediction

gobiidae <- readPNG("data/fish_vignette/gobiidae.png")
muraenidae <- readPNG("data/fish_vignette/muraenidae.png")
pomacentridae <- readPNG("data/fish_vignette/pomacentridae.png")
labridae <- readPNG("data/fish_vignette/labridae.png")
mugilidae <- readPNG("data/fish_vignette/mugilidae.png")
kyphosidae <- readPNG("data/fish_vignette/kyphosidae.png")

e <- ggplot(fam_summary, aes(log_motu, log_checklist))+
  annotation_custom(rasterGrob(gobiidae), xmin = 5.1, xmax = 6, ymin = 4.2, ymax = 5)+
  annotation_custom(rasterGrob(muraenidae), xmin = 4.3, xmax = 5.2, ymin = 3, ymax = 3.5)+
  annotation_custom(rasterGrob(pomacentridae), xmin = 3.8, xmax = 4.7, ymin = 5.4, ymax = 6)+
  annotation_custom(rasterGrob(labridae), xmin = 4.9, xmax = 5.6, ymin = 5.5, ymax = 6.1)+
  annotation_custom(rasterGrob(mugilidae), xmin = 3.3, xmax = 4.3, ymin = 0.8, ymax = 1.4)+
  annotation_custom(rasterGrob(kyphosidae), xmin = 1.2, xmax = 2, ymin = 3.3, ymax = 4)+
  geom_point(size=2)+
  geom_abline(slope = 1, intercept = 0, color="red", size=0.8)+
  geom_abline(slope = 0.8, intercept = 0.49, size=0.8)+
  annotate(geom="text", x=6, y=1, label="y = 0.8x+0.49\nR² = 0.6\np < 0.001", hjust=1, size=3.2) +
  xlim(0,6)+
  ylim(0,6)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size=9),
        axis.title.x = element_text(size=9),
        plot.title = element_text(size = 12, face="bold"),
        plot.margin=unit(c(0,0.1,0.1,0.2), "cm")) + 
  labs(x="log(1+Nb of MOTUs per family in eDNA)", y="log(1+Nb of species per family in VC)")+
  ggtitle("e")




## For panel f : Run scripts "scripts/03_Analysis/03a_family_presence_proportions.R"
#                             "scripts/03_Analysis/04b_family_proportion_coverage_resolution.R"
## Or load Rdata :
load("Rdata/family_proportion_region_main.rdata")

# plot panel f : Proportion of each family globally and in each region

f <- ggplot(families_prop_global_main, aes(x=reorder(family, prop), y = prop, fill = Province)) + 
  geom_bar(stat="identity", show.legend = TRUE) + 
  theme_bw() +
  scale_fill_manual(values =c("#a6611a", "#E5A729",  "#015462", "#b2182b", "#80cdc1"))+ 
  labs(x="", y="Percentage of MOTUs \nat global scale")+
  ggtitle("f")+
  scale_y_continuous(labels = percent, breaks = c(0, 0.05, 0.10))+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA), 
        plot.title = element_text(size = 12, face="bold"),
        axis.title.x = element_text(size=10),
        plot.margin=unit(c(0.1,0.4,0.2,0), "cm"))+
  coord_flip()


# plot all together

a_e <- ggarrange(a, b, c, d, e, nrow = 3, ncol=2)
Fig2 <- ggarrange(a_e, f, ncol=2, widths = c(2,1))
Fig2
ggsave("outputs/00_Figures_for_paper/Figure2.png", width = 7.5, height = 7.5)


## legend added on powerpoint afterward
