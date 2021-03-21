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


## For panel a : Run scripts "scripts/03_Analysis/05_richess_motu_family.R"  
## Or load Rdata :
load("Rdata/richness_station_site_region.rdata")
load("Rdata/richness_motu_region.rdata")

# plot panel a left : motu richness per region

all_motus <- ggplot(rich_site, aes(col=province))+
  geom_point(aes(x=dist_to_CT, y=motu), shape=17, size=2, alpha=0.7, show.legend = FALSE) +
  geom_errorbar(aes(x=dist_to_CT, ymin=mean_motu-sd_motu, ymax=mean_motu+sd_motu), show.legend = FALSE, alpha=0.7)+
  geom_point(aes(x=dist_to_CT, y=mean_motu), shape=21, size=2, fill="white", alpha=0.7, show.legend = FALSE) +
  geom_point(data=all_province, aes(x=dist_to_CT, y=n_motus), shape=23, size=2.5, show.legend = FALSE) +
  geom_vline(xintercept = 14000, linetype="dashed", color="grey")+
  scale_color_manual(values=c("#a6611a", "#E5A729","#b2182b", "#80cdc1", "#015462"))+ 
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size = 10, face = "bold"), 
        plot.margin=unit(c(0.2,0.4,0,0.1), "cm"), 
        text = element_text(size=12))+
  labs(x="",y="MOTU richness")

# plot panel a right : family richness per region

load("Rdata/richness_family_region.rdata")
all_family <- ggplot(rich_site, aes(col=province))+
  geom_jitter(aes(x=dist_to_CT, y=family), shape=17, size=2, alpha=0.7, show.legend = FALSE) +
  geom_errorbar(aes(x=dist_to_CT, ymin=mean_family-sd_family, ymax=mean_family+sd_family), show.legend = FALSE, alpha=0.7)+
  geom_jitter(aes(x=dist_to_CT, y=mean_family), shape=21, size=2, alpha=0.7, fill="white", show.legend = FALSE) +
  geom_point(data=all_province, aes(x=dist_to_CT, y=n_family), shape=23, size=2.5, show.legend = FALSE) +
  geom_vline(xintercept = 14000, linetype="dashed", color="grey")+
  scale_color_manual(values=c("#a6611a", "#E5A729","#b2182b", "#80cdc1", "#015462"))+ 
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size = 10, face = "bold"), 
        plot.margin=unit(c(0.2,0.4,0,0.1), "cm"), 
        text = element_text(size=12))+
  labs(x="",y="Family richness")

# plot panel a 

plot <- ggarrange(all_motus, all_family, ncol = 2, nrow=1)
x.grob <- textGrob("Distance to Coral Triangle (km, W-E)", 
                   gp=gpar(fontface="bold", col="black", fontsize=10), vjust = -0.5)

a <- grid.arrange(plot, bottom=x.grob)

ggsave(a, file="outputs/00_Figures_for_paper/Extended_Data/ED_Figure4.png", width=8, height = 4)
