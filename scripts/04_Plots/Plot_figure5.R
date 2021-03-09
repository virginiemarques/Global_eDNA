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

all_motus <- ggplot(rich_site, aes(col=region))+
  geom_jitter(aes(x=dist_to_CT, y=motu), shape=17, size=2, alpha=0.7, show.legend = FALSE) +
  geom_errorbar(aes(x=dist_to_CT, ymin=mean_motu-sd_motu, ymax=mean_motu+sd_motu), show.legend = FALSE, alpha=0.7)+
  geom_jitter(aes(x=dist_to_CT, y=mean_motu), shape=21, size=2, fill="white", alpha=0.7, show.legend = FALSE) +
  geom_point(data=all_region, aes(x=dist_to_CT, y=n_motus), shape=23, size=2.5, show.legend = FALSE) +
  geom_vline(xintercept = 14000, linetype="dashed", color="grey")+
  scale_color_manual(values=c("#E5A729", "#80cdc1", "#a6611a", "#b2182b", "#015462"))+ 
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
all_family <- ggplot(rich_site, aes(col=region))+
  geom_jitter(aes(x=dist_to_CT, y=family), shape=17, size=2, alpha=0.7, show.legend = FALSE) +
  geom_errorbar(aes(x=dist_to_CT, ymin=mean_family-sd_family, ymax=mean_family+sd_family), show.legend = FALSE, alpha=0.7)+
  geom_jitter(aes(x=dist_to_CT, y=mean_family), shape=21, size=2, alpha=0.7, fill="white", show.legend = FALSE) +
  geom_point(data=all_region, aes(x=dist_to_CT, y=n_family), shape=23, size=2.5, show.legend = FALSE) +
  geom_vline(xintercept = 14000, linetype="dashed", color="grey")+
  scale_color_manual(values=c("#E5A729", "#80cdc1", "#a6611a", "#b2182b", "#015462"))+ 
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

## For panel b : Run scripts "scripts/03_Analysis/03a_family_presence_proportions.R"
#                             "scripts/03_Analysis/03b_random_proportion_families.R"
## Or load Rdata :
load("Rdata/family_proportion_per_site.rdata")
load("Rdata/CI_null_model_family_proportions.rdata")



# plot panel b

family <- c("Acanthuridae", "Labridae", "Serranidae", "Carangidae", "Pomacentridae","Gobiidae")

prop <- vector("list")
for (i in 1:length(family)) {
  fam <- df_all_site[df_all_site$family == family[i],]
  fam_CI <- CI_family[CI_family$family == family[i],]
  prop[[i]] <- ggplot(fam, aes(n_motus_total, prop))+
    geom_smooth(data=fam_CI, aes(n_total, upper), col="black", size=0.5, show.legend = FALSE)+
    geom_smooth(data=fam_CI, aes(n_total, lower), col="black", size=0.5, show.legend = FALSE)+
    geom_point(size=2, aes(colour=region))+
    xlim(0, 600)+
    ylim(0,0.3)+
    theme_bw()+
    theme(legend.position = "none")+
    scale_color_manual(values =c("#E5A729", "#80cdc1", "#a6611a", "#b2182b", "#015462"))+ 
    labs(title=family[i], x="", y="")+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_rect(fill = NA),
          axis.title.y = element_text(size = 10, face = "bold"),
          plot.title = element_text(size=10, face = "bold"),
          plot.margin=unit(c(0,0.1,0,0), "cm"))
}

plot <- ggarrange(plotlist = prop, ncol=2, nrow = 3, common.legend = TRUE, legend = "top")

x.grob <- textGrob("Total number of MOTUs per site", 
                   gp=gpar(fontface="bold", col="black", fontsize=10))
y.grob <- textGrob("Proportion of MOTUs assigned to the family in each site", 
                   gp=gpar(fontface="bold", col="black", fontsize=10), rot = 90)
b <- grid.arrange(plot, bottom=x.grob, left=y.grob)


# plot all together

Fig2 <- ggarrange(a, b, nrow = 2, ncol = 1, labels = c("A", "B"), heights = c(1,2.5))
Fig2
ggsave("outputs/00_Figures_for_paper/Figure5.png", width = 7.5, height = 8, dpi = 300)
