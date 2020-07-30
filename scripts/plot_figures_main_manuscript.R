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

#################################################################################################################################•
#### Plot figure 1

# a
load("Rdata/accumulation_asymptote_motus_all.rdata")
a <- ggplot(df_motus) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.8, fill="#d2981a") +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "solid", size = 1, col="#d2981a") +
  annotate(geom="text", x=250, y=2789+150, label="Asymptote : 2789", hjust=1, color="#d2981a", size=3.2) +
  annotate(geom="text", x=250, y=2160+150, label="eDNA MOTUs : 2160", hjust=1, color="#d2981a", size=3.2) +
  annotate(geom="text", x=250, y=600, label="Slope = 2.3", hjust=1, alpha=0.7, size=3.2) +
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
        plot.margin=unit(c(0,0.1,0,0.1), "cm")) + 
  ggtitle("a")

# b
load("Rdata/accumulation_species_RLS.rdata")
b <- ggplot(all_accumulation_species_RLS_df) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.5) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 1, alpha=0.7) +
  annotate(geom="text", x=2990, y=2403+150, label="Asymptote : 2403",hjust=1, alpha=0.7, size=3.2) +
  annotate(geom="text", x=2990, y=1887+150, label="RLS Species : 1887",hjust=1, alpha=0.7, size=3.2) +
  annotate(geom="text", x=2990, y=600, label="Slope = 1.69",hjust=1, alpha=0.7, size=3.2) +
  ylim(0,3000)+
  xlab("Number of transects") +
  ylab("")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        panel.border = element_rect(fill = NA),
        plot.margin=unit(c(0,0.15,0,0), "cm")) + 
  ggtitle("b")

# c
load("Rdata/accumulation_asymptote_families_all.rdata")
c <- ggplot(df_fam) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.8, fill="#457277") +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "solid", size = 1, col="#457277") +
  annotate(geom="text", x=250, y=177+10, label="Asymptote : 177",hjust=1,color="#457277", size=3.2) +
  annotate(geom="text", x=250, y=143+10, label="eDNA Family : 143",hjust=1,color="#457277", size=3.2) +
  annotate(geom="text", x=250, y=30, label="Slope = 1.85",hjust=1, alpha=0.7, size=3.2)+
  ylim(0,190)+
  labs(y="", x="Number of samples")+
  ylab("Number of families") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        panel.border = element_rect(fill = NA),
        plot.margin=unit(c(0,0.1,0,0.1), "cm")) + 
  ggtitle("c")

# d 
load("Rdata/accumulation_families_RLS.rdata")
d <- ggplot(all_accumulation_RLS_df) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.5) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 1, alpha=0.7) +
  annotate(geom="text", x=2990, y=118+10, label="Asymptote : 118",hjust=1, alpha=0.7, size=3.2) +
  annotate(geom="text", x=2990, y=96+10, label="RLS Family : 96",hjust=1, alpha=0.7, size=3.2)+
  annotate(geom="text", x=2990, y=30, label="Slope = 1.44",hjust=1, alpha=0.7, size=3.2)+
  ylim(0,190)+
  xlab("Number of transects") +
  ylab("")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        panel.border = element_rect(fill = NA),
        plot.margin=unit(c(0,0.15,0,0), "cm")) + 
  ggtitle("d")

# e 
load("Rdata/all_predictions_motus_species.rdata")

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
  geom_abline(slope = 0.85, intercept = 0.5, size=0.8)+
  annotate(geom="text", x=6, y=1, label="y = 0.83x+0.6\nR² = 0.56\np < 0.001", hjust=1, size=3.2) +
  xlim(0,6)+
  ylim(0,6)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        plot.margin=unit(c(0,0.1,0.1,0.2), "cm")) + 
  labs(x="log(1+Nb of MOTUs in eDNA)", y="log(1+Nb of species in RLS)")+
  ggtitle("e")

# f
load("Rdata/family_proportion_region_main.rdata")
f <- ggplot(families_prop_global_main, aes(x=reorder(family, prop), y = prop, fill = Region)) + 
  geom_bar(stat="identity", show.legend = TRUE) + 
  theme_bw() +
  scale_fill_manual(values =c("#4F4D1D", "#E5A729",  "#C67052", "#863b34", "#8AAE8A"))+ 
  labs(x="", y="Proportion of MOTUs \nat global scale")+
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


a_e <- ggarrange(a, b, c, d, e, nrow = 3, ncol=2)
all <- ggarrange(a_e, f, ncol=2, widths = c(2,1))
all
ggsave("outputs/Figures papier/Figure1.png", width = 7.5, height = 7.5)


###########################################################################################################################
### Plot figure 2

# a
load("Rdata/richness_station_site_region.rdata")
load("Rdata/richness_motu_region.rdata")
all_motus <- ggplot(rich_site, aes(col=region))+
  geom_jitter(aes(x=dist_to_CT, y=motu), shape=17, size=2, alpha=0.7, show.legend = FALSE) +
  geom_errorbar(aes(x=dist_to_CT, ymin=mean_motu-sd_motu, ymax=mean_motu+sd_motu), show.legend = FALSE, alpha=0.7)+
  geom_jitter(aes(x=dist_to_CT, y=mean_motu), shape=21, size=2, fill="white", alpha=0.7, show.legend = FALSE) +
  geom_point(data=all_region, aes(x=dist_to_CT, y=n_motus), shape=23, size=2.5, show.legend = FALSE) +
  geom_vline(xintercept = 14000, linetype="dashed", color="grey")+
  scale_color_manual(values=c("#E5A729", "#8AAE8A", "#4F4D1D", "#863b34", "#C67052"))+ 
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size = 10, face = "bold"), 
        plot.margin=unit(c(0.2,0.4,0,0.1), "cm"), 
        text = element_text(size=12))+
  labs(x="",y="MOTU richness")

load("Rdata/richness_family_region.rdata")
all_family <- ggplot(rich_site, aes(col=region))+
  geom_jitter(aes(x=dist_to_CT, y=family), shape=17, size=2, alpha=0.7, show.legend = FALSE) +
  geom_errorbar(aes(x=dist_to_CT, ymin=mean_family-sd_family, ymax=mean_family+sd_family), show.legend = FALSE, alpha=0.7)+
  geom_jitter(aes(x=dist_to_CT, y=mean_family), shape=21, size=2, alpha=0.7, fill="white", show.legend = FALSE) +
  geom_point(data=all_region, aes(x=dist_to_CT, y=n_family), shape=23, size=2.5, show.legend = FALSE) +
  geom_vline(xintercept = 14000, linetype="dashed", color="grey")+
  scale_color_manual(values=c("#E5A729", "#8AAE8A", "#4F4D1D", "#863b34", "#C67052"))+ 
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size = 10, face = "bold"), 
        plot.margin=unit(c(0.2,0.4,0,0.1), "cm"), 
        text = element_text(size=12))+
  labs(x="",y="Family richness")

plot <- ggarrange(all_motus, all_family, ncol = 2, nrow=1)
x.grob <- textGrob("Distance to Coral Triangle (km, W-E)", 
                   gp=gpar(fontface="bold", col="black", fontsize=10), vjust = -0.5)

a <- grid.arrange(plot, bottom=x.grob)


# b
load("Rdata/family_proportion_per_site.rdata")
load("Rdata/CI_null_model_family_proportions.rdata")
family <- c("Acanthuridae", "Labridae", "Serranidae", "Carangidae", "Pomacentridae","Gobiidae")

prop <- vector("list")
for (i in 1:length(family)) {
  fam <- df_all_site[df_all_site$family == family[i],]
  fam_CI <- CI_family[CI_family$family == family[i],]
  prop[[i]] <- ggplot(fam, aes(n_motus_total, prop))+
    geom_smooth(data=fam_CI, aes(n_total, upper), col="black", size=0.5, show.legend = FALSE)+
    geom_smooth(data=fam_CI, aes(n_total, lower), col="black", size=0.5, show.legend = FALSE)+
    geom_point(size=2, aes(colour=region))+
    #scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3))+
    xlim(0, 600)+
    ylim(0,0.3)+
    theme_bw()+
    theme(legend.position = "none")+
    scale_color_manual(values =c("#E5A729", "#8AAE8A", "#4F4D1D", "#863b34", "#C67052"))+ 
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


ggarrange(a, b, nrow = 2, ncol = 1, labels = c("a", "b"), heights = c(1,2.5))

ggsave("outputs/Figures papier/Figure2.png", width = 7.5, height = 8)


#############################################################################################################################
### Plot figure 3

# a = upset plot

load("Rdata/upset_plot_motus_region.rdata")

# b = plot log10(occurence species/motus)

# fit models on eDNA
load("Rdata/rarete_motu_station.rdata")
tab=as.data.frame(motu_station)

ls=fitsad(tab[,2], "ls")
ln=fitsad(tab[,2],"lnorm")
po=fitsad(tab[,2], "pareto")

AICtab(ls, po, ln,weights=TRUE)


# fit models on RLS
load("Rdata/rarete_species_transects.rdata")
tab2 <- species_transects

ls2=fitsad(tab2[,2], "ls")
ln2=fitsad(tab2[,2],"lnorm")
po2=fitsad(tab2[,2], "pareto")

AICtab(ls2, po2, ln2, weights=TRUE)

# transform log10
tab$log10_nstation <- log10(tab$n)
tab$log10_nmotus <- log10(tab$n_motus)

tab2$log10_ntransect <- log10(tab2[,1])
tab2$log10_nspecies <- log10(tab2[,2])

# linear regressions log-log
lm(tab$log10_nmotus ~ tab$log10_nstation)
lm(tab2$log10_nspecies ~ tab2$log10_ntransect)

# plot figure 3b (with linear regression slopes)
ggplot(tab, aes(x=log10_nstation, y=log10_nmotus))+
  geom_point(colour="#d2981a", size=2, show.legend = TRUE)+
  geom_abline(intercept = 3.259, slope = -1.707, colour="#d2981a")+
  geom_point(data=tab2, aes(x=log10_ntransect, y=log10_nspecies), size=2, show.legend = TRUE)+
  geom_abline(intercept = 2.457, slope = -1.004)+
  xlim(0,3)+
  ylim(0,3)+
  annotate(geom="text", x=3, y=3, label="eDNA MOTUs ~ stations, beta=0.45", hjust=1, size=4, colour="#d2981a") +
  annotate(geom="text", x=3, y=2.8, label="RLS species ~ transects, beta=1.21", hjust=1, size=4) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        plot.title = element_text(size=12, face = "bold"))+
  labs(x="log10(Number of sampling units)",y="log10(Number of species/MOTUs)")

ggsave("outputs/Figures papier/Figure4b.png")





# Histogrammes rarete supp 

load("Rdata/02_clean_all.Rdata")
df_all_filters <- df_all_filters %>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")) %>%
  filter(sample_method !="niskin" & region!="East_Pacific" & comment %ni% c("Distance decay 600m", "Distance decay 300m") & station!="glorieuse_distance_300m")%>%
  filter(project != "SEAMOUNTS") %>% 
  filter(habitat_type %ni% c("BAIE", "Sommet"))


occu_edna <- df_all_filters %>%
  group_by(sequence) %>%
  summarise(n = n_distinct(station))
  
occu_edna$perc <- (occu_edna$n/145)*100

a <- ggplot(occu_edna, aes(x=reorder(sequence, 1-perc), y=perc))+
  geom_bar(stat = "identity", color="#d2981a")+
  ylim(0,65)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.text.x = element_blank())+
  labs(x="eDNA MOTUs", y="Percentage of stations where MOTU detected")
ggsave("outputs/Figures papier/Figure4b.png")

# c percentage of transects where each species detected
occu_RLS <- as.data.frame(occ_RLS)
occu_RLS$species <- rownames(occu_RLS)
occu_RLS$perc <- (occu_RLS$occ_RLS/2813)*100

b <- ggplot(occu_RLS, aes(x=reorder(species, 1-perc), y=perc))+
  geom_bar(stat = "identity")+
  ylim(0,65)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.text.x = element_blank())+
  labs(x="RLS species", y="Percentage of transects where species detected")
ggsave("outputs/Figures papier/Figure4c.png")

ggarrange(a, b, nrow=2, labels = c("a", "b"))
ggsave("outputs/Figures papier/ED_Figurexx.png")

## plot les deux courbes rarete sur le meme graphe
occu_edna <- dplyr::arrange(occu_edna, desc(n))
occu_edna$rank <- c(1:2160)
occu_RLS <- dplyr::arrange(occu_RLS, desc(perc))
occu_RLS$rank <- c(1:1786)


ggplot(occu_edna, aes(x=rank, y=perc))+
  geom_point(color="#d2981a")+
  geom_point(data=occu_RLS, aes(x=rank, y=perc))+
  geom_segment(aes(x=1786, y=0, xend=1786, yend=20), linetype="dashed", size=1)+
  geom_segment(aes(x=2160, y=0, xend=2160, yend=20), linetype="dashed", size=1, color="#d2981a")+
  annotate(geom="text", x=1900, y=22, label="1786 species", hjust=1, size=4) +
  annotate(geom="text", x=2300, y=22, label="2160 species", hjust=1, size=4, color="#d2981a") +
  annotate(geom="text", x=2200, y=55, label="ADN environnemental", hjust=1, size=6, color="#d2981a") +
  annotate(geom="text", x=2200, y=60, label="Plongée", hjust=1, size=6) +
  ylim(0,65)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.text.x = element_blank())+
  labs(x="", y="")
ggsave("c:/Users/mathon/Desktop/rareté.png")
