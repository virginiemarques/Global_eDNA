library(tidyverse)
library(reshape2)
library(gridExtra)
library(grid)
library(cowplot)
library(ggplot2)
library(ggpubr)

load("Rdata/02_clean_all.Rdata")
'%ni%' <- Negate("%in%")
#Remove estuary stations and deep niskin station
df_all_filters <- df_all_filters %>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")) %>%
  filter(sample_method !="niskin" & region!="East_Pacific" & comment %ni% c("Distance decay 600m", "Distance decay 300m") & station!="glorieuse_distance_300m")%>%
  filter(project != "SEAMOUNTS") %>% 
  filter(habitat_type %ni% c("BAIE", "Sommet"))

df_all_filters <- df_all_filters %>%
  filter(!is.na(new_family_name))

global_motu <- df_all_filters %>%
  distinct(sequence, .keep_all = TRUE)
global_family<- as.data.frame(unique(global_motu$new_family_name))
colnames(global_family) <- "family"

# sample random family assignations for site with 1-800 species, 1000 times (long! load Rdata)
random_sp <- vector("list", 800)
random_prop_tot <- data.frame()
rep <- 800

for (j in seq(1:1000)) {
  random_prop <- vector("list" )
  for (i in 1:rep) {
    random_sp[[i]] <- as.data.frame(global_motu[sample(nrow(global_motu), i, replace = TRUE), ])
    random_sp[[i]] <- data.frame(table(random_sp[[i]]$new_family_name))
    colnames(random_sp[[i]]) <- c("family", "n_motus")
    random_sp[[i]] <- full_join(random_sp[[i]], global_family, by="family")
    random_sp[[i]][is.na(random_sp[[i]])] <- 0
    random_sp[[i]]$n_total <- sum(random_sp[[i]]$n_motus)
    random_sp[[i]]$prop <- random_sp[[i]]$n_motus / random_sp[[i]]$n_total
    random_prop <- rbind(random_prop, random_sp[[i]])
  }
  random_prop_tot <- rbind(random_prop_tot, random_prop)
}

save(random_prop_tot, file = "Rdata/random_family_proportions.rdata")

load("Rdata/random_family_proportions.rdata")
# calculate 2.5 and 97.5 quantile for each sample size and each family
family <- family <- c("Acanthuridae", "Chaetodontidae", "Labridae", "Lutjanidae", "Serranidae", "Carangidae", "Pomacentridae", "Apogonidae", "Gobiidae")

CI_family <- data.frame()
quant <- data.frame(upper=numeric(800), lower=numeric(800))
for (j in 1:length(family)) {
  df <- random_prop_tot[random_prop_tot$family==family[j],]
  for (i in seq(1:800)) {
    df2 <- df[df$n_total==i,]
    quant[i,] <- quantile(df2$prop, c(.975, .025))
    quant$n_total <- seq(1:800)
    quant$family <- family[j]
  }
  CI_family <- rbind(CI_family, quant)
}


load("Rdata/family_proportion_per_site.rdata")

prop <- vector("list")
for (i in 1:length(family)) {
  fam <- df_all_site[df_all_site$family == family[i],]
  fam_CI <- CI_family[CI_family$family == family[i],]
  prop[[i]] <- ggplot(fam, aes(n_motus_total, prop))+
    geom_point(size=2, aes(colour=region))+
    geom_line(data=fam_CI, aes(n_total, upper))+
    geom_line(data=fam_CI, aes(n_total, lower))+
    scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3))+
    xlim(0, 800)+
    ylim(0,0.3)+
    theme_bw()+
    theme(legend.position = "none")+
    scale_color_manual(values =c("#E5A729", "#8AAE8A", "#4F4D1D", "#C67052"))+ #863b34
    labs(title=family[i], x="", y="")+
    theme(plot.title = element_text(size = 10, face="bold"), plot.margin=unit(c(0,0.1,0,0), "cm"))
}



plot <- ggarrange(plotlist = prop, ncol=3, nrow = 3, common.legend = TRUE, legend = "top")
plot
x.grob <- textGrob("Total number of MOTUs per site", 
                   gp=gpar(fontface="bold", col="black", fontsize=12))
y.grob <- textGrob("Proportion of MOTUs assigned to the family in each site", 
                   gp=gpar(fontface="bold", col="black", fontsize=12), rot = 90)
plot_grid <- grid.arrange(plot, bottom=x.grob, left=y.grob)



load("Rdata/plot_richness_site~dist_CT.rdata")

ggarrange(plot_all_rich_site, plot_grid, nrow = 2, ncol = 1, labels = c("A", "B"), heights = c(1,3))

