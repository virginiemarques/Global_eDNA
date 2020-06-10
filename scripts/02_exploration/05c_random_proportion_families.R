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

unique_motus <- unique(df_all_filters$sequence)
unique_family <- as.data.frame(unique(df_all_filters$new_family_name))
colnames(unique_family) <- "family"

global_family <- df_all_filters %>%
  distinct(sequence, new_family_name)
colnames(global_family) <- c("motu", "family")

#calculate frequency of motus
site <- unique(df_all_filters$site)

df_site <- data.frame(motu=character(), stringsAsFactors = FALSE)
for (i in 1:length(site)) {
  df <- df_all_filters[df_all_filters$site==site[i],] %>%
    distinct(sequence, site)
  colnames(df) <- c("motu", site[i])
  df_site <- full_join(df_site, df, by="motu")
}
rownames(df_site) <- df_site[,1]
df_site <- decostand(df_site[,c(-1)], "pa",na.rm = TRUE)
df_site[is.na(df_site)] <- 0
df_site <- as.data.frame(t(df_site))

prob_motu <- colSums(df_site)/25

# sample random family assignations for site with 1-800 species, 1000 times (long! load Rdata)
random_sp <- vector("list", 800)
random_prop_tot <- data.frame()
rep <- 800

for (j in seq(1:1000)) {
  random_prop <- vector("list" )
  for (i in 1:rep) {
    random_sp[[i]] <- as.data.frame(sample(unique_motus, i, replace = TRUE, prob = prob_motu))
    colnames(random_sp[[i]]) <- "motu"
    random_sp[[i]] <- left_join(random_sp[[i]], global_family, by="motu")
    random_sp[[i]] <- data.frame(table(random_sp[[i]]$family))
    colnames(random_sp[[i]]) <- c("family", "n_motus")
    random_sp[[i]] <- full_join(random_sp[[i]], unique_family, by="family")
    random_sp[[i]][is.na(random_sp[[i]])] <- 0
    random_sp[[i]]$n_total <- sum(random_sp[[i]]$n_motus)
    random_sp[[i]]$prop <- random_sp[[i]]$n_motus / random_sp[[i]]$n_total
    random_prop <- rbind(random_prop, random_sp[[i]])
  }
  random_prop_tot <- rbind(random_prop_tot, random_prop)
}

save(random_prop_tot, file = "c:/Users/mathon/Desktop/PhD/Projets/Megafauna/Global_eDNA/Rdata/random_family_proportions.rdata")

load("c:/Users/mathon/Desktop/PhD/Projets/Megafauna/Global_eDNA/Rdata/random_family_proportions.rdata")



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

ggarrange(plot_all_rich_site, plot_grid, nrow = 2, ncol = 1, labels = c("A", "B"), heights = c(1,2.5))

ggsave("outputs/Figures papier/Figure3.png", width = 7.5, height = 9)
