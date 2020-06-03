library(tidyverse)
library(reshape2)
library(gridExtra)
library(grid)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(vegan)

load("Rdata/02_clean_all.Rdata")
'%ni%' <- Negate("%in%")
#Remove estuary stations and deep niskin station
df_all_filters <- df_all_filters %>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")) %>%
  filter(sample_method !="niskin" & region!="East_Pacific" & comment %ni% c("Distance decay 600m", "Distance decay 300m") & station!="glorieuse_distance_300m")%>%
  filter(project != "SEAMOUNTS") %>% 
  filter(habitat_type %ni% c("BAIE", "Sommet"))

# remove motus not assigned to family
df_all_filters <- df_all_filters %>%
  filter(!is.na(new_family_name))


motu_fam <- df_all_filters[,c("sequence", "new_family_name")]%>%
  distinct(sequence, .keep_all = TRUE)


site <- unique(df_all_filters$site)
family <- unique(df_all_filters$new_family_name)

# calculate motus total per site
site_motu <- data.frame(site=character(25), total_motus=numeric(25), stringsAsFactors = FALSE)
for (i in 1:length(site)) {
  site_motu[i,"site"] <- site[[i]]
  site_motu[i,"total_motus"] <- df_all_filters[df_all_filters$site==site[i],]%>%summarise(n=n_distinct(sequence))
}


# create matrice site-motus
df_site=data.frame(motu=character())

for (i in 1:length(site)) {
  df <- df_all_filters[df_all_filters$site == site[i],] %>%
    distinct(sequence, site)
  colnames(df) <- c("motu", site[i], family)
  df_site <- full_join(df_site, df, by="motu")
}
rownames(df_site) <- df_site[,1]
df_site <- decostand(df_site[,c(-1)], "pa",na.rm = TRUE)
df_site[is.na(df_site)] <- 0
df_site <- as.data.frame(t(df_site))


# random matrices for 1000 iterations
rep = 1000
null.algo <- nullmodel(df_site, "curveball")
random_mat <- simulate(null.algo, nsim=rep)

# get the family proportion for each family, in each site, in each random simulation
df_tot <- data.frame()
for (k in 1:1000) {
  rd_mat <- as.data.frame(random_mat[,,k])
  rd_mat$site <- rownames(rd_mat)
  rd_mat<- melt(rd_mat, id="site")
  colnames(rd_mat) <- c("site", "sequence", "presence")
  rd_mat <- left_join(rd_mat, motu_fam, by="sequence")
  rd_mat <- left_join(rd_mat, site_motu, by="site")
  
  df_prop_site <- data.frame(site=character(), family=character(), total_motus=numeric(), prop=numeric(), stringsAsFactors = FALSE)
  df_prop <- data.frame()
  for (i in 1:length(site)) {
    df <- rd_mat %>%
      subset(site==site[i])
    total_motus <- sum(df$presence)
    for (j in 1:length(family)) {
      df2 <- df%>%
        subset(new_family_name==family[j])
      df_prop_site[j, "site"] <- site[i]
      df_prop_site[j, "family"] <- family[j]
      df_prop_site[j, "total_motus"] <- total_motus
      df_prop_site[j, "prop"] <- sum(df2$presence)/total_motus
    }
    df_prop <- rbind(df_prop, df_prop_site)
  }
  df_tot <- rbind(df_tot, df_prop)
}

random_prop_tot <- df_tot

save(random_prop_tot, file = "c:/Users/mathon/Desktop/PhD/Projets/Megafauna/Global_eDNA/Rdata/random_family_proportions.rdata")



# calculate 2.5 and 97.5 quantile for each sample size and each family
load("c:/Users/mathon/Desktop/PhD/Projets/Megafauna/Global_eDNA/Rdata/random_family_proportions.rdata")

family <- c("Acanthuridae", "Chaetodontidae", "Labridae", "Lutjanidae", "Serranidae", "Carangidae", "Pomacentridae", "Apogonidae", "Gobiidae")

CI_family <- data.frame()
quant <- data.frame(upper=numeric(), lower=numeric())
for (j in 1:length(family)) {
  df <- random_prop_tot[random_prop_tot$family==family[j],]
  for (i in 1:length(site)) {
    df2 <- df[df$site==site[i],]
    quant[i,] <- quantile(df2$prop, c(.975, .025))
    quant[i, "n_total"] <- unique(df2$total_motus)
    quant$family <- family[j]
  }
  CI_family <- rbind(CI_family, quant)
}


# plot on top of family proportions in our data
load("Rdata/family_proportion_per_site.rdata")

prop <- vector("list")
for (i in 1:length(family)) {
  fam <- df_all_site[df_all_site$family == family[i],]
  fam_CI <- CI_family[CI_family$family == family[i],]
  prop[[i]] <- ggplot(fam, aes(n_motus_total, prop))+
    #geom_smooth(data=fam_CI, aes(n_total, upper), col="black", size=0.5, show.legend = FALSE)+
    #geom_smooth(data=fam_CI, aes(n_total, lower), col="black", size=0.5, show.legend = FALSE)+
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


# assemble with part A of the figure 3
load("Rdata/plot_richness_site~dist_CT.rdata")

ggarrange(plot_all_rich_site, plot_grid, nrow = 2, ncol = 1, labels = c("A", "B"), heights = c(1,2.5))

ggsave("outputs/Figures papier/Figure3.png", width = 7.5, height = 9)
