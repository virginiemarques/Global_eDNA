library(tidyverse)
library(reshape2)
library(gridExtra)
library(grid)
library(cowplot)
library(ggplot2)
library(ggpubr)

setwd("c:/Users/mathon/Desktop/linux/Global_eDNA/")
load("Rdata/02_clean_all.Rdata")

#Remove estuary stations and deep niskin station
df_all_filters <- subset(df_all_filters, !(station %in% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")))
df_all_filters <- subset(df_all_filters, sample_method!="niskin")
df_all_filters <- subset(df_all_filters, region!="East_Pacific")
df_all_filters <- subset(df_all_filters, !(comment %in% c("Distance decay 600m", "Distance decay 300m")))
df_all_filters <- subset(df_all_filters, station!="glorieuse_distance_300m")


# Proportion is calculated as (nb of unique motus per family)/(total nb of unique motus per site/region), using only amplicons assigned to family level minimum



## Proportion of families in Caribbean

  ## Caribbean total
caribbean <- df_all_filters %>%
  filter(region=="Caribbean") 

cari_motu <- caribbean %>%
  distinct(sequence, .keep_all = TRUE)

count_families_caribbean <- data.frame(table(cari_motu$new_family_name))
colnames(count_families_caribbean) <- c("family", "n_motus")
count_families_caribbean$n_motus_total <- nrow(cari_motu)
count_families_caribbean$prop <- count_families_caribbean$n_motus / count_families_caribbean$n_motus_total

ggplot(count_families_caribbean, aes(x=reorder(family, prop), y = prop)) + 
  geom_bar(stat="identity") + 
  theme_bw() +
  labs(x="Family", y="Proportion")+
  theme(axis.text.x=element_text(angle = 0, hjust = 0)) +
  theme(legend.position = "none")+
  coord_flip()

ggsave("outputs/05_family_proportion/02_based_on_species_presence/per region/family_proportion_caribbean.png", width=6, height=16)
write.csv(count_families_caribbean, "outputs/05_family_proportion/02_based_on_species_presence/per region/family_proportion_caribbean.csv")

  ## Caribbean sites
site <- c(unique(caribbean$site))

count_families_site_caribbean=NULL 

for (i in 1:length(site)) {
  s <- site[i]
  cari_site <- caribbean[caribbean$site == site[i],]
  cari_motu_site <- cari_site%>%
    distinct(sequence, .keep_all = TRUE)
  count_families <- data.frame(table(cari_motu_site$new_family_name))
  colnames(count_families) <- c("family", "n_motus")
  count_families$n_motus_total <- nrow(cari_motu_site)
  count_families$prop <- count_families$n_motus / count_families$n_motus_total
  count_families <- count_families[order(count_families$prop, decreasing = TRUE),]
  count_families$site <- s
  count_families$region <- "Caribbean"
  count_families_site_caribbean <- rbind(count_families_site_caribbean, count_families)
}

write.csv(count_families_site_caribbean, "outputs/05_family_proportion/02_based_on_species_presence/per site/family_proportion_site_caribbean.csv")

count_families_site_caribbean <- count_families_site_caribbean %>%
  ungroup() %>%
  arrange(site, prop) %>%
  mutate(order = row_number())

ggplot(count_families_site_caribbean, aes(order, prop)) + 
  geom_bar(stat="identity") +
  facet_wrap(~site, scales ="free_y", as.table = FALSE)+
  theme_bw() +
  labs(x="Family", y="Proportion")+
  theme(axis.text.x=element_text(angle = 0, hjust = 0)) +
  theme(legend.position = "none")+
  coord_flip()+
  scale_x_continuous(breaks=count_families_site_caribbean$order, labels=count_families_site_caribbean$family, expand = c(0,0))

ggsave("outputs/05_family_proportion/02_based_on_species_presence/per site/family_proportion_site_caribbean.png", width=20, height=16)


  ## Caribbean station

station <- c(unique(caribbean$station))

count_families_station_caribbean=NULL 

for (i in 1:length(station)) {
  cari_station <- caribbean[caribbean$station == station[i],]
  st <- station[i]
  s <- unique(cari_station$site)
  cari_motu_station <- cari_station%>%
    distinct(sequence, .keep_all = TRUE)
  count_families <- data.frame(table(cari_motu_station$new_family_name))
  colnames(count_families) <- c("family", "n_motus")
  count_families$n_motus_total <- nrow(cari_motu_station)
  count_families$prop <- count_families$n_motus / count_families$n_motus_total
  count_families <- count_families[order(count_families$prop, decreasing = TRUE),]
  count_families$station <- st
  count_families$site <- s
  count_families$region <- "Caribbean"
  count_families_station_caribbean <- rbind(count_families_station_caribbean, count_families)
}

write.csv(count_families_station_caribbean, "outputs/05_family_proportion/02_based_on_species_presence/per station/family_proportion_station_caribbean.csv")


## Frequency of families in Lengguru
  ## Lengguru total

lengguru <- df_all_filters %>%
  filter(region=="Central_IndoPacific") 

leng_motu <- lengguru %>%
  distinct(sequence, .keep_all = TRUE)

count_families_lengguru <- data.frame(table(leng_motu$new_family_name))
colnames(count_families_lengguru) <- c("family", "n_motus")
count_families_lengguru$n_motus_total <- nrow(leng_motu)
count_families_lengguru$prop <- count_families_lengguru$n_motus / count_families_lengguru$n_motus_total


ggplot(count_families_lengguru, aes(x=reorder(family, prop), y = prop)) + 
  geom_bar(stat="identity") + 
  theme_bw() +
  labs(x="Family", y="Proportion")+
  theme(axis.text.x=element_text(angle = 0, hjust = 0)) +
  theme(legend.position = "none")+
  coord_flip()

ggsave("outputs/05_family_proportion/02_based_on_species_presence/per region/family_proportion_lengguru.png", width=6, height=16)
write.csv(count_families_lengguru, "outputs/05_family_proportion/02_based_on_species_presence/per region/family_proportion_lengguru.csv")

  ## Lengguru sites

site <- c(unique(lengguru$site))

count_families_site_lengguru=NULL 

for (i in 1:length(site)) {
  s <- site[i]
  leng_site <- lengguru[lengguru$site == site[i],]
  leng_motu_site <- leng_site%>%
    distinct(sequence, .keep_all = TRUE)
  count_families <- data.frame(table(leng_motu_site$new_family_name))
  colnames(count_families) <- c("family", "n_motus")
  count_families$n_motus_total <- nrow(leng_motu_site)
  count_families$prop <- count_families$n_motus / count_families$n_motus_total
  count_families <- count_families[order(count_families$prop, decreasing = TRUE),]
  count_families$site <- s
  count_families$region <- "Central_IndoPacific"
  count_families_site_lengguru <- rbind(count_families_site_lengguru, count_families)
}

write.csv(count_families_site_lengguru, "outputs/05_family_proportion/02_based_on_species_presence/per site/family_proportion_site_lengguru.csv")

count_families_site_lengguru <- count_families_site_lengguru %>%
  ungroup() %>%
  arrange(site, prop) %>%
  mutate(order = row_number())

  ## plot by 3 sites, because plot too small otherwise
site_sub <- c("pulau_pisang", "lobo", "tanjung_boi")
subset1 <- count_families_site_lengguru %>%
  filter(site%in%site_sub)
ggplot(subset1, aes(order, prop)) + 
  geom_bar(stat="identity") +
  facet_wrap(~site, scales ="free_y")+
  theme_bw() +
  labs(x="Family", y="Proportion")+
  theme(axis.text.x=element_text(angle = 0, hjust = 0)) +
  theme(legend.position = "none")+
  coord_flip()+
  scale_x_continuous(breaks=subset1$order, labels=subset1$family, expand = c(0,0))


ggsave("outputs/05_family_proportion/02_based_on_species_presence/per site/family_proportion_lengguru_site7-9.png", width=20, height=16)

  ## Lengguru station

station <- c(unique(lengguru$station))

count_families_station_lengguru=NULL 

for (i in 1:length(station)) {
  leng_station <- lengguru[lengguru$station == station[i],]
  st <- station[i]
  s <- unique(leng_station$site)
  leng_motu_station <- leng_station%>%
    distinct(sequence, .keep_all = TRUE)
  count_families <- data.frame(table(leng_motu_station$new_family_name))
  colnames(count_families) <- c("family", "n_motus")
  count_families$n_motus_total <- nrow(leng_motu_station)
  count_families$prop <- count_families$n_motus / count_families$n_motus_total
  count_families <- count_families[order(count_families$prop, decreasing = TRUE),]
  count_families$station <- st
  count_families$site <- s
  count_families$region <- "Central_IndoPacific"
  count_families_station_lengguru <- rbind(count_families_station_lengguru, count_families)
}

write.csv(count_families_station_lengguru, "outputs/05_family_proportion/02_based_on_species_presence/per station/family_proportion_station_lengguru.csv")


## proportion of families in Fakarava
  ## par site = region

fakarava <- df_all_filters %>%
  filter(region=="Central_Pacific")

faka_motu <- fakarava %>%
  distinct(sequence, .keep_all = TRUE)

count_families_site_fakarava <- data.frame(table(faka_motu$new_family_name))
colnames(count_families_site_fakarava) <- c("family", "n_motus")
count_families_site_fakarava$n_motus_total <- nrow(faka_motu)
count_families_site_fakarava$site <- "fakarava"
count_families_site_fakarava$region <- "Central_Pacific"
count_families_site_fakarava$prop <- count_families_site_fakarava$n_motus / count_families_site_fakarava$n_motus_total


ggplot(count_families_site_fakarava, aes(x=reorder(family, prop), y = prop)) + 
  geom_bar(stat="identity") + 
  theme_bw() +
  labs(x="Family", y="Proportion")+
  theme(axis.text.x=element_text(angle = 0, hjust = 0)) +
  theme(legend.position = "none")+
  coord_flip()

ggsave("outputs/05_family_proportion/02_based_on_species_presence/per region/family_proportion_fakarava.png", width=6, height=16)
write.csv(count_families_site_fakarava, "outputs/05_family_proportion/02_based_on_species_presence/per region/family_proportion_fakarava.csv")

  ## per station
station <- c(unique(fakarava$station))

count_families_station_fakarava=NULL 

for (i in 1:length(station)) {
  faka_station <- fakarava[fakarava$station == station[i],]
  st <- station[i]
  faka_motu_station <- faka_station%>%
    distinct(sequence, .keep_all = TRUE)
  count_families <- data.frame(table(faka_motu_station$new_family_name))
  colnames(count_families) <- c("family", "n_motus")
  count_families$n_motus_total <- nrow(faka_motu_station)
  count_families$prop <- count_families$n_motus / count_families$n_motus_total
  count_families <- count_families[order(count_families$prop, decreasing = TRUE),]
  count_families$station <- st
  count_families$site <- "fakarava"
  count_families$region <- "Central_Pacific"
  count_families_station_fakarava <- rbind(count_families_station_fakarava, count_families)
}

write.csv(count_families_station_fakarava, "outputs/05_family_proportion/02_based_on_species_presence/per station/family_proportion_station_fakarava.csv")


## frequency of families in Eparses
  ## eparse total

eparse <- df_all_filters %>%
  filter(region=="West_Indian") 

eparse_motu <- eparse %>%
  distinct(sequence, .keep_all = TRUE)

count_families_eparse <- data.frame(table(eparse_motu$new_family_name))
colnames(count_families_eparse) <- c("family", "n_motus")
count_families_eparse$n_motus_total <- nrow(eparse_motu)
count_families_eparse$prop <- count_families_eparse$n_motus / count_families_eparse$n_motus_total


ggplot(count_families_eparse, aes(x=reorder(family, prop), y = prop)) + 
  geom_bar(stat="identity") + 
  theme_bw() +
  labs(x="Family", y="Proportion")+
  theme(axis.text.x=element_text(angle = 0, hjust = 0)) +
  theme(legend.position = "none")+
  coord_flip()

ggsave("outputs/05_family_proportion/02_based_on_species_presence/per region/family_proportion_eparse.png", width=6, height=16)
write.csv(count_families_eparse, "outputs/05_family_proportion/02_based_on_species_presence/per region/family_proportion_eparse.csv")

  ## eparse sites

site <- c(unique(eparse$site))

count_families_site_eparse=NULL 

for (i in 1:length(site)) {
  s <- site[i]
  eparse_site <- eparse[eparse$site == site[i],]
  eparse_motu_site <- eparse_site%>%
    distinct(sequence, .keep_all = TRUE)
  count_families <- data.frame(table(eparse_motu_site$new_family_name))
  colnames(count_families) <- c("family", "n_motus")
  count_families$n_motus_total <- nrow(eparse_motu_site)
  count_families$prop <- count_families$n_motus / count_families$n_motus_total
  count_families <- count_families[order(count_families$prop, decreasing = TRUE),]
  count_families$site <- s
  count_families$region <- "West_Indian"
  count_families_site_eparse <- rbind(count_families_site_eparse, count_families)
}

write.csv(count_families_site_eparse, "outputs/05_family_proportion/02_based_on_species_presence/per site/family_proportion_site_eparse.csv")

count_families_site_eparse <- count_families_site_eparse %>%
  ungroup() %>%
  arrange(site, prop) %>%
  mutate(order = row_number())


ggplot(count_families_site_eparse, aes(order, prop)) + 
  geom_bar(stat="identity") +
  facet_wrap(~site, scales ="free_y", nrow = 1)+
  theme_bw() +
  labs(x="Family", y="Proportion")+
  theme(axis.text.x=element_text(angle = 0, hjust = 0)) +
  theme(legend.position = "none")+
  coord_flip()+
  scale_x_continuous(breaks=count_families_site_eparse$order, labels=count_families_site_eparse$family, expand = c(0,0))


ggsave("outputs/05_family_proportion/02_based_on_species_presence/per site/family_proportion_site_eparse.png", width=20, height=16)

  ## eparse station

station <- c(unique(eparse$station))

count_families_station_eparse=NULL 

for (i in 1:length(station)) {
  eparse_station <- eparse[eparse$station == station[i],]
  st <- station[i]
  s <- unique(eparse_station$site)
  eparse_motu_station <- eparse_station%>%
    distinct(sequence, .keep_all = TRUE)
  count_families <- data.frame(table(eparse_motu_station$new_family_name))
  colnames(count_families) <- c("family", "n_motus")
  count_families$n_motus_total <- nrow(eparse_motu_station)
  count_families$prop <- count_families$n_motus / count_families$n_motus_total
  count_families <- count_families[order(count_families$prop, decreasing = TRUE),]
  count_families$station <- st
  count_families$site <- s
  count_families$region <- "West_Indian"
  count_families_station_eparse <- rbind(count_families_station_eparse, count_families)
}

write.csv(count_families_station_eparse, "outputs/05_family_proportion/02_based_on_species_presence/per station/family_proportion_station_eparse.csv")


## Proportion of families global


global_motu <- df_all_filters %>%
  distinct(sequence, .keep_all = TRUE)


count_families_global <- data.frame(table(global_motu$new_family_name))
colnames(count_families_global) <- c("family", "n_motus")
save(count_families_global, file = "Rdata/nb_motus_per_family_global.Rdata")

count_families_global$n_motus_total <- nrow(global_motu)
count_families_global$prop <- count_families_global$n_motus / count_families_global$n_motus_total

write.csv(count_families_global, "outputs/05_family_proportion/02_based_on_species_presence/family_proportion_global.csv")

ggplot(count_families_global, aes(x=reorder(family, prop), y = prop)) + 
  geom_bar(stat="identity") + 
  theme_bw() +
  labs(x="Family", y="Proportion")+
  theme(legend.position = "none")+
  theme(axis.text.x=element_text(angle = 0, hjust = 0)) + 
  coord_flip()

ggsave("outputs/05_family_proportion/02_based_on_species_presence/family_proportion_global.png", width=6, height=16)



## plot proportion of each family in each region on global scale
families_prop_global <- left_join(count_families_global, count_families_lengguru[,1:2], by="family" )
families_prop_global <- left_join(families_prop_global, count_families_caribbean[,1:2], by="family" )
families_prop_global <- left_join(families_prop_global, count_families_site_fakarava[,1:2], by="family" )
families_prop_global <- left_join(families_prop_global, count_families_eparse[,1:2], by="family" )
families_prop_global <- families_prop_global[, c(-3)]
colnames(families_prop_global) <- c("family", "n_global", "prop_global", "n_lengguru", "n_caribbean", "n_fakarava", "n_eparse")
families_prop_global[is.na(families_prop_global)] <- 0

for (i in 1:dim(families_prop_global)[1]) {
    families_prop_global[i,"new_n_leng"] <- (families_prop_global[i,"n_lengguru"]*families_prop_global[i,"n_global"])/sum(families_prop_global[i,4:7])
    families_prop_global[i,"new_n_car"] <- (families_prop_global[i,"n_caribbean"]*families_prop_global[i,"n_global"])/sum(families_prop_global[i,4:7])
    families_prop_global[i,"new_n_faka"] <- (families_prop_global[i,"n_fakarava"]*families_prop_global[i,"n_global"])/sum(families_prop_global[i,4:7])
    families_prop_global[i,"new_n_eparse"] <- (families_prop_global[i,"n_eparse"]*families_prop_global[i,"n_global"])/sum(families_prop_global[i,4:7])
}

families_prop_global$Lengguru <- (families_prop_global$new_n_leng*families_prop_global$prop_global)/families_prop_global$n_global
families_prop_global$Caribbean <- (families_prop_global$new_n_car*families_prop_global$prop_global)/families_prop_global$n_global
families_prop_global$Fakarava <- (families_prop_global$new_n_faka*families_prop_global$prop_global)/families_prop_global$n_global
families_prop_global$Eparse <- (families_prop_global$new_n_eparse*families_prop_global$prop_global)/families_prop_global$n_global

families_prop_global <- families_prop_global[,c(1,12:15)]
families_prop_global2 <- melt(families_prop_global)
colnames(families_prop_global2) <- c("family", "Region", "prop")

save(families_prop_global2, file = "Rdata/family_proportion_global.Rdata")

ggplot(families_prop_global2, aes(x=reorder(family, prop), y = prop, fill = Region)) + 
  geom_bar(stat="identity", show.legend = TRUE) + 
  theme_bw() +
  scale_fill_manual(values =c("#8AAE8A", "#E5A729", "#4F4D1D", "#C67052"))+
  labs(x="Family", y="Proportion")+
  coord_flip()

ggsave("outputs/05_family_proportion/02_based_on_species_presence/family_proportion_global_region.png", width=6, height=16)


## Bellwood figures : proportion of families per site

#family_leng <- unique(lengguru$new_family_name)
#family_car <- unique(caribbean$new_family_name)
#family_faka <- unique(fakarava$new_family_name)

#fam <- intersect(family_car, family_faka)
#family <- intersect(fam, family_leng)

df_all_site <- rbind(count_families_site_lengguru, count_families_site_caribbean, count_families_site_eparse, count_families_site_fakarava)


family <- c("Acanthuridae", "Chaetodontidae", "Labridae", "Lutjanidae", "Serranidae", "Carangidae", "Pomacentridae", "Apogonidae", "Gobiidae")

prop <- vector("list")
for (i in 1:length(family)) {
  fam <- df_all_site[df_all_site$family == family[i],]
  prop[[i]] <- ggplot(fam, aes(n_motus_total, prop, ymin=0, ymax=0.2, colour=region))+
    geom_point(size=2)+
    scale_y_continuous(breaks = c(0, 0.1, 0.2))+
    xlim(0, 800)+
    theme(legend.position = "none")+
    scale_color_manual(values =c("#E5A729", "#8AAE8A", "#4F4D1D", "#C67052"))+
    labs(title=family[i], x="", y="")+
    theme(plot.title = element_text(size = 10, face="bold"), plot.margin=unit(c(0,0.1,0,0), "cm"))
}


plot <- ggarrange(plotlist = prop, ncol=3, nrow = 3, common.legend = TRUE, legend = "top")
x.grob <- textGrob("Total number of MOTUs per site", 
                   gp=gpar(fontface="bold", col="black", fontsize=12))
y.grob <- textGrob("Proportion of MOTUs assigned to the family in each site", 
                   gp=gpar(fontface="bold", col="black", fontsize=12), rot = 90)
plot_grid <- grid.arrange(plot, bottom=x.grob, left=y.grob)



load("Rdata/plot_richness_site~dist_CT.rdata")

ggarrange(plot_all_rich_site, plot_grid, nrow = 2, ncol = 1, labels = c("A", "B"), heights = c(1,3))



## same for stations
df_all_station <- rbind(count_families_station_lengguru, count_families_station_caribbean, count_families_station_eparse, count_families_station_fakarava)


family <- c("Acanthuridae", "Chaetodontidae", "Labridae", "Lutjanidae", "Serranidae", "Carangidae", "Pomacentridae", "Apogonidae", "Gobiidae")

prop <- vector("list")
for (i in 1:length(family)) {
  fam <- df_all_station[df_all_station$family == family[i],]
  prop[[i]] <- ggplot(fam, aes(n_motus_total, prop, ymin=0, ymax=0.5, colour=region))+
    geom_point(size=2)+
    scale_y_continuous(breaks = c(0, 0.2, 0.4))+
    xlim(0, 310)+
    theme(legend.position = "none")+
    scale_color_manual(values =c("#E5A729", "#8AAE8A", "#4F4D1D", "#C67052"))+
    labs(title=family[i], x="", y="")+
    theme(plot.title = element_text(size = 10, face="bold"), plot.margin=unit(c(0,0.1,0,0), "cm"))
}


plot <- ggarrange(plotlist = prop, ncol=3, nrow = 3, common.legend = TRUE, legend = "top")
x.grob <- textGrob("Total number of MOTUs per station", 
                   gp=gpar(fontface="bold", col="black", fontsize=12))
y.grob <- textGrob("Proportion of MOTUs assigned to the family in each station", 
                   gp=gpar(fontface="bold", col="black", fontsize=12), rot = 90)
plot_grid <- grid.arrange(plot, bottom=x.grob, left=y.grob)



load("Rdata/plot_richness~dist_CT.rdata")

ggarrange(plot_all_rich_station, plot_grid, nrow = 2, ncol = 1, labels = c("A", "B"), heights = c(1,3))

## test chi²

  ## on all families, all sites
test1 <- chisq.test(df_all_site$prop)

  ## selected families, all sites
subset_fam <- df_all_site[df_all_site$family%in%family,]
test2 <- chisq.test(subset_fam$prop)

## selected families, by sites
subset_fam <- subset_fam[, -c(2,3,6)]
subset_fam_spread <- spread(subset_fam, family, prop)
subset_fam_spread[is.na(subset_fam_spread)] <- 0

test3 <- chisq.test(subset_fam_spread[,c(-1)])
