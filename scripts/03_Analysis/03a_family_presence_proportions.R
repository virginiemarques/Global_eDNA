library(tidyverse)
library(reshape2)
library(gridExtra)
library(grid)
library(cowplot)
library(ggplot2)
library(ggpubr)

load("Rdata/02-clean-data.Rdata")
'%ni%' <- Negate("%in%")

#Remove estuary stations and deep niskin station
df_all_filters <- df_all_filters %>%
  filter(province %in% c("Western_Indian_Ocean", "Southeast_Polynesia", "Tropical_Northwestern_Atlantic", "Western_Coral_Triangle", "Tropical_Southwestern_Pacific"))%>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3", "glorieuse_distance_300m")) %>%
  filter(sample_method !="niskin" & comment %ni% c("Distance decay 600m", "Distance decay 300m"))%>%
  filter(project != "Curacao") %>%
  filter(habitat=="marine")%>%
  filter(habitat_type %ni% c("BAIE"))%>%
  filter(depth<40) %>%
  filter(family_name_corrected %ni% "Salmonidae")


df_all_filters <- df_all_filters %>%
  filter(!is.na(family_name_corrected))

## Proportion of families in Caribbean

  ## Caribbean total
caribbean <- df_all_filters %>%
  filter(province=="Tropical_Northwestern_Atlantic")%>%
  filter(site35!="")

cari_motu <- caribbean %>%
  distinct(sequence, .keep_all = TRUE)

count_families_caribbean <- data.frame(table(cari_motu$family_name_corrected))
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

ggsave("outputs/04_family_proportion/Region/family_proportion_caribbean.png", width=6, height=16)
write.csv(count_families_caribbean, "outputs/04_family_proportion/Region/family_proportion_caribbean.csv")

  ## Caribbean sites
site <- c(unique(caribbean$site35))

count_families_site_caribbean=NULL 

for (i in 1:length(site)) {
  s <- site[i]
  cari_site <- caribbean[caribbean$site35 == site[i],]
  cari_motu_site <- cari_site%>%
    distinct(sequence, .keep_all = TRUE)
  count_families <- data.frame(table(cari_motu_site$family_name_corrected))
  colnames(count_families) <- c("family", "n_motus")
  count_families$n_motus_total <- nrow(cari_motu_site)
  count_families$prop <- count_families$n_motus / count_families$n_motus_total
  count_families <- count_families[order(count_families$prop, decreasing = TRUE),]
  count_families$site <- s
  count_families$province <- "Tropical_Northwestern_Atlantic"
  count_families_site_caribbean <- rbind(count_families_site_caribbean, count_families)
}

write.csv(count_families_site_caribbean, "outputs/04_family_proportion/Site/family_proportion_site_caribbean.csv")

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

ggsave("outputs/04_family_proportion/Site/family_proportion_site_caribbean.png", width=20, height=16)


  ## Caribbean station

station <- c(unique(caribbean$station))

count_families_station_caribbean=NULL 

for (i in 1:length(station)) {
  cari_station <- caribbean[caribbean$station == station[i],]
  st <- station[i]
  s <- unique(cari_station$site35)
  cari_motu_station <- cari_station%>%
    distinct(sequence, .keep_all = TRUE)
  count_families <- data.frame(table(cari_motu_station$family_name_corrected))
  colnames(count_families) <- c("family", "n_motus")
  count_families$n_motus_total <- nrow(cari_motu_station)
  count_families$prop <- count_families$n_motus / count_families$n_motus_total
  count_families <- count_families[order(count_families$prop, decreasing = TRUE),]
  count_families$station <- st
  count_families$site <- s
  count_families$province <- "Tropical_Northwestern_Atlantic"
  count_families_station_caribbean <- rbind(count_families_station_caribbean, count_families)
}

write.csv(count_families_station_caribbean, "outputs/04_family_proportion/Station/family_proportion_station_caribbean.csv")


## Frequency of families in Lengguru
  ## Lengguru total

lengguru <- df_all_filters %>%
  filter(province=="Western_Coral_Triangle") 

leng_motu <- lengguru %>%
  distinct(sequence, .keep_all = TRUE)

count_families_lengguru <- data.frame(table(leng_motu$family_name_corrected))
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

ggsave("outputs/04_family_proportion/Region/family_proportion_lengguru.png", width=6, height=16)
write.csv(count_families_lengguru, "outputs/04_family_proportion/Region/family_proportion_lengguru.csv")

  ## Lengguru sites

site <- c(unique(lengguru$site35))

count_families_site_lengguru=NULL 

for (i in 1:length(site)) {
  s <- site[i]
  leng_site <- lengguru[lengguru$site35 == site[i],]
  leng_motu_site <- leng_site%>%
    distinct(sequence, .keep_all = TRUE)
  count_families <- data.frame(table(leng_motu_site$family_name_corrected))
  colnames(count_families) <- c("family", "n_motus")
  count_families$n_motus_total <- nrow(leng_motu_site)
  count_families$prop <- count_families$n_motus / count_families$n_motus_total
  count_families <- count_families[order(count_families$prop, decreasing = TRUE),]
  count_families$site <- s
  count_families$province <- "Western_Coral_Triangle"
  count_families_site_lengguru <- rbind(count_families_site_lengguru, count_families)
}

write.csv(count_families_site_lengguru, "outputs/04_family_proportion/Site/family_proportion_site_lengguru.csv")

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


ggsave("outputs/04_family_proportion/Site/family_proportion_lengguru_site7-9.png", width=20, height=16)

  ## Lengguru station

station <- c(unique(lengguru$station))

count_families_station_lengguru=NULL 

for (i in 1:length(station)) {
  leng_station <- lengguru[lengguru$station == station[i],]
  st <- station[i]
  s <- unique(leng_station$site35)
  leng_motu_station <- leng_station%>%
    distinct(sequence, .keep_all = TRUE)
  count_families <- data.frame(table(leng_motu_station$family_name_corrected))
  colnames(count_families) <- c("family", "n_motus")
  count_families$n_motus_total <- nrow(leng_motu_station)
  count_families$prop <- count_families$n_motus / count_families$n_motus_total
  count_families <- count_families[order(count_families$prop, decreasing = TRUE),]
  count_families$station <- st
  count_families$site <- s
  count_families$province <- "Western_Coral_Triangle"
  count_families_station_lengguru <- rbind(count_families_station_lengguru, count_families)
}

write.csv(count_families_station_lengguru, "outputs/04_family_proportion/Station/family_proportion_station_lengguru.csv")


## proportion of families in Fakarava
  ## par site = region

fakarava <- df_all_filters %>%
  filter(province=="Southeast_Polynesia")

faka_motu <- fakarava %>%
  distinct(sequence, .keep_all = TRUE)

count_families_site_fakarava <- data.frame(table(faka_motu$family_name_corrected))
colnames(count_families_site_fakarava) <- c("family", "n_motus")
count_families_site_fakarava$n_motus_total <- nrow(faka_motu)
count_families_site_fakarava$site <- "fakarava"
count_families_site_fakarava$province <- "Southeast_Polynesia"
count_families_site_fakarava$prop <- count_families_site_fakarava$n_motus / count_families_site_fakarava$n_motus_total


ggplot(count_families_site_fakarava, aes(x=reorder(family, prop), y = prop)) + 
  geom_bar(stat="identity") + 
  theme_bw() +
  labs(x="Family", y="Proportion")+
  theme(axis.text.x=element_text(angle = 0, hjust = 0)) +
  theme(legend.position = "none")+
  coord_flip()

ggsave("outputs/04_family_proportion/Region/family_proportion_fakarava.png", width=6, height=16)
write.csv(count_families_site_fakarava, "outputs/04_family_proportion/Region/family_proportion_fakarava.csv")

  ## per station
station <- c(unique(fakarava$station))

count_families_station_fakarava=NULL 

for (i in 1:length(station)) {
  faka_station <- fakarava[fakarava$station == station[i],]
  st <- station[i]
  faka_motu_station <- faka_station%>%
    distinct(sequence, .keep_all = TRUE)
  count_families <- data.frame(table(faka_motu_station$family_name_corrected))
  colnames(count_families) <- c("family", "n_motus")
  count_families$n_motus_total <- nrow(faka_motu_station)
  count_families$prop <- count_families$n_motus / count_families$n_motus_total
  count_families <- count_families[order(count_families$prop, decreasing = TRUE),]
  count_families$station <- st
  count_families$site <- "fakarava"
  count_families$province <- "Southeast_Polynesia"
  count_families_station_fakarava <- rbind(count_families_station_fakarava, count_families)
}

write.csv(count_families_station_fakarava, "outputs/04_family_proportion/Station/family_proportion_station_fakarava.csv")


## frequency of families in Eparses
  ## eparse total

eparse <- df_all_filters %>%
  filter(province=="Western_Indian_Ocean") 

eparse_motu <- eparse %>%
  distinct(sequence, .keep_all = TRUE)

count_families_eparse <- data.frame(table(eparse_motu$family_name_corrected))
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

ggsave("outputs/04_family_proportion/Region/family_proportion_eparse.png", width=6, height=16)
write.csv(count_families_eparse, "outputs/04_family_proportion/Region/family_proportion_eparse.csv")

  ## eparse sites

site <- c(unique(eparse$site35))

count_families_site_eparse=NULL 

for (i in 1:length(site)) {
  s <- site[i]
  eparse_site <- eparse[eparse$site35 == site[i],]
  eparse_motu_site <- eparse_site%>%
    distinct(sequence, .keep_all = TRUE)
  count_families <- data.frame(table(eparse_motu_site$family_name_corrected))
  colnames(count_families) <- c("family", "n_motus")
  count_families$n_motus_total <- nrow(eparse_motu_site)
  count_families$prop <- count_families$n_motus / count_families$n_motus_total
  count_families <- count_families[order(count_families$prop, decreasing = TRUE),]
  count_families$site <- s
  count_families$province <- "Western_Indian_Ocean"
  count_families_site_eparse <- rbind(count_families_site_eparse, count_families)
}

write.csv(count_families_site_eparse, "outputs/04_family_proportion/Site/family_proportion_site_eparse.csv")

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


ggsave("outputs/04_family_proportion/Site/family_proportion_site_eparse.png", width=20, height=16)

  ## eparse station

station <- c(unique(eparse$station))

count_families_station_eparse=NULL 

for (i in 1:length(station)) {
  eparse_station <- eparse[eparse$station == station[i],]
  st <- station[i]
  s <- unique(eparse_station$site35)
  eparse_motu_station <- eparse_station%>%
    distinct(sequence, .keep_all = TRUE)
  count_families <- data.frame(table(eparse_motu_station$family_name_corrected))
  colnames(count_families) <- c("family", "n_motus")
  count_families$n_motus_total <- nrow(eparse_motu_station)
  count_families$prop <- count_families$n_motus / count_families$n_motus_total
  count_families <- count_families[order(count_families$prop, decreasing = TRUE),]
  count_families$station <- st
  count_families$site <- s
  count_families$province <- "Western_Indian_Ocean"
  count_families_station_eparse <- rbind(count_families_station_eparse, count_families)
}

write.csv(count_families_station_eparse, "outputs/04_family_proportion/Station/family_proportion_station_eparse.csv")


## Frequency of families in caledonia
## caledonia total

caledonia <- df_all_filters %>%
  filter(province=="Tropical_Southwestern_Pacific") 

cal_motu <- caledonia %>%
  distinct(sequence, .keep_all = TRUE)

count_families_caledonia <- data.frame(table(cal_motu$family_name_corrected))
colnames(count_families_caledonia) <- c("family", "n_motus")
count_families_caledonia$n_motus_total <- nrow(cal_motu)
count_families_caledonia$prop <- count_families_caledonia$n_motus / count_families_caledonia$n_motus_total


ggplot(count_families_caledonia, aes(x=reorder(family, prop), y = prop)) + 
  geom_bar(stat="identity") + 
  theme_bw() +
  labs(x="Family", y="Proportion")+
  theme(axis.text.x=element_text(angle = 0, hjust = 0)) +
  theme(legend.position = "none")+
  coord_flip()

ggsave("outputs/04_family_proportion/Region/family_proportion_caledonia.png", width=6, height=16)
write.csv(count_families_caledonia, "outputs/04_family_proportion/Region/family_proportion_caledonia.csv")

## caledonia sites

site <- c(unique(caledonia$site35))

count_families_site_caledonia=NULL 

for (i in 1:length(site)) {
  s <- site[i]
  cal_site <- caledonia[caledonia$site35 == site[i],]
  cal_motu_site <- cal_site%>%
    distinct(sequence, .keep_all = TRUE)
  count_families <- data.frame(table(cal_motu_site$family_name_corrected))
  colnames(count_families) <- c("family", "n_motus")
  count_families$n_motus_total <- nrow(cal_motu_site)
  count_families$prop <- count_families$n_motus / count_families$n_motus_total
  count_families <- count_families[order(count_families$prop, decreasing = TRUE),]
  count_families$site <- s
  count_families$province <- "Tropical_Southwestern_Pacific"
  count_families_site_caledonia <- rbind(count_families_site_caledonia, count_families)
}

write.csv(count_families_site_caledonia, "outputs/04_family_proportion/Site/family_proportion_site_caledonia.csv")

count_families_site_caledonia <- count_families_site_caledonia %>%
  ungroup() %>%
  arrange(site, prop) %>%
  mutate(order = row_number())

## plot by 3 sites, because plot too small otherwise
site_sub <- c("STVINCENT", "BOURAIL", "NOUMEA")
subset1 <- count_families_site_caledonia %>%
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


ggsave("outputs/04_family_proportion/Site/family_proportion_caledonia_site4-6.png", width=20, height=16)

## caledonia station

station <- c(unique(caledonia$station))

count_families_station_caledonia=NULL 

for (i in 1:length(station)) {
  cal_station <- caledonia[caledonia$station == station[i],]
  st <- station[i]
  s <- unique(cal_station$site35)
  cal_motu_station <- cal_station%>%
    distinct(sequence, .keep_all = TRUE)
  count_families <- data.frame(table(cal_motu_station$family_name_corrected))
  colnames(count_families) <- c("family", "n_motus")
  count_families$n_motus_total <- nrow(cal_motu_station)
  count_families$prop <- count_families$n_motus / count_families$n_motus_total
  count_families <- count_families[order(count_families$prop, decreasing = TRUE),]
  count_families$station <- st
  count_families$site <- s
  count_families$province <- "Tropical_Southwestern_Pacific"
  count_families_station_caledonia <- rbind(count_families_station_caledonia, count_families)
}

write.csv(count_families_station_caledonia, "outputs/04_family_proportion/Station/family_proportion_station_caledonia.csv")



## Proportion of families global


global_motu <- df_all_filters %>%
  distinct(sequence, .keep_all = TRUE)


count_families_global <- data.frame(table(global_motu$family_name_corrected))
colnames(count_families_global) <- c("family", "n_motus")
save(count_families_global, file = "Rdata/nb_motus_per_family_global.Rdata")

count_families_global$n_motus_total <- nrow(global_motu)
count_families_global$prop <- count_families_global$n_motus / count_families_global$n_motus_total

write.csv(count_families_global, "outputs/04_family_proportion/family_proportion_global.csv")

ggplot(count_families_global, aes(x=reorder(family, prop), y = prop)) + 
  geom_bar(stat="identity") + 
  theme_bw() +
  labs(x="Family", y="Proportion")+
  theme(legend.position = "none")+
  theme(axis.text.x=element_text(angle = 0, hjust = 0)) + 
  coord_flip()

ggsave("outputs/04_family_proportion/family_proportion_global.png", width=6, height=16)



## plot proportion of each family in each region on global scale
families_prop_global <- left_join(count_families_global, count_families_site_fakarava[,1:2], by="family" )
families_prop_global <- left_join(families_prop_global, count_families_caribbean[,1:2], by="family" )
families_prop_global <- left_join(families_prop_global, count_families_eparse[,1:2], by="family" )
families_prop_global <- left_join(families_prop_global, count_families_caledonia[,1:2], by="family" )
families_prop_global <- left_join(families_prop_global, count_families_lengguru[,1:2], by="family" )

families_prop_global <- families_prop_global[, c(-3)]
colnames(families_prop_global) <- c("family", "n_global", "prop_global", "n_fakarava", "n_caribbean", "n_eparse", "n_caledonia", "n_lengguru")
families_prop_global[is.na(families_prop_global)] <- 0

for (i in 1:dim(families_prop_global)[1]) {
  families_prop_global[i,"new_n_faka"] <- (families_prop_global[i,"n_fakarava"]*families_prop_global[i,"n_global"])/sum(families_prop_global[i,4:8])
  families_prop_global[i,"new_n_car"] <- (families_prop_global[i,"n_caribbean"]*families_prop_global[i,"n_global"])/sum(families_prop_global[i,4:8])
  families_prop_global[i,"new_n_eparse"] <- (families_prop_global[i,"n_eparse"]*families_prop_global[i,"n_global"])/sum(families_prop_global[i,4:8])
  families_prop_global[i,"new_n_caledonia"] <- (families_prop_global[i,"n_caledonia"]*families_prop_global[i,"n_global"])/sum(families_prop_global[i,4:8])
  families_prop_global[i,"new_n_leng"] <- (families_prop_global[i,"n_lengguru"]*families_prop_global[i,"n_global"])/sum(families_prop_global[i,4:8])
    
}

families_prop_global$Fakarava <- (families_prop_global$new_n_faka*families_prop_global$prop_global)/families_prop_global$n_global
families_prop_global$Caribbean <- (families_prop_global$new_n_car*families_prop_global$prop_global)/families_prop_global$n_global
families_prop_global$Eparse <- (families_prop_global$new_n_eparse*families_prop_global$prop_global)/families_prop_global$n_global
families_prop_global$Caledonia <- (families_prop_global$new_n_caledonia*families_prop_global$prop_global)/families_prop_global$n_global
families_prop_global$Lengguru <- (families_prop_global$new_n_leng*families_prop_global$prop_global)/families_prop_global$n_global

families_prop_global <- families_prop_global[,c(1,14:18)]
families_prop_global2 <- reshape2::melt(families_prop_global)
colnames(families_prop_global2) <- c("family", "Province", "prop")

save(families_prop_global, file = "Rdata/family_proportion_global.Rdata")

ggplot(families_prop_global2, aes(x=reorder(family, prop), y = prop, fill = Province)) + 
  geom_bar(stat="identity", show.legend = TRUE) + 
  theme_bw() +
  scale_fill_manual(values =c("#a6611a", "#E5A729", "#015462", "#b2182b", "#80cdc1"))+ 
  labs(x="Family", y="Proportion")+
  coord_flip()

ggsave("outputs/04_family_proportion/family_proportion_global_region.png", width=6, height=16)




## Bellwood figures : proportion of families per site

df_all_site <- rbind(count_families_site_lengguru[,-7], count_families_site_caribbean[,-7], count_families_site_eparse[,-7], count_families_site_fakarava, count_families_site_caledonia[,-7])
save(df_all_site, file = "Rdata/family_proportion_per_site.rdata")

family <- c("Acanthuridae", "Chaetodontidae", "Labridae", "Lutjanidae", "Serranidae", "Carangidae", "Pomacentridae", "Apogonidae", "Gobiidae")

prop <- vector("list")
for (i in 1:length(family)) {
  fam <- df_all_site[df_all_site$family == family[i],]
  prop[[i]] <- ggplot(fam, aes(n_motus_total, prop, ymin=0, ymax=0.2, colour=province))+
    geom_point(size=2)+
    scale_y_continuous(breaks = c(0, 0.1, 0.2))+
    xlim(0, 800)+
    theme(legend.position = "none")+
    theme_bw()+
    scale_color_manual(values =c("#E5A729", "#80cdc1", "#a6611a", "#b2182b", "#015462"))+ 
    labs(title=family[i], x="", y="")+
    theme(plot.title = element_text(size = 10, face="bold"), plot.margin=unit(c(0,0.1,0,0), "cm"))
}


plot <- ggarrange(plotlist = prop, ncol=3, nrow = 3, common.legend = TRUE, legend = "top")
x.grob <- textGrob("Total number of MOTUs per site", 
                   gp=gpar(fontface="bold", col="black", fontsize=12))
y.grob <- textGrob("Proportion of MOTUs assigned to the family in each site", 
                   gp=gpar(fontface="bold", col="black", fontsize=12), rot = 90)
plot_grid <- grid.arrange(plot, bottom=x.grob, left=y.grob)





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

