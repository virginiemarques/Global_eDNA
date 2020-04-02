library(tidyverse)
library(reshape2)

setwd("c:/Users/mathon/Desktop/linux/Global_eDNA/")
load("Rdata/02_clean_all.Rdata")

df_all_filters <- subset(df_all_filters, !(station %in% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")))


# Proportion is calculated as (nb of unique motus per family)/(total nb of unique motus per site/region), using only amplicons assigned to family level minimum

## frequency of families in Caribbean

  ## Caribbean total
caribbean <- df_all_filters %>%
  filter(region=="Caribbean" & !is.na(new_family_name))

cari_motu <- caribbean %>%
  distinct(sequence, .keep_all = TRUE)

count_families_caribbean <- data.frame(table(cari_motu$new_family_name))
colnames(count_families_caribbean) <- c("family", "n_motus")
count_families_caribbean$n_motus_total <- nrow(cari_motu)
count_families_caribbean$prop <- count_families_caribbean$n_motus / count_families_caribbean$n_motus_total

ggplot(count_families_caribbean, aes(x=reorder(family, prop), y = prop, fill = prop)) + 
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
  count_families_site_caribbean <- rbind(count_families_site_caribbean, count_families)
}

write.csv(count_families_site_caribbean, "outputs/05_family_proportion/02_based_on_species_presence/per site/family_proportion_site_caribbean.csv")

count_families_site_caribbean <- count_families_site_caribbean %>%
  ungroup() %>%
  arrange(site, prop) %>%
  mutate(order = row_number())

ggplot(count_families_site_caribbean, aes(order, prop, fill = prop)) + 
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
  count_families_station_caribbean <- rbind(count_families_station_caribbean, count_families)
}

write.csv(count_families_station_caribbean, "outputs/05_family_proportion/02_based_on_species_presence/per station/family_proportion_station_caribbean.csv")


## Frequency of families in Lengguru
  ## Lengguru total

lengguru <- df_all_filters %>%
  filter(region=="West_Papua" & !is.na(new_family_name))

leng_motu <- lengguru %>%
  distinct(sequence, .keep_all = TRUE)

count_families_lengguru <- data.frame(table(leng_motu$new_family_name))
colnames(count_families_lengguru) <- c("family", "n_motus")
count_families_lengguru$n_motus_total <- nrow(leng_motu)
count_families_lengguru$prop <- count_families_lengguru$n_motus / count_families_lengguru$n_motus_total


ggplot(count_families_lengguru, aes(x=reorder(family, prop), y = prop, fill = prop)) + 
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
ggplot(subset1, aes(order, prop, fill = prop)) + 
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
  count_families_station_lengguru <- rbind(count_families_station_lengguru, count_families)
}

write.csv(count_families_station_lengguru, "outputs/05_family_proportion/02_based_on_species_presence/per station/family_proportion_station_lengguru.csv")



## frequency of families in Mediterranean
## frequency of families in Eparses


## Bellwood figures : proportion of families per site

df_all_site <- rbind(count_families_site_caribbean[,c(-6)], count_families_site_lengguru[,c(-6)])


family <- c("Serranidae", "Lutjanidae", "Acanthuridae", "Pomacentridae", "Balistidae", "Lethrinidae", "Scombridae", "Exocoetidae", "Myctophidae", "Apogonidae", "Carangidae", "Dasyatidae", "Haemulidae", "Clupeidae", "Gobiidae" )
prop <- vector("list")
for (i in 1:length(family)) {
  fam <- df_all_site[df_all_site$family == family[i],]
  prop[[i]] <- ggplot(fam, aes(n_motus_total, prop, ymin=0, ymax=1))+
    geom_point()+
    labs(title=family[i], x="total number of motus per site", y="Proportion")
  ggsave(filename=paste("outputs/05_family_proportion/02_based_on_species_presence/per site/prop_",family[i],".png", sep = ""))
}


## Bellwood figures : proportion of families per station

df_all_station <- rbind(count_families_station_caribbean, count_families_station_lengguru)


family <- c("Serranidae", "Lutjanidae", "Acanthuridae", "Pomacentridae", "Balistidae", "Lethrinidae", "Scombridae", "Exocoetidae", "Myctophidae", "Apogonidae", "Carangidae", "Dasyatidae", "Haemulidae", "Clupeidae", "Gobiidae" )
prop <- vector("list")
for (i in 1:length(family)) {
  fam <- df_all_station[df_all_station$family == family[i],]
  prop[[i]] <- ggplot(fam, aes(n_motus_total, prop, ymin=0, ymax=1, colour=site))+
    geom_point()+
    labs(title=family[i], x="total number of motus per station", y="Proportion")+
    theme(legend.text = element_text(size = 6))+
    scale_y_continuous(breaks=c(0, 0.2, 0.6, 1))
  ggsave(filename=paste("outputs/05_family_proportion/02_based_on_species_presence/per station/prop_",family[i],".png", sep = ""), width=5, height=4)
}
