library(tidyverse)
library(reshape2)

setwd("c:/Users/mathon/Desktop/linux/Global_eDNA/")
load("Rdata/02_clean_all.Rdata")

df_all_filters <- subset(df_all_filters, !(station %in% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")))
df_all_filters <- subset(df_all_filters, sample_method!="niskin")
df_all_filters <- subset(df_all_filters, region!="East_Pacific")



# Proportion is calculated as (nb of amplicons per family)/(total nb of amplicons per site/region), using only amplicons assigned to family level minimum

## frequency of families in Caribbean

  ## Caribbean total
caribbean <- df_all_filters %>%
  filter(region=="Caribbean" & !is.na(new_family_name))

count_families_caribbean <- data.frame(table(caribbean$new_family_name))
colnames(count_families_caribbean) <- c("family", "n_motus")
count_families_caribbean$n_motus_total <- nrow(caribbean)
count_families_caribbean$prop <- count_families_caribbean$n_motus / count_families_caribbean$n_motus_total

ggplot(count_families_caribbean, aes(x=reorder(family, prop), y = prop, fill = prop)) + 
  geom_bar(stat="identity") + 
  theme_bw() +
  labs(x="Family", y="Frequency")+
  theme(axis.text.x=element_text(angle = 0, hjust = 0)) +
  theme(legend.position = "none")+
  coord_flip()

ggsave("outputs/05_family_proportion/01_based_on_species_occurences/family_frequency_caribbean.png", width=6, height=16)
write.csv(count_families_caribbean, "outputs/05_family_proportion/01_based_on_species_occurences/family_frequency_caribbean.csv")

  ## Caribbean sites
site <- c(unique(caribbean$site))

count_families_site_caribbean=NULL 

for (i in 1:length(site)) {
  s <- site[i]
  cari_site <- caribbean[caribbean$site == site[i],]
  count_families <- data.frame(table(cari_site$new_family_name))
  colnames(count_families) <- c("family", "n_motus")
  count_families$n_motus_total <- nrow(cari_site)
  count_families$prop <- count_families$n_motus / count_families$n_motus_total
  count_families <- count_families[order(count_families$prop, decreasing = TRUE),]
  count_families$site <- s
  count_families_site_caribbean <- rbind(count_families_site_caribbean, count_families)
}

write.csv(count_families_site_caribbean, "outputs/05_family_proportion/01_based_on_species_occurences/family_frequency_site_caribbean.csv")

count_families_site_caribbean <- count_families_site_caribbean %>%
  ungroup() %>%
  arrange(site, prop) %>%
  mutate(order = row_number())

ggplot(count_families_site_caribbean, aes(order, prop, fill = prop)) + 
  geom_bar(stat="identity") +
  facet_wrap(~site, scales ="free_y", as.table = FALSE)+
  theme_bw() +
  labs(x="Family", y="Frequency")+
  theme(axis.text.x=element_text(angle = 0, hjust = 0)) +
  theme(legend.position = "none")+
  coord_flip()+
  scale_x_continuous(breaks=count_families_site_caribbean$order, labels=count_families_site_caribbean$family, expand = c(0,0))

ggsave("outputs/05_family_proportion/01_based_on_species_occurences/family_frequency_site_caribbean.png", width=20, height=16)


## Frequency of families in Lengguru
  ## Lengguru total

lengguru <- df_all_filters %>%
  filter(region=="West_Papua" & !is.na(new_family_name))

count_families_lengguru <- data.frame(table(lengguru$new_family_name))
colnames(count_families_lengguru) <- c("family", "n_motus")
count_families_lengguru$n_motus_total <- nrow(lengguru)
count_families_lengguru$prop <- count_families_lengguru$n_motus / count_families_lengguru$n_motus_total


ggplot(count_families_lengguru, aes(x=reorder(family, prop), y = prop, fill = prop)) + 
  geom_bar(stat="identity") + 
  theme_bw() +
  labs(x="Family", y="Frequency")+
  theme(axis.text.x=element_text(angle = 0, hjust = 0)) +
  theme(legend.position = "none")+
  coord_flip()

ggsave("outputs/05_family_proportion/01_based_on_species_occurences/family_frequency_lengguru.png", width=6, height=16)
write.csv(count_families_lengguru, "outputs/05_family_proportion/01_based_on_species_occurences/family_frequency_lengguru.csv")

  ## Lengguru sites

site <- c(unique(lengguru$site))

count_families_site_lengguru=NULL 

for (i in 1:length(site)) {
  s <- site[i]
  leng_site <- lengguru[lengguru$site == site[i],]
  count_families <- data.frame(table(leng_site$new_family_name))
  colnames(count_families) <- c("family", "n_motus")
  count_families$n_motus_total <- nrow(leng_site)
  count_families$prop <- count_families$n_motus / count_families$n_motus_total
  count_families <- count_families[order(count_families$prop, decreasing = TRUE),]
  count_families$site <- s
  count_families_site_lengguru <- rbind(count_families_site_lengguru, count_families)
}

write.csv(count_families_site_lengguru, "outputs/05_family_proportion/01_based_on_species_occurences/family_frequency_site_lengguru.csv")

count_families_site_lengguru <- count_families_site_lengguru %>%
  ungroup() %>%
  arrange(site, prop) %>%
  mutate(order = row_number())

  ## plot by 3 sites, because plot too small otherwise
site_sub <- c("pulau_aiduma", "pulau_aiduma_ext")
subset1 <- count_families_site_lengguru %>%
  filter(site%in%site_sub)
ggplot(subset1, aes(order, prop, fill = prop)) + 
  geom_bar(stat="identity") +
  facet_wrap(~site, scales ="free_y")+
  theme_bw() +
  labs(x="Family", y="Frequency")+
  theme(axis.text.x=element_text(angle = 0, hjust = 0)) +
  theme(legend.position = "none")+
  coord_flip()+
  scale_x_continuous(breaks=subset1$order, labels=subset1$family, expand = c(0,0))


ggsave("outputs/05_family_proportion/01_based_on_species_occurences/family_frequency_lengguru_site10-11.png", width=20, height=16)



## frequency of families in Eparses

  ## eparse total

eparse <- df_all_filters %>%
  filter(region=="" & !is.na(new_family_name))

count_families_eparse <- data.frame(table(eparse$new_family_name))
colnames(count_families_eparse) <- c("family", "n_motus")
count_families_eparse$n_motus_total <- nrow(eparse)
count_families_eparse$prop <- count_families_eparse$n_motus / count_families_eparse$n_motus_total


ggplot(count_families_eparse, aes(x=reorder(family, prop), y = prop, fill = prop)) + 
  geom_bar(stat="identity") + 
  theme_bw() +
  labs(x="Family", y="Frequency")+
  theme(axis.text.x=element_text(angle = 0, hjust = 0)) +
  theme(legend.position = "none")+
  coord_flip()

ggsave("outputs/05_family_proportion/01_based_on_species_occurences/family_frequency_eparse.png", width=6, height=16)
write.csv(count_families_eparse, "outputs/05_family_proportion/01_based_on_species_occurences/family_frequency_eparse.csv")

  ## eparse sites

site <- c(unique(eparse$site))

count_families_site_eparse=NULL 

for (i in 1:length(site)) {
  s <- site[i]
  eparse_site <- eparse[eparse$site == site[i],]
  count_families <- data.frame(table(eparse_site$new_family_name))
  colnames(count_families) <- c("family", "n_motus")
  count_families$n_motus_total <- nrow(eparse_site)
  count_families$prop <- count_families$n_motus / count_families$n_motus_total
  count_families <- count_families[order(count_families$prop, decreasing = TRUE),]
  count_families$site <- s
  count_families_site_eparse <- rbind(count_families_site_eparse, count_families)
}

write.csv(count_families_site_eparse, "outputs/05_family_proportion/01_based_on_species_occurences/family_frequency_site_eparse.csv")

count_families_site_eparse <- count_families_site_eparse %>%
  ungroup() %>%
  arrange(site, prop) %>%
  mutate(order = row_number())

## plot by 3 sites, because plot too small otherwise
site_sub <- c("", "")
subset1 <- count_families_site_eparse %>%
  filter(site%in%site_sub)
ggplot(subset1, aes(order, prop, fill = prop)) + 
  geom_bar(stat="identity") +
  facet_wrap(~site, scales ="free_y")+
  theme_bw() +
  labs(x="Family", y="Frequency")+
  theme(axis.text.x=element_text(angle = 0, hjust = 0)) +
  theme(legend.position = "none")+
  coord_flip()+
  scale_x_continuous(breaks=subset1$order, labels=subset1$family, expand = c(0,0))


ggsave("outputs/05_family_proportion/01_based_on_species_occurences/family_frequency_eparse_site...png", width=20, height=16)

