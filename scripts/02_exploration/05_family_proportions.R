library(tidyverse)
library(reshape2)

setwd("c:/Users/mathon/Desktop/linux/Global_eDNA/")
load("Rdata/02_clean_all.Rdata")


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

ggsave("outputs/05_family_proportion/family_frequency_caribbean.png", width=6, height=16)
write.csv(count_families_caribbean, "outputs/05_family_proportion/family_frequency_caribbean.csv")

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

write.csv(count_families_site_caribbean, "outputs/05_family_proportion/family_frequency_site_caribbean.csv")

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

ggsave("outputs/05_family_proportion/family_frequency_site_caribbean.png", width=20, height=16)


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

ggsave("outputs/05_family_proportion/family_frequency_lengguru.png", width=6, height=16)
write.csv(count_families_lengguru, "outputs/05_family_proportion/family_frequency_lengguru.csv")

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

write.csv(count_families_site_lengguru, "outputs/05_family_proportion/family_frequency_site_lengguru.csv")

count_families_site_lengguru <- count_families_site_lengguru %>%
  ungroup() %>%
  arrange(site, prop) %>%
  mutate(order = row_number())

  ## plot by 3 sites, because plot too small otherwise
site_sub <- c("lobo", "pulau_aiduma", "pulau_aiduma_ext")
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


ggsave("outputs/05_family_proportion/family_frequency_lengguru_site7-9.png", width=20, height=16)


## frequency of families in Mediterranean
## frequency of families in Eparses


## Bellwood figures : proportion of families per site

df_all_site <- rbind(count_families_site_caribbean[,c(-6)], count_families_site_lengguru[,c(-6)])


family <- c("Serranidae", "Lutjanidae", "Acanthuridae", "Pomacentridae", "Balistidae", "Lethrinidae", "Scombridae", "Exocoetidae", "Myctophidae", "Apogonidae", "Carangidae", "Dasyatidae", "Haemulidae", "Clupeidae", "Gobiidae" )
prop <- vector("list")
for (i in 1:length(family)) {
  fam <- df_all_site[df_all_site$family == family[i],]
  prop[[i]] <- ggplot(fam, aes(n_motus_total, prop))+
    geom_point()+
    labs(title=family[i], x="total number of amplicons per site", y="Proportion")
  ggsave(filename=paste("outputs/05_family_proportion/prop_",family[i],".png", sep = ""))
}
