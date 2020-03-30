library(tidyverse)
library(reshape2)

setwd("c:/Users/mathon/Desktop/linux/Global_eDNA/")
load("Rdata/02_clean_all.Rdata")


## frequency of families in Caribbean

  ## Caribbean total
caribbean <- df_all_filters %>%
  filter(region=="Caribbean")

otu_names_family <- caribbean %>%
  distinct(sequence, new_rank_ncbi, new_family_name, new_scientific_name_ncbi, best_identity_database)

count_families_caribbean <- data.frame(table(otu_names_family$new_family_name))

ggplot(count_families_caribbean, aes(x=reorder(Var1, Freq), y = Freq, fill = Freq)) + 
  geom_bar(stat="identity") + 
  theme_bw() +
  labs(x="Family", y="Frequency")+
  theme(axis.text.x=element_text(angle = 0, hjust = 0)) +
  theme(legend.position = "none")+
  coord_flip()

ggsave("outputs/05_family_proportion/family_frequency_caribbean.png", width=6, height=16)
write.csv(count_families_caribbean, "outputs/05_family_proportion/family_frequency_caribbean.csv")

  ## Caribbean sites (je sais pas comment trier les fréquences par site dans le graph)
site <- c(unique(caribbean$site))

count_families_site_caribbean=NULL 

for (i in 1:length(site)) {
  s <- site[i]
  otu_names_family <- caribbean[caribbean$site == site[i],] %>%
    distinct(sequence, new_rank_ncbi, new_family_name, new_scientific_name_ncbi, best_identity_database)
  count_families <- data.frame(table(otu_names_family$new_family_name))
  count_families <- count_families[order(count_families$Freq, decreasing = TRUE),]
  count_families$site <- s
  count_families_site_caribbean <- rbind(count_families_site_caribbean, count_families)
}

write.csv(count_families_site_caribbean, "outputs/05_family_proportion/family_frequency_site_caribbean.csv")

count_families_site_caribbean <- count_families_site_caribbean %>%
  ungroup() %>%
  arrange(site, Freq) %>%
  mutate(order = row_number())

ggplot(count_families_site_caribbean, aes(order, Freq, fill = Freq)) + 
  geom_bar(stat="identity") +
  facet_wrap(~site, scales ="free_y", as.table = FALSE)+
  theme_bw() +
  labs(x="Family", y="Frequency")+
  theme(axis.text.x=element_text(angle = 0, hjust = 0)) +
  theme(legend.position = "none")+
  coord_flip()+
  scale_x_continuous(breaks=count_families_site_caribbean$order, labels=count_families_site_caribbean$Var1, expand = c(0,0))

ggsave("outputs/05_family_proportion/family_frequency_site_caribbean.png", width=20, height=16)


## Frequency of families in Lengguru
  ## Lengguru total

lengguru <- df_all_filters %>%
  filter(region=="West_Papua")

otu_names_family <- lengguru %>%
  distinct(sequence, new_rank_ncbi, new_family_name, new_scientific_name_ncbi, best_identity_database)


count_families_lengguru <- data.frame(table(otu_names_family$new_family_name))

ggplot(count_families_lengguru, aes(x=reorder(Var1, Freq), y = Freq, fill = Freq)) + 
  geom_bar(stat="identity") + 
  theme_bw() +
  labs(x="Family", y="Frequency")+
  theme(axis.text.x=element_text(angle = 0, hjust = 0)) +
  theme(legend.position = "none")+
  coord_flip()

ggsave("outputs/05_family_proportion/family_frequency_lengguru.png", width=6, height=16)
write.csv(count_families_lengguru, "outputs/05_family_proportion/family_frequency_lengguru.csv")

  ## Lengguru sites (je sais pas comment trier les fréquences par site, dans le graph)

site <- c(unique(lengguru$site))

count_families_site_lengguru=NULL 

for (i in 1:length(site)) {
  s <- site[i]
  otu_names_family <- lengguru[lengguru$site == site[i],] %>%
    distinct(sequence, new_rank_ncbi, new_family_name, new_scientific_name_ncbi, best_identity_database)
  count_families <- data.frame(table(otu_names_family$new_family_name))
  count_families <- count_families[order(count_families$Freq, decreasing = TRUE),]
  count_families$site <- s
  count_families_site_lengguru <- rbind(count_families_site_lengguru, count_families)
}

write.csv(count_families_site_lengguru, "outputs/05_family_proportion/family_frequency_site_lengguru.csv")

count_families_site_lengguru <- count_families_site_lengguru %>%
  ungroup() %>%
  arrange(site, Freq) %>%
  mutate(order = row_number())

  ## plot by 3 sites, because too big otherwise
site_sub <- c("tanjung_boi", "triton_bay")
subset1 <- count_families_site_lengguru %>%
  filter(site%in%site_sub)
ggplot(subset1, aes(order, Freq, fill = Freq)) + 
  geom_bar(stat="identity") +
  facet_wrap(~site, scales ="free_y")+
  theme_bw() +
  labs(x="Family", y="Frequency")+
  theme(axis.text.x=element_text(angle = 0, hjust = 0)) +
  theme(legend.position = "none")+
  coord_flip()+
  scale_x_continuous(breaks=subset1$order, labels=subset1$Var1, expand = c(0,0))


ggsave("outputs/05_family_proportion/family_frequency_lengguru_site10-11.png", width=20, height=16)
