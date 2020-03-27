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

  ## Caribbean sites (work in progress)
site <- c(unique(caribbean$site))

count_families_site_caribbean <- data.frame(site=character(), family=character(), count=numeric(), stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  s <- site[i]
  otu_names_family <- caribbean[caribbean$site == site[i],] %>%
    distinct(sequence, new_rank_ncbi, new_family_name, new_scientific_name_ncbi, best_identity_database)
  count_families <- table(otu_names_family$new_family_name)
  
 
}


## Frequency of families in Lengguru
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


  