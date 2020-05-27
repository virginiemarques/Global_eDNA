load("Rdata/02_clean_all.Rdata")

'%ni%' <- Negate("%in%")

df_all_filters <- df_all_filters %>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")) %>%
  filter(sample_method !="niskin" & region!="East_Pacific" & comment %ni% c("Distance decay 600m", "Distance decay 300m") & station!="glorieuse_distance_300m")%>%
  filter(project != "SEAMOUNTS") %>% 
  filter(habitat_type %ni% c("BAIE", "Sommet"))



## total number of reads after cleaning
sum(df_all_filters$count_reads)

## mean read count per sample
count_reads <- data.frame()
s <- unique(df_all_filters$sample_name_all_pcr)
for (i in 1:length(s)) {
  df <- df_all_filters %>%
    filter(sample_name_all_pcr == s[i]) %>%
    distinct(amplicon, .keep_all=TRUE)
  count_reads[i,1] <- sum(df$count_reads_all_pcr)
}
mean(count_reads$V1)

# number of samples, replicates and station
df_all_filters %>% summarise(n=n_distinct(sample_name_all_pcr))
df_all_filters %>% summarise(n=n_distinct(sample_name))
df_all_filters %>% summarise(n=n_distinct(station))

## nb of unique Motus
df_all_filters %>%
  summarise(n = n_distinct(sequence))


## number of species
df_all_filters %>%
  filter(!is.na(new_species_name))%>%
  summarise(n = n_distinct(new_species_name))

## unassigned reads
df <- df_all_filters %>%
  subset(new_rank_ncbi == "higher")
sum(df$count_reads)

## unassigned motus
df_all_filters %>%
  filter(new_rank_ncbi == "higher")%>%
  summarise(n = n_distinct(amplicon)) %>% pull()

## assigned motus
df_all_filters %>%
  filter(new_rank_ncbi != "higher")%>%
  summarise(n = n_distinct(amplicon)) %>% pull()



library(reshape2)
# count RLS
coral_fishes <- read.csv("data/coral_fishes.csv", sep = ";", header = FALSE)
colnames(coral_fishes) <- "Species"
coral_fishes2 <- colsplit(coral_fishes$Species, " ", c("Genus", "Species"))
coral_fishes$Genus <- coral_fishes2$Genus


library(rfishbase)
all_fishbase <- load_taxa()
coral_fishes <- left_join(coral_fishes, all_fishbase[,c("Species", "Family")], by="Species")


n_distinct(coral_fishes$Species)
n_distinct(coral_fishes$Family)
n_distinct(coral_fishes$Genus)

write.csv(coral_fishes, file = "data/coral_fishes2.csv")
