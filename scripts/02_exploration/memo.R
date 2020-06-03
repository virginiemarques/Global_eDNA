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
  summarise(n = n_distinct(sequence))

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


## nb motus Labridae, gobiidae et pomacentridae
df_all_filters %>%
  filter(new_family_name == "Poridae")%>%
  summarise(n = n_distinct(sequence))

family_coverage %>%
  filter(coef_resolution==1)%>%
  summarise(n=n_distinct(Family))


library(reshape2)
# count RLS
RLS_species <- read.csv("data/RLS/RLS_species.csv", sep = ";", stringsAsFactors = FALSE, check.names = FALSE)
RLS_species <- RLS_species %>%
  filter(realm%ni%c("Temperate Australasia", "Temperate Northern Atlantic", "Temperate Southern Africa", "Temperate Northern Pacific"))
RLS_sp <- RLS_species[,c(12:2051)]
RLS_sp <- RLS_sp[,colSums(RLS_sp)>0]
coral_fishes <- as.data.frame(colnames(RLS_sp))
colnames(coral_fishes) <- "Species"
coral_fishes2 <- colsplit(coral_fishes$Species, " ", c("Genus", "Species"))
coral_fishes$Genus <- coral_fishes2$Genus


library(rfishbase)
all_fishbase <- load_taxa()
coral_fishes <- left_join(coral_fishes, all_fishbase[,c("Species", "Family")], by="Species")


n_distinct(coral_fishes$Species)
n_distinct(coral_fishes$Family)
n_distinct(coral_fishes$Genus)

write.csv(coral_fishes, file = "data/RLS/coral_fishes2.csv")
