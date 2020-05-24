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


## nb of unique Motus
df_all_filters %>%
  summarise(n = n_distinct(sequence))


## unassigned reads
df <- df_all_filters %>%
  subset(new_rank_ncbi == "higher")
sum(df$count_reads)


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
