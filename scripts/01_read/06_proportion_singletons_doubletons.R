library(tidyverse)
library(reshape2)

setwd("c:/Users/mathon/Desktop/linux/Global_eDNA/")
load("Rdata/02_clean_all.Rdata")

df_all_filters <- subset(df_all_filters, !(station %in% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")))
df_all_filters <- subset(df_all_filters, sample_method!="niskin")
df_all_filters <- subset(df_all_filters, region!="East_Pacific")


gamma_global <- as.numeric(df_all_filters %>%
                             summarise(n = n_distinct(sequence)))

# Calculate proportion of singletons in samples_all_pcr 


sample <- unique(df_all_filters$sample_name_all_pcr)

df_sample=data.frame(motu=character())

for (i in 1:length(sample)) {
  df <- df_all_filters[df_all_filters$sample_name_all_pcr == sample[i],] %>%
    distinct(sequence, sample_name_all_pcr)
  colnames(df) <- c("motu", sample[i])
  df_sample <- full_join(df_sample, df, by="motu")
}
rownames(df_sample) <- df_sample[,1]
df_sample <- decostand(df_sample[,2:173], "pa",na.rm = TRUE)
df_sample[is.na(df_sample)] <- 0


singleton_sample <- nrow(subset(df_sample, rowSums(df_sample[,1:172]) == 1))
doubleton_sample <- nrow(subset(df_sample, rowSums(df_sample[,1:172]) == 2))

perc_sing_sample <- (singleton_sample*100)/gamma_global
perc_doub_sample <- (doubleton_sample*100)/gamma_global


# Calculate proportion of singletons in replicates (sample_name) 


replicate <- unique(df_all_filters$sample_name)

df_replicate=data.frame(motu=character())

for (i in 1:length(replicate)) {
  df <- df_all_filters[df_all_filters$sample_name == replicate[i],] %>%
    distinct(sequence, sample_name)
  colnames(df) <- c("motu", replicate[i])
  df_replicate <- full_join(df_replicate, df, by="motu")
}
rownames(df_replicate) <- df_replicate[,1]
df_replicate<- decostand(df_replicate[,2:1820], "pa",na.rm = TRUE)
df_replicate[is.na(df_replicate)] <- 0


singleton_replicate <- nrow(subset(df_replicate, rowSums(df_replicate[,1:1819]) == 1))
doubleton_replicate <- nrow(subset(df_replicate, rowSums(df_replicate[,1:1819]) == 2))

perc_sing_replicate <- (singleton_replicate*100)/gamma_global
perc_doub_replicate <- (doubleton_replicate*100)/gamma_global
