library(tidyverse)

load("Rdata/02_clean_all.Rdata")

#Remove estuary stations and deep niskin station
df_all_filters <- subset(df_all_filters, !(station %in% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")))
df_all_filters <- subset(df_all_filters, sample_method!="niskin")
df_all_filters <- subset(df_all_filters, region!="East_Pacific")

df_all_filters <- df_all_filters %>%
  filter(!is.na(new_family_name))


family <- unique(df_all_filters$new_family_name)

prop_similarity <- data.frame(stringsAsFactors = FALSE)

for (i in 1:length(family)) {
  df <- subset(df_all_filters, new_family_name==family[i])
  ntot <- nrow(df)
  prop_similarity[i,"family"] <- family[i]
  prop_similarity[i,"p80_85"] <- (nrow(subset(df, best_identity_database > 0.80 & best_identity_database <= 0.85))*100)/ntot
  prop_similarity[i,"p85_90"] <- (nrow(subset(df, best_identity_database > 0.85 & best_identity_database <= 0.90))*100)/ntot
  prop_similarity[i,"p90_95"] <- (nrow(subset(df, best_identity_database > 0.90 & best_identity_database <= 0.95))*100)/ntot
  prop_similarity[i,"p95_100"] <- (nrow(subset(df, best_identity_database > 0.95 & best_identity_database <= 1))*100)/ntot
}

save(prop_similarity, file = "Rdata/proportion_motus_similarity_threshold.Rdata")
