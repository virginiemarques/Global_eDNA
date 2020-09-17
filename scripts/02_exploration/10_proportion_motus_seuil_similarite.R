library(tidyverse)

load("Rdata/02_clean_all.Rdata")
'%ni%' <- Negate("%in%")
#Remove estuary stations and deep niskin station
df_all_filters <- df_all_filters %>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")) %>%
  filter(sample_method !="niskin" & region!="East_Pacific" & comment %ni% c("Distance decay 600m", "Distance decay 300m") & station!="glorieuse_distance_300m")%>%
  filter(project != "SEAMOUNTS") %>% 
  filter(habitat_type %ni% c("BAIE", "Sommet"))


df_all_filters <- df_all_filters %>%
  filter(!is.na(new_family_name))

# Calculate proportion of MOTUs in each family, within each similarity intervals
family <- unique(df_all_filters$new_family_name)

prop_similarity <- data.frame(stringsAsFactors = FALSE)

for (i in 1:length(family)) {
  df <- subset(df_all_filters, new_family_name==family[i])
  ntot <- nrow(df)
  prop_similarity[i,"family"] <- family[i]
  prop_similarity[i,"85-90%"] <- (nrow(subset(df, best_identity_database > 0.85 & best_identity_database <= 0.90))*100)/ntot
  prop_similarity[i,"90-95%"] <- (nrow(subset(df, best_identity_database > 0.90 & best_identity_database <= 0.95))*100)/ntot
  prop_similarity[i,"95-<100%"] <- (nrow(subset(df, best_identity_database > 0.95 & best_identity_database < 1))*100)/ntot
  prop_similarity[i,"100%"] <- (nrow(subset(df, best_identity_database == 1))*100)/ntot
  
}

save(prop_similarity, file = "Rdata/proportion_motus_similarity_threshold.Rdata")
