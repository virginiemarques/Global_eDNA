library(reshape2)
# count RLS
RLS_species <- read.csv("data/RLS/RLS_species_NEW.csv", sep = ";", stringsAsFactors = FALSE, check.names = FALSE)
RLS_species <- RLS_species %>%
  filter(Realm%ni%c("Temperate Australasia", "Temperate Northern Atlantic", "Temperate Southern Africa", "Temperate Northern Pacific"))
RLS_sp <- RLS_species[,c(10:2165)]
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
