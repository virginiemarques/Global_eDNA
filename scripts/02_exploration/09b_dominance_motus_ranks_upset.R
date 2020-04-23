library(UpSetR)
library(devtools)
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

# ---------------------------------------------------------------- # 

# Lib
library(tidyverse)
library(reshape2)
library(vegan)
library(betapart)
library(ggpubr)

"%ni%" <- Negate("%in%")

# data
load("Rdata/02_clean_all.Rdata")

#Remove estuary stations and deep niskin station
df_all_filters <- subset(df_all_filters, !(station %in% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")))
df_all_filters <- subset(df_all_filters, sample_method!="niskin")
df_all_filters <- subset(df_all_filters, region!="East_Pacific")

# on est d'accord qu'on enlève toujours les MOTUs non assignés au dessus de la famille dans ce MS? 
df_all_filters <- df_all_filters %>%
  filter(!is.na(new_family_name))

# the motus 
df_regions <- split(df_all_filters, df_all_filters$region)

# Construct the inital matrix - add some important infos as a side
sampling_effort <-  df_all_filters %>%
  group_by(region) %>%
  summarize(n_samples = n_distinct(sample_name_all_pcr))

matrix_motus <- df_all_filters %>%
  distinct(sequence, new_scientific_name_ncbi) %>%
  mutate(West_papua = ifelse(sequence %in% df_regions$West_Papua$sequence, 1, 0), 
         French_polynesia = ifelse(sequence %in% df_regions$French_Polynesia$sequence, 1, 0),
         Caribbean = ifelse(sequence %in% df_regions$Caribbean$sequence, 1, 0)) %>%
  left_join(., sampling_effort)




# contruct matrix
mm = make_comb_mat(matrix_motus)
mm

# Remove the intersections with no match
mm = mm[comb_degree(mm) > 0]

# Color panel 

library(nationalparkcolors)
pal <- park_palette("Everglades", 3)
pal

# the upset plot

p <- UpSet(mm, 
      comb_order = order(-comb_size(mm)), 
      comb_col = pal[comb_degree(mm)])
p

# Save
png('outputs/09_dominance_motus_ranks/upset_plot_motus.png', width = 6, height=3, units = "in", res=150)
p
dev.off()


# Supp Settings

# Ideas: rajouter le nombre de familles/genres sur le cote, par region? en boxplot. Le nombre de samples peut être ? 

# remove the sole intersections
UpSet(mm[comb_degree(mm) == 2])


























# -------------------------------------- # 
# Tuto 

movies = read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), 
                  header = TRUE, sep = ";")
head(movies) # `make_comb_mat()` automatically ignores the first two columns

m = make_comb_mat(movies, top_n_sets = 6)
m

m = m[comb_degree(m) > 0]
UpSet(m)


# --------------------- # 
# Add graphs on the side 

m = make_comb_mat(movies[, genre])

comb_elements = lapply(comb_name(m), function(nm) extract_comb(m, nm))
years = lapply(comb_elements, function(ind) movies$ReleaseDate[ind])
rating = lapply(comb_elements, function(ind) movies$AvgRating[ind])
watches = lapply(comb_elements, function(ind) movies$Watches[ind])


comb_elements = lapply(comb_name(mm), function(nm) extract_comb(mm, nm))












  



