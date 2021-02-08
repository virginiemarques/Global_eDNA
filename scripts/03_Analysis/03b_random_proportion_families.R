library(tidyverse)
library(reshape2)
library(gridExtra)
library(grid)
library(cowplot)
library(ggplot2)
library(ggpubr)

load("Rdata/02_clean_all.Rdata")
'%ni%' <- Negate("%in%")
#Remove estuary stations and deep niskin station
df_all_filters <- df_all_filters %>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")) %>%
  filter(sample_method !="niskin" & province!="Tropical_East_Pacific" & comment %ni% c("Distance decay 600m", "Distance decay 300m") & station!="glorieuse_distance_300m")%>%
  filter(project != "SEAMOUNTS") %>% 
  filter(habitat_type %ni% c("BAIE", "Sommet"))

df_all_filters <- df_all_filters %>%
  filter(!is.na(new_family_name))

unique_motus <- unique(df_all_filters$sequence)
unique_family <- as.data.frame(unique(df_all_filters$new_family_name))
colnames(unique_family) <- "family"

global_family <- df_all_filters %>%
  distinct(sequence, new_family_name)
colnames(global_family) <- c("motu", "family")

#calculate frequency of motus
site <- unique(df_all_filters$site)

df_site <- data.frame(motu=character(), stringsAsFactors = FALSE)
for (i in 1:length(site)) {
  df <- df_all_filters[df_all_filters$site==site[i],] %>%
    distinct(sequence, site)
  colnames(df) <- c("motu", site[i])
  df_site <- full_join(df_site, df, by="motu")
}
rownames(df_site) <- df_site[,1]
df_site <- df_site[order(rownames(df_site)),]
df_site <- decostand(df_site[,c(-1)], "pa",na.rm = TRUE)
df_site[is.na(df_site)] <- 0
df_site <- as.data.frame(t(df_site))

prob_motu <- colSums(df_site)/25

# sample random family assignations for site with 1-800 species, 1000 times (long! load Rdata)
random_sp <- vector("list", 800)
random_prop_tot <- data.frame()
rep <- 800

for (j in seq(1:1000)) {
  random_prop <- vector("list" )
  for (i in 1:rep) {
    random_sp[[i]] <- as.data.frame(sample(unique_motus, i, replace = TRUE, prob = prob_motu))
    colnames(random_sp[[i]]) <- "motu"
    random_sp[[i]] <- left_join(random_sp[[i]], global_family, by="motu")
    random_sp[[i]] <- data.frame(table(random_sp[[i]]$family))
    colnames(random_sp[[i]]) <- c("family", "n_motus")
    random_sp[[i]] <- full_join(random_sp[[i]], unique_family, by="family")
    random_sp[[i]][is.na(random_sp[[i]])] <- 0
    random_sp[[i]]$n_total <- sum(random_sp[[i]]$n_motus)
    random_sp[[i]]$prop <- random_sp[[i]]$n_motus / random_sp[[i]]$n_total
    random_prop <- rbind(random_prop, random_sp[[i]])
  }
  random_prop_tot <- rbind(random_prop_tot, random_prop)
}

save(random_prop_tot, file = "c:/Users/mathon/Desktop/PhD/Projets/Megafauna/Global_eDNA/Rdata/random_family_proportions.rdata")

load("c:/Users/mathon/Desktop/PhD/Projets/Megafauna/Global_eDNA/Rdata/random_family_proportions.rdata")



# calculate 2.5 and 97.5 quantile for each sample size and each family
family <- c("Acanthuridae", "Chaetodontidae", "Labridae", "Lutjanidae", "Serranidae", "Carangidae", "Pomacentridae", "Apogonidae", "Gobiidae")

CI_family <- data.frame()
quant <- data.frame(upper=numeric(800), lower=numeric(800))
for (j in 1:length(family)) {
  df <- random_prop_tot[random_prop_tot$family==family[j],]
  for (i in seq(1:800)) {
    df2 <- df[df$n_total==i,]
    quant[i,] <- quantile(df2$prop, c(.975, .025))
    quant$n_total <- seq(1:800)
    quant$family <- family[j]
  }
  CI_family <- rbind(CI_family, quant)
}

save(CI_family, file = "Rdata/CI_null_model_family_proportions.rdata")

