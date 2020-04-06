library(tidyverse)
library(reshape2)
library(vegan)
library(betapart)

setwd("c:/Users/mathon/Desktop/linux/Global_eDNA/")
load("Rdata/02_clean_all.Rdata")

#Remove estuary stations and deep niskin station
df_all_filters <- subset(df_all_filters, !(station %in% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")))
df_all_filters <- subset(df_all_filters, sample_method!="niskin")


# gamma global =1603

gamma_global <- as.numeric(df_all_filters %>%
  summarise(n = n_distinct(sequence)))
  

region <- unique(df_all_filters$region)
site <- unique(df_all_filters$site)
station <- unique(df_all_filters$station)

# calculate alpha diversity per station
alpha_station=data.frame(site=character(94), station=character(94), motu=numeric(94), stringsAsFactors = FALSE)

for (i in 1:length(station)) {
  st <- station[i]
  s <- unique(df_all_filters[df_all_filters$station == station[i],]$site)
  motu <- df_all_filters[df_all_filters$station == station[i],] %>%
    summarise(n = n_distinct(sequence))
  alpha_station[i,1] <- s
  alpha_station[i,2] <- st
  alpha_station[i,3] <- motu
}

mean_alpha_station <- mean(alpha_station$motu)

# calculate alpha diversity per site

alpha_site=data.frame(region=character(16), site=character(16), motu=numeric(16), stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  s <- site[i]
  r <- unique(df_all_filters[df_all_filters$site == site[i],]$region)
  motu <- df_all_filters[df_all_filters$site == site[i],] %>%
    summarise(n = n_distinct(sequence))
  alpha_site[i,1] <- r
  alpha_site[i,2] <- s
  alpha_site[i,3] <- motu
}

# calculate beta inter-station

beta_station <- data.frame(site=character(16), beta=numeric(16), stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  s <- site[i]
  gamma <- alpha_site[alpha_site$site==site[i],]$motu
  alpha <- mean(alpha_station[alpha_station$site==site[i],]$motu)
  beta_station[i,1] <- s
  beta_station[i,2] <- gamma - alpha
}

mean_beta_station <- mean(beta_station$beta)


# calculate alpha diversity per region

alpha_region=data.frame(region=character(4), motu=numeric(4), stringsAsFactors = FALSE)

for (i in 1:length(region)) {
  r <- region[i]
  motu <- df_all_filters[df_all_filters$region == region[i],] %>%
    summarise(n = n_distinct(sequence))
  alpha_region[i,1] <- r
  alpha_region[i,2] <- motu
}

# calculate beta inter-site

beta_site <- data.frame(region=character(4), beta=numeric(4), stringsAsFactors = FALSE)

for (i in 1:length(region)) {
  r <- region[i]
  gamma <- alpha_region[alpha_region$region==region[i],]$motu
  alpha <- mean(alpha_site[alpha_site$region==region[i],]$motu)
  beta_site[i,1] <- r
  beta_site[i,2] <- gamma - alpha
  
}

mean_beta_site <- mean(beta_site$beta)

# calculate beta inter-region

beta_region <- gamma_global - mean_alpha_station - mean_beta_site - mean_beta_station

beta_region+mean_beta_site+mean_beta_station+mean_alpha_station

div_partition <- data.frame(component=c("mean_alpha_station", "mean_beta_station", "mean_beta_site", "beta_region"), value= c(mean_alpha_station, mean_beta_station, mean_beta_site, beta_region), percent=numeric(4))
div_partition$percent <- (div_partition$value*100)/gamma_global

write.csv(div_partition, "outputs/06_diversity_partitioning/diversity_partitioning.csv")



# check with betapart
df_region=data.frame(motu=character())

for (i in 1:length(region)) {
  df <- df_all_filters[df_all_filters$region == region[i],] %>%
    distinct(sequence, region)
  colnames(df) <- c("motu", region[i])
  df_region <- full_join(df_region, df, by="motu")
}
rownames(df_region) <- df_region[,1]
df_region <- decostand(df_region[,2:5], "pa",na.rm = TRUE)
df_region[is.na(df_region)] <- 0
df_region <- as.data.frame(t(df_region))

beta<- betapart.core(df_region)
beta.multi(beta, "jaccard")
