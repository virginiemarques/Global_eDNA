library(tidyverse)
library(reshape2)
library(vegan)
library(betapart)

setwd("c:/Users/mathon/Desktop/linux/Global_eDNA/")
load("Rdata/02_clean_all.Rdata")

#Remove estuary stations and deep niskin station
df_all_filters <- subset(df_all_filters, !(station %in% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")))
df_all_filters <- subset(df_all_filters, sample_method!="niskin")

# different tests
  ## 1. all MOTUs
  ## 2. all assigned MOTUs
df_all_filters <- df_all_filters %>%
  filter(!is.na(new_family_name))
  ## all assigned MOTUs without cryptobenthics
df_all_filters <- subset(df_all_filters, !(new_family_name %in% c("Tripterygiidae", "Grammatidae", "Aploactinidae", "Creediidae", "Gobiidae", "Chaenopsidae", "Gobiesocidae", "Labrisomidae", "Pseudochromidae", "Bythitidae", "Plesiopidae", "Dactyloscopidae", "Blenniidae", "Apogonidae", "Callionymidae", "Opistognathidae", "Syngnathidae")))
  ## crypto uniquement
df_all_filters <- subset(df_all_filters, new_family_name %in% c("Tripterygiidae", "Grammatidae", "Aploactinidae", "Creediidae", "Gobiidae", "Chaenopsidae", "Gobiesocidae", "Labrisomidae", "Pseudochromidae", "Bythitidae", "Plesiopidae", "Dactyloscopidae", "Blenniidae", "Apogonidae", "Callionymidae", "Opistognathidae", "Syngnathidae"))


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
sd_alpha_station <- sd(alpha_station$motu)

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
sd_beta_station <- sd(beta_station$beta)

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
sd_beta_site <- sd(beta_site$beta)

# calculate beta inter-region

beta_region <- gamma_global - mean_alpha_station - mean_beta_site - mean_beta_station

beta_region+mean_beta_site+mean_beta_station+mean_alpha_station

div_partition <- data.frame(component=c("mean_alpha_station", "mean_beta_station", "mean_beta_site", "beta_region"), 
                            value= c(mean_alpha_station, mean_beta_station, mean_beta_site, beta_region), 
                            sd=c(sd_alpha_station, sd_beta_station, sd_beta_site, 0),
                            percent=numeric(4))
div_partition$percent <- (div_partition$value*100)/gamma_global

write.csv(div_partition, "outputs/06_diversity_partitioning/diversity_partitioning_only_crypto.csv")



# calculate beta, turnover and nestedness with betapart
  ## beta inter-regions

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

b <- betapart.core(df_region)
beta <- beta.multi(b, "jaccard")
betaregion <- data.frame(scale="region", total=beta$beta.JAC, turnover=beta$beta.JTU, nestedness=beta$beta.JNE)

  ## beta inter-site

df_site=vector("list", 4)
betasite <- data.frame(scale="site", total=numeric(4), turnover=numeric(4), nestedness=numeric(4))


for (i in 1:length(region)) {
  df <- df_all_filters[df_all_filters$region == region[i],]
  site <- unique(df$site)
  df_site[[i]] <- data.frame(motu=character(), stringsAsFactors = FALSE)
    for (j in 1:length(site)) {
      df2 <- df[df$site==site[j],] %>%
        distinct(sequence, site)
      colnames(df2) <- c("motu", site[j])
      df_site[[i]] <- full_join(df_site[[i]], df2, by="motu")
    }
  rownames(df_site[[i]]) <- df_site[[i]][,1]
  df_site[[i]] <- decostand(df_site[[i]][,c(-1)], "pa",na.rm = TRUE)
  df_site[[i]][is.na(df_site[[i]])] <- 0
  df_site[[i]] <- as.data.frame(t(df_site[[i]]))
  
  b <- betapart.core(df_site[[i]])
  beta <- beta.multi(b, "jaccard")
  betasite[i,2] <- beta$beta.JAC
  betasite[i,3] <- beta$beta.JTU
  betasite[i,4] <- beta$beta.JNE
}


## beta inter-station

site <- unique(df_all_filters$site)
df_station=vector("list", 16)
betastation <- data.frame(scale="station", total=numeric(16), turnover=numeric(16), nestedness=numeric(16))


for (i in 1:length(site)) {
  df <- df_all_filters[df_all_filters$site == site[i],]
  station <- unique(df$station)
  df_station[[i]] <- data.frame(motu=character(), stringsAsFactors = FALSE)
  for (j in 1:length(station)) {
    df2 <- df[df$station==station[j],] %>%
      distinct(sequence, station)
    colnames(df2) <- c("motu", station[j])
    df_station[[i]] <- full_join(df_station[[i]], df2, by="motu")
  }
  rownames(df_station[[i]]) <- df_station[[i]][,1]
  df_station[[i]] <- decostand(df_station[[i]][,c(-1)], "pa",na.rm = TRUE)
  df_station[[i]][is.na(df_station[[i]])] <- 0
  df_station[[i]] <- as.data.frame(t(df_station[[i]]))
  
  b <- betapart.core(df_station[[i]])
  beta <- beta.multi(b, "jaccard")
  betastation[i,2] <- beta$beta.JAC
  betastation[i,3] <- beta$beta.JTU
  betastation[i,4] <- beta$beta.JNE
}


betatotal <- rbind(betaregion, betasite, betastation)

beta_melt <- melt(betatotal)

ggplot(beta_melt, aes(variable, value, colour=scale))+
  geom_boxplot()+
  labs(x=" Beta component", y="Jaccard dissimilarity")

ggsave("outputs/06_diversity_partitioning/diversity_partitioning_only_crypto.png")  
