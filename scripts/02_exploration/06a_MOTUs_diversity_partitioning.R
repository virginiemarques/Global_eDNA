library(tidyverse)
library(reshape2)
library(vegan)
library(betapart)


load("Rdata/02_clean_all.Rdata")
'%ni%' <- Negate("%in%")
#Remove estuary stations and deep niskin station
df_all_filters <- df_all_filters %>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")) %>%
  filter(sample_method !="niskin" & region!="East_Pacific" & comment %ni% c("Distance decay 600m", "Distance decay 300m") & station!="glorieuse_distance_300m")%>%
  filter(project != "SEAMOUNTS") %>% 
  filter(habitat_type %ni% c("BAIE", "Sommet"))


# different tests
  ## 1. all MOTUs
  ## 2. all assigned MOTUs
df_all_filters <- df_all_filters %>%
  filter(!is.na(new_species_name))
  ## 4. all assigned MOTUs without cryptobenthics
df_all_filters <- subset(df_all_filters, !(new_family_name %in% c("Tripterygiidae", "Grammatidae", "Aploactinidae", "Creediidae", "Gobiidae", "Chaenopsidae", "Gobiesocidae", "Labrisomidae", "Pseudochromidae", "Bythitidae", "Plesiopidae", "Dactyloscopidae", "Blenniidae", "Apogonidae", "Callionymidae", "Opistognathidae", "Syngnathidae")))
  ## 5. crypto only
df_all_filters <- subset(df_all_filters, new_family_name %in% c("Tripterygiidae", "Grammatidae", "Aploactinidae", "Creediidae", "Gobiidae", "Chaenopsidae", "Gobiesocidae", "Labrisomidae", "Pseudochromidae", "Bythitidae", "Plesiopidae", "Dactyloscopidae", "Blenniidae", "Apogonidae", "Callionymidae", "Opistognathidae", "Syngnathidae"))
  ## 6. pelagic only
load("Rdata/pelagic_family.Rdata")
df_all_filters <- subset(df_all_filters, new_family_name %in% pelagic_family$family_name)



# gamma global =2160

gamma_global <- as.numeric(df_all_filters %>%
  summarise(n = n_distinct(sequence)))
  

Region <- unique(df_all_filters$region)
Site <- unique(df_all_filters$site)
Station <- unique(df_all_filters$station)

# calculate alpha region

alpha_region=data.frame(region=character(5), motu=numeric(5), stringsAsFactors = FALSE)

for (i in 1:length(Region)) {
  r <- Region[i]
  motu <- df_all_filters[df_all_filters$region == Region[i],] %>%
    summarise(n = n_distinct(sequence))
  alpha_region[i,1] <- r
  alpha_region[i,2] <- motu
}

# Calculate beta interregion
beta_region <- data.frame(alpha=mean(alpha_region$motu), gamma=gamma_global, beta=numeric(1), scale="inter-region")
beta_region$beta <- beta_region$gamma - beta_region$alpha


# calculate alpha site

alpha_site=data.frame(region=character(), site=character(), motu=numeric(), stringsAsFactors = FALSE)

for (i in 1:length(Site)) {
  s <- Site[i]
  r <- unique(df_all_filters[df_all_filters$site == Site[i],]$region)
  motu <- df_all_filters[df_all_filters$site == Site[i],] %>%
    summarise(n = n_distinct(sequence))
  alpha_site[i,1] <- r
  alpha_site[i,2] <- s
  alpha_site[i,3] <- motu
}


# calculate beta inter-site

beta_site <- data.frame(region=character(5), alpha=numeric(5), gamma=numeric(5), beta=numeric(5), scale="inter-site", stringsAsFactors = FALSE)

for (i in 1:length(Region)) {
  r <- Region[i]
  gamma <- alpha_region[alpha_region$region==Region[i],]$motu
  alpha <- mean(alpha_site[alpha_site$region==Region[i],]$motu)
  beta_site[i,1] <- r
  beta_site[i,2] <- alpha
  beta_site[i,3] <- gamma
  beta_site[i,4] <- gamma - alpha
  
}

mean_beta_site <- mean(beta_site$beta) 
sd_beta_site <- sd(beta_site$beta)



# calculate alpha station
alpha_station=data.frame(region=character(), site=character(), station=character(), motu=numeric(), stringsAsFactors = FALSE)

for (i in 1:length(Station)) {
  st <- Station[i]
  r <- unique(df_all_filters[df_all_filters$station == Station[i],]$region)
  s <- unique(df_all_filters[df_all_filters$station == Station[i],]$site)
  motu <- df_all_filters[df_all_filters$station == Station[i],] %>%
    summarise(n = n_distinct(sequence))
  alpha_station[i,1] <- r
  alpha_station[i,2] <- s
  alpha_station[i,3] <- st
  alpha_station[i,4] <- motu
}

mean_a_station <- data.frame(stringsAsFactors = FALSE)

for (i in 1:length(Site)) {
  df <- alpha_station %>%
    filter(site==Site[i])
  mean_a_station[i,1] <- mean(df$motu)
  mean_a_station[i,2] <- Site[i]
  
}

mean_alpha_station <- mean(mean_a_station$V1)
sd_alpha_station <- sd(mean_a_station$V1)

# calculate beta inter-station

beta_station <- data.frame(region=character(25), site=character(25), alpha=numeric(25), gamma=numeric(25), beta=numeric(25), scale="inter-station", stringsAsFactors = FALSE)

for (i in 1:length(Site)) {
  s <- Site[i]
  r <- unique(alpha_station[alpha_station$site==Site[i],]$region)
  gamma <- alpha_site[alpha_site$site==Site[i],]$motu
  alpha <- mean(alpha_station[alpha_station$site==Site[i],]$motu)
  beta_station[i,1] <- r
  beta_station[i,2] <- s
  beta_station[i,3] <- alpha
  beta_station[i,4] <- gamma
  beta_station[i,5] <- gamma - alpha
}

mean_b_station <- data.frame(stringsAsFactors = FALSE)

for (i in 1:length(Region)) {
  df <- beta_station %>%
    subset(region==Region[i])
  mean_b_station[i,1] <- mean(df$beta)
  
}

mean_beta_station <- mean(mean_b_station$V1)
sd_beta_station <- sd(mean_b_station$V1)

beta_region$beta+mean_beta_site+mean_beta_station+mean_alpha_station

# calculate diversity partitioning

div_partition <- data.frame(component=c("mean_alpha_station", "mean_beta_station", "mean_beta_site", "beta_region"), 
                            value=c(mean_alpha_station, mean_beta_station, mean_beta_site, beta_region$beta),
                            sd=c(sd_alpha_station, sd_beta_station, sd_beta_site, 0),
                            percent=numeric(4))
div_partition$percent <- (div_partition$value*100)/gamma_global

write.csv(div_partition, "outputs/06_diversity_partitioning/all_motus_diversity_partitioning.csv")




# calculate beta, turnover and nestedness with betapart
  ## beta inter-regions

df_region=data.frame(motu=character())

for (i in 1:length(Region)) {
  df <- df_all_filters[df_all_filters$region == Region[i],] %>%
    distinct(sequence, region)
  colnames(df) <- c("motu", Region[i])
  df_region <- full_join(df_region, df, by="motu")
}
rownames(df_region) <- df_region[,1]
df_region <- decostand(df_region[,2:6], "pa",na.rm = TRUE)
df_region[is.na(df_region)] <- 0
df_region <- as.data.frame(t(df_region))

b <- betapart.core(df_region)
beta <- beta.multi(b, "jaccard")
betaregion <- data.frame(scale="inter-region", total=beta$beta.JAC, turnover=beta$beta.JTU, nestedness=beta$beta.JNE)

  ## beta inter-site

df_site=vector("list", 5)
betasite <- data.frame(scale="inter-site", total=numeric(5), turnover=numeric(5), nestedness=numeric(5))


for (i in 1:length(Region)) {
  df <- df_all_filters[df_all_filters$region == Region[i],]
  Site <- unique(df$site)
  df_site[[i]] <- data.frame(motu=character(), stringsAsFactors = FALSE)
    for (j in 1:length(Site)) {
      df2 <- df[df$site==Site[j],] %>%
        distinct(sequence, site)
      colnames(df2) <- c("motu", Site[j])
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

Site <- unique(df_all_filters$site)
df_station=vector("list", 25)
betastation <- data.frame(scale="inter-station", total=numeric(25), turnover=numeric(25), nestedness=numeric(25))


for (i in 1:length(Site)) {
  df <- df_all_filters[df_all_filters$site == Site[i],]
  Station <- unique(df$station)
  df_station[[i]] <- data.frame(motu=character(), stringsAsFactors = FALSE)
  for (j in 1:length(Station)) {
    df2 <- df[df$station==Station[j],] %>%
      distinct(sequence, station)
    colnames(df2) <- c("motu", Station[j])
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
write.csv(betatotal, "outputs/06_diversity_partitioning/beta_motus.csv")

beta_melt <- reshape2::melt(betatotal)

beta_motus_eDNA <- ggplot(beta_melt, aes(variable, value, colour=scale))+
  geom_boxplot()+
  theme_bw()+
  ylim(0,1)+
  theme(legend.position = "none")+
  labs(x="", y="Jaccard dissimilarity (motus)")

ggsave("outputs/06_diversity_partitioning/all_motus_diversity_partitioning.png")  
save(beta_motus_eDNA, file = "Rdata/beta_motus_eDNA.rdata")




