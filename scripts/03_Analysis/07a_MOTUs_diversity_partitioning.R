library(tidyverse)
library(reshape2)
library(vegan)
library(betapart)
library(ggpubr)


load("Rdata/02-clean-data.Rdata")
'%ni%' <- Negate("%in%")
#Remove estuary stations and deep niskin station
df_all_filters <- df_all_filters %>%
  filter(province %in% c("Western_Indian_Ocean", "Southeast_Polynesia", "Tropical_Northwestern_Atlantic", "Western_Coral_Triangle", "Tropical_Southwestern_Pacific"))%>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3", "glorieuse_distance_300m")) %>%
  filter(sample_method !="niskin" & comment %ni% c("Distance decay 600m", "Distance decay 300m"))%>%
  filter(project.y != "SEAMOUNTS") %>%
  filter(habitat=="marine")%>%
  filter(habitat_type %ni% c("BAIE"))


# different tests
  ## 1. all MOTUs
  ## 2. all assigned MOTUs
df_all_filters <- df_all_filters %>%
  filter(!is.na(new_species_name))
  ## 3. all assigned MOTUs without cryptobenthics
df_all_filters <- subset(df_all_filters, !(family_name_corrected %in% c("Tripterygiidae", "Grammatidae", "Aploactinidae", "Creediidae", "Gobiidae", "Chaenopsidae", "Gobiesocidae", "Labrisomidae", "Pseudochromidae", "Bythitidae", "Plesiopidae", "Dactyloscopidae", "Blenniidae", "Apogonidae", "Callionymidae", "Opistognathidae", "Syngnathidae")))
  ## 4. crypto only
df_all_filters <- subset(df_all_filters, family_name_corrected %in% c("Tripterygiidae", "Grammatidae", "Aploactinidae", "Creediidae", "Gobiidae", "Chaenopsidae", "Gobiesocidae", "Labrisomidae", "Pseudochromidae", "Bythitidae", "Plesiopidae", "Dactyloscopidae", "Blenniidae", "Apogonidae", "Callionymidae", "Opistognathidae", "Syngnathidae"))
  ## 5. pelagic only
load("Rdata/pelagic_family.Rdata")
df_all_filters <- subset(df_all_filters, family_name_corrected %in% pelagic_family$family_name)



# gamma global =2160

gamma_global <- as.numeric(df_all_filters %>%
  summarise(n = n_distinct(sequence)))
  

Province <- unique(df_all_filters$province)
Site <- unique(df_all_filters$site35)
Station <- unique(df_all_filters$station)

# calculate alpha region

alpha_province=data.frame(province=character(), motu=numeric(), stringsAsFactors = FALSE)

for (i in 1:length(Province)) {
  r <- Province[i]
  motu <- df_all_filters[df_all_filters$province == Province[i],] %>%
    summarise(n = n_distinct(sequence))
  alpha_province[i,1] <- r
  alpha_province[i,2] <- motu
}

# Calculate beta interregion
beta_province <- data.frame(alpha=mean(alpha_province$motu), gamma=gamma_global, beta=numeric(1), scale="inter-province")
beta_province$beta <- beta_province$gamma - beta_province$alpha
save(beta_province, file="Rdata/beta_province.rdata")

# calculate alpha site

alpha_site=data.frame(province=character(), site=character(), motu=numeric(), stringsAsFactors = FALSE)

for (i in 1:length(Site)) {
  s <- Site[i]
  r <- unique(df_all_filters[df_all_filters$site35 == Site[i],]$province)
  motu <- df_all_filters[df_all_filters$site35 == Site[i],] %>%
    summarise(n = n_distinct(sequence))
  alpha_site[i,1] <- r
  alpha_site[i,2] <- s
  alpha_site[i,3] <- motu
}


# calculate beta inter-site

beta_site <- data.frame(province=character(length(Province)), alpha=numeric(length(Province)), gamma=numeric(length(Province)), beta=numeric(length(Province)), scale="inter-site", stringsAsFactors = FALSE)

for (i in 1:length(Province)) {
  r <- Province[i]
  gamma <- alpha_province[alpha_province$province==Province[i],]$motu
  alpha <- mean(alpha_site[alpha_site$province==Province[i],]$motu)
  beta_site[i,1] <- r
  beta_site[i,2] <- alpha
  beta_site[i,3] <- gamma
  beta_site[i,4] <- gamma - alpha
  
}

mean_beta_site <- mean(beta_site$beta) 
sd_beta_site <- sd(beta_site$beta)
save(beta_site, file="Rdata/beta_site.rdata")


# calculate alpha station
alpha_station=data.frame(province=character(), site=character(), station=character(), motu=numeric(), stringsAsFactors = FALSE)

for (i in 1:length(Station)) {
  st <- Station[i]
  r <- unique(df_all_filters[df_all_filters$station == Station[i],]$province)
  s <- unique(df_all_filters[df_all_filters$station == Station[i],]$site35)
  motu <- df_all_filters[df_all_filters$station == Station[i],] %>%
    summarise(n = n_distinct(sequence))
  alpha_station[i,1] <- r
  alpha_station[i,2] <- s
  alpha_station[i,3] <- st
  alpha_station[i,4] <- motu
}

mean_a_station <- data.frame(stringsAsFactors = FALSE)
save(alpha_station, file="Rdata/alpha_station.rdata")
for (i in 1:length(Site)) {
  df <- alpha_station %>%
    filter(site==Site[i])
  mean_a_station[i,1] <- mean(df$motu)
  mean_a_station[i,2] <- Site[i]
  
}

mean_alpha_station <- mean(mean_a_station$V1)
sd_alpha_station <- sd(mean_a_station$V1)

# calculate beta inter-station

beta_station <- data.frame(province=character(length(Site)), site=character(length(Site)), alpha=numeric(length(Site)), gamma=numeric(length(Site)), beta=numeric(length(Site)), scale="inter-station", stringsAsFactors = FALSE)

for (i in 1:length(Site)) {
  s <- Site[i]
  r <- unique(alpha_station[alpha_station$site==Site[i],]$province)
  gamma <- alpha_site[alpha_site$site==Site[i],]$motu
  alpha <- mean(alpha_station[alpha_station$site==Site[i],]$motu)
  beta_station[i,1] <- r
  beta_station[i,2] <- s
  beta_station[i,3] <- alpha
  beta_station[i,4] <- gamma
  beta_station[i,5] <- gamma - alpha
}

mean_b_station <- data.frame(stringsAsFactors = FALSE)
save(beta_station, file="Rdata/beta_station.rdata")

for (i in 1:length(Province)) {
  df <- beta_station %>%
    subset(province==Province[i])
  mean_b_station[i,1] <- mean(df$beta)
  
}

mean_beta_station <- mean(mean_b_station$V1)
sd_beta_station <- sd(mean_b_station$V1)

beta_province$beta+mean_beta_site+mean_beta_station+mean_alpha_station

# calculate diversity partitioning

div_partition <- data.frame(component=c("mean_alpha_station", "mean_beta_station", "mean_beta_site", "beta_province"), 
                            value=c(mean_alpha_station, mean_beta_station, mean_beta_site, beta_province$beta),
                            sd=c(sd_alpha_station, sd_beta_station, sd_beta_site, 0),
                            percent=numeric(4))
div_partition$percent <- (div_partition$value*100)/gamma_global

write.csv(div_partition, "outputs/07_diversity_partitioning/all_motus_diversity_partitioning.csv")




# calculate beta, turnover and nestedness with betapart
  ## beta inter-regions

df_province=data.frame(motu=character())

for (i in 1:length(Province)) {
  df <- df_all_filters[df_all_filters$province == Province[i],] %>%
    distinct(sequence, province)
  colnames(df) <- c("motu", Province[i])
  df_province <- full_join(df_province, df, by="motu")
}
rownames(df_province) <- df_province[,1]
df_province <- decostand(df_province[,2:6], "pa",na.rm = TRUE)
df_province[is.na(df_province)] <- 0
df_province <- as.data.frame(t(df_province))

b <- betapart.core(df_province)
beta <- beta.multi(b, "jaccard")
betaprovince <- data.frame(scale="inter-province", total=beta$beta.JAC, turnover=beta$beta.JTU, nestedness=beta$beta.JNE)

  ## beta inter-site

df_site=vector("list", length(Province))
betasite <- data.frame(scale="inter-site", total=numeric(length(Province)), turnover=numeric(length(Province)), nestedness=numeric(length(Province)))


for (i in 1:length(Province)) {
  df <- df_all_filters[df_all_filters$province == Province[i],]
  Site <- unique(df$site35)
  df_site[[i]] <- data.frame(motu=character(), stringsAsFactors = FALSE)
    for (j in 1:length(Site)) {
      df2 <- df[df$site35==Site[j],] %>%
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

Site <- unique(df_all_filters$site35)
df_station=vector("list", 25)
betastation <- data.frame(scale="inter-station", total=numeric(length(Site)), turnover=numeric(length(Site)), nestedness=numeric(length(Site)))


for (i in 1:length(Site)) {
  df <- df_all_filters[df_all_filters$site35 == Site[i],]
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


betatotal <- rbind(betaprovince, betasite, betastation)
write.csv(betatotal, "outputs/07_diversity_partitioning/beta_motus.csv")

beta_melt <- reshape2::melt(betatotal)

beta_motus_eDNA <- ggplot(beta_melt, aes(variable, value, colour=scale))+
  geom_boxplot()+
  theme_bw()+
  ylim(0,1)+
  theme(legend.position = "none")+
  labs(x="", y="Jaccard dissimilarity (motus)")

ggsave("outputs/07_diversity_partitioning/all_motus_diversity_partitioning.png")  
save(beta_motus_eDNA, file = "Rdata/beta_motus_eDNA.rdata")

