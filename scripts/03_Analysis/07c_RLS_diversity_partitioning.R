library(tidyverse)
library(reshape2)
library(vegan)
library(betapart)
'%ni%' <- Negate("%in%")

# gamma global 

RLS_species <- read.csv("data/RLS/RLS_species_NEW.csv", sep = ";", stringsAsFactors = FALSE)
RLS_species <- RLS_species[, c(1,2,11,7,18:2173)]
RLS_species <- reshape2::melt(RLS_species, id=c("SurveyID", "Station", "site35", "Province"))
RLS_species <- RLS_species%>%
  filter(value!=0)
RLS_species <- RLS_species[,-6]
colnames(RLS_species) <- c("SurveyID", "Station", "site35", "Province", "Species")

gamma_global <- as.numeric(RLS_species %>%
                             summarise(n = n_distinct(Species)))


province <- unique(RLS_species$Province)
site <- unique(RLS_species$site35)
station <- unique(RLS_species$Station)
transect <- unique(RLS_species$SurveyID)

# calculate alpha province

alpha_province=data.frame(province=character(), species=numeric(), stringsAsFactors = FALSE)

for (i in 1:length(province)) {
  r <- province[i]
  species <- RLS_species[RLS_species$Province == province[i],] %>%
    summarise(n = n_distinct(Species))
  alpha_province[i,1] <- r
  alpha_province[i,2] <- species
}

# Calculate beta interprovince
beta_province <- data.frame(alpha=mean(alpha_province$species), gamma=gamma_global, beta=numeric(1), scale="inter-province")
beta_province$beta <- beta_province$gamma - beta_province$alpha


# calculate alpha site

alpha_site=data.frame(province=character(), site=character(), species=numeric(), stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  s <- site[i]
  r <- unique(RLS_species[RLS_species$site35 == site[i],]$Province)
  species <- RLS_species[RLS_species$site35 == site[i],] %>%
    summarise(n = n_distinct(Species))
  alpha_site[i,1] <- r
  alpha_site[i,2] <- s
  alpha_site[i,3] <- species
}


# calculate beta inter-site

beta_site <- data.frame(province=character(length(province)), alpha=numeric(length(province)), gamma=numeric(length(province)), beta=numeric(length(province)), scale="inter-site", stringsAsFactors = FALSE)

for (i in 1:length(province)) {
  r <- province[i]
  gamma <- alpha_province[alpha_province$province==province[i],]$species
  alpha <- mean(alpha_site[alpha_site$province==province[i],]$species)
  beta_site[i,1] <- r
  beta_site[i,2] <- alpha
  beta_site[i,3] <- gamma
  beta_site[i,4] <- gamma - alpha
  
}

mean_beta_site <- mean(beta_site$beta) 
sd_beta_site <- sd(beta_site$beta)



# calculate alpha station
alpha_station=data.frame(province=character(), site=character(), station=character(), species=numeric(), stringsAsFactors = FALSE)

for (i in 1:length(station)) {
  st <- station[i]
  r <- unique(RLS_species[RLS_species$Station == station[i],]$Province)
  s <- unique(RLS_species[RLS_species$Station == station[i],]$site35)
  species <- RLS_species[RLS_species$Station == station[i],] %>%
    summarise(n = n_distinct(Species))
  alpha_station[i,1] <- r
  alpha_station[i,2] <- s
  alpha_station[i,3] <- st
  alpha_station[i,4] <- species
}

alpha_station$species <- as.numeric(alpha_station$species)
mean_a_station <- data.frame(stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  df <- alpha_station %>%
    subset(site==site[i])
  mean_a_station[i,1] <- mean(df$species)
  
}


mean_alpha_station <- mean(mean_a_station$V1)
sd_alpha_station <- sd(mean_a_station$V1)

# calculate beta inter-station

beta_station <- data.frame(province=character(length(station)), site=character(length(station)), alpha=numeric(length(station)), gamma=numeric(length(station)), beta=numeric(length(station)), scale="inter-station", stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  s <- site[i]
  r <- unique(alpha_station[alpha_station$site==site[i],]$province)
  gamma <- alpha_site[alpha_site$site==site[i],]$species
  alpha <- mean(alpha_station[alpha_station$site==site[i],]$species)
  beta_station[i,1] <- r
  beta_station[i,2] <- s
  beta_station[i,3] <- alpha
  beta_station[i,4] <- gamma
  beta_station[i,5] <- gamma - alpha
}


mean_b_station <- data.frame(stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  df <- beta_station %>%
    subset(site==site[i])
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


write.csv(div_partition, "outputs/07_diversity_partitioning/RLS_diversity_partitioning.csv", row.names = F)





#__________________________________________________________________________
### Beta On species ###
#__________________________________________________________________________
RLS_species <- read.csv("data/RLS/RLS_species_NEW.csv", sep = ";", stringsAsFactors = FALSE)
RLS_sp <- RLS_species[,c(18:2173)]
RLS_sp <- RLS_sp[,colSums(RLS_sp)>0]
RLS_species <- cbind(RLS_species[,c("SurveyID", "Station", "site35", "Province")], RLS_sp)
save(RLS_species, file = "Rdata/RLS_species_clean.rdata")

province <- unique(RLS_species$Province)
site <- unique(RLS_species$site35)
station <- unique(RLS_species$Station)
transect <- unique(RLS_species$SurveyID)

# calculate beta, turnover and nestedness with betapart
## beta inter-provinces

df_province=data.frame()


for (i in 1:length(province)) {
  df <- RLS_species[RLS_species$Province == province[i],]
  df <- df[,c(5:ncol(df))]
  df_province[i,"province"] <- province[i]
  df_province[i, 2:(ncol(df)+1)] <- colSums(df)
}
rownames(df_province) <- df_province[,1]
df_province <- decostand(df_province[,2:ncol(df_province)], "pa",na.rm = TRUE)

b <- betapart.core(df_province)
beta <- beta.multi(b, "jaccard")
betaprovince <- data.frame(scale="inter-province", total=beta$beta.JAC, turnover=beta$beta.JTU, nestedness=beta$beta.JNE)

  ## beta inter-site

df_site=vector("list", length(province))
betasite <- data.frame(scale="inter-site", total=numeric(length(province)), turnover=numeric(length(province)), nestedness=numeric(length(province)))


for (i in 1:length(province)) {
  df <- RLS_species[RLS_species$Province == province[i],]
  site <- unique(df$site35)
  df_site[[i]] <- data.frame(stringsAsFactors = FALSE)
    for (j in 1:length(site)) {
      df2 <- df[df$site35==site[j],]
      df2 <- df2[,c(5:ncol(df2))]
      df_site[[i]][j,"site"]<- site[j]
      df_site[[i]][j,2:(ncol(df2)+1)] <- colSums(df2)
    }
  rownames(df_site[[i]]) <- df_site[[i]][,1]
  df_site[[i]] <- decostand(df_site[[i]][,c(-1)], "pa",na.rm = TRUE)
  
  
  b <- betapart.core(df_site[[i]])
  beta <- beta.multi(b, "jaccard")
  betasite[i,2] <- beta$beta.JAC
  betasite[i,3] <- beta$beta.JTU
  betasite[i,4] <- beta$beta.JNE
}


## beta inter-station

site <- unique(RLS_species$site35)
df_station=vector("list", length(site))
betastation <- data.frame(scale="inter-transect", total=numeric(length(site)), turnover=numeric(length(site)), nestedness=numeric(length(site)))


for (i in 1:length(site)) {
  df <- RLS_species[RLS_species$site35 == site[i],]
  station <- unique(df$Station)
  df_station[[i]] <- data.frame(stringsAsFactors = FALSE)
  for (j in 1:length(station)) {
    df2 <- df[df$Station==station[j],]
    df2 <- df2[,c(5:ncol(df2))]
    df_station[[i]][j, "station"] <- station[j]
    df_station[[i]][j, 2:(ncol(df2)+1)] <- colSums(df2)
  }
  rownames(df_station[[i]]) <- df_station[[i]][,1]
  df_station[[i]] <- decostand(df_station[[i]][,c(-1)], "pa",na.rm = TRUE)
  
  
  b <- betapart.core(df_station[[i]])
  beta <- beta.multi(b, "jaccard")
  betastation[i,2] <- beta$beta.JAC
  betastation[i,3] <- beta$beta.JTU
  betastation[i,4] <- beta$beta.JNE
}
betastation <- betastation %>%
  filter(!is.na(total))

betatotal <- rbind(betaprovince, betasite, betastation)
write.csv(betatotal, "outputs/07_diversity_partitioning/beta_RLS_species.csv")
beta_melt <- reshape2::melt(betatotal)

beta_species_RLS <- ggplot(beta_melt, aes(variable, value, colour=scale))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x=" Beta component", y="Jaccard dissimilarity (RLS Species)")

ggsave("outputs/07_diversity_partitioning/beta_RLS_species_Province.png")  
save(beta_species_RLS, file="Rdata/beta_species_RLS.rdata")


#__________________________________________________________________________
### Beta On families ###
#__________________________________________________________________________
RLS_families <- read.csv("data/RLS/RLS_families_NEW.csv", sep = ";", stringsAsFactors = FALSE)
RLS_fam <- RLS_families[,c(17:137)]
RLS_fam <- RLS_fam[,colSums(RLS_fam)>0]
RLS_families <- cbind(RLS_families[,c("SurveyID", "Station", "site35", "Province")], RLS_fam)


province <- unique(RLS_families$Province)
site <- unique(RLS_families$site35)
station <- unique(RLS_families$Station)
transect <- unique(RLS_families$SurveyID)

# calculate beta, turnover and nestedness with betapart
## beta inter-provinces

df_province=data.frame()


for (i in 1:length(province)) {
  df <- RLS_families[RLS_families$Province == province[i],]
  df <- df[,c(5:ncol(df))]
  df_province[i,"province"] <- province[i]
  df_province[i, 2:(ncol(df)+1)] <- colSums(df)
}
rownames(df_province) <- df_province[,1]
df_province <- decostand(df_province[,2:ncol(df_province)], "pa",na.rm = TRUE)

b <- betapart.core(df_province)
beta <- beta.multi(b, "jaccard")
betaprovince <- data.frame(scale="inter-province", total=beta$beta.JAC, turnover=beta$beta.JTU, nestedness=beta$beta.JNE)

## beta inter-site

df_site=vector("list", length(province))
betasite <- data.frame(scale="inter-site", total=numeric(length(province)), turnover=numeric(length(province)), nestedness=numeric(length(province)))


for (i in 1:length(province)) {
  df <- RLS_families[RLS_families$Province == province[i],]
  site <- unique(df$site35)
  df_site[[i]] <- data.frame(stringsAsFactors = FALSE)
  for (j in 1:length(site)) {
    df2 <- df[df$site35==site[j],]
    df2 <- df2[,c(5:ncol(df2))]
    df_site[[i]][j,"site"]<- site[j]
    df_site[[i]][j,2:(ncol(df2)+1)] <- colSums(df2)
  }
  rownames(df_site[[i]]) <- df_site[[i]][,1]
  df_site[[i]] <- decostand(df_site[[i]][,c(-1)], "pa",na.rm = TRUE)
  
  
  b <- betapart.core(df_site[[i]])
  beta <- beta.multi(b, "jaccard")
  betasite[i,2] <- beta$beta.JAC
  betasite[i,3] <- beta$beta.JTU
  betasite[i,4] <- beta$beta.JNE
}


## beta inter-station

site <- unique(RLS_families$site35)
df_station=vector("list", length(site))
betastation <- data.frame(scale="inter-station", total=numeric(length(site)), turnover=numeric(length(site)), nestedness=numeric(length(site)))


for (i in 1:length(site)) {
  df <- RLS_families[RLS_families$site35 == site[i],]
  station <- unique(df$Station)
  df_station[[i]] <- data.frame(stringsAsFactors = FALSE)
  for (j in 1:length(station)) {
    df2 <- df[df$Station==station[j],]
    df2 <- df2[,c(5:ncol(df2))]
    df_station[[i]][j, "station"] <- station[j]
    df_station[[i]][j, 2:(ncol(df2)+1)] <- colSums(df2)
  }
  rownames(df_station[[i]]) <- df_station[[i]][,1]
  df_station[[i]] <- decostand(df_station[[i]][,c(-1)], "pa",na.rm = TRUE)
  
  
  b <- betapart.core(df_station[[i]])
  beta <- beta.multi(b, "jaccard")
  betastation[i,2] <- beta$beta.JAC
  betastation[i,3] <- beta$beta.JTU
  betastation[i,4] <- beta$beta.JNE
}


betatotal <- rbind(betaprovince, betasite, betastation)
write.csv(betatotal, "outputs/07_diversity_partitioning/beta_RLS_families.csv")
beta_melt <- reshape2::melt(betatotal)

beta_family_RLS <- ggplot(beta_melt, aes(variable, value, colour=scale))+
  geom_boxplot()+
  theme_bw()+
  labs(x=" Beta component", y="Jaccard dissimilarity (RLS families)")

ggsave("outputs/07_diversity_partitioning/beta_RLS_families_Province.png")  
save(beta_family_RLS, file="Rdata/beta_family_RLS.rdata")


load("Rdata/beta_family_RLS.rdata")
load("Rdata/beta_species_RLS.rdata")
load("Rdata/beta_motus_eDNA.rdata")
load("Rdata/beta_family_eDNA.rdata")


ggarrange(beta_motus_eDNA, beta_family_eDNA, beta_species_RLS, beta_family_RLS, nrow=2, ncol=2, labels = c("a", "b", "c", "d"), label.y = c(1, 1, 1.1,1.1), widths = c(1,1.4))

ggsave("outputs/00_Figures_for_paper/Extended_Data/ED_Figure8.png")
