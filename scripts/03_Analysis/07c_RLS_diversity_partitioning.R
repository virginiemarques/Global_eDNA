library(tidyverse)
library(reshape2)
library(vegan)
library(betapart)
'%ni%' <- Negate("%in%")

# gamma global 

RLS_species <- read.csv("data/RLS/RLS_species_NEW.csv", sep = ";", stringsAsFactors = FALSE)
RLS_species <- RLS_species %>%
  filter(Realm%ni%c("Temperate Australasia", "Temperate Northern Atlantic", "Temperate Southern Africa", "Temperate Northern Pacific"))
RLS_species <- RLS_species[, c(1,2,6,10:2165)]
RLS_species <- reshape2::melt(RLS_species, id=c("SurveyID", "SiteCode", "Ecoregion"))
RLS_species <- RLS_species%>%
  filter(value!=0)
RLS_species <- RLS_species[,-5]
colnames(RLS_species) <- c("SurveyID", "SiteCode", "Ecoregion", "Species")

gamma_global <- as.numeric(RLS_species %>%
                             summarise(n = n_distinct(Species)))


region <- unique(RLS_species$Ecoregion)
site <- unique(RLS_species$SiteCode)
transect <- unique(RLS_species$SurveyID)

# calculate alpha region

alpha_region=data.frame(region=character(48), species=numeric(48), stringsAsFactors = FALSE)

for (i in 1:length(region)) {
  r <- region[i]
  species <- RLS_species[RLS_species$Ecoregion == region[i],] %>%
    summarise(n = n_distinct(Species))
  alpha_region[i,1] <- r
  alpha_region[i,2] <- species
}

# Calculate beta interregion
beta_region <- data.frame(alpha=mean(alpha_region$species), gamma=gamma_global, beta=numeric(1), scale="inter-region")
beta_region$beta <- beta_region$gamma - beta_region$alpha


# calculate alpha site

alpha_site=data.frame(region=character(), site=character(), species=numeric(), stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  s <- site[i]
  r <- unique(RLS_species[RLS_species$SiteCode == site[i],]$Ecoregion)
  species <- RLS_species[RLS_species$SiteCode == site[i],] %>%
    summarise(n = n_distinct(Species))
  alpha_site[i,1] <- r
  alpha_site[i,2] <- s
  alpha_site[i,3] <- species
}


# calculate beta inter-site

beta_site <- data.frame(region=character(48), alpha=numeric(48), gamma=numeric(48), beta=numeric(48), scale="inter-site", stringsAsFactors = FALSE)

for (i in 1:length(region)) {
  r <- region[i]
  gamma <- alpha_region[alpha_region$region==region[i],]$species
  alpha <- mean(alpha_site[alpha_site$region==region[i],]$species)
  beta_site[i,1] <- r
  beta_site[i,2] <- alpha
  beta_site[i,3] <- gamma
  beta_site[i,4] <- gamma - alpha
  
}

mean_beta_site <- mean(beta_site$beta) 
sd_beta_site <- sd(beta_site$beta)



# calculate alpha transect
alpha_transect=data.frame(region=character(), site=character(), transect=character(), species=numeric(), stringsAsFactors = FALSE)

for (i in 1:length(transect)) {
  st <- transect[i]
  r <- unique(RLS_species[RLS_species$SurveyID == transect[i],]$Ecoregion)
  s <- unique(RLS_species[RLS_species$SurveyID == transect[i],]$SiteCode)
  species <- RLS_species[RLS_species$SurveyID == transect[i],] %>%
    summarise(n = n_distinct(Species))
  alpha_transect[i,1] <- r
  alpha_transect[i,2] <- s
  alpha_transect[i,3] <- st
  alpha_transect[i,4] <- species
}

alpha_transect$species <- as.numeric(alpha_transect$species)
mean_a_transect <- data.frame(stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  df <- alpha_transect %>%
    subset(site==site[i])
  mean_a_transect[i,1] <- mean(df$species)
  
}


mean_alpha_transect <- mean(mean_a_transect$V1)
sd_alpha_transect <- sd(mean_a_transect$V1)

# calculate beta inter-station

beta_transect <- data.frame(region=character(1311), site=character(1311), alpha=numeric(1311), gamma=numeric(1311), beta=numeric(1311), scale="inter-transect", stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  s <- site[i]
  r <- unique(alpha_transect[alpha_transect$site==site[i],]$region)
  gamma <- alpha_site[alpha_site$site==site[i],]$species
  alpha <- mean(alpha_transect[alpha_transect$site==site[i],]$species)
  beta_transect[i,1] <- r
  beta_transect[i,2] <- s
  beta_transect[i,3] <- alpha
  beta_transect[i,4] <- gamma
  beta_transect[i,5] <- gamma - alpha
}


mean_b_transect <- data.frame(stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  df <- beta_transect %>%
    subset(site==site[i])
  mean_b_transect[i,1] <- mean(df$beta)
  
}


mean_beta_transect <- mean(mean_b_transect$V1)
sd_beta_transect <- sd(mean_b_transect$V1)

beta_region$beta+mean_beta_site+mean_beta_transect+mean_alpha_transect

# calculate diversity partitioning

div_partition <- data.frame(component=c("mean_alpha_transect", "mean_beta_transect", "mean_beta_site", "beta_region"), 
                            value=c(mean_alpha_transect, mean_beta_transect, mean_beta_site, beta_region$beta),
                            sd=c(sd_alpha_transect, sd_beta_transect, sd_beta_site, 0),
                            percent=numeric(4))
div_partition$percent <- (div_partition$value*100)/gamma_global


write.csv(div_partition, "outputs/07_diversity_partitioning/RLS_diversity_partitioning.csv")





#__________________________________________________________________________
### Beta On species ###
#__________________________________________________________________________
RLS_species <- read.csv("data/RLS/RLS_species_NEW.csv", sep = ";", stringsAsFactors = FALSE)
RLS_species <- RLS_species %>%
  subset(Realm%ni%c("Temperate Australasia", "Temperate Northern Atlantic", "Temperate Southern Africa", "Temperate Northern Pacific"))
RLS_sp <- RLS_species[,c(10:2165)]
RLS_sp <- RLS_sp[,colSums(RLS_sp)>0]
RLS_species <- cbind(RLS_species[,c(1,2,6)], RLS_sp)
save(RLS_species, file = "Rdata/RLS_species_clean.rdata")

region <- unique(RLS_species$Ecoregion)
site <- unique(RLS_species$SiteCode)
transect <- unique(RLS_species$SurveyID)

# calculate beta, turnover and nestedness with betapart
## beta inter-regions

df_region=data.frame()


for (i in 1:length(region)) {
  df <- RLS_species[RLS_species$Ecoregion == region[i],]
  df <- df[,c(4:ncol(df))]
  df_region[i,"region"] <- region[i]
  df_region[i, 2:1888] <- colSums(df)
}
rownames(df_region) <- df_region[,1]
df_region <- decostand(df_region[,2:ncol(df_region)], "pa",na.rm = TRUE)

b <- betapart.core(df_region)
beta <- beta.multi(b, "jaccard")
betaregion <- data.frame(scale="inter-region", total=beta$beta.JAC, turnover=beta$beta.JTU, nestedness=beta$beta.JNE)

  ## beta inter-site

df_site=vector("list", 48)
betasite <- data.frame(scale="inter-site", total=numeric(48), turnover=numeric(48), nestedness=numeric(48))


for (i in 1:length(region)) {
  df <- RLS_species[RLS_species$Ecoregion == region[i],]
  site <- unique(df$SiteCode)
  df_site[[i]] <- data.frame(stringsAsFactors = FALSE)
    for (j in 1:length(site)) {
      df2 <- df[df$SiteCode==site[j],]
      df2 <- df2[,c(4:ncol(df2))]
      df_site[[i]][j,"site"]<- site[j]
      df_site[[i]][j,2:1888] <- colSums(df2)
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

site <- unique(RLS_species$SiteCode)
df_transect=vector("list", 1311)
betatransect <- data.frame(scale="inter-transect", total=numeric(1311), turnover=numeric(1311), nestedness=numeric(1311))


for (i in 1:length(site)) {
  df <- RLS_species[RLS_species$SiteCode == site[i],]
  transect <- unique(df$SurveyID)
  df_transect[[i]] <- data.frame(stringsAsFactors = FALSE)
  for (j in 1:length(transect)) {
    df2 <- df[df$SurveyID==transect[j],]
    df2 <- df2[,c(4:ncol(df2))]
    df_transect[[i]][j, "transect"] <- transect[j]
    df_transect[[i]][j, 2:1888] <- colSums(df2)
  }
  rownames(df_transect[[i]]) <- df_transect[[i]][,1]
  df_transect[[i]] <- decostand(df_transect[[i]][,c(-1)], "pa",na.rm = TRUE)
  
  
  b <- betapart.core(df_transect[[i]])
  beta <- beta.multi(b, "jaccard")
  betatransect[i,2] <- beta$beta.JAC
  betatransect[i,3] <- beta$beta.JTU
  betatransect[i,4] <- beta$beta.JNE
}


betatotal <- rbind(betaregion, betasite, betatransect)
write.csv(betatotal, "outputs/07_diversity_partitioning/beta_RLS_species.csv")
beta_melt <- reshape2::melt(betatotal)

beta_species_RLS <- ggplot(beta_melt, aes(variable, value, colour=scale))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x=" Beta component", y="Jaccard dissimilarity (RLS Species)")

ggsave("outputs/07_diversity_partitioning/beta_RLS_species_Ecoregion.png")  
save(beta_species_RLS, file="Rdata/beta_species_RLS.rdata")


#__________________________________________________________________________
### Beta On families ###
#__________________________________________________________________________
RLS_families <- read.csv("data/RLS/RLS_families_NEW.csv", sep = ",", stringsAsFactors = FALSE)
RLS_families <- RLS_families %>%
  filter(Realm%ni%c("Temperate Australasia", "Temperate Northern Atlantic", "Temperate Southern Africa", "Temperate Northern Pacific"))
RLS_fam <- RLS_families[,c(10:130)]
RLS_fam <- RLS_fam[,colSums(RLS_fam)>0]
RLS_families <- cbind(RLS_families[,c(1,2,6)], RLS_fam)


region <- unique(RLS_families$Ecoregion)
site <- unique(RLS_families$SiteCode)
transect <- unique(RLS_families$SurveyID)

# calculate beta, turnover and nestedness with betapart
## beta inter-regions

df_region=data.frame()


for (i in 1:length(region)) {
  df <- RLS_families[RLS_families$Ecoregion == region[i],]
  df <- df[,c(4:ncol(df))]
  df_region[i,"region"] <- region[i]
  df_region[i, 2:97] <- colSums(df)
}
rownames(df_region) <- df_region[,1]
df_region <- decostand(df_region[,2:ncol(df_region)], "pa",na.rm = TRUE)

b <- betapart.core(df_region)
beta <- beta.multi(b, "jaccard")
betaregion <- data.frame(scale="inter-region", total=beta$beta.JAC, turnover=beta$beta.JTU, nestedness=beta$beta.JNE)

## beta inter-site

df_site=vector("list", 48)
betasite <- data.frame(scale="inter-site", total=numeric(48), turnover=numeric(48), nestedness=numeric(48))


for (i in 1:length(region)) {
  df <- RLS_families[RLS_families$Ecoregion == region[i],]
  site <- unique(df$SiteCode)
  df_site[[i]] <- data.frame(stringsAsFactors = FALSE)
  for (j in 1:length(site)) {
    df2 <- df[df$SiteCode==site[j],]
    df2 <- df2[,c(4:ncol(df2))]
    df_site[[i]][j,"site"]<- site[j]
    df_site[[i]][j,2:97] <- colSums(df2)
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

site <- unique(RLS_families$SiteCode)
df_transect=vector("list", 1311)
betatransect <- data.frame(scale="inter-transect", total=numeric(1311), turnover=numeric(1311), nestedness=numeric(1311))


for (i in 1:length(site)) {
  df <- RLS_families[RLS_families$SiteCode == site[i],]
  transect <- unique(df$SurveyID)
  df_transect[[i]] <- data.frame(stringsAsFactors = FALSE)
  for (j in 1:length(transect)) {
    df2 <- df[df$SurveyID==transect[j],]
    df2 <- df2[,c(4:ncol(df2))]
    df_transect[[i]][j, "transect"] <- transect[j]
    df_transect[[i]][j, 2:97] <- colSums(df2)
  }
  rownames(df_transect[[i]]) <- df_transect[[i]][,1]
  df_transect[[i]] <- decostand(df_transect[[i]][,c(-1)], "pa",na.rm = TRUE)
  
  
  b <- betapart.core(df_transect[[i]])
  beta <- beta.multi(b, "jaccard")
  betatransect[i,2] <- beta$beta.JAC
  betatransect[i,3] <- beta$beta.JTU
  betatransect[i,4] <- beta$beta.JNE
}


betatotal <- rbind(betaregion, betasite, betatransect)
write.csv(betatotal, "outputs/07_diversity_partitioning/beta_RLS_families.csv")
beta_melt <- reshape2::melt(betatotal)

beta_family_RLS <- ggplot(beta_melt, aes(variable, value, colour=scale))+
  geom_boxplot()+
  theme_bw()+
  labs(x=" Beta component", y="Jaccard dissimilarity (RLS families)")

ggsave("outputs/07_diversity_partitioning/beta_RLS_families_Ecoregion.png")  
save(beta_family_RLS, file="Rdata/beta_family_RLS.rdata")


load("Rdata/beta_family_RLS.rdata")
load("Rdata/beta_species_RLS.rdata")
load("Rdata/beta_motus_eDNA.rdata")
load("Rdata/beta_family_eDNA.rdata")


ggarrange(beta_motus_eDNA, beta_family_eDNA, beta_species_RLS, beta_family_RLS, nrow=2, ncol=2, labels = c("A", "B", "C", "D"), label.y = c(1, 1, 1.1,1.1), widths = c(1,1.4))

ggsave("outputs/00_Figures_for_paper/Extended_Data/ED_Figure8.png")
