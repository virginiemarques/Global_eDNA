library(tidyverse)
library(reshape2)
library(vegan)
library(betapart)
'%ni%' <- Negate("%in%")

#__________________________________________________________________________
### On species ###
#__________________________________________________________________________
RLS_species <- read.csv("data/RLS/RLS_species.csv", sep = ";", stringsAsFactors = FALSE)
RLS_species <- RLS_species %>%
  filter(realm%ni%c("Temperate Australasia", "Temperate Northern Atlantic", "Temperate Southern Africa", "Temperate Northern Pacific"))
RLS_sp <- RLS_species[,c(12:2051)]
RLS_sp <- RLS_sp[,colSums(RLS_sp)>0]
RLS_species <- cbind(RLS_species[,c(1,2,9)], RLS_sp)


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
  df_region[i, 2:1787] <- colSums(df)
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
      df_site[[i]][j,2:1787] <- colSums(df2)
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
df_transect=vector("list", 1310)
betatransect <- data.frame(scale="inter-transect", total=numeric(1310), turnover=numeric(1310), nestedness=numeric(1310))


for (i in 1:length(site)) {
  df <- RLS_species[RLS_species$SiteCode == site[i],]
  transect <- unique(df$SurveyID)
  df_transect[[i]] <- data.frame(stringsAsFactors = FALSE)
  for (j in 1:length(transect)) {
    df2 <- df[df$SurveyID==transect[j],]
    df2 <- df2[,c(4:ncol(df2))]
    df_transect[[i]][j, "transect"] <- transect[j]
    df_transect[[i]][j, 2:1787] <- colSums(df2)
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

beta_melt <- reshape2::melt(betatotal)

ggplot(beta_melt, aes(variable, value, colour=scale))+
  geom_boxplot()+
  theme_bw()+
  labs(x=" Beta component", y="Jaccard dissimilarity (Species)")

ggsave("outputs/06_diversity_partitioning/beta_RLS_species_Ecoregion.png")  



#__________________________________________________________________________
### On families ###
#__________________________________________________________________________
RLS_families <- read.csv("data/RLS/RLS_families.csv", sep = ";", stringsAsFactors = FALSE)
RLS_families <- RLS_families %>%
  filter(realm%ni%c("Temperate Australasia", "Temperate Northern Atlantic", "Temperate Southern Africa", "Temperate Northern Pacific"))
RLS_fam <- RLS_families[,c(12:128)]
RLS_fam <- RLS_fam[,colSums(RLS_fam)>0]
RLS_families <- cbind(RLS_families[,c(1,2,9)], RLS_fam)


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
  df_region[i, 2:95] <- colSums(df)
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
    df_site[[i]][j,2:95] <- colSums(df2)
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
df_transect=vector("list", 1310)
betatransect <- data.frame(scale="inter-transect", total=numeric(1310), turnover=numeric(1310), nestedness=numeric(1310))


for (i in 1:length(site)) {
  df <- RLS_families[RLS_families$SiteCode == site[i],]
  transect <- unique(df$SurveyID)
  df_transect[[i]] <- data.frame(stringsAsFactors = FALSE)
  for (j in 1:length(transect)) {
    df2 <- df[df$SurveyID==transect[j],]
    df2 <- df2[,c(4:ncol(df2))]
    df_transect[[i]][j, "transect"] <- transect[j]
    df_transect[[i]][j, 2:95] <- colSums(df2)
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

beta_melt <- reshape2::melt(betatotal)

ggplot(beta_melt, aes(variable, value, colour=scale))+
  geom_boxplot()+
  theme_bw()+
  labs(x=" Beta component", y="Jaccard dissimilarity (families)")

ggsave("outputs/06_diversity_partitioning/beta_RLS_families_Ecoregion.png")  
