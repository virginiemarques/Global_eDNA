library(tidyverse)
library(reshape2)
library(vegan)
library(betapart)

setwd("c:/Users/mathon/Desktop/linux/Global_eDNA/")
load("Rdata/02_clean_all.Rdata")
'%ni%' <- Negate("%in%")
#Remove estuary stations and deep niskin station
df_all_filters <- df_all_filters %>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")) %>%
  filter(sample_method !="niskin" & region!="East_Pacific" & comment %ni% c("Distance decay 600m", "Distance decay 300m") & station!="glorieuse_distance_300m")%>%
  filter(project != "SEAMOUNTS") %>% 
  filter(habitat_type %ni% c("BAIE", "Sommet"))


df_all_filters <- df_all_filters %>%
  filter(!is.na(new_family_name))


# gamma global =143

gamma_global <- as.numeric(df_all_filters %>%
  summarise(n = n_distinct(new_family_name)))
  

Region <- unique(df_all_filters$region)
Site <- unique(df_all_filters$site)
Station <- unique(df_all_filters$station)


# calculate alpha diversity per region

alpha_region=data.frame(region=character(5), family=numeric(5), stringsAsFactors = FALSE)

for (i in 1:length(Region)) {
  r <- Region[i]
  family <- df_all_filters[df_all_filters$region == Region[i],] %>%
    summarise(n = n_distinct(new_family_name))
  alpha_region[i,1] <- r
  alpha_region[i,2] <- family
}


# calculate beta inter-region
beta_region <- data.frame(alpha=mean(alpha_region$family), gamma=gamma_global, beta=numeric(1), scale="inter-region")
beta_region$beta <- beta_region$gamma - beta_region$alpha


# calculate alpha diversity per site

alpha_site=data.frame(region=character(), site=character(), family=numeric(), stringsAsFactors = FALSE)

for (i in 1:length(Site)) {
  s <- Site[i]
  r <- unique(df_all_filters[df_all_filters$site == Site[i],]$region)
  family <- df_all_filters[df_all_filters$site == Site[i],] %>%
    summarise(n = n_distinct(new_family_name))
  alpha_site[i,1] <- r
  alpha_site[i,2] <- s
  alpha_site[i,3] <- family
}


# calculate beta inter-site

beta_site <- data.frame(region=character(5), alpha=numeric(5), gamma=numeric(5), beta=numeric(5), scale="inter-site", stringsAsFactors = FALSE)

for (i in 1:length(Region)) {
  r <- Region[i]
  gamma <- alpha_region[alpha_region$region==Region[i],]$family
  alpha <- mean(alpha_site[alpha_site$region==Region[i],]$family)
  beta_site[i,1] <- r
  beta_site[i,2] <- alpha
  beta_site[i,3] <- gamma
  beta_site[i,4] <- gamma - alpha
  
}

mean_beta_site <- mean(beta_site$beta)
sd_beta_site <- sd(beta_site$beta)



# calculate alpha diversity per station
alpha_station=data.frame(region=character(), site=character(), station=character(), family=numeric(), stringsAsFactors = FALSE)

for (i in 1:length(Station)) {
  st <- station[i]
  r <- unique(df_all_filters[df_all_filters$station == station[i],]$region)
  s <- unique(df_all_filters[df_all_filters$station == station[i],]$site)
  family <- df_all_filters[df_all_filters$station == station[i],] %>%
    summarise(n = n_distinct(new_family_name))
  alpha_station[i,1] <- r
  alpha_station[i,2] <- s
  alpha_station[i,3] <- st
  alpha_station[i,4] <- family
}

mean_a_station <- data.frame(stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  df <- alpha_station %>%
    subset(site==site[i])
  mean_a_station[i,1] <- mean(df$motu)
  
}

mean_alpha_station <- mean(mean_a_station)
sd_alpha_station <- sd(mean_a_station)


# calculate beta inter-station

beta_station <- data.frame(region=character(25), site=character(25), alpha=numeric(25), gamma=numeric(25), beta=numeric(25), scale="inter-station", stringsAsFactors = FALSE)

for (i in 1:length(Site)) {
  s <- Site[i]
  r <- unique(alpha_station[alpha_station$site==Site[i],]$region)
  gamma <- alpha_site[alpha_site$site==Site[i],]$family
  alpha <- mean(alpha_station[alpha_station$site==Site[i],]$family)
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
                            value= c(mean_alpha_station, mean_beta_station, mean_beta_site, beta_region$beta), 
                            sd=c(sd_alpha_station, sd_beta_station, sd_beta_site, 0),
                            percent=numeric(4))
div_partition$percent <- (div_partition$value*100)/gamma_global

write.csv(div_partition, "outputs/06_diversity_partitioning/Family_diversity_partitioning.csv")


# plot alpha and beta ~ gamma at each spatial scale

alpha_beta <- rbind(beta_station[,c(-1)], beta_site[,c(-1)], beta_region)


ggplot(alpha_beta, aes(colour=scale))+
  geom_point(data=beta_region, aes(gamma, alpha))+
  geom_point(data=beta_region, aes(gamma, beta))+
  geom_smooth(method=lm, se=FALSE, aes(x=gamma, y=alpha))+
  geom_smooth(method=lm,  linetype="longdash", se=FALSE, aes(x=gamma, y=beta))+
  geom_abline(intercept = 0, slope=1)+
  geom_abline(intercept = 0, slope=0.5, linetype="dotted")+
  ylim(0,130)+
  xlim(0,130)+
  labs(x="Regional family richness", y="Alpha and Beta richness")

ggsave("outputs/06_diversity_partitioning/family_alpha_beta~gamma.png")

# calculate beta, turnover and nestedness with betapart
  ## beta inter-regions

df_region=data.frame(family=character())

for (i in 1:length(Region)) {
  df <- df_all_filters[df_all_filters$region == Region[i],] %>%
    distinct(new_family_name, region)
  colnames(df) <- c("family", Region[i])
  df_region <- full_join(df_region, df, by="family")
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
  df_site[[i]] <- data.frame(family=character(), stringsAsFactors = FALSE)
    for (j in 1:length(Site)) {
      df2 <- df[df$site==Site[j],] %>%
        distinct(new_family_name, site)
      colnames(df2) <- c("family", Site[j])
      df_site[[i]] <- full_join(df_site[[i]], df2, by="family")
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
  df_station[[i]] <- data.frame(family=character(), stringsAsFactors = FALSE)
  for (j in 1:length(Station)) {
    df2 <- df[df$station==Station[j],] %>%
      distinct(new_family_name, station)
    colnames(df2) <- c("family", Station[j])
    df_station[[i]] <- full_join(df_station[[i]], df2, by="family")
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
write.csv(betatotal, "outputs/06_diversity_partitioning/beta_families.csv")

beta_melt <- reshape2::melt(betatotal)

beta_family_eDNA <- ggplot(beta_melt, aes(variable, value, colour=scale))+
  geom_boxplot()+
  ylim(0,1)+
  theme_bw()+
  labs(x=" Beta component", y="Jaccard dissimilarity (family)")

ggsave("outputs/06_diversity_partitioning/family_diversity_partitioning.png")  
save(beta_family_eDNA, file = "Rdata/beta_family_eDNA.rdata")
