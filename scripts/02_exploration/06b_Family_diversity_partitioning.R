library(tidyverse)
library(reshape2)
library(vegan)
library(betapart)

setwd("c:/Users/mathon/Desktop/linux/Global_eDNA/")
load("Rdata/02_clean_all.Rdata")

#Remove estuary stations and deep niskin station
df_all_filters <- subset(df_all_filters, !(station %in% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")))
df_all_filters <- subset(df_all_filters, sample_method!="niskin")
df_all_filters <- subset(df_all_filters, region!="East_Pacific")
df_all_filters <- subset(df_all_filters, !(comment %in% c("Distance decay 600m", "Distance decay 300m")))
df_all_filters <- subset(df_all_filters, station!="glorieuse_distance_300m")


df_all_filters <- df_all_filters %>%
  filter(!is.na(new_family_name))


# gamma global =130

gamma_global <- as.numeric(df_all_filters %>%
  summarise(n = n_distinct(new_family_name)))
  

region <- unique(df_all_filters$region)
site <- unique(df_all_filters$site)
station <- unique(df_all_filters$station)


# calculate alpha diversity per region

alpha_region=data.frame(region=character(5), family=numeric(5), stringsAsFactors = FALSE)

for (i in 1:length(region)) {
  r <- region[i]
  family <- df_all_filters[df_all_filters$region == region[i],] %>%
    summarise(n = n_distinct(new_family_name))
  alpha_region[i,1] <- r
  alpha_region[i,2] <- family
}


# calculate beta inter-region
beta_region <- data.frame(alpha=mean(alpha_region$family), gamma=gamma_global, beta=numeric(1), scale="inter-region")
beta_region$beta <- beta_region$gamma - beta_region$alpha


# calculate alpha diversity per site

alpha_site=data.frame(region=character(), site=character(), family=numeric(), stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  s <- site[i]
  r <- unique(df_all_filters[df_all_filters$site == site[i],]$region)
  family <- df_all_filters[df_all_filters$site == site[i],] %>%
    summarise(n = n_distinct(new_family_name))
  alpha_site[i,1] <- r
  alpha_site[i,2] <- s
  alpha_site[i,3] <- family
}


# calculate beta inter-site

beta_site <- data.frame(region=character(5), alpha=numeric(5), gamma=numeric(5), beta=numeric(5), scale="inter-site", stringsAsFactors = FALSE)

for (i in 1:length(region)) {
  r <- region[i]
  gamma <- alpha_region[alpha_region$region==region[i],]$family
  alpha <- mean(alpha_site[alpha_site$region==region[i],]$family)
  beta_site[i,1] <- r
  beta_site[i,2] <- alpha
  beta_site[i,3] <- gamma
  beta_site[i,4] <- gamma - alpha
  
}

mean_beta_site <- mean(beta_site$beta)
sd_beta_site <- sd(beta_site$beta)



# calculate alpha diversity per station
alpha_station=data.frame(region=character(), site=character(), station=character(), family=numeric(), stringsAsFactors = FALSE)

for (i in 1:length(station)) {
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

alpha_station_car <- alpha_station %>%
  subset(region == "Caribbean")
mean_alpha_station_car <- mean(alpha_station_car$family)

alpha_station_leng <- alpha_station %>%
  subset(region == "Central_IndoPacific")
mean_alpha_station_leng <- mean(alpha_station_leng$family)

alpha_station_faka <- alpha_station %>%
  subset(region == "Central_Pacific")
mean_alpha_station_faka <- mean(alpha_station_faka$family)

alpha_station_eparse <- alpha_station %>%
  subset(region == "West_Indian")
mean_alpha_station_eparse <- mean(alpha_station_eparse$family)

alpha_station_cal <- alpha_station %>%
  subset(region == "South_West_Pacific")
mean_alpha_station_cal <- mean(alpha_station_cal$family)

a_station <- c(mean_alpha_station_car, mean_alpha_station_faka, mean_alpha_station_leng, mean_alpha_station_eparse, mean_alpha_station_cal)
mean_alpha_station <- mean(a_station)
sd_alpha_station <- sd(a_station)


# calculate beta inter-station

beta_station <- data.frame(region=character(25), site=character(25), alpha=numeric(25), gamma=numeric(25), beta=numeric(25), scale="inter-station", stringsAsFactors = FALSE)

for (i in 1:length(site)) {
  s <- site[i]
  r <- unique(alpha_station[alpha_station$site==site[i],]$region)
  gamma <- alpha_site[alpha_site$site==site[i],]$family
  alpha <- mean(alpha_station[alpha_station$site==site[i],]$family)
  beta_station[i,1] <- r
  beta_station[i,2] <- s
  beta_station[i,3] <- alpha
  beta_station[i,4] <- gamma
  beta_station[i,5] <- gamma - alpha
}

beta_station_car <- beta_station %>%
  subset(region == "Caribbean")
mean_beta_station_car <- mean(beta_station_car$beta)

beta_station_leng <- beta_station %>%
  subset(region == "Central_IndoPacific")
mean_beta_station_leng <- mean(beta_station_leng$beta)

beta_station_eparse <- beta_station %>%
  subset(region == "West_Indian")
mean_beta_station_eparse <- mean(beta_station_eparse$beta)

beta_station_cal <- beta_station %>%
  subset(region == "South_West_Pacific")
mean_beta_station_cal <- mean(beta_station_cal$beta)

beta_station_faka <- beta_station[beta_station$region=="Central_Pacific",]$beta

b_station <- c(mean_beta_station_car, beta_station_faka, mean_beta_station_leng, mean_beta_station_eparse, mean_beta_station_cal)
mean_beta_station <- mean(b_station)
sd_beta_station <- sd(b_station)

beta_region$beta+mean_beta_site+mean_beta_station+mean_alpha_station

# calculate diversity partitioning

div_partition <- data.frame(component=c("mean_alpha_station", "mean_beta_station", "mean_beta_site", "beta_region"), 
                            value= c(mean_alpha_station, mean_beta_station, mean_beta_site, beta_region$beta), 
                            sd=c(sd_alpha_station, sd_beta_station, sd_beta_site, 0),
                            percent=numeric(5))
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

for (i in 1:length(region)) {
  df <- df_all_filters[df_all_filters$region == region[i],] %>%
    distinct(new_family_name, region)
  colnames(df) <- c("family", region[i])
  df_region <- full_join(df_region, df, by="family")
}
rownames(df_region) <- df_region[,1]
df_region <- decostand(df_region[,2:5], "pa",na.rm = TRUE)
df_region[is.na(df_region)] <- 0
df_region <- as.data.frame(t(df_region))

b <- betapart.core(df_region)
beta <- beta.multi(b, "jaccard")
betaregion <- data.frame(scale="inter-region", total=beta$beta.JAC, turnover=beta$beta.JTU, nestedness=beta$beta.JNE)

  ## beta inter-site

df_site=vector("list", 5)
betasite <- data.frame(scale="inter-site", total=numeric(5), turnover=numeric(5), nestedness=numeric(5))


for (i in 1:length(region)) {
  df <- df_all_filters[df_all_filters$region == region[i],]
  site <- unique(df$site)
  df_site[[i]] <- data.frame(family=character(), stringsAsFactors = FALSE)
    for (j in 1:length(site)) {
      df2 <- df[df$site==site[j],] %>%
        distinct(new_family_name, site)
      colnames(df2) <- c("family", site[j])
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

site <- unique(df_all_filters$site)
df_station=vector("list", 25)
betastation <- data.frame(scale="inter-station", total=numeric(25), turnover=numeric(25), nestedness=numeric(25))


for (i in 1:length(site)) {
  df <- df_all_filters[df_all_filters$site == site[i],]
  station <- unique(df$station)
  df_station[[i]] <- data.frame(family=character(), stringsAsFactors = FALSE)
  for (j in 1:length(station)) {
    df2 <- df[df$station==station[j],] %>%
      distinct(new_family_name, station)
    colnames(df2) <- c("family", station[j])
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

beta_melt <- melt(betatotal)

ggplot(beta_melt, aes(variable, value, colour=scale))+
  geom_boxplot()+
  ylim(0,1)+
  labs(x=" Beta component", y="Jaccard dissimilarity (family)")

ggsave("outputs/06_diversity_partitioning/family_diversity_partitioning.png")  
