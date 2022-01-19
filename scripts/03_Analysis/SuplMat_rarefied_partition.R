library(tidyverse)
library(reshape2)
library(vegan)
library(betapart)
library(ggpubr)

# Functions
"%ni%" <- Negate("%in%")

# Data - load rarefied data 
load("Rdata/rarefied_regions.rdata")
load("Rdata/rarefied_sites.rdata")


#-----------------------------------------------------------------------------------------------------------
# Partition on list rarefied_sites
#-----------------------------------------------------------------------------------------------------------

list_partition_rarefied_sites <- lapply(rarefied_sites, function(x){


# gamma global 

gamma_global <- as.numeric(x %>%
                             summarise(n = n_distinct(sequence)))


Province <- unique(x$province)
Site <- unique(x$site35)
Station <- unique(x$station)

# calculate alpha region

alpha_province=data.frame(province=character(), motu=numeric(), stringsAsFactors = FALSE)

for (i in 1:length(Province)) {
  r <- Province[i]
  motu <- x[x$province == Province[i],] %>%
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
  r <- unique(x[x$site35 == Site[i],]$province)
  motu <- x[x$site35 == Site[i],] %>%
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
  r <- unique(x[x$station == Station[i],]$province)
  s <- unique(x[x$station == Station[i],]$site35)
  motu <- x[x$station == Station[i],] %>%
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

return(div_partition)

})

partition_rarefied_sites <- bind_rows(list_partition_rarefied_sites)


mean_partition_site <- data.frame(component=unique(partition_rarefied_sites$component), mean=numeric(4), sd=numeric(4))
component <- unique(partition_rarefied_sites$component)
for (i in 1:length(component)) {
  df <- partition_rarefied_sites %>%
    filter(component==component[i])
  mean_partition_site[i,2] <- mean(df$percent)
  mean_partition_site[i,3] <- sd(df$percent)
  
}


#-----------------------------------------------------------------------------------------------------------
# Partition on list rarefied_regions
#-----------------------------------------------------------------------------------------------------------

list_partition_rarefied_regions <- lapply(rarefied_regions, function(x){
  
  
  # gamma global 
  
  gamma_global <- as.numeric(x %>%
                               summarise(n = n_distinct(sequence)))
  
  
  Province <- unique(x$province)
  Site <- unique(x$site35)
  Station <- unique(x$station)
  
  # calculate alpha region
  
  alpha_province=data.frame(province=character(), motu=numeric(), stringsAsFactors = FALSE)
  
  for (i in 1:length(Province)) {
    r <- Province[i]
    motu <- x[x$province == Province[i],] %>%
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
    r <- unique(x[x$site35 == Site[i],]$province)
    motu <- x[x$site35 == Site[i],] %>%
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
    r <- unique(x[x$station == Station[i],]$province)
    s <- unique(x[x$station == Station[i],]$site35)
    motu <- x[x$station == Station[i],] %>%
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
  
  return(div_partition)
  
})

partition_rarefied_regions <- bind_rows(list_partition_rarefied_regions)


mean_partition_region <- data.frame(component=unique(partition_rarefied_regions$component), mean=numeric(4), sd=numeric(4))
component <- unique(partition_rarefied_regions$component)
for (i in 1:length(component)) {
  df <- partition_rarefied_regions %>%
    filter(component==component[i])
  mean_partition_region[i,2] <- mean(df$percent)
  mean_partition_region[i,3] <- sd(df$percent)
  
}


write.csv(mean_partition_site, "outputs/07_diversity_partitioning/SUP_diversity_partition_rarefied_sites.csv", row.names = F)
write.csv(mean_partition_region, "outputs/07_diversity_partitioning/SUP_diversity_partition_rarefied_regions.csv", row.names = F)
