# Accumulation curves on genus, family, order on RLS data

# Lib 
library(tidyverse)
library(reshape)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(vegan)
library(ggthemes)


# 
'%ni%' <- Negate("%in%")

# Functions
source('scripts/03_Analysis/00_functions.R')


# Functions 
akaike.weights <- function(x)
{
  x <- x[!is.na(x)]
  delta.aic <- x - min(x, na.rm = TRUE)
  rel.LL <- exp(-0.5 * delta.aic)
  sum.LL <- sum(rel.LL, na.rm = TRUE)
  weights.aic <- rel.LL/sum.LL
  return(list(deltaAIC = delta.aic, rel.LL = rel.LL, weights = weights.aic))
}

# Load data RLS

RLS_species <- read.csv("data/RLS/RLS_species_NEW.csv", sep = ";", stringsAsFactors = FALSE)

list_region_RLS <- split(RLS_species, RLS_species$Province)
lapply(list_region_RLS, function(x){
  length(unique(x$SurveyID))
})




# ------------------------------------------------------------------------------- # 
#### On family ----
# ------------------------------------------------------------------------------- # 

RLS_families <- read.csv("data/RLS/RLS_families_NEW.csv", sep = ",", stringsAsFactors = FALSE)
RLS_fam <- RLS_families[,c(11:131)]
RLS_fam <- RLS_fam[,colSums(RLS_fam)>0]
RLS_families <- cbind(RLS_families$SurveyID, RLS_fam)
# method
method_accumulation = "exact"

# Add a global saturation curve for RLS
all_accumulation_RLS <- specaccum(RLS_families[,2:ncol(RLS_families)], 
                                  method = method_accumulation, permutations = 100,
                                  conditioned =TRUE, gamma = "jack1")

# Accumulation for plot 
all_accumulation_RLS_df <- data.frame(richness = all_accumulation_RLS$richness, sd = all_accumulation_RLS$sd, sites = all_accumulation_RLS$sites)

# Asymptote of RLS  


# models asymptotes
lomo_edna <- fitspecaccum(all_accumulation_RLS, "lomolino")
aic_lomo_edna <-AIC(lomo_edna)
mm_edna <- fitspecaccum(all_accumulation_RLS, "michaelis-menten")
aic_mm_edna <- AIC(mm_edna)
gom_edna <- fitspecaccum(all_accumulation_RLS, "gompertz")
aic_gom_edna <- AIC(gom_edna)
asy_edna <- fitspecaccum(all_accumulation_RLS, "asymp")
aic_asy_edna <- AIC(asy_edna)
gis_edna <- fitspecaccum(all_accumulation_RLS, "logis")
aic_gis_edna <- AIC(gis_edna)

# compute results
res_edna <- matrix(NA,nrow = 5, ncol = 3)
rownames(res_edna) <- c("lomolino", "michaelis-menten", "gompertz", "asymp", "logis")
colnames(res_edna) <- c("AIC", "Asymptote", "Weigth")
res_edna[,"AIC"] <- c(aic_lomo_edna, aic_mm_edna, aic_gom_edna, aic_asy_edna, aic_gis_edna)
res_edna[,"Weigth"] <- akaike.weights(c(aic_lomo_edna, aic_mm_edna, aic_gom_edna, aic_asy_edna, aic_gis_edna))$weights
res_edna[,"Asymptote"] <- c(coef(lomo_edna)[[1]], coef(mm_edna)[[1]], coef(gom_edna)[[1]], coef(asy_edna)[[1]], coef(gis_edna)[[1]])

# Calculation asymptote value
asymp_RLS <- weighted.mean(res_edna[,"Asymptote"], res_edna[,"Weigth"])

# Add to DF
all_accumulation_RLS_df$asymptote <- asymp_RLS
save(all_accumulation_RLS_df, file = "Rdata/accumulation_families_RLS.rdata")

# plot RLS
family_RLS <- ggplot(all_accumulation_RLS_df) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.5) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 1, alpha=0.7) +
  annotate(geom="text", x=2990, y=118+10, label="Asymptote : 118",hjust=1, alpha=0.7, size=3.5) +
  annotate(geom="text", x=2990, y=96+10, label="RLS Family : 96",hjust=1, alpha=0.7, size=3.5)+
  annotate(geom="text", x=2990, y=30, label="Slope = 1.44",hjust=1, alpha=0.7)+
  ylim(0,190)+
  xlab("Number of transects") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill = NA)) + 
  ggtitle("d")

family_RLS


#combined plots


# ------------------------------------------------------------------------------- # 
#### On species ----
# ------------------------------------------------------------------------------- # 

RLS_species <- read.csv("data/RLS/RLS_species_NEW.csv", sep = ";", stringsAsFactors = FALSE, check.names = FALSE)
RLS_sp <- RLS_species[,c(12:2167)]
RLS_sp <- RLS_sp[,colSums(RLS_sp)>0]
RLS_species <- cbind(RLS_species$SurveyID, RLS_sp)

# method
method_accumulation = "exact"

# Add a global saturation curve for RLS
all_accumulation_species_RLS <- specaccum(RLS_species[,2:ncol(RLS_species)], 
                                          method = method_accumulation, permutations = 100,
                                          conditioned =TRUE, gamma = "jack1")

# Accumulation for plot 
all_accumulation_species_RLS_df <- data.frame(richness = all_accumulation_species_RLS$richness, sd = all_accumulation_species_RLS$sd, sites = all_accumulation_species_RLS$sites)


# Asymptote of RLS  


# models asymptotes
lomo_edna <- fitspecaccum(all_accumulation_species_RLS, "lomolino")
aic_lomo_edna <-AIC(lomo_edna)
mm_edna <- fitspecaccum(all_accumulation_species_RLS, "michaelis-menten")
aic_mm_edna <- AIC(mm_edna)
gom_edna <- fitspecaccum(all_accumulation_species_RLS, "gompertz")
aic_gom_edna <- AIC(gom_edna)
asy_edna <- fitspecaccum(all_accumulation_species_RLS, "asymp")
aic_asy_edna <- AIC(asy_edna)
gis_edna <- fitspecaccum(all_accumulation_species_RLS, "logis")
aic_gis_edna <- AIC(gis_edna)

# compute results
res_edna <- matrix(NA,nrow = 5, ncol = 3)
rownames(res_edna) <- c("lomolino", "michaelis-menten", "gompertz", "asymp", "logis")
colnames(res_edna) <- c("AIC", "Asymptote", "Weigth")
res_edna[,"AIC"] <- c(aic_lomo_edna, aic_mm_edna, aic_gom_edna, aic_asy_edna, aic_gis_edna)
res_edna[,"Weigth"] <- akaike.weights(c(aic_lomo_edna, aic_mm_edna, aic_gom_edna, aic_asy_edna, aic_gis_edna))$weights
res_edna[,"Asymptote"] <- c(coef(lomo_edna)[[1]], coef(mm_edna)[[1]], coef(gom_edna)[[1]], coef(asy_edna)[[1]], coef(gis_edna)[[1]])

# Calculation asymptote value
asymp_species_RLS <- weighted.mean(res_edna[,"Asymptote"], res_edna[,"Weigth"])

# Add to DF
all_accumulation_species_RLS_df$asymptote <- asymp_species_RLS

# Save
save(all_accumulation_species_RLS_df, file = "Rdata/accumulation_species_RLS.rdata")

# plot
species_RLS <- ggplot(all_accumulation_species_RLS_df) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.5) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 1, alpha=0.7) +
  annotate(geom="text", x=2990, y=2322+150, label="Asymptote : 2322",hjust=1, alpha=0.7, size=3.5) +
  annotate(geom="text", x=2990, y=1786+150, label="RLS Species : 1786",hjust=1, alpha=0.7, size=3.5) +
  annotate(geom="text", x=2990, y=600, label="Slope = 1.68",hjust=1, alpha=0.7) +
  ylim(0,3000)+
  xlab("Number of transects") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill = NA)) + 
  ggtitle("b")

species_RLS


