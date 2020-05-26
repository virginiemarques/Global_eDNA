# Accumulation curves on genus, family, order 
# Lib 
library(tidyverse)
library(reshape)
library(ggplot2)
library(ggpubr)
library(dplyr)

# data
load("Rdata/02_clean_all.Rdata")

# 
'%ni%' <- Negate("%in%")

# Functions
source('scripts/02_exploration/00_functions.R')


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

# ------------------------------------------------------------------------------- # 
#### On family ----
# ------------------------------------------------------------------------------- # 

RLS_families <- read.csv("data/RLS_families.csv", sep = ";")
RLS_families <- RLS_families[,c(2,12:128)]

# method
method_accumulation = "exact"

# Add a global saturation curve for RLS
all_accumulation_RLS <- specaccum(RLS_families[,2:ncol(RLS_families)], 
                      method = method_accumulation, permutations = 100,
                      conditioned =TRUE, gamma = "jack1")

# Accumulation for plot 
all_accumulation_RLS_df <- data.frame(richness = all_accumulation_RLS$richness, sd = all_accumulation_RLS$sd, sites = all_accumulation_RLS$sites)

# Save
save(all_accumulation_RLS_df, file = "Rdata/accumulation_families_RLS.rdata")

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

# plot RLS
family_RLS <- ggplot(all_accumulation_RLS_df) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.5) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 1, alpha=0.7) +
  annotate(geom="text", x=2990, y=160+10, label="RLS Family : 160",hjust=1, alpha=0.7) +
  ylim(0,190)+
  ylab("Number of families") +
  xlab("Number of transects") +
  theme_bw() + 
  ggtitle("C")
  
family_RLS


# plot combined with our data
load("Rdata/accumulation_asymptote_families_all.rdata")

fam_edna <- ggplot(df_fam) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.8, fill="#457277") +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "solid", size = 1, col="#457277") +
  annotate(geom="text", x=250, y=178+10, label="eDNA Family : 178",hjust=1,color="#457277") +
  ylim(0,190)+
  labs(y="", x="Number of samples")+
  theme_bw() + 
  ggtitle("D")

fam_edna

#combined plots

ggarrange(family_RLS, fam_edna +rremove("ylab"), nrow = 1, ncol=2)
ggsave("outputs/03_accumulation_curves/accumulation_families_panels.png")
# ------------------------------------------------------------------------------- # 
#### On species ----
# ------------------------------------------------------------------------------- # 

RLS_species <- read.csv("data/RLS_species.csv", sep = ";")
RLS_species <- RLS_species[,c(2,12:2051)]

# method
method_accumulation = "exact"

# Add a global saturation curve for RLS
all_accumulation_species_RLS <- specaccum(RLS_species[,2:ncol(RLS_species)], 
                                  method = method_accumulation, permutations = 100,
                                  conditioned =TRUE, gamma = "jack1")

# Accumulation for plot 
all_accumulation_species_RLS_df <- data.frame(richness = all_accumulation_species_RLS$richness, sd = all_accumulation_species_RLS$sd, sites = all_accumulation_species_RLS$sites)

# Save
save(all_accumulation_species_RLS_df, file = "Rdata/accumulation_species_RLS.rdata")

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
res_edna <- matrix(NA,nrow = 4, ncol = 3)
rownames(res_edna) <- c("lomolino", "michaelis-menten", "asymp", "logis")
colnames(res_edna) <- c("AIC", "Asymptote", "Weigth")
res_edna[,"AIC"] <- c(aic_lomo_edna, aic_mm_edna, aic_asy_edna, aic_gis_edna)
res_edna[,"Weigth"] <- akaike.weights(c(aic_lomo_edna, aic_mm_edna, aic_asy_edna, aic_gis_edna))$weights
res_edna[,"Asymptote"] <- c(coef(lomo_edna)[[1]], coef(mm_edna)[[1]], coef(asy_edna)[[1]], coef(gis_edna)[[1]])

# Calculation asymptote value
asymp_species_RLS <- weighted.mean(res_edna[,"Asymptote"], res_edna[,"Weigth"])

# Add to DF
all_accumulation_species_RLS_df$asymptote <- asymp_species_RLS

# plot
species_RLS <- ggplot(all_accumulation_species_RLS_df) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.5) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 1, alpha=0.7) +
  annotate(geom="text", x=2990, y=2733+150, label="RLS Species : 2733",hjust=1, alpha=0.7) +
  ylim(0,3000)+
  ylab("Number of Species / MOTUs") +
  xlab("Number of transects") +
  theme_bw() + 
  ggtitle("A")

species_RLS


# plot combined with our data
load("Rdata/accumulation_asymptote_motus_all.rdata")

species_edna <- ggplot(df_motus) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.8, fill="#d2981a") +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "solid", size = 1, col="#d2981a") +
  annotate(geom="text", x=250, y=2818+150, label="eDNA MOTUs : 2818",hjust=1,
           color="#d2981a") +
  ylim(0,3000)+
  ylab("") +
  xlab("Number of samples")+
  theme_bw() + 
  ggtitle("B")

species_edna

ggarrange(species_RLS, species_edna, nrow = 1, ncol=2)
ggsave("outputs/03_accumulation_curves/accumulation_species_motus_panels.png")
