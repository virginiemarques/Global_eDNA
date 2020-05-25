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

RLS_families <- read.csv("data/RLS_families.csv", sep = ";")
RLS_families <- RLS_families[,c(2,12:128)]

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

# method
method_accumulation = "random"

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

# plot
family_RLS <- ggplot(all_accumulation_RLS_df) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.5) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 1) +
  ylab("Number of families") +
  xlab("Number of transects") +
  theme_bw() + 
  ggtitle("Families") + 
  
family
