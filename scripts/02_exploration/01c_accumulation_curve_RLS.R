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


# ------------------------------------------------------------------------------- # 
#### On family ----
# ------------------------------------------------------------------------------- # 

# rank_specify
rank_choice = 'family'
sample = "SiteCode"

# Add a global saturation curve, i.e. all samples together?
all_accumulation_RLS <- accumulation_curve_df(RLS_families, species_unit = rank_choice, column_station = sample) %>%
  mutate(project_name = "All") %>%
  select(project_name, richness, sd, sites)

save(all_accumulation_RLS, file = "Rdata/accumulation_families_RLS.rdata")

# Asymptote of all plots 
all_asymptote_RLS <- asymptote_mm(RLS_families, species_unit = rank_choice, column_station = sample) %>%
  mutate(project_name = "All") %>%
  select(project_name, asymptote)

save(all_asymptote_RLS, file = "Rdata/asymptote_families_RLS.rdata")

all_accumulation_RLS$asymptote <- all_asymptote_RLS$asymptote

# plot
family_RLS <- ggplot(all_accumulation_RLS) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.5) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 1) +
  ylab("Number of families") +
  xlab("Number of transects") +
  theme_bw() + 
  ggtitle("Families") + 
  
family
