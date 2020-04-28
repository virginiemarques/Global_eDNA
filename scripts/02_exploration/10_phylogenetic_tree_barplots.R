library(tidyverse)
library(reshape2)
library(vegan)
library(ggplot2)

# code for the figure like de Vargas 2015
load("Rdata/02_clean_all.Rdata")
df_all_filters <- subset(df_all_filters, !(station %in% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")))
df_all_filters <- subset(df_all_filters, sample_method!="niskin")
df_all_filters <- subset(df_all_filters, region!="East_Pacific")


# load data to plot by the phylogenetic tree

load("Rdata/nb_motus_per_family_global.Rdata")
load("Rdata/nb_reads_per_family_global.Rdata")
family_coverage <- read.csv("outputs/01_read_data_stats/family_resolution_coefs.csv")

family <- unique(df_all_filters$new_family_name)
family_coverage <- subset(family_coverage, Family%in%family)
