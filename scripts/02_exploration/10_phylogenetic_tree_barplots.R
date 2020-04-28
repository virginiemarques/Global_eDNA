library(tidyverse)
library(reshape2)
library(vegan)
library(ggplot2)

# code for the figure like de Vargas 2015

# load data to plot by the phylogenetic tree
load("Rdata/nb_motus_per_family_global.Rdata")
load("Rdata/nb_reads_per_family_global.Rdata")
family_coverage <- read.csv("outputs/01_read_data_stats/family_resolution_coefs.csv")

