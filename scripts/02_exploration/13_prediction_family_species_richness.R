library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(rfishbase)
library(grid)

## build figure prediction of number of species (Fig 1C)
'%ni%' <- Negate("%in%")

# Functions
source('scripts/02_exploration/00_functions.R')


# family asymptote  + n_motus data
asymptote_fam <- read.csv("outputs/10_acculation_curve_family/table_asymptote_family20.csv")
asymptote_fam$fam <- c("Acnt", "Apgn", "Blst", "Crng", "Gobd", "Hlcn", "Lbrd", "Lthr", "Ltjn", "Mrnd", "Mctp", "Pmcnth", "Pmcntr", "Srrn")

# coral fish species data
coral_fishes <- read.csv("data/Coral_fishes2.csv", sep=";")

# our data
df_all_filters <- df_all_filters %>%
  #filter(new_rank_ncbi != "higher") %>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")) %>%
  filter(sample_method !="niskin" & region!="East_Pacific" & comment %ni% c("Distance decay 600m", "Distance decay 300m") & station!="glorieuse_distance_300m")


## select abundant families
fam_summary <- df_all_filters %>%
  filter(!is.na(new_family_name)) %>%
  group_by(new_family_name) %>%
  summarise(n_motus = n_distinct(sequence), 
            n_samples = n_distinct(sample_name_all_pcr),
            n_station = n_distinct(station))

families20 <- fam_summary %>%
  filter(n_motus > 20 & n_station > 20) %>%
  distinct(new_family_name) %>% pull()

# calculate number of species per family on the checklist
for (i in 1:length(families20)) {
  asymptote_fam[i,"n_species_checklist"] <- coral_fishes %>%
    subset(Family==families20[i]) %>%
    summarize(n_distinct(Species))
}

# calculate number of taxa per family in our data

for (i in 1:length(families20)) {
  asymptote_fam[i,"n_species"] <- df_all_filters %>%
    subset(new_family_name==families20[i]) %>%
    filter(!is.na(new_scientific_name_ncbi)) %>%
    summarize(n_distinct(new_scientific_name_ncbi))
}


for (i in 1:length(families20)) {
  asymptote_fam[i, "delta_taxa"] <- ((asymptote_fam[i, "n_species_checklist"] - asymptote_fam[i, "n_species"])*100)/asymptote_fam[i, "n_species_checklist"]
  asymptote_fam[i, "delta_motu"] <- ((asymptote_fam[i, "n_species_checklist"] - asymptote_fam[i, "n_motus"])*100)/asymptote_fam[i, "n_species_checklist"]
  asymptote_fam[i, "delta_asym"] <- ((asymptote_fam[i, "n_species_checklist"] - asymptote_fam[i, "n_asymtote"])*100)/asymptote_fam[i, "n_species_checklist"]
  
}

asymptote_fam[asymptote_fam$family=="Myctophidae", 7:9] <- 0


# linear regression and plot for taxa
lm_taxa <- lm(n_species_checklist~n_species, data=asymptote_fam)
summary(lm_taxa)

plot_taxa <- ggplot(asymptote_fam, aes(n_species, n_species_checklist))+
  geom_point(size=2)+
  geom_text_repel(aes(label=fam), size=3.5, point.padding = 0.05)+
  geom_abline(slope = 1, intercept = 0, color="red", size=0.8)+
  geom_abline(slope = 3.8, intercept = -1.7, size=0.8)+
  xlim(0,270)+
  ylim(0,270)+
  theme_bw()+
  labs(x="Number of taxa",
       y="Number of species in the checklist")

grob <- grobTree(textGrob("y = 3.85x-1.7\nR² = 0.74\np < 0.001", x=0.8, y=0.2))
plot_taxa2 <- plot_taxa+annotation_custom(grob)
plot_taxa2


# linear regression and plot for motus
lm_motu <- lm(n_species_checklist~n_motus, data=asymptote_fam)
summary(lm_motu)

plot_motu <- ggplot(asymptote_fam, aes(n_motus, n_species_checklist))+
  geom_point(size=2)+
  geom_text_repel(aes(label=fam), size=3.5, point.padding = 0.05)+
  geom_abline(slope = 1, intercept = 0, color="red", size=0.8)+
  geom_abline(slope = 2.365, intercept = -36.3, size=0.8)+
  xlim(0,270)+
  ylim(0,270)+
  theme_bw()+
  labs(x="Number of MOTUs",
       y="Number of species in the checklist")
grob <- grobTree(textGrob("y = 2.36x-36.3\nR² = 0.65\np < 0.001", x=0.8, y=0.2))
plot_motu2 <- plot_motu+annotation_custom(grob)
plot_motu2


# linear regression and plot for asymptotes
lm_asym <- lm(n_species_checklist~n_asymtote, data=asymptote_fam)
summary(lm_asym)

plot_asym <- ggplot(asymptote_fam, aes(n_asymtote, n_species_checklist))+
  geom_point(size=2)+
  geom_text_repel(aes(label=fam), size=3.5, point.padding = 0.05)+
  geom_abline(slope = 1, intercept = 0, color="red", size=0.8)+
  geom_abline(slope = 1.54, intercept = -10.15, size=0.8)+
  xlim(0,270)+
  ylim(0,270)+
  theme_bw()+
  labs(x="MOTUs asymptote",
       y="Number of species in the checklist")
grob <- grobTree(textGrob("y = 1.54x-10.1\nR² = 0.48\np < 0.005", x=0.8, y=0.2))
plot_asym2 <- plot_asym+annotation_custom(grob)
plot_asym2


# delta checklist - taxa
lm_taxa <- lm(delta_taxa~n_species_checklist, data=asymptote_fam)
summary(lm_taxa)

delta_taxa <- ggplot(asymptote_fam, aes(n_species_checklist, delta_taxa))+
  geom_point(size=2)+
  geom_text_repel(aes(label=fam), size=3.5, point.padding = 0.05)+
  ylim(0,100)+
  theme_bw()+
  labs(x="Number of species in the checklist",
       y="Δ checklist - taxa (%)")
grob <- grobTree(textGrob("R² = 0.12\np = 0.012", x=0.8, y=0.1))
delta_taxa2 <- delta_taxa+annotation_custom(grob)
delta_taxa2

# delta checklist - motu
lm_motu <- lm(delta_motu~n_species_checklist, data=asymptote_fam)
summary(lm_motu)

delta_motu <- ggplot(asymptote_fam, aes(n_species_checklist, delta_motu))+
  geom_point(size=2)+
  geom_text_repel(aes(label=fam), size=3.5, point.padding = 0.05)+
  ylim(-70,100)+
  theme_bw()+
  labs(x="Number of species in the checklist",
       y="Δ checklist - MOTUs (%)")
grob <- grobTree(textGrob("R² = 0.29\np = 0.028", x=0.8, y=0.1))
delta_motu2 <- delta_motu+annotation_custom(grob)
delta_motu2

# delta checklist - asym
lm_asym <- lm(delta_asym~n_species_checklist, data=asymptote_fam)
summary(lm_asym)

delta_asym <- ggplot(asymptote_fam, aes(n_species_checklist, delta_asym))+
  geom_point(size=2)+
  geom_text_repel(aes(label=fam), size=3.5, point.padding = 0.05)+
  ylim(-135,100)+
  theme_bw()+
  labs(x="Number of species in the checklist",
       y="Δ checklist - MOTUs asymp (%)")
grob <- grobTree(textGrob("R² = 0.14\np = 0.10", x=0.8, y=0.1))
delta_asym2 <- delta_asym+annotation_custom(grob)
delta_asym2


# plot all together

plot_all <- ggarrange(plot_taxa2, plot_motu2, plot_asym2, delta_taxa2, delta_motu2, delta_asym2, nrow = 2, ncol = 3, labels = c("A", "B","C", "D", "E", "F"), label.x = 0.17, label.y = 1)
plot_all
