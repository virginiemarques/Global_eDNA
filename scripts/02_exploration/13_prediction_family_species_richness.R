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
load("Rdata/02_clean_all.Rdata")

# family asymptote  + n_motus data
asymptote_fam <- read.csv("outputs/10_acculation_curve_family/table_asymptote_family10.csv")

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

families10 <- fam_summary %>%
  filter(n_motus > 10 & n_station > 10) %>%
  distinct(new_family_name) %>% pull()

# calculate number of species per family on the checklist
for (i in 1:length(families10)) {
  asymptote_fam[i,"n_species_checklist"] <- coral_fishes %>%
    subset(Family==families10[i]) %>%
    summarize(n_distinct(Species))
}

# calculate number of taxa per family in our data

for (i in 1:length(families10)) {
  asymptote_fam[i,"n_species"] <- df_all_filters %>%
    subset(new_family_name==families10[i]) %>%
    filter(!is.na(new_scientific_name_ncbi)) %>%
    summarize(n_distinct(new_scientific_name_ncbi))
}

# calculate delta checklist - us (%)
for (i in 1:length(families10)) {
  asymptote_fam[i, "delta_taxa"] <- ((asymptote_fam[i, "n_species_checklist"] - asymptote_fam[i, "n_species"])*100)/asymptote_fam[i, "n_species_checklist"]
  asymptote_fam[i, "delta_motu"] <- ((asymptote_fam[i, "n_species_checklist"] - asymptote_fam[i, "n_motus"])*100)/asymptote_fam[i, "n_species_checklist"]
  asymptote_fam[i, "delta_asym"] <- ((asymptote_fam[i, "n_species_checklist"] - asymptote_fam[i, "n_asymtote"])*100)/asymptote_fam[i, "n_species_checklist"]
  
}

# calculate percentage us ~ checklist
for (i in 1:length(families10)) {
  asymptote_fam[i, "perc_taxa"] <- (asymptote_fam[i, "n_species"]*100)/asymptote_fam[i, "n_species_checklist"]
  asymptote_fam[i, "perc_motu"] <- (asymptote_fam[i, "n_motus"]*100)/asymptote_fam[i, "n_species_checklist"]
  asymptote_fam[i, "perc_asym"] <- (asymptote_fam[i, "n_asymtote"]*100)/asymptote_fam[i, "n_species_checklist"]
  
}

# remove families not in RLS
asymptote_fam <- asymptote_fam %>%
  filter(n_species_checklist!=0)
asymptote_fam$fam <- c("Acnt", "Apgn", "Blst", "Crng", "Crch", "Chtd", "Clpd","Dstd", "Gobd", "Hlcn", "Lbrd", "Lthr", "Ltjn", "Mgld", "Mlld", "Mrnd", "Pmcnth", "Pmcntr", "Scbrd", "Srrn", "Sndt", "Ttrdt")


log(asymptote_fam$n_species)



# linear regression and plot for taxa
lm_taxa <- lm(log(n_species_checklist)~log(n_species), data=asymptote_fam)
summary(lm_taxa)

plot_taxa <- ggplot(asymptote_fam, aes(log(n_species), log(n_species_checklist)))+
  geom_point(size=2)+
  geom_text(aes(label=fam), size=3, position = position_jitter(width=0.5, height = 0.5))+
  geom_abline(slope = 1, intercept = 0, color="red", size=0.8)+
  geom_abline(slope = 0.9, intercept = 1.16, size=0.8)+
  xlim(0,6)+
  ylim(0,6)+
  theme_bw()+
  labs(x="log(Number of taxa)",
       y="log(Number of species in the checklist)")

grob <- grobTree(textGrob("y = 0.9x+1.16\nR² = 0.32\np < 0.005", x=0.8, y=0.2))
plot_taxa2 <- plot_taxa+annotation_custom(grob)
plot_taxa2


# linear regression and plot for motus
lm_motu <- lm(log(n_species_checklist)~log(n_motus), data=asymptote_fam)
summary(lm_motu)

plot_motu <- ggplot(asymptote_fam, aes(log(n_motus), log(n_species_checklist)))+
  geom_point(size=2)+
  geom_text(aes(label=fam), size=3, position = position_jitter(width=0.5, height = 0.5))+
  geom_abline(slope = 1, intercept = 0, color="red", size=0.8)+
  geom_abline(slope = 1.27, intercept = -0.68, size=0.8)+
  xlim(0,6)+
  ylim(0,6)+
  theme_bw()+
  labs(x="log(Number of MOTUs)",
       y="log(Number of species in the checklist)")
grob <- grobTree(textGrob("y = 1.27x-0.68\nR² = 0.59\np < 0.001", x=0.8, y=0.2))
plot_motu2 <- plot_motu+annotation_custom(grob)
plot_motu2


# linear regression and plot for asymptotes
lm_asym <- lm(log(n_species_checklist)~log(n_asymtote), data=asymptote_fam)
summary(lm_asym)

plot_asym <- ggplot(asymptote_fam, aes(log(n_asymtote), log(n_species_checklist)))+
  geom_point(size=2)+
  geom_text(aes(label=fam), size=3, position = position_jitter(width=0.5, height = 0.5))+
  geom_abline(slope = 1, intercept = 0, color="red", size=0.8)+
  geom_abline(slope = 1.17, intercept = -0.5, size=0.8)+
  xlim(0,6)+
  ylim(0,6)+
  theme_bw()+
  labs(x="log(MOTUs asymptote)",
       y="log(Number of species in the checklist)")
grob <- grobTree(textGrob("y = 1.17x-0.5\nR² = 0.56\np < 0.001", x=0.8, y=0.2))
plot_asym2 <- plot_asym+annotation_custom(grob)
plot_asym2


# delta checklist - taxa
lm_taxa <- lm(delta_taxa~n_species_checklist, data=asymptote_fam)
summary(lm_taxa)

delta_taxa <- ggplot(asymptote_fam, aes(n_species_checklist, delta_taxa))+
  geom_point(size=2)+
  geom_text(aes(label=fam), size=3, position = position_jitter(width=20, height = 30))+
  ylim(-300,110)+
  theme_bw()+
  labs(x="Number of species in the checklist",
       y="Δ checklist - taxa (%)")
grob <- grobTree(textGrob("R² = 0.08\np = 0.11", x=0.8, y=0.1))
delta_taxa2 <- delta_taxa+annotation_custom(grob)
delta_taxa2

# delta checklist - motu
lm_motu <- lm(delta_motu~n_species_checklist, data=asymptote_fam)
summary(lm_motu)

delta_motu <- ggplot(asymptote_fam, aes(n_species_checklist, delta_motu))+
  geom_point(size=2)+
  geom_text(aes(label=fam), size=3, position = position_jitter(width=20, height = 30))+
  ylim(-500,110)+
  theme_bw()+
  labs(x="Number of species in the checklist",
       y="Δ checklist - MOTUs (%)")
grob <- grobTree(textGrob("R² = 0.10\np = 0.078", x=0.8, y=0.1))
delta_motu2 <- delta_motu+annotation_custom(grob)
delta_motu2

# delta checklist - asym
lm_asym <- lm(delta_asym~n_species_checklist, data=asymptote_fam)
summary(lm_asym)

delta_asym <- ggplot(asymptote_fam, aes(n_species_checklist, delta_asym))+
  geom_point(size=2)+
  geom_text(aes(label=fam), size=3, position = position_jitter(width=20, height = 30))+
  ylim(-650,110)+
  theme_bw()+
  labs(x="Number of species in the checklist",
       y="Δ checklist - MOTUs asymp (%)")
grob <- grobTree(textGrob("R² = 0.08\np = 0.10", x=0.8, y=0.1))
delta_asym2 <- delta_asym+annotation_custom(grob)
delta_asym2


# perc checklist - taxa
lm_taxa <- lm(perc_taxa~n_species_checklist, data=asymptote_fam)
summary(lm_taxa)

perc_taxa <- ggplot(asymptote_fam, aes(n_species_checklist, perc_taxa))+
  geom_point(size=2)+
  geom_text(aes(label=fam), size=3, position = position_jitter(width=20, height = 30))+
  geom_hline(yintercept=100, color="red")+
  theme_bw()+
  labs(x="Number of species in the checklist",
       y="%(taxa ~ checklist)")
grob <- grobTree(textGrob("R² = 0.08\np = 0.11", x=0.8, y=0.9))
perc_taxa2 <- perc_taxa+annotation_custom(grob)
perc_taxa2

# perc checklist - motu
lm_motu <- lm(perc_motu~n_species_checklist, data=asymptote_fam)
summary(lm_motu)

perc_motu <- ggplot(asymptote_fam, aes(n_species_checklist, perc_motu))+
  geom_point(size=2)+
  geom_text(aes(label=fam), size=3, position = position_jitter(width=20, height = 30))+
  geom_hline(yintercept=100, color="red")+
  theme_bw()+
  labs(x="Number of species in the checklist",
       y="%(MOTUs ~ checklist)")
grob <- grobTree(textGrob("R² = 0.10\np = 0.078", x=0.8, y=0.9))
perc_motu2 <- perc_motu+annotation_custom(grob)
perc_motu2

# perc checklist - asym
lm_asym <- lm(perc_asym~n_species_checklist, data=asymptote_fam)
summary(lm_asym)

perc_asym <- ggplot(asymptote_fam, aes(n_species_checklist, perc_asym))+
  geom_point(size=2)+
  geom_text(aes(label=fam), size=3, position = position_jitter(width=20, height = 30))+
  geom_hline(yintercept=100, color="red")+
  theme_bw()+
  labs(x="Number of species in the checklist",
       y="%(MOTUs asymp ~ checklist)")
grob <- grobTree(textGrob("R² = 0.08\np = 0.10", x=0.8, y=0.9))
perc_asym2 <- perc_asym+annotation_custom(grob)
perc_asym2




# plot all together (les labels c'est degueulasse)
plot_rich_perc <- ggarrange(plot_taxa2, plot_motu2, plot_asym2, perc_taxa2, perc_motu2, perc_asym2, nrow = 2, ncol = 3, labels = c("A", "B","C", "D", "E", "F"), label.x = -0.007, label.y = 1)
plot_rich_perc



plot_rich_delta <- ggarrange(plot_taxa2, plot_motu2, plot_asym2, delta_taxa2, delta_motu2, delta_asym2, nrow = 2, ncol = 3, labels = c("A", "B","C", "D", "E", "F"), label.x = -0.007, label.y = 1)
plot_rich_delta
