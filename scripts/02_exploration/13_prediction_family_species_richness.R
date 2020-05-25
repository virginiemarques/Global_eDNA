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
fam_coral  <- coral_fishes %>%
  group_by(Family) %>%
  summarise(n_species = n_distinct(Species))

# our data
df_all_filters <- df_all_filters %>%
  #filter(new_rank_ncbi != "higher") %>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")) %>%
  filter(sample_method !="niskin" & region!="East_Pacific" & comment %ni% c("Distance decay 600m", "Distance decay 300m") & station!="glorieuse_distance_300m")

fam_summary <- df_all_filters %>%
  filter(!is.na(new_family_name)) %>%
  group_by(new_family_name) %>%
  summarise(n_motus = n_distinct(sequence))
            
fam_summary <- fam_summary%>%
  filter(new_family_name%in%(fam_coral %>%filter(n_species > 5) %>%distinct(Family) %>% pull()))

families10 <- unique(fam_summary$new_family_name)


# calculate number of species per family on the checklist
for (i in 1:length(families10)) {
  fam_summary[i,"n_species_checklist"] <- coral_fishes %>%
    subset(Family==families10[i]) %>%
    summarize(n_distinct(Species))
}

# calculate number of taxa per family in our data

for (i in 1:length(families10)) {
  fam_summary[i,"n_taxa"] <- df_all_filters %>%
    subset(new_family_name==families10[i]) %>%
    filter(!is.na(new_scientific_name_ncbi)) %>%
    summarize(n_distinct(new_scientific_name_ncbi))
}

# calculate delta checklist - us (%)
for (i in 1:length(families10)) {
  fam_summary[i, "delta_taxa"] <- ((fam_summary[i, "n_species_checklist"] - fam_summary[i, "n_taxa"])*100)/fam_summary[i, "n_species_checklist"]
  fam_summary[i, "delta_motu"] <- ((fam_summary[i, "n_species_checklist"] - fam_summary[i, "n_motus"])*100)/fam_summary[i, "n_species_checklist"]
  #fam_summary[i, "delta_asym"] <- ((fam_summary[i, "n_species_checklist"] - fam_summary[i, "n_asymtote"])*100)/fam_summary[i, "n_species_checklist"]
  
}

# calculate percentage us ~ checklist
for (i in 1:length(families10)) {
  fam_summary[i, "perc_taxa"] <- (fam_summary[i, "n_taxa"]*100)/fam_summary[i, "n_species_checklist"]
  fam_summary[i, "perc_motu"] <- (fam_summary[i, "n_motus"]*100)/fam_summary[i, "n_species_checklist"]
  #fam_summary[i, "perc_asym"] <- (fam_summary[i, "n_asymtote"]*100)/fam_summary[i, "n_species_checklist"]
  
}


# tranform log richness
for (i in 1:length(families10)) {
  fam_summary[i, "log_taxa"] <- log10(fam_summary[i, "n_taxa"])
  fam_summary[i, "log_motu"] <- log10(fam_summary[i, "n_motus"])
  #fam_summary[i, "log_asym"] <- log10(fam_summary[i, "n_asymtote"])
  fam_summary[i, "log_checklist"] <- log10(fam_summary[i, "n_species_checklist"])
  
}

# add short family name

#fam_summary$fam <- c("Acnt", "Apgn", "Blst", "Blnd", "Crng", "Chtd", "Dstd", "Gobd", "Hmld", "Hlcn", "Lbrd", "Lthr", "Ltjn", "Mnct", "Mlld", "Mrnd", "Pmcnth", "Pmcntr", "Scbr", "Srrn", "Sgnd", "Synd", "Ttrdt")



# linear regression and plot for taxa
lm_taxa <- lm(log_checklist~log_taxa, data=fam_summary)
summary(lm_taxa)

plot_taxa <- ggplot(fam_summary, aes(log_taxa, log_checklist))+
  geom_point(size=2)+
  #geom_text(aes(label=fam), size=3, position = position_jitter(width=0.5, height = 0.5))+
  geom_abline(slope = 1, intercept = 0, color="red", size=0.8)+
  geom_abline(slope = 0.58, intercept = 1, size=0.8)+
  xlim(0,2.5)+
  ylim(0,2.5)+
  theme_bw()+
  labs(x="log(Number of taxa)",
       y="log(Number of species in the checklist)")

grob <- grobTree(textGrob("y = 0.58x+1\nR² = 0.48\np < 0.001", x=0.8, y=0.2))
plot_taxa2 <- plot_taxa+annotation_custom(grob)
plot_taxa2


# linear regression and plot for motus
lm_motu <- lm(log_checklist~log_motu, data=fam_summary)
summary(lm_motu)

plot_motu <- ggplot(fam_summary, aes(log_motu, log_checklist))+
  geom_point(size=2)+
  #geom_text(aes(label=fam), size=3, position = position_jitter(width=0.5, height = 0.5))+
  geom_abline(slope = 1, intercept = 0, color="red", size=0.8)+
  geom_abline(slope = 0.55, intercept = 0.86, size=0.8)+
  xlim(0,2.5)+
  ylim(0,2.5)+
  theme_bw()+
  labs(x="log(Number of MOTUs)",
       y="log(Number of species in the checklist)")
grob <- grobTree(textGrob("y = 0.55x+0.86\nR² = 0.54\np < 0.001", x=0.8, y=0.2))
plot_motu2 <- plot_motu+annotation_custom(grob)
plot_motu2


# linear regression and plot for asymptotes
lm_asym <- lm(log_checklist~log_asym, data=fam_summary)
summary(lm_asym)

plot_asym <- ggplot(fam_summary, aes(log_asym, log_checklist))+
  geom_point(size=2)+
  geom_text(aes(label=fam), size=3, position = position_jitter(width=0.5, height = 0.5))+
  geom_abline(slope = 1, intercept = 0, color="red", size=0.8)+
  geom_abline(slope = 0.65, intercept = 1.6, size=0.8)+
  xlim(2,6)+
  ylim(2,6)+
  theme_bw()+
  labs(x="log(MOTUs asymptote)",
       y="log(Number of species in the checklist)")
grob <- grobTree(textGrob("y = 0.65x+1.6\nR² = 0.4\np < 0.001", x=0.8, y=0.2))
plot_asym2 <- plot_asym+annotation_custom(grob)
plot_asym2


# delta checklist - taxa
lm_taxa <- lm(delta_taxa~n_species_checklist, data=fam_summary)
summary(lm_taxa)

delta_taxa <- ggplot(fam_summary, aes(n_species_checklist, delta_taxa))+
  geom_point(size=2)+
  geom_text(aes(label=fam), size=3, position = position_jitter(width=10, height = 10))+
  ylim(0,110)+
  theme_bw()+
  labs(x="Number of species in the checklist",
       y="Δ checklist - taxa (%)")
grob <- grobTree(textGrob("R² = -0.01\np = 0.39", x=0.8, y=0.1))
delta_taxa2 <- delta_taxa+annotation_custom(grob)
delta_taxa2

# delta checklist - motu
lm_motu <- lm(delta_motu~n_species_checklist, data=fam_summary)
summary(lm_motu)

delta_motu <- ggplot(fam_summary, aes(n_species_checklist, delta_motu))+
  geom_point(size=2)+
  geom_text(aes(label=fam), size=3, position = position_jitter(width=10, height = 10))+
  ylim(-100,110)+
  theme_bw()+
  labs(x="Number of species in the checklist",
       y="Δ checklist - MOTUs (%)")
grob <- grobTree(textGrob("R² = 0.04\np = 0.18", x=0.8, y=0.1))
delta_motu2 <- delta_motu+annotation_custom(grob)
delta_motu2


# delta checklist - asym
lm_asym <- lm(delta_asym~n_species_checklist, data=fam_summary)
summary(lm_asym)

delta_asym <- ggplot(fam_summary, aes(n_species_checklist, delta_asym))+
  geom_point(size=2)+
  geom_text(aes(label=fam), size=3, position = position_jitter(width=10, height = 10))+
  ylim(-150,110)+
  theme_bw()+
  labs(x="Number of species in the checklist",
       y="Δ checklist - MOTUs asymp (%)")
grob <- grobTree(textGrob("R² = 0.01\np = 0.27", x=0.8, y=0.1))
delta_asym2 <- delta_asym+annotation_custom(grob)
delta_asym2


# perc checklist - taxa
lm_taxa <- lm(perc_taxa~n_species_checklist, data=fam_summary)
summary(lm_taxa)

perc_taxa <- ggplot(fam_summary, aes(n_species_checklist, perc_taxa))+
  geom_point(size=2)+
  #geom_text(aes(label=fam), size=3, position = position_jitter(width=10, height = 10))+
  geom_hline(yintercept=100, color="red")+
  theme_bw()+
  labs(x="Number of species in the checklist",
       y="%(taxa ~ checklist)")
grob <- grobTree(textGrob("R² = 0.02\np = 0.41", x=0.8, y=0.9))
perc_taxa2 <- perc_taxa+annotation_custom(grob)
perc_taxa2

# perc checklist - motu
lm_motu <- lm(perc_motu~n_species_checklist, data=fam_summary)
summary(lm_motu)

perc_motu <- ggplot(fam_summary, aes(n_species_checklist, perc_motu))+
  geom_point(size=2)+
  #geom_text(aes(label=fam), size=3, position = position_jitter(width=10, height = 10))+
  geom_hline(yintercept=100, color="red")+
  theme_bw()+
  labs(x="Number of species in the checklist",
       y="%(MOTUs ~ checklist)")
grob <- grobTree(textGrob("R² = 0.006\np = 0.62", x=0.8, y=0.9))
perc_motu2 <- perc_motu+annotation_custom(grob)
perc_motu2

# perc checklist - asym
lm_asym <- lm(perc_asym~n_species_checklist, data=fam_summary)
summary(lm_asym)

perc_asym <- ggplot(fam_summary, aes(n_species_checklist, perc_asym))+
  geom_point(size=2)+
  geom_text(aes(label=fam), size=3, position = position_jitter(width=15, height = 15))+
  geom_hline(yintercept=100, color="red")+
  theme_bw()+
  labs(x="Number of species in the checklist",
       y="%(MOTUs asymp ~ checklist)")
grob <- grobTree(textGrob("R² = 0.05\np = 0.14", x=0.8, y=0.9))
perc_asym2 <- perc_asym+annotation_custom(grob)
perc_asym2




# plot all together (les labels c'est degueulasse)
plot_rich_perc <- ggarrange(plot_taxa2, plot_motu2, perc_taxa2, perc_motu2, nrow = 2, ncol = 2, labels = c("A", "B","C", "D"), label.x = -0.007, label.y = 1)
plot_rich_perc



plot_rich_delta <- ggarrange(plot_taxa2, plot_motu2, plot_asym2, delta_taxa2, delta_motu2, delta_asym2, nrow = 2, ncol = 3, labels = c("A", "B","C", "D", "E", "F"), label.x = -0.007, label.y = 1)
plot_rich_delta
