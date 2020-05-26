library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(rfishbase)
library(grid)

## build figure prediction of number of species (Fig 1C)
'%ni%' <- Negate("%in%")

load("Rdata/02_clean_all.Rdata")


# coral fish species data
coral_fishes <- read.csv("data/Coral_fishes2.csv", sep=";")
fam_coral  <- coral_fishes %>%
  group_by(Family) %>%
  summarise(n_species = n_distinct(Species))

# our data
df_all_filters <- df_all_filters %>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")) %>%
  filter(sample_method !="niskin" & region!="East_Pacific" & comment %ni% c("Distance decay 600m", "Distance decay 300m") & station!="glorieuse_distance_300m")%>%
  filter(project != "SEAMOUNTS") %>% 
  filter(habitat_type %ni% c("BAIE", "Sommet"))

fam_summary <- df_all_filters %>%
  filter(!is.na(new_family_name)) %>%
  group_by(new_family_name) %>%
  summarise(n_motus = n_distinct(sequence))%>%
  filter(new_family_name%in%(fam_coral %>%filter(n_species > 3) %>%distinct(Family) %>% pull()))
            
families3 <- unique(fam_summary$new_family_name)


# calculate number of species per family on the checklist
for (i in 1:length(families3)) {
  fam_summary[i,"n_species_checklist"] <- coral_fishes %>%
    subset(Family==families3[i]) %>%
    summarize(n_distinct(Species))
}

# calculate number of taxa per family in our data

for (i in 1:length(families3)) {
  fam_summary[i,"n_taxa"] <- df_all_filters %>%
    subset(new_family_name==families3[i]) %>%
    filter(!is.na(new_scientific_name_ncbi)) %>%
    summarize(n_distinct(new_scientific_name_ncbi))
}

# calculate delta checklist - us (%)
for (i in 1:length(families3)) {
  fam_summary[i, "delta_taxa"] <- ((fam_summary[i, "n_species_checklist"] - fam_summary[i, "n_taxa"])*100)/fam_summary[i, "n_species_checklist"]
  fam_summary[i, "delta_motu"] <- ((fam_summary[i, "n_species_checklist"] - fam_summary[i, "n_motus"])*100)/fam_summary[i, "n_species_checklist"]
  
}

# calculate percentage us ~ checklist
for (i in 1:length(families3)) {
  fam_summary[i, "perc_taxa"] <- (fam_summary[i, "n_taxa"]*100)/fam_summary[i, "n_species_checklist"]
  fam_summary[i, "perc_motu"] <- (fam_summary[i, "n_motus"]*100)/fam_summary[i, "n_species_checklist"]
  
}


# tranform log10 richness
for (i in 1:length(families3)) {
  fam_summary[i, "log_taxa"] <- log10(fam_summary[i, "n_taxa"])
  fam_summary[i, "log_motu"] <- log10(fam_summary[i, "n_motus"])
  fam_summary[i, "log_checklist"] <- log10(fam_summary[i, "n_species_checklist"])
  
}

# add short family name
#fam_summary$fam <- c("Acnt", "Apgn", "Blst", "Blnd", "Crng", "Chtd", "Dstd", "Gobd", "Hmld", "Hlcn", "Lbrd", "Lthr", "Ltjn", "Mnct", "Mlld", "Mrnd", "Pmcnth", "Pmcntr", "Scbr", "Srrn", "Sgnd", "Synd", "Ttrdt")



## plot log(number of species in checklist) ~ log(our data)
#-----------------------------------------------------------------------------------------------------------


# linear regression and plot for taxa
lm_taxa <- lm(log_checklist~log_taxa, data=fam_summary)
summary(lm_taxa)

plot_taxa <- ggplot(fam_summary, aes(log_taxa, log_checklist))+
  geom_point(size=2)+
  #geom_text(aes(label=fam), size=3, position = position_jitter(width=0.5, height = 0.5))+
  geom_abline(slope = 1, intercept = 0, color="red", size=0.8)+
  geom_abline(slope = 0.6, intercept = 0.9, size=0.8)+
  xlim(0,2.5)+
  ylim(0,2.5)+
  theme_bw()+
  labs(x="log(Number of taxa)",
       y="log(Number of species in the checklist)")

grob <- grobTree(textGrob("y = 0.6x+0.9\nR² = 0.42\np < 0.001", x=0.8, y=0.2))
plot_taxa2 <- plot_taxa+annotation_custom(grob)
plot_taxa2


# linear regression and plot for motus
lm_motu <- lm(log_checklist~log_motu, data=fam_summary)
summary(lm_motu)

plot_motu <- ggplot(fam_summary, aes(log_motu, log_checklist))+
  geom_point(size=2)+
  #geom_text(aes(label=fam), size=3, position = position_jitter(width=0.5, height = 0.5))+
  geom_abline(slope = 1, intercept = 0, color="red", size=0.8)+
  geom_abline(slope = 0.57, intercept = 0.78, size=0.8)+
  xlim(0,2.5)+
  ylim(0,2.5)+
  theme_bw()+
  labs(x="log(Number of MOTUs)",
       y="log(Number of species in the checklist)")
grob <- grobTree(textGrob("y = 0.55x+0.78\nR² = 0.49\np < 0.001", x=0.8, y=0.2))
plot_motu2 <- plot_motu+annotation_custom(grob)
plot_motu2



## plot percentage of our data ~ checklist
#-----------------------------------------------------------------------------------------------------------


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
grob <- grobTree(textGrob("R² = 0.01\np = 0.24", x=0.8, y=0.9))
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
grob <- grobTree(textGrob("R² = 0.01\np = 0.39", x=0.8, y=0.9))
perc_motu2 <- perc_motu+annotation_custom(grob)
perc_motu2




# plot all together (les labels c'est degueulasse)
plot_rich_perc <- ggarrange(plot_taxa2, plot_motu2, perc_taxa2, perc_motu2, nrow = 2, ncol = 2, labels = c("A", "B","C", "D"), label.x = -0.007, label.y = 1)
plot_rich_perc







## plots not used ##
#-----------------------------------------------------------------------------------------------------------


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

#-------------------------------------------------------------------------------------------------