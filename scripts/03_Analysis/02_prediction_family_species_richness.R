library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(rfishbase)
library(grid)
library(png)

## build figure prediction of number of species (Fig 1e)
'%ni%' <- Negate("%in%")

load("Rdata/02_clean_all.Rdata")


# coral fish species data
coral_fishes <- read.csv("data/RLS/Coral_fishes2.csv", sep=";", stringsAsFactors = FALSE)
fam_coral  <- coral_fishes %>%
  group_by(Family) %>%
  summarise(n_species = n_distinct(Species))

fam_coral2 <- unique(fam_coral$Family)

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
  #filter(new_family_name%in%(fam_coral %>%filter(n_species > 3) %>%distinct(Family) %>% pull()))
  filter(new_family_name%in%fam_coral2)
            
families3 <- unique(fam_summary$new_family_name)


# calculate number of species per family on the checklist
for (i in 1:length(families3)) {
  fam_summary[i,"n_species_checklist"] <- coral_fishes %>%
    subset(Family==families3[i]) %>%
    dplyr::summarize(n_distinct(Species))
}

# calculate number of species per family in our data

for (i in 1:length(families3)) {
  fam_summary[i,"n_species"] <- df_all_filters %>%
    subset(new_family_name==families3[i]) %>%
    filter(!is.na(new_species_name)) %>%
    dplyr::summarize(n_distinct(new_species_name))
}

# calculate delta checklist - us (%)
for (i in 1:length(families3)) {
  fam_summary[i, "delta_species"] <- ((fam_summary[i, "n_species_checklist"] - fam_summary[i, "n_species"])*100)/fam_summary[i, "n_species_checklist"]
  fam_summary[i, "delta_motu"] <- ((fam_summary[i, "n_species_checklist"] - fam_summary[i, "n_motus"])*100)/fam_summary[i, "n_species_checklist"]
  
}

# calculate percentage us ~ checklist
for (i in 1:length(families3)) {
  fam_summary[i, "perc_species"] <- (fam_summary[i, "n_species"]*100)/fam_summary[i, "n_species_checklist"]
  fam_summary[i, "perc_motu"] <- (fam_summary[i, "n_motus"]*100)/fam_summary[i, "n_species_checklist"]
  
}


# tranform log10 richness
for (i in 1:length(families3)) {
  fam_summary[i, "log_species"] <- log1p(fam_summary[i, "n_species"])
  fam_summary[i, "log_motu"] <- log1p(fam_summary[i, "n_motus"])
  fam_summary[i, "log_checklist"] <- log1p(fam_summary[i, "n_species_checklist"])
  
}

# add short family name
fam_summary$fam <- gsub("(?<=\\S)[aeiouy]", "", fam_summary$new_family_name, perl = TRUE)
fam_summary$fam <- substring(fam_summary$fam, 1, 4)

save(fam_summary, file="Rdata/all_predictions_motus_species.rdata")

## plot log(number of species in checklist) ~ log(our data)
#-----------------------------------------------------------------------------------------------------------


# linear regression and plot for species
lm_species <- lm(log_checklist~log_species, data=fam_summary)
summary(lm_species)

plot_species <- ggplot(fam_summary, aes(log_species, log_checklist))+
  geom_point(size=2)+
  geom_abline(slope = 1, intercept = 0, color="red", size=0.8)+
  geom_abline(slope = 0.76, intercept = 1.54, size=0.8)+
  xlim(0,6)+
  ylim(0,6)+
  theme_bw()+
  labs(x="log(1+Number of species in eDNA)",
       y="log(1+Number of species in RLS)")+
  annotate(geom="text", x=6, y=1, label="y = 0.72x+1.61\nR² = 0.35\np < 0.001", hjust=1, size=3.5) +
  ggtitle("e")

save(plot_species, file="Rdata/plot_prediction_species.rdata")

## pearson correlation test
library(Hmisc)
mcor <- cor(fam_summary[,2:3], method = c("pearson"))
rcorr(as.matrix(fam_summary[,10:11]), type = c("pearson"))

cor.test(fam_summary$n_motus, fam_summary$n_species_checklist, method = c("pearson"))

library(corrplot)
corrplot(mcor, type="upper", order="hclust", tl.col="black", tl.srt=45)

# linear regression and plot for motus
lm_motu <- lm(log_checklist~log_motu, data=fam_summary)
summary(lm_motu)
gobiidae <- readPNG("data/fish_vignette/gobiidae.png")
muraenidae <- readPNG("data/fish_vignette/muraenidae.png")
pomacentridae <- readPNG("data/fish_vignette/pomacentridae.png")
labridae <- readPNG("data/fish_vignette/labridae.png")
mugilidae <- readPNG("data/fish_vignette/mugilidae.png")
kyphosidae <- readPNG("data/fish_vignette/kyphosidae.png")


plot_motu <- ggplot(fam_summary, aes(log_motu, log_checklist))+
  annotation_custom(rasterGrob(gobiidae), xmin = 5.1, xmax = 6, ymin = 4.2, ymax = 5)+
  annotation_custom(rasterGrob(muraenidae), xmin = 4.3, xmax = 5.2, ymin = 3, ymax = 3.5)+
  annotation_custom(rasterGrob(pomacentridae), xmin = 3.8, xmax = 4.7, ymin = 5.4, ymax = 6)+
  annotation_custom(rasterGrob(labridae), xmin = 4.9, xmax = 5.6, ymin = 5.5, ymax = 6.1)+
  annotation_custom(rasterGrob(mugilidae), xmin = 3.1, xmax = 4.1, ymin = 0.8, ymax = 1.4)+
  annotation_custom(rasterGrob(kyphosidae), xmin = 1.2, xmax = 2, ymin = 3.3, ymax = 4)+
  geom_point(size=2)+
  geom_abline(slope = 1, intercept = 0, color="red", size=0.8)+
  geom_abline(slope = 0.85, intercept = 0.5, size=0.8)+
  annotate(geom="text", x=6, y=1, label="y = 0.83x+0.57\nR² = 0.56\np < 0.001", hjust=1, size=3.5) +
  xlim(0,6)+
  ylim(0,6)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill = NA)) + 
  labs(x="log(1+Number of MOTUs in eDNA)", y="log(1+Number of species in RLS)")+
  ggtitle("f")

lm_motu <- lm(log_checklist~log_motu +offset(log_motu), data=fam_summary)
summary(lm_motu)

save(plot_motu, file="Rdata/plot_prediction_motus.rdata")


