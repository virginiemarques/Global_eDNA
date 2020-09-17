# Accumulation curves on genus, family, order 
# Lib 
library(tidyverse)
library(ggplot2)
library(ggpubr)

# data
load("Rdata/02_clean_all.Rdata")

# 
'%ni%' <- Negate("%in%")

# Functions
source('scripts/02_exploration/00_functions.R')

# Filter
df_all_filters <- df_all_filters %>%
  #filter(new_rank_ncbi != "higher") %>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")) %>%
  filter(sample_method !="niskin" & region!="East_Pacific" & comment %ni% c("Distance decay 600m", "Distance decay 300m") & station!="glorieuse_distance_300m")

# Summary of the families
fam_summary <- df_all_filters %>%
  filter(!is.na(new_family_name)) %>%
  group_by(new_family_name) %>%
  summarise(n_motus = n_distinct(sequence), 
            n_samples = n_distinct(sample_name_all_pcr),
            n_station = n_distinct(station))

# extract abundant families from RLS data
coral_fishes <- read.csv("data/Coral_fishes2.csv", sep=";")
fam_coral  <- coral_fishes %>%
  group_by(Family) %>%
  summarise(n_species = n_distinct(Species))

#select the same families in our data
families10 <- fam_summary%>%
  filter(new_family_name%in%(fam_coral %>%filter(n_species > 10) %>%distinct(Family) %>% pull()))
families10 <- families10 %>%
  filter(n_motus > 5)
families10 <- unique(families10$new_family_name)

# Set df to fill 
asymptote_fam <- data.frame(family=character(),
                          n_motus=character(), 
                          n_asymtote=numeric(), 
                          stringsAsFactors = FALSE)

# Complete samples
complete_samples <- unique(df_all_filters$sample_name_all_pcr)

# ------------------------------------------------------------------------------- # 
####       ASYMPTOTES  ----
# ------------------------------------------------------------------------------- # 

# loop over
for (i in 1:length(families10)){
  
  # Debug
  #i=1
  fam <- families10[[i]]
  
  # n_motus
  df_temp <- df_all_filters %>% 
    filter(new_family_name %in% fam)
  
  n_motus <- length(unique(df_temp$sequence))
  n_motus
  
  # Asymptote
  liste_accumulation <- accumulation_curve_all_samples(df_temp, species_unit = "sequence", 
                                 column_station = "sample_name_all_pcr", 
                                 complete_samples = unique(df_all_filters$sample_name_all_pcr))
  
  n_asymptote <- round(liste_accumulation$asymptote,0)
  
  # fill 
  asymptote_fam[i,"family"] <- fam
  asymptote_fam[i,"n_motus"] <- n_motus
  asymptote_fam[i,"n_asymtote"] <- n_asymptote
  
  #
  print(paste(i ,"-", fam))
}

# ------------------------------------------------------------------------------- # 
####       PLOTS ACCUMULATION  ----
# ------------------------------------------------------------------------------- # 

# Save all the families to make a giant plot 
liste_plots_family <- list()
liste_ggarrange <- list()

# Print the curves for each family 
for (i in 1:length(families20)){
  
  # Debug
  #i=1
  fam <- families20[[i]]
  
  # n_motus
  df_temp <- df_all_filters %>% 
    filter(new_family_name %in% fam) %>%
    select(sequence, sample_name_all_pcr, count_reads, new_family_name)
  
  liste_accumulation <- accumulation_curve_all_samples(df_temp, species_unit = "sequence", 
                                                       column_station = "sample_name_all_pcr", 
                                                       complete_samples = unique(df_all_filters$sample_name_all_pcr))
  
  n_asymptote <- round(liste_accumulation$asymptote,0)
  liste_accumulation$accumulation_plot$asymptote <- n_asymptote
  
  # Plot
  plot_acc_family <- ggplot(liste_accumulation$accumulation_plot) + 
    geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.7) +
    geom_line(aes(x = sites, y = richness)) +
    geom_hline(yintercept = n_asymptote, size = 0.7) +
    ylab("Number of MOTUs") +
    xlab("Samples (filter)") +
    theme_bw() + 
    ggtitle(fam)
  
  plot_acc_family
  # Save
  liste_plots_family[[i]] <- liste_accumulation$accumulation_plot
  liste_ggarrange[[i]] <- plot_acc_family
    
  #
  print(paste(i ,"-", fam, " - ", n_asymptote))
}

# Names
names(liste_plots_family) <- families20

# Plot all 
ggarrange(plotlist=liste_ggarrange)

# Unlist 
df_family_accumulation <- bind_rows(liste_plots_family, .id = "family")

# Plot
ggplot(df_family_accumulation, aes(fill = family)) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.5, fill= "black") +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), size = 0.7) +
  facet_wrap(~family, scales = 'free') + 
  ylab("Number of MOTUs") +
  xlab("Samples (filter)") +
  theme_bw() + 
  theme(legend.position = "none")

ggsave("outputs/10_acculation_curve_family/all_family_20_accumulation_curve_no_scale.png", width=8, height=8)

# Plot 2
ggplot(df_family_accumulation, aes(fill = family)) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd),  alpha = 0.5, fill= "black") +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), size = 0.7) +
  facet_wrap(~family) + 
  ylab("Number of MOTUs") +
  xlab("Samples (filter)") +
  theme_bw() + 
  theme(legend.position = "none")

ggsave("outputs/10_acculation_curve_family/all_family_20_accumulation_curve_scaled.png", width=8, height=8)



# ------------------------------------------------------------------------------- # 
####       EXPORT   ----
# ------------------------------------------------------------------------------- # 

# Export 
write.csv(asymptote_fam, "outputs/10_acculation_curve_family/table_asymptote_family10.csv", row.names = F)



