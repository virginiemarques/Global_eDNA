library(tidyverse)
library(reshape2)
library(vegan)
library(betapart)
library(ggpubr)

# Data
load("Rdata/02-clean-data.Rdata")
family_reference_coef <- read.csv("outputs/01_read_data_stats/family_resolution_coefs.csv")
'%ni%' <- Negate("%in%")
#Remove estuary stations and deep niskin station
df_all_filters <- df_all_filters %>%
  filter(province %in% c("Western_Indian_Ocean", "Southeast_Polynesia", "Tropical_Northwestern_Atlantic", "Western_Coral_Triangle", "Tropical_Southwestern_Pacific"))%>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3", "glorieuse_distance_300m")) %>%
  filter(sample_method !="niskin" & comment %ni% c("Distance decay 600m", "Distance decay 300m"))%>%
  filter(project.y != "SEAMOUNTS") %>%
  filter(habitat=="marine")%>%
  filter(habitat_type %ni% c("BAIE"))

df_all_filters <- df_all_filters %>%
  filter(!is.na(family_name_corrected))

# Count the abondance of MOTUs (global)
family_motus <- df_all_filters %>%
  group_by(family_name_corrected) %>%
  summarise(n = n_distinct(sequence))

# See the coefs in the eDNA dataset global
family_global_coef <- family_reference_coef %>%
  inner_join(., family_motus, by = c("Family" = "family_name_corrected")) %>%
  mutate(plot_low_coef_seq = ifelse(coef_sequencing <= 0.1, "low", "fine"), 
         plot_low_coef_reso = ifelse(coef_resolution <= 0.7, "low", "fine"), 
         plot_low_n_motus = ifelse(n >5, ">5 MOTUs", "<= 5 MOTUs")) 

# --------------------------------------------------------------------------- # 
# Plot for each coef - > 5 MOTUs
# Colors
cols <- c("low" = "red", "fine" = "blue")

# Sequencing biais
plot_sequencing_bias <- ggplot(family_global_coef %>% filter(n >5), aes(x=reorder(Family, n), y=coef_sequencing)) +
  geom_segment( aes(x=reorder(Family, n), xend=reorder(Family, n), y=0, yend=coef_sequencing), color="skyblue") +
  geom_point( aes(color=plot_low_coef_seq), size=6, alpha=0.6, shape = 21, show.legend = FALSE) +
  geom_text(aes(label = n)) +
  theme_bw() + coord_flip() +
  scale_y_continuous(limits= c(0,1), breaks = seq(0,1,0.1)) + 
  scale_color_manual(values = cols) + 
  ggtitle("Sequencing biais (families > 5 MOTUs)", 
          subtitle = "Number in the bubble represents the number of MOTUs")

plot_sequencing_bias

# Resolution biais
plot_resolution_bias <- ggplot(family_global_coef %>% filter(n >5), aes(x=reorder(Family, n), y=coef_resolution)) +
  geom_segment( aes(x=reorder(Family, n), xend=reorder(Family, n), y=0, yend=coef_resolution), color="skyblue") +
  geom_point( aes(color=plot_low_coef_reso), size=6, alpha=0.6, shape = 21, show.legend = FALSE) +
  geom_text(aes(label = n)) +
  theme_bw() + coord_flip() +
  scale_y_continuous(limits= c(0,1), breaks = seq(0,1,0.1)) + 
  scale_color_manual(values = cols) + 
  ggtitle("Resolution biais (families > 5 MOTUs)", 
          subtitle = "Number in the bubble represents the number of MOTUs")

plot_resolution_bias

# Arrange both plots 
ggarrange(plot_sequencing_bias, plot_resolution_bias + rremove("ylab"), ncol = 2)

# Save
ggsave("outputs/08_bias_reference_sequencing/biais_reference_families_5.png", width = 10, height=10)

# --------------------------------------------------------------------------- # 
# Plot for each coef - <6 MOTUs

# Sequencing biais_2
plot_sequencing_bias_low <- ggplot(family_global_coef %>% filter(n < 6), aes(x=reorder(Family, n), y=coef_sequencing)) +
  geom_segment( aes(x=reorder(Family, n), xend=reorder(Family, n), y=0, yend=coef_sequencing), color="skyblue") +
  geom_point( aes(color=plot_low_coef_seq), size=6, alpha=0.6, shape = 21, show.legend = FALSE) +
  geom_text(aes(label = n)) +
  theme_bw() + coord_flip() +
  scale_y_continuous(limits= c(0,1), breaks = seq(0,1,0.1)) + 
  scale_color_manual(values = cols) + 
  ggtitle("Sequencing biais (families < 6 MOTUs)", 
          subtitle = "Number in the bubble represents the number of MOTUs")

plot_sequencing_bias_low

# Resolution biais
plot_resolution_bias_low <- ggplot(family_global_coef %>% filter(n < 6), aes(x=reorder(Family, n), y=coef_resolution)) +
  geom_segment( aes(x=reorder(Family, n), xend=reorder(Family, n), y=0, yend=coef_resolution), color="skyblue") +
  geom_point( aes(color=plot_low_coef_reso), size=6, alpha=0.6, shape = 21, show.legend = FALSE) +
  geom_text(aes(label = n)) +
  theme_bw() + coord_flip() +
  scale_y_continuous(limits= c(0,1), breaks = seq(0,1,0.1)) + 
  scale_color_manual(values = cols) + 
  ggtitle("Resolution biais (families <6 MOTUs)", 
          subtitle = "Number in the bubble represents the number of MOTUs")

plot_resolution_bias_low

# Arrange both plots 
ggarrange(plot_sequencing_bias_low, plot_resolution_bias_low + rremove("ylab"), ncol = 2)

# Save
ggsave("outputs/08_bias_reference_sequencing/biais_reference_families_less_5.png", width = 10, height=12)

# Plot 2 - combine all information
library(ggrepel)

ggplot(family_global_coef, aes(y = coef_sequencing, x = coef_resolution)) + 
  #geom_text_repel(aes(label = Family)) + 
  geom_text(aes(label = Family), check_overlap = T) + 
  geom_point(alpha=0.3) + 
  theme_bw() + 
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) + 
  geom_hline(yintercept=0.5, col = "red") + 
  geom_vline(xintercept=0.5, col = "red") + 
  facet_wrap(~plot_low_n_motus)

# Save
ggsave("outputs/08_bias_reference_sequencing/biais_reference_resolution_sequencing.png", width = 15, height=12)






