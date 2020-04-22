library(tidyverse)
library(rfishbase)
# Load all data
all_fishbase <- load_taxa()
# Complete infos - depth and pelagic
all_fishbase2 <- all_fishbase %>%
  left_join(., rfishbase::species(.$Species,  fields=c("Species", "DemersPelag", "Length")))
# Compute for family level
pelagic_family <- all_fishbase2 %>%
  # Pelagic status
  mutate(pelagic_status = case_when(
    DemersPelag %in% c("pelagic-neritic", "pelagic-oceanic", "bathypelagic", "pelagic") ~ "pelagic",
    TRUE ~ "not_pelagic"
  )) %>%
  group_by(Family) %>%
  # Deep + pelagic count
  summarise(n_tot = n_distinct(Species),
            n_pelagic = sum(pelagic_status == "pelagic"),
            percent_pelagic = round(n_pelagic/n_tot, 2)
  ) %>%
  mutate(pelagic_categ_family = ifelse(percent_pelagic >= 0.70, "pelagic", "not_pelagic")) %>%
  rename(family_name = Family) %>%
  select(-n_tot, -n_pelagic, -percent_pelagic) %>%
  subset(pelagic_categ_family=="pelagic")

save(pelagic_family, file="Rdata/pelagic_family.Rdata")

