# Accumulation curves on genus, family, order on eDNA and RLS data

# Lib 
library(tidyverse)
library(ggplot2)
library(ggpubr)

# data
load("Rdata/02-clean-data.Rdata")

# 
'%ni%' <- Negate("%in%")

# Functions
source('scripts/03_Analysis/00_functions.R')

# On the df as well
df_all_filters <- df_all_filters %>%
  filter(province %in% c("Western_Indian_Ocean", "Southeast_Polynesia", "Tropical_Northwestern_Atlantic", "Western_Coral_Triangle", "Tropical_Southwestern_Pacific"))%>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3", "glorieuse_distance_300m")) %>%
  filter(sample_method !="niskin" & comment %ni% c("Distance decay 600m", "Distance decay 300m"))%>%
  filter(project != "Curacao") %>%
  filter(habitat=="marine")%>%
  filter(habitat_type %ni% c("BAIE"))


# ------------------------------------------------------------------------------- # 
#### On eDNA MOTUs ----
# ------------------------------------------------------------------------------- # 

# rank_specify
rank_choice = 'sequence'

## select cryptic species
cryptic_family <- c("Tripterygiidae", "Grammatidae", "Aploactinidae", "Creediidae", "Gobiidae", "Chaenopsidae", "Gobiesocidae", "Labrisomidae", "Pseudochromidae", "Bythitidae", "Plesiopidae", "Dactyloscopidae", "Blenniidae", "Apogonidae", "Callionymidae", "Opistognathidae", "Syngnathidae", "Kurtidae")
cryptic_order <- c("Kurtiformes", "Gobiiformes", "Blenniiformes", "Syngnathiformes")
df_all_filters <- filter(df_all_filters, order_name %in% cryptic_order | family_name_corrected %in% cryptic_family)

# Add a global saturation curve
all_accumulation_cryp <- accumulation_curve_df(df_all_filters, species_unit = rank_choice) %>%
  dplyr::mutate(project_name = "All") %>%
  select(project_name, richness, sd, sites)


# Asymptote of all plots 
all_asymptote_cryp <- asymptote_mm(df_all_filters, species_unit = rank_choice) %>%
  dplyr::mutate(project_name = "All") %>%
  select(project_name, asymptote, slope)


# 
df_join_all_cryp <- all_accumulation_cryp %>%
  left_join(., all_asymptote_cryp, by = "project_name") %>%
  dplyr::mutate(position_asymptote_y = 1.05*asymptote, 
         position_asymptote_x = max(sites),
         position_slope_y = 0.20 * max(asymptote))
df_join_all_cryp$family <- "cryptobenthic"

## select pelagic species
load("Rdata/02-clean-data.Rdata")
df_all_filters <- df_all_filters %>%
  filter(province %in% c("Western_Indian_Ocean", "Southeast_Polynesia", "Tropical_Northwestern_Atlantic", "Western_Coral_Triangle", "Tropical_Southwestern_Pacific"))%>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3", "glorieuse_distance_300m")) %>%
  filter(sample_method !="niskin" & comment %ni% c("Distance decay 600m", "Distance decay 300m"))%>%
  filter(project != "Curacao") %>%
  filter(habitat=="marine")%>%
  filter(habitat_type %ni% c("BAIE"))

load("Rdata/pelagic_family.Rdata")
df_all_filters <- subset(df_all_filters, family_name_corrected %in% pelagic_family$family_name)

# Add a global saturation curve
all_accumulation_pel <- accumulation_curve_df(df_all_filters, species_unit = rank_choice) %>%
  dplyr::mutate(project_name = "All") %>%
  select(project_name, richness, sd, sites)


# Asymptote of all plots 
all_asymptote_pel <- asymptote_mm(df_all_filters, species_unit = rank_choice) %>%
  dplyr::mutate(project_name = "All") %>%
  select(project_name, asymptote, slope)

# 
df_join_all_pel <- all_accumulation_pel %>%
  left_join(., all_asymptote_pel, by = "project_name") %>%
  dplyr::mutate(position_asymptote_y = 1.05*asymptote, 
         position_asymptote_x = max(sites),
         position_slope_y = 0.20 * max(asymptote))
df_join_all_pel$family <- "pelagic"


## select demersal species
load("Rdata/02-clean-data.Rdata")
df_all_filters <- df_all_filters %>%
  filter(province %in% c("Western_Indian_Ocean", "Southeast_Polynesia", "Tropical_Northwestern_Atlantic", "Western_Coral_Triangle", "Tropical_Southwestern_Pacific"))%>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3", "glorieuse_distance_300m")) %>%
  filter(sample_method !="niskin" & comment %ni% c("Distance decay 600m", "Distance decay 300m"))%>%
  filter(project != "Curacao") %>%
  filter(habitat=="marine")%>%
  filter(habitat_type %ni% c("BAIE"))


df_all_filters <- df_all_filters%>%
  subset(family_name_corrected %ni% pelagic_family$family_name & family_name_corrected %ni% cryptic_family & order_name %ni% cryptic_order)

# Add a global saturation curve
all_accumulation_dem <- accumulation_curve_df(df_all_filters, species_unit = rank_choice) %>%
  dplyr::mutate(project_name = "All") %>%
  select(project_name, richness, sd, sites)


# Asymptote of all plots 
all_asymptote_dem <- asymptote_mm(df_all_filters, species_unit = rank_choice) %>%
  dplyr::mutate(project_name = "All") %>%
  select(project_name, asymptote, slope)

# 
df_join_all_dem <- all_accumulation_dem %>%
  left_join(., all_asymptote_dem, by = "project_name") %>%
  dplyr::mutate(position_asymptote_y = 1.05*asymptote, 
         position_asymptote_x = max(sites),
         position_slope_y = 0.20 * max(asymptote))
df_join_all_dem$family <- "demersal"


#df_join_all <- rbind(df_join_all_cryp, df_join_all_pel, df_join_all_dem)

# Plots
edna_crypto <- ggplot(df_join_all_cryp, aes(xmax=250)) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd), alpha = 0.5, fill="#d2981a") +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 0.5, col="#d2981a") +
  ylab("Number of MOTUs") +
  xlab("") +
  ggtitle("eDNA")+
  theme(axis.text.x = element_blank())+
  facet_wrap(~family, scales = "free") +
  geom_text(aes(x = position_asymptote_x, y =position_asymptote_y, hjust = 1, label = paste(round(asymptote, 0), "MOTUs")), col = "black", size=3.5)+
  geom_text(aes(x = position_asymptote_x, y =position_slope_y, hjust = 1, label = paste("slope=",round(slope, 1))), col = "black", size=3.5)+
  theme_bw()

edna_pelagic <- ggplot(df_join_all_pel) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd), alpha = 0.5, fill="#d2981a") +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 0.5, col="#d2981a") +
  ylab("Number of MOTUs") +
  xlab("") +
  theme(axis.text.x = element_blank())+
  facet_wrap(~family, scales = "free") +
  geom_text(aes(x = position_asymptote_x, y =position_asymptote_y, hjust = 1, label = paste(round(asymptote, 0), "MOTUs")), col = "black", size=3.5)+
  geom_text(aes(x = position_asymptote_x, y =position_slope_y, hjust = 1, label = paste("slope=",round(slope, 1))), col = "black", size=3.5)+
  theme_bw()

edna_dem <- ggplot(df_join_all_dem) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd), alpha = 0.5, fill="#d2981a") +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 0.5, col="#d2981a") +
  ylab("Number of MOTUs") +
  xlab("Samples (filter)") +
  facet_wrap(~family, scales = "free") +
  geom_text(aes(x = position_asymptote_x, y =position_asymptote_y, hjust = 1, label = paste(round(asymptote, 0), "MOTUs")), col = "black", size=3.5)+
  geom_text(aes(x = position_asymptote_x, y =position_slope_y, hjust = 1, label = paste("slope=",round(slope, 1))), col = "black", size=3.5)+
  theme_bw()







# ------------------------------------------------------------------------------- # 
#### On RLS species ----
# ------------------------------------------------------------------------------- # 
# Functions 
akaike.weights <- function(x)
{
  x <- x[!is.na(x)]
  delta.aic <- x - min(x, na.rm = TRUE)
  rel.LL <- exp(-0.5 * delta.aic)
  sum.LL <- sum(rel.LL, na.rm = TRUE)
  weights.aic <- rel.LL/sum.LL
  return(list(deltaAIC = delta.aic, rel.LL = rel.LL, weights = weights.aic))
}


RLS_species <- read.csv("data/RLS/RLS_species_NEW.csv", sep = ";", stringsAsFactors = FALSE, check.names = FALSE)
RLS_sp <- RLS_species[,c(18:2173)]
RLS_sp <- RLS_sp[,colSums(RLS_sp)>0]
RLS_species <- cbind(RLS_species$SurveyID, RLS_sp)
colnames(RLS_species)[1] <- "SurveyID"

library(rfishbase)
all_fishbase <- load_taxa()


# select cryptobenthic species

fishbase_crypto <- all_fishbase %>%
  filter(Family %in% cryptic_family)
RLS_crypto <- RLS_species %>%
  select(one_of(fishbase_crypto$Species))
RLS_crypto <- cbind(RLS_species$SurveyID, RLS_crypto)
colnames(RLS_crypto)[1] <- "SurveyID"
# method
method_accumulation = "exact"

# Add a global saturation curve for RLS
all_accumulation_crypto_RLS <- specaccum(RLS_crypto[,2:ncol(RLS_crypto)], 
                                          method = method_accumulation, permutations = 100,
                                          conditioned =TRUE, gamma = "jack1")

# Accumulation for plot 
all_accumulation_crypto_RLS_df <- data.frame(richness = all_accumulation_crypto_RLS$richness, sd = all_accumulation_crypto_RLS$sd, sites = all_accumulation_crypto_RLS$sites)


# Asymptote of RLS  
# models asymptotes
lomo_edna <- fitspecaccum(all_accumulation_crypto_RLS, "lomolino")
aic_lomo_edna <-AIC(lomo_edna)
mm_edna <- fitspecaccum(all_accumulation_crypto_RLS, "michaelis-menten")
aic_mm_edna <- AIC(mm_edna)
gom_edna <- fitspecaccum(all_accumulation_crypto_RLS, "gompertz")
aic_gom_edna <- AIC(gom_edna)
asy_edna <- fitspecaccum(all_accumulation_crypto_RLS, "asymp")
aic_asy_edna <- AIC(asy_edna)
gis_edna <- fitspecaccum(all_accumulation_crypto_RLS, "logis")
aic_gis_edna <- AIC(gis_edna)

# compute results
res_edna <- matrix(NA,nrow = 4, ncol = 3)
rownames(res_edna) <- c("lomolino", "michaelis-menten", "asymp", "logis")
colnames(res_edna) <- c("AIC", "Asymptote", "Weigth")
res_edna[,"AIC"] <- c(aic_lomo_edna, aic_mm_edna, aic_asy_edna, aic_gis_edna)
res_edna[,"Weigth"] <- akaike.weights(c(aic_lomo_edna, aic_mm_edna, aic_asy_edna, aic_gis_edna))$weights
res_edna[,"Asymptote"] <- c(coef(lomo_edna)[[1]], coef(mm_edna)[[1]], coef(asy_edna)[[1]], coef(gis_edna)[[1]])

# Calculation asymptote value
asymp_crypto_RLS <- weighted.mean(res_edna[,"Asymptote"], res_edna[,"Weigth"])

# Add to DF
all_accumulation_crypto_RLS_df$asymptote <- asymp_crypto_RLS
all_accumulation_crypto_RLS_df$slope <- coef(lomo_edna)[[3]]
all_accumulation_crypto_RLS_df$family <- "cryptobenthic"
all_accumulation_crypto_RLS_df$position_asymptote_y = 1.05*all_accumulation_crypto_RLS_df$asymptote
all_accumulation_crypto_RLS_df$position_asymptote_x = max(all_accumulation_crypto_RLS_df$sites)
all_accumulation_crypto_RLS_df$position_slope_y = 0.20 * max(all_accumulation_crypto_RLS_df$asymptote)

# select pelagic species

fishbase_pelagic <- all_fishbase %>%
  filter(Family %in% pelagic_family$family_name)
RLS_pelagic <- RLS_species %>%
  select(one_of(fishbase_pelagic$Species))
RLS_pelagic <- cbind(RLS_species$SurveyID, RLS_pelagic)
colnames(RLS_pelagic)[1] <- "SurveyID"
# method
method_accumulation = "exact"

# Add a global saturation curve for RLS
all_accumulation_pelagic_RLS <- specaccum(RLS_pelagic[,2:ncol(RLS_pelagic)], 
                                         method = method_accumulation, permutations = 100,
                                         conditioned =TRUE, gamma = "jack1")

# Accumulation for plot 
all_accumulation_pelagic_RLS_df <- data.frame(richness = all_accumulation_pelagic_RLS$richness, sd = all_accumulation_pelagic_RLS$sd, sites = all_accumulation_pelagic_RLS$sites)


# Asymptote of RLS  
# models asymptotes
lomo_edna <- fitspecaccum(all_accumulation_pelagic_RLS, "lomolino")
aic_lomo_edna <-AIC(lomo_edna)
mm_edna <- fitspecaccum(all_accumulation_pelagic_RLS, "michaelis-menten")
aic_mm_edna <- AIC(mm_edna)
gom_edna <- fitspecaccum(all_accumulation_pelagic_RLS, "gompertz")
aic_gom_edna <- AIC(gom_edna)
asy_edna <- fitspecaccum(all_accumulation_pelagic_RLS, "asymp")
aic_asy_edna <- AIC(asy_edna)
gis_edna <- fitspecaccum(all_accumulation_pelagic_RLS, "logis")
aic_gis_edna <- AIC(gis_edna)

# compute results
res_edna <- matrix(NA,nrow = 4, ncol = 3)
rownames(res_edna) <- c("lomolino", "michaelis-menten", "asymp", "logis")
colnames(res_edna) <- c("AIC", "Asymptote", "Weigth")
res_edna[,"AIC"] <- c(aic_lomo_edna, aic_mm_edna, aic_asy_edna, aic_gis_edna)
res_edna[,"Weigth"] <- akaike.weights(c(aic_lomo_edna, aic_mm_edna, aic_asy_edna, aic_gis_edna))$weights
res_edna[,"Asymptote"] <- c(coef(lomo_edna)[[1]], coef(mm_edna)[[1]], coef(asy_edna)[[1]], coef(gis_edna)[[1]])

# Calculation asymptote value
asymp_pelagic_RLS <- weighted.mean(res_edna[,"Asymptote"], res_edna[,"Weigth"])

# Add to DF
all_accumulation_pelagic_RLS_df$asymptote <- asymp_pelagic_RLS
all_accumulation_pelagic_RLS_df$slope <- coef(lomo_edna)[[3]]
all_accumulation_pelagic_RLS_df$family <- "pelagic"
all_accumulation_pelagic_RLS_df$position_asymptote_y = 1.05*all_accumulation_pelagic_RLS_df$asymptote
all_accumulation_pelagic_RLS_df$position_asymptote_x = max(all_accumulation_pelagic_RLS_df$sites)
all_accumulation_pelagic_RLS_df$position_slope_y = 0.20 * max(all_accumulation_pelagic_RLS_df$asymptote)

# select demersal species

fishbase_demersal <- all_fishbase %>%
  filter(Family %ni% cryptic_family & Family %ni% pelagic_family$family_name)
RLS_demersal <- RLS_species %>%
  select(one_of(fishbase_demersal$Species))
RLS_demersal <- cbind(RLS_species$SurveyID, RLS_demersal)
colnames(RLS_demersal)[1] <- "SurveyID"
# method
method_accumulation = "exact"

# Add a global saturation curve for RLS
all_accumulation_demersal_RLS <- specaccum(RLS_demersal[,2:ncol(RLS_demersal)], 
                                         method = method_accumulation, permutations = 100,
                                         conditioned =TRUE, gamma = "jack1")

# Accumulation for plot 
all_accumulation_demersal_RLS_df <- data.frame(richness = all_accumulation_demersal_RLS$richness, sd = all_accumulation_demersal_RLS$sd, sites = all_accumulation_demersal_RLS$sites)


# Asymptote of RLS  
# models asymptotes
lomo_edna <- fitspecaccum(all_accumulation_demersal_RLS, "lomolino")
aic_lomo_edna <-AIC(lomo_edna)
mm_edna <- fitspecaccum(all_accumulation_demersal_RLS, "michaelis-menten")
aic_mm_edna <- AIC(mm_edna)
gom_edna <- fitspecaccum(all_accumulation_demersal_RLS, "gompertz")
aic_gom_edna <- AIC(gom_edna)
asy_edna <- fitspecaccum(all_accumulation_demersal_RLS, "asymp")
aic_asy_edna <- AIC(asy_edna)
gis_edna <- fitspecaccum(all_accumulation_demersal_RLS, "logis")
aic_gis_edna <- AIC(gis_edna)

# compute results
res_edna <- matrix(NA,nrow = 5, ncol = 3)
rownames(res_edna) <- c("lomolino", "michaelis-menten", "gompertz", "asymp", "logis")
colnames(res_edna) <- c("AIC", "Asymptote", "Weigth")
res_edna[,"AIC"] <- c(aic_lomo_edna, aic_mm_edna, aic_gom_edna, aic_asy_edna, aic_gis_edna)
res_edna[,"Weigth"] <- akaike.weights(c(aic_lomo_edna, aic_mm_edna, aic_gom_edna, aic_asy_edna, aic_gis_edna))$weights
res_edna[,"Asymptote"] <- c(coef(lomo_edna)[[1]], coef(mm_edna)[[1]], coef(gom_edna)[[1]], coef(asy_edna)[[1]], coef(gis_edna)[[1]])

# Calculation asymptote value
asymp_demersal_RLS <- weighted.mean(res_edna[,"Asymptote"], res_edna[,"Weigth"])

# Add to DF
all_accumulation_demersal_RLS_df$asymptote <- asymp_demersal_RLS
all_accumulation_demersal_RLS_df$slope <- coef(lomo_edna)[[3]]
all_accumulation_demersal_RLS_df$family <- "demersal"
all_accumulation_demersal_RLS_df$position_asymptote_y = 1.05*all_accumulation_demersal_RLS_df$asymptote
all_accumulation_demersal_RLS_df$position_asymptote_x = max(all_accumulation_demersal_RLS_df$sites)
all_accumulation_demersal_RLS_df$position_slope_y = 0.20 * max(all_accumulation_demersal_RLS_df$asymptote)


##join all

#df_join_all_RLS <- rbind(all_accumulation_demersal_RLS_df, all_accumulation_crypto_RLS_df, all_accumulation_pelagic_RLS_df)



RLS_crypto <- ggplot(all_accumulation_crypto_RLS_df, aes(ymax=600)) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd), alpha = 0.5) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 0.5) +
  ylab("Number of species") +
  xlab("") +
  ggtitle("Visual Census")+
  theme(axis.text.x = element_blank())+
  facet_wrap(~family, scales = "free") +
  geom_text(aes(x = position_asymptote_x, y =position_asymptote_y, hjust = 1, label = paste(round(asymptote, 0), "species")), col = "black", size=3.5)+
  geom_text(aes(x = position_asymptote_x, y =position_slope_y, hjust = 1, label = paste("slope=",round(slope, 1))), col = "black", size=3.5)+
  theme_bw()

RLS_pel <- ggplot(all_accumulation_pelagic_RLS_df, aes(ymax=40)) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd), alpha = 0.5) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 0.5) +
  ylab("Number of species") +
  xlab("") +
  theme(axis.text.x = element_blank())+
  facet_wrap(~family, scales = "free") +
  geom_text(aes(x = position_asymptote_x, y =position_asymptote_y+1, hjust = 1, label = paste(round(asymptote, 0), "species")), col = "black", size=3.5)+
  geom_text(aes(x = position_asymptote_x, y =position_slope_y, hjust = 1, label = paste("slope=",round(slope, 1))), col = "black", size=3.5)+
  theme_bw()

RLS_dem <- ggplot(all_accumulation_demersal_RLS_df, aes(ymax=2000)) + 
  geom_ribbon(aes(x = sites, ymin = richness-sd, ymax = richness+sd), alpha = 0.5) +
  geom_line(aes(x = sites, y = richness)) +
  geom_hline(aes(yintercept = asymptote), linetype = "dashed", size = 0.5) +
  ylab("Number of species") +
  xlab("Transects") +
  facet_wrap(~family, scales = "free") +
  geom_text(aes(x = position_asymptote_x, y =position_asymptote_y+10, hjust = 1, label = paste(round(asymptote, 0), "species")), col = "black", size=3.5)+
  geom_text(aes(x = position_asymptote_x, y =position_slope_y, hjust = 1, label = paste("slope=",round(slope, 1))), col = "black", size=3.5)+
  theme_bw()


# --------------------------------------------------------------------- # 
#### Final figure - combine all levels  ----
# --------------------------------------------------------------------- # 


ggarrange(edna_crypto, RLS_crypto, edna_pelagic, RLS_pel, edna_dem, RLS_dem, nrow=3, ncol=2, labels=c("a", "", "b", "", "c", ""), heights = c(1.1, 1,1))
ggsave("outputs/00_Figures_for_paper/Figure3.png", width=6, height = 8)





