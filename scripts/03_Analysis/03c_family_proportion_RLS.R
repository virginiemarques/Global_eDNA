library(tidyverse)
library(reshape2)
library(vegan)
library(betapart)
library(gridExtra)
library(grid)
library(cowplot)
library(ggplot2)
library(ggpubr)
'%ni%' <- Negate("%in%")

# load and format data

RLS_species <- read.csv("data/RLS/RLS_species_NEW.csv", sep = ";", stringsAsFactors = F, check.names = F)
RLS_species <- RLS_species[, c(1,11,7,18:2173)]
RLS_species <- reshape2::melt(RLS_species, id=c("SurveyID", "site35", "Province"))
RLS_species <- RLS_species%>%
  filter(value!=0)
RLS_species <- RLS_species[,-5]
colnames(RLS_species) <- c("SurveyID", "site35", "Province", "Species")


# add family name
library(rfishbase)
all_fishbase <- load_taxa()

RLS_species <- left_join(RLS_species, all_fishbase[,c(2,5)], by="Species")

# count family proportion at each site
count_families_site_RLS=NULL
site <- unique(RLS_species$site35)

for (i in 1:length(site)) {
  s <- site[i]
  RLS_site <- RLS_species[RLS_species$site35 == site[i],]
  p <- unique(RLS_site$Province)
  RLS_sp_site <- RLS_site%>%
    distinct(Species, .keep_all = TRUE)
  count_families <- data.frame(table(RLS_sp_site$Family))
  colnames(count_families) <- c("family", "n_species")
  count_families$n_sp_total <- sum(count_families$n_species)
  count_families$prop <- count_families$n_species / count_families$n_sp_total
  count_families <- count_families[order(count_families$prop, decreasing = TRUE),]
  count_families$site35 <- s
  count_families$province <- p
  count_families_site_RLS <- rbind(count_families_site_RLS, count_families)
}

write.csv(count_families_site_RLS, "outputs/04_family_proportion/Site/family_proportion_site_RLS.csv")
save(count_families_site_RLS, file="Rdata/family_proportion_site_RLS.rdata")


#---------------------------------------------------------------------------------------------------------------------
# random family proportion
  # prepare data
  
global_family <- RLS_species %>%
  distinct(Species, Family)
colnames(global_family) <- c("species", "family")

unique_family <- as.data.frame(unique(RLS_species$Family))
colnames(unique_family) <- "family"
unique_family<- unique_family %>%
  filter(!is.na(family))


unique_sp <- unique(RLS_species$Species)

# calculate species frequency

df_site <- data.frame(species=character(), stringsAsFactors = FALSE)
for (i in 1:length(site)) {
  df <- RLS_species[RLS_species$site35==site[i],] %>%
    distinct(Species, site35)
  colnames(df) <- c("species", site[i])
  df_site <- full_join(df_site, df, by="species")
}
rownames(df_site) <- df_site[,1]
df_site <- df_site[order(rownames(df_site)),]
df_site <- decostand(df_site[,c(-1)], "pa",na.rm = TRUE)
df_site[is.na(df_site)] <- 0
df_site <- as.data.frame(t(df_site))

prob_species <- colSums(df_site)/length(site)


# sample random family assignations for site with 1-300 species, 1000 times (long! load Rdata)
random_sp <- vector("list", 300)
random_prop_tot4 <- data.frame()
rep <- 300


for (j in seq(1:100)) {
  random_prop <- vector("list", 100)
  for (i in 1:rep) {
    random_sp[[i]] <- as.data.frame(sample(unique_sp, i, replace = TRUE, prob = prob_species))
    colnames(random_sp[[i]]) <- "species"
    random_sp[[i]] <- left_join(random_sp[[i]], global_family, by="species")
    random_sp[[i]] <- data.frame(table(random_sp[[i]]$family))
    #colnames(random_sp[[i]]) <- c("family", "n_species")
    random_sp[[i]] <- full_join(random_sp[[i]], unique_family, by=c("Var1"="family"))
    random_sp[[i]][is.na(random_sp[[i]])] <- 0
    random_sp[[i]]$n_total <- sum(random_sp[[i]]$Freq)
    random_sp[[i]]$prop <- random_sp[[i]]$Freq / random_sp[[i]]$n_total
    random_prop[[j]] <- rbind(random_prop[[j]], random_sp[[i]])
  }
  random_prop_tot4 <- rbind(random_prop_tot4, random_prop[[j]])
}

random_prop_tot <- rbind(random_prop_tot1, random_prop_tot2, random_prop_tot3, random_prop_tot4, random_prop_tot5, random_prop_tot6, random_prop_tot7, random_prop_tot8, random_prop_tot9, random_prop_tot10)
colnames(random_prop_tot) <- c("family", "n_species", "n_total", "prop")
save(random_prop_tot, file = "c:/Users/mathon/Desktop/PhD/Projets/Megafauna/Global_eDNA/Rdata/random_family_proportions_RLS.rdata")

load("c:/Users/mathon/Desktop/PhD/Projets/Megafauna/Global_eDNA/Rdata/random_family_proportions_RLS.rdata")

# calculate 2.5 and 97.5 quantile for each sample size and each family
family <- c("Acanthuridae", "Chaetodontidae", "Labridae", "Lutjanidae", "Serranidae", "Carangidae", "Pomacentridae", "Apogonidae", "Gobiidae")

CI_family <- data.frame()
quant <- data.frame(upper=numeric(300), lower=numeric(300))
for (j in 1:length(family)) {
  df <- random_prop_tot[random_prop_tot$family==family[j],]
  for (i in seq(1:300)) {
    df2 <- df[df$n_total==i,]
    quant[i,] <- quantile(df2$prop, c(.975, .025))
    quant$n_total <- seq(1:300)
    quant$family <- family[j]
  }
  CI_family <- rbind(CI_family, quant)
}

save(CI_family, file = "Rdata/CI_null_model_family_proportions_RLS.rdata")



# plot for important families

family <- c("Acanthuridae", "Chaetodontidae", "Labridae", "Lutjanidae", "Serranidae", "Carangidae", "Pomacentridae", "Apogonidae", "Gobiidae")

prop <- vector("list")
for (i in 1:length(family)) {
  fam <- count_families_site_RLS[count_families_site_RLS$family == family[i],]
  fam_CI <- CI_family[CI_family$family == family[i],]
  prop[[i]] <- ggplot(fam, aes(n_sp_total, prop))+
    geom_smooth(data=fam_CI, aes(n_total, upper), col="black", size=0.5, show.legend = FALSE)+
    geom_smooth(data=fam_CI, aes(n_total, lower), col="black", size=0.5, show.legend = FALSE)+
    geom_point(size=1)+
    scale_y_continuous(breaks = c(0, 0.1, 0.2))+
    xlim(0, 300)+
    theme(legend.position = "none")+
    theme_bw()+
    labs(title=family[i], x="", y="")+
    theme(plot.title = element_text(size = 10, face="bold"), plot.margin=unit(c(0,0.1,0,0), "cm"))
}


plot <- ggarrange(plotlist = prop, ncol=3, nrow = 3, common.legend = TRUE, legend = "top")
x.grob <- textGrob("Total number of species per site", 
                   gp=gpar(fontface="bold", col="black", fontsize=12))
y.grob <- textGrob("Proportion of species in the family in each site", 
                   gp=gpar(fontface="bold", col="black", fontsize=12), rot = 90)
plot_grid <- grid.arrange(plot, bottom=x.grob, left=y.grob)

ggsave(plot_grid, file="outputs/00_Figures_for_paper/Extended_Data/ED_family_stability_RLS.png")


