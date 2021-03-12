library(tidyverse)
library(reshape2)
library(vegan)
library(betapart)
library(ggpubr)


load("Rdata/02-clean-data.Rdata")
'%ni%' <- Negate("%in%")
#Remove estuary stations and deep niskin station
df_all_filters <- df_all_filters %>%
  filter(province %in% c("Western_Indian_Ocean", "Southeast_Polynesia", "Tropical_Northwestern_Atlantic", "Western_Coral_Triangle", "Tropical_Southwestern_Pacific"))%>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3", "glorieuse_distance_300m")) %>%
  filter(sample_method !="niskin" & comment %ni% c("Distance decay 600m", "Distance decay 300m"))%>%
  filter(project != "Curacao") %>%
  filter(habitat=="marine")%>%
  filter(habitat_type %ni% c("BAIE"))




Station <- unique(df_all_filters$station)
df_sample=vector("list", length(Station))
betastation <- data.frame()


for (i in 1:length(Station)) {
  df <- df_all_filters[df_all_filters$station==Station[i],]
  sample <- unique(df$sample_name_all_pcr)
  region <- unique(df$province)
  df_sample[[i]] <- data.frame(motu=character(), stringsAsFactors = FALSE)
  for (j in 1:length(sample)) {
    df2 <- df[df$sample_name_all_pcr==sample[j],] %>%
      distinct(sequence, sample_name_all_pcr)
    colnames(df2) <- c("motu", sample[j])
    df_sample[[i]] <- full_join(df_sample[[i]], df2, by="motu")
  }
  rownames(df_sample[[i]]) <- df_sample[[i]][,1]
  df_sample[[i]] <- decostand(df_sample[[i]][,c(-1)], "pa",na.rm = TRUE)
  df_sample[[i]][is.na(df_sample[[i]])] <- 0
  df_sample[[i]] <- as.data.frame(t(df_sample[[i]]))
  
  b <- betapart.core(df_sample[[i]])
  beta <- beta.multi(b, "jaccard")
  betastation[i,1]<- beta$beta.JAC
  betastation[i,2] <- Station[i]
  betastation[i,3] <- region
  colnames(betastation) <- c("beta", "station", "region")
}


ggplot(betastation, aes(region, beta))+
  geom_boxplot()+
  theme_bw()+
  ylim(0,1)+
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x="", y="Jaccard dissimilarity (motus)")

ggsave("outputs/00_Figures_for_paper/Extended_Data/ED_beta_inter-replicate.png")
