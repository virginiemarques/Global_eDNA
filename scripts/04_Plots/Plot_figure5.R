library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggpubr)
library(scales)
library(plyr)
library(dplyr)
library(conflicted)
library(gambin)
library(sads)
library(vegan)
library(bbmle)
library(nlreg)
library(MASS)
library(fitdistrplus)

# 
conflict_prefer("filter", "dplyr")
conflict_prefer("summarise", "dplyr")
'%ni%' <- Negate("%in%")

#________________________________________________________________________________________________________________________________
## Plot panel a (or run script "scripts/03_Analysis/06a_Upset_plot_motu_family_distribution.R")
#________________________________________________________________________________________________________________________________

load("Rdata/upset_plot_motus_region.rdata")
p1




#_______________________________________________________________________________________________________________________________
## Plot panel b
#_______________________________________________________________________________________________________________________________

## Run script "scripts/03_Analysis/07a_MOTUs_diversity_partitioning.R"
## Or load Rdata :
load("Rdata/beta_region.rdata")
load("Rdata/beta_site.rdata")
load("Rdata/beta_station.rdata")
load("Rdata/alpha_station.rdata")

Region <- unique(beta_site$province)

gamma_global = 2023

all_data <- data.frame()

all_data[1,1] <- "gamma_global"
all_data[1,2] <- "all"
all_data[1,3] <- 100
all_data[1,4] <- 0
colnames(all_data) <- c("scale", "region", "mean", "sd")

all_data[2,1] <- "Beta_inter-region"
all_data[2,2] <- "all"
all_data[2,3] <- 74
all_data[2,4] <- 0

for (i in 1:length(Region)) {
  all_data[i+2,1] <- "mean_Beta_inter-site"
  all_data[i+2,2] <- Region[i]
  all_data[i+2,3] <- beta_site %>%
    filter(province==Region[i])%>%
    summarize(n=beta*100/gamma_global)
  all_data[i+2,4] <- NA
}

for (i in 1:length(Region)) {
  all_data[i+7,1] <- "mean_Beta_inter-station"
  all_data[i+7,2] <- Region[i]
  all_data[i+7,3] <- beta_station %>%
    filter(province==Region[i])%>%
    summarize(n=mean(beta)*100/gamma_global)
  all_data[i+7,4] <- beta_station %>%
    filter(province==Region[i])%>%
    summarize(n=sd(beta)*100/gamma_global)
}

for (i in 1:length(Region)) {
  all_data[i+12,1] <- "mean_alpha_station"
  all_data[i+12,2] <- Region[i]
  all_data[i+12,3] <- alpha_station %>%
    filter(province==Region[i])%>%
    summarize(n=mean(motu)*100/gamma_global)
  all_data[i+12,4] <- alpha_station %>%
    filter(province==Region[i])%>%
    summarize(n=sd(motu)*100/gamma_global)
}
all_data[9,4] <- 0
all_data[4,3] <- 0.1
all_data[2,4] <- NA


all_data <- all_data %>%
  group_by(scale)%>%
  mutate(position= rank(mean))

all_data <- all_data%>%
  group_by(scale)%>%
  mutate(label= case_when(
    scale=="Beta_inter-region" ~ 74,
    scale=="mean_Beta_inter-site" ~ 14.5,
    scale=="mean_Beta_inter-station" ~ 6.2,
    scale=="mean_alpha_station" ~ 5.3,
    scale=="gamma_global" ~ 100
  ))



p2 <-ggplot(all_data[3:17,], aes(x=reorder(scale,mean), y=mean, fill=region, group= position))+
  geom_col(stat = "identity", position=position_dodge(), width=0.5)+
  scale_fill_manual(values=c("#a6611a", "#E5A729","#b2182b", "#80cdc1", "#015462"))+
  geom_errorbar(data=all_data[3:17,], aes(ymin=mean-sd, ymax=mean+sd), width=.3,
                position=position_dodge(0.5), size=0.5)+
  geom_errorbar(aes(ymin=label, ymax=label), width=0.5, size=1)+
  ylim(0,100)+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        plot.title = element_text(size = 12, face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin=unit(c(0,0,0,0), "cm"))+
  coord_flip()

p1 <-ggplot(all_data[1:2,], aes(x=reorder(scale,mean), y=mean, fill=region, group= position))+
  geom_bar(stat = "identity", position=position_dodge(), width=0.2, fill="grey")+
  geom_errorbar(aes(ymin=label, ymax=label), width=0.4, size=1)+
  ylim(0,100)+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        plot.title = element_text(size = 12, face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin=unit(c(0,0,0,0), "cm"))+
  coord_flip()

Fig3 <- ggarrange(p1, p2, nrow=2, heights = c(1, 2.5))

ggsave(Fig3, file="outputs/00_Figures_for_paper/Figure5b.png")


### Complete figure 3 is assembled on powerpoint "outputs/00Figures_fot_paper/Figure3.ppt"

