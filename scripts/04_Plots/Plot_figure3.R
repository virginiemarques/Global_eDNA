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


load("Rdata/02_clean_all.Rdata")

# Remove estuary stations and deep niskin station
df_all_filters <- df_all_filters %>%
  filter(station %ni% c("estuaire_rio_don_diego_1", "estuaire_rio_don_diego_2", "estuaire_rio_don_diego_3")) %>%
  filter(sample_method !="niskin" & region!="East_Pacific" & comment %ni% c("Distance decay 600m", "Distance decay 300m") & station!="glorieuse_distance_300m")%>%
  filter(project != "SEAMOUNTS") %>% 
  filter(habitat_type %ni% c("BAIE", "Sommet"))

unique(df_all_filters$region)

# the motus 
df_regions <- split(df_all_filters, df_all_filters$region)
df_site <- split(df_all_filters, df_all_filters$site)


pal <- c("#8AAE8A", "#E5A729", "#C67052", "#4F4D1D", "#863b34") 
pal

# Construct the inital matrix - add some important infos as a side
motus_fam_region <-  df_all_filters %>%
  group_by(region) %>%
  dplyr::summarize(n_motus = n_distinct(sequence), 
                   n_family = n_distinct(new_family_name)) %>%
  as.data.frame() %>%
  dplyr::arrange(n_motus)

# metadata 
metadata1 <- df_all_filters %>%
  distinct(region) %>% 
  dplyr::mutate(sets = region) %>%
  select(sets, region) %>%
  dplyr::mutate(color = case_when(
    region == "South_West_Pacific" ~ pal[5],
    region == "Central_Pacific" ~ pal[4],
    region == "West_Indian" ~ pal[3],
    region == "Caribbean" ~ pal[2], 
    region == "Central_Indo_Pacific" ~ pal[1]
  )) %>%
  as.data.frame() %>%
  left_join(., motus_fam_region)

# MOTUs
matrix_motus <- df_all_filters %>%
  distinct(sequence, new_scientific_name_ncbi) %>%
  dplyr::mutate(`Central_Indo_Pacific` = ifelse(sequence %in% df_regions$`Central_Indo_Pacific`$sequence, 1, 0), 
                Central_Pacific = ifelse(sequence %in% df_regions$Central_Pacific$sequence, 1, 0),
                Caribbean = ifelse(sequence %in% df_regions$Caribbean$sequence, 1, 0),
                West_Indian = ifelse(sequence %in% df_regions$West_Indian$sequence, 1, 0),
                South_West_Pacific = ifelse(sequence %in% df_regions$South_West_Pacific$sequence, 1, 0)) %>%
  as.data.frame()


# Plot MOTUs distribution in regions
p1 <- upset(matrix_motus, 
            mb.ratio = c(0.7, 0.3),
            order.by = c("freq"),
            mainbar.y.label = "Number of MOTUs", 
            sets.x.label = "Number of MOTUs", 
            text.scale = c(1.2, 1.2, 1.2,1.2,1.2,1.2), 
            # Color bar 
            sets.bar.color=c("#8AAE8A", "#863b34", "#E5A729", "#C67052", "#4F4D1D"), 
            # Color matrix
            set.metadata = list(
              data = metadata1,
              plots = list(list(
                type = "matrix_rows",
                column = "region", 
                colors = c(`Central_Indo_Pacific` = pal[1], Caribbean =  pal[2], West_Indian = pal[3], Central_Pacific =  pal[4], South_West_Pacific = pal[5]),
                alpha = 0.3
              ))
            ))
p1

save(p1, file = "Rdata/upset_plot_motus_region.rdata")
png('outputs/00_Figures_for_paper/Figure3a.png', width = 7, height=3.6, units = "in", res=300)
p1
grid.text("Regions - MOTUs",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()


#_______________________________________________________________________________________________________________________________
## Plot panel b
#_______________________________________________________________________________________________________________________________

## Run script "scripts/03_Analysis/07a_MOTUs_diversity_partitioning.R"
## Or load Rdata :
load("Rdata/beta_region.rdata")
load("Rdata/beta_site.rdata")
load("Rdata/beta_station.rdata")
load("Rdata/alpha_station.rdata")

gamma_global = 2160

all_data <- data.frame(scale=character(), region=character(), mean=numeric(), sd=numeric())

all_data[1,1] <- "gamma_global"
all_data[1,2] <- "all"
all_data[1,3] <- 100
all_data[1,4] <- 0

all_data[2,1] <- "Beta_inter-region"
all_data[2,2] <- "all"
all_data[2,3] <- (beta_region$beta)*100/gamma_global
all_data[2,4] <- 0

for (i in 1:length(Region)) {
  all_data[i+2,1] <- "mean_Beta_inter-site"
  all_data[i+2,2] <- Region[i]
  all_data[i+2,3] <- beta_site %>%
    filter(region==Region[i])%>%
    summarize(n=beta*100/gamma_global)
  all_data[i+2,4] <- NA
}

for (i in 1:length(Region)) {
  all_data[i+7,1] <- "mean_Beta_inter-station"
  all_data[i+7,2] <- Region[i]
  all_data[i+7,3] <- beta_station %>%
    filter(region==Region[i])%>%
    summarize(n=mean(beta)*100/gamma_global)
  all_data[i+7,4] <- beta_station %>%
    filter(region==Region[i])%>%
    summarize(n=sd(beta)*100/gamma_global)
}

for (i in 1:length(Region)) {
  all_data[i+12,1] <- "mean_alpha_station"
  all_data[i+12,2] <- Region[i]
  all_data[i+12,3] <- alpha_station %>%
    filter(region==Region[i])%>%
    summarize(n=mean(motu)*100/gamma_global)
  all_data[i+12,4] <- alpha_station %>%
    filter(region==Region[i])%>%
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
    scale=="Beta_inter-region" ~ 74.5,
    scale=="mean_Beta_inter-site" ~ 14,
    scale=="mean_Beta_inter-station" ~ 7,
    scale=="mean_alpha_station" ~ 4.5,
    scale=="gamma_global" ~ 100
  ))

p2 <-ggplot(all_data[3:17,], aes(x=reorder(scale,mean), y=mean, fill=region, group= position))+
  geom_col(stat = "identity", position=position_dodge(), width=0.5)+
  scale_fill_manual(values=c("#E5A729", "#8AAE8A", "#4F4D1D",  "#863b34", "#C67052"))+
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

ggsave(Fig3, "outputs/00_Figures_for_paper/Figure3b.png", dpi=500)


### Complete figure 3 is assembled on powerpoint "outputs/00Figures_fot_paper/Figure3.ppt"

