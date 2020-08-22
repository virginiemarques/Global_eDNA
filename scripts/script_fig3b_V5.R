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

library(nlreg)

# repertoire David

setwd("/Users/davidmouillot/Documents/articles/en cours/Global eDNA")


# fit log-log models

load("Rdata/rarete_motu_station.rdata")
tab=as.data.frame(motu_station)

head(tab)

load("Rdata/rarete_species_transects.rdata")
tab2 <- species_transects

head(tab2)

# power law

edna.po <- nlreg(n_motus ~ b0*n^b1, weights = ~ ( 1+n^g )^2,start = c(b0 = 426, b1 = -1, g = 1),
                     data = tab)

AIC(edna.po)

rls.po <- nlreg(Freq ~ b0*occ_RLS^b1, weights = ~ ( 1+occ_RLS^g )^2,start = c(b0 = 256, b1 = -1, g = 1),
                 data = tab2)

AIC(rls.po)

# log series

edna.ls <- nlreg(n_motus ~ b0*(1/n)*exp(-b2*n), weights = ~ ( 1+n^g )^2,start = c(b0 = 426, b2 = 1, g = 1),
                 data = tab)

AIC(edna.ls)

rls.ls <- nlreg(Freq ~ b0*(1/occ_RLS)*exp(-b2*occ_RLS), weights = ~ ( 1+occ_RLS^g )^2,start = c(b0 = 256, b2 = 1, g = 1),
                data = tab2)

AIC(rls.ls)

# power bended

edna.pb <- nlreg(n_motus ~ b0*(n^b1)*exp(-b2*n), weights = ~ ( 1+n^g )^2,start = c(b0 = 426, b1=-1,b2 = 1, g = 1),
                 data = tab)

AIC(edna.pb)

rls.pb <- nlreg(Freq ~ b0*(occ_RLS^b1)*exp(-b2*occ_RLS), weights = ~ ( 1+occ_RLS^g )^2,start = c(b0 = 256, b1=-1,b2 = 1, g = 1),
                data = tab2)

AIC(rls.pb)


AICtab(edna.po, edna.ls, edna.pb, weights=TRUE)

AICtab(rls.po, rls.ls, rls.pb, weights=TRUE)

# predict for each model
tab$pb <- predict(edna.pb)
tab$po <- predict(edna.po)
tab$ls <- predict(edna.ls)

tab2$pb <- predict(rls.pb)
tab2$po <- predict(rls.po)
tab2$ls <- predict(rls.ls)



# plot figure 4a edna
ggplot(tab, aes(x=log10(n), y=log10(n_motus)))+
  geom_point(colour="#d2981a", size=2, show.legend = TRUE)+
  geom_line(aes(x=log10(n), y=log10(po)), linetype = "dashed", size = 0.8)+
  geom_line(aes(x=log10(n), y=log10(pb)), linetype = "solid", size = 0.8)+
  geom_line(aes(x=log10(n), y=log10(ls)), linetype = "dotted", size = 0.8)+
  xlim(0,2)+
  ylim(0,3)+
  annotate(geom="text", x=2, y=3, label="eDNA MOTUs ~ stations", hjust=1, size=4, colour="#d2981a") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        plot.title = element_text(size=12, face = "bold"))+
  labs(x="log10(Number of station)",y="log10(Number of MOTUs)")

ggsave("outputs/Figures papier/Figure4a.png")



# plot figure 4b rls
ggplot(tab2, aes(x=log10(occ_RLS), y=log10(Freq)))+
  geom_point(colour = "darkgrey", size=2, show.legend = TRUE)+
  geom_line(aes(x=log10(occ_RLS), y=log10(po)), linetype = "dashed", size = 0.8)+
  geom_line(aes(x=log10(occ_RLS), y=log10(pb)), linetype = "solid", size = 0.8)+
  geom_line(aes(x=log10(occ_RLS), y=log10(ls)), linetype = "dotted", size = 0.8)+
  xlim(0,3)+
  ylim(0,3)+
  annotate(geom="text", x=3, y=3, label="RLS species ~ transects", hjust=1, size=4) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        plot.title = element_text(size=12, face = "bold"))+
  labs(x="log10(Number of transect)",y="log10(Number of species)")

ggsave("outputs/Figures papier/Figure4b.png")
