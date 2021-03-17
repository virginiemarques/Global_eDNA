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
library(broom)

conflict_prefer("summarise", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("ggsave", "ggplot2")

## For panel a : Run script "scripts/03_Analysis/06b_Histograms_motus_families_ranks.R
## Or load Rdata :
load("Rdata/rarete_motu_station.rdata")



## For panel b load Rdata :


RLS_species <- read.csv("data/RLS/RLS_species_NEW.csv", sep = ";", stringsAsFactors = F, check.names = F)
RLS_species <- RLS_species[, c(1,11,7,18:2173)]
RLS_species <- reshape2::melt(RLS_species, id=c("SurveyID", "site35", "Province"))
RLS_species <- RLS_species%>%
  filter(value!=0)
RLS_species <- RLS_species[,-5]
colnames(RLS_species) <- c("SurveyID", "site35", "Province", "Species")

species_transects <- RLS_species %>%
  group_by(Species) %>%
  summarise(n = n_distinct(SurveyID)) %>%
  ungroup() %>%
  group_by(n) %>%
  summarise(Freq = n_distinct(Species))
colnames(species_transects) <- c("occ_RLS", "Freq")
save(species_transects, file="Rdata/rarete_species_transects.rdata")


# fit log-log models

tab=as.data.frame(motu_station)
tab$logn=log10(tab$n)
tab$logmotus=log10(tab$n_motus)
tab$mlogn=-log10(tab$n)

tab2 <- species_transects
tab2$logocc <- log10(tab2$occ_RLS)
tab2$logfreq <- log10(tab2$Freq)
tab2$mlogocc <- -log10(tab2$occ_RLS)



# fit non linear regression by maximum likelihood

  # power law (pareto)

    # eDNA

linmod.motus <- lm(log10(n_motus) ~ log10(n), data = tab)

edna.po <- mle2(n_motus ~ dnbinom(mu=b0*n^b1, size=exp(logdisp)),data=tab,
                start=list(b0 = 10^(coef(linmod.motus)[1]), b1 = coef(linmod.motus)[2],logdisp=0))

confint(edna.po,method="quad")
edna.po_tidy <- (as.data.frame(tidy(edna.po)))
# recuperer parametre
edna.po.AIC <- AIC(edna.po)

edna.po.int <- edna.po_tidy %>%
  filter(term=="b0")%>%
  select(estimate) 
  
edna.po.int_ci <- edna.po_tidy %>%
  filter(term=="b0")%>%
  select(std.error)

edna.po.sl <- edna.po_tidy %>%
  filter(term=="b1")%>%
  select(estimate)

edna.po.sl_ci <- edna.po_tidy %>%
  filter(term=="b1")%>%
  select(std.error)

      # RLS

linmod.rls <- lm(log10(Freq) ~ log10(occ_RLS), data = tab2)

rls.po <- mle2(Freq ~ dnbinom(mu=b0*occ_RLS^b1, size=exp(logdisp)),data=tab2,
               start=list(b0 = 10^(coef(linmod.rls)[1]), b1 = coef(linmod.rls)[2],logdisp=0))

confint(rls.po,method="quad")

rls.po_tidy <- (as.data.frame(tidy(rls.po)))
# recuperer parametre
rls.po.AIC <- AIC(rls.po)

rls.po.int <- rls.po_tidy %>%
  filter(term=="b0")%>%
  select(estimate) 

rls.po.int_ci <- rls.po_tidy %>%
  filter(term=="b0")%>%
  select(std.error)

rls.po.sl <- rls.po_tidy %>%
  filter(term=="b1")%>%
  select(estimate)

rls.po.sl_ci <- rls.po_tidy %>%
  filter(term=="b1")%>%
  select(std.error)

# log series
    #eDNA
edna.ls <- mle2(n_motus ~ dnbinom(mu=b0*(1/n)*exp(-b2*n), size=exp(logdisp)),data=tab,control=list(maxit=1E5,trace=0),
                start=list(b0 = 426, b2 = 0,logdisp=0))

confint(edna.ls,method="quad")

edna.ls_tidy <- (as.data.frame(tidy(edna.ls)))
# recuperer parametre
edna.ls.AIC <- AIC(edna.ls)

edna.ls.int <- edna.ls_tidy %>%
  filter(term=="b0")%>%
  select(estimate) 

edna.ls.int_ci <- edna.ls_tidy %>%
  filter(term=="b0")%>%
  select(std.error)

edna.ls.sl <- -1

edna.ls.bend <- edna.ls_tidy %>%
  filter(term=="b2")%>%
  select(estimate) 

edna.ls.bend_ci <- edna.ls_tidy %>%
  filter(term=="b2")%>%
  select(std.error)

    # RLS
rls.ls <- mle2(Freq ~ dnbinom(mu=b0*(1/occ_RLS)*exp(-b2*occ_RLS), size=exp(logdisp)),data=tab2,
               start=list(b0 = 256, b2 = 0,logdisp=0))

confint(rls.ls,method="quad")
rls.ls_tidy <- (as.data.frame(tidy(rls.ls)))
# recuperer parametre
rls.ls.AIC <- AIC(rls.ls)

rls.ls.int <- rls.ls_tidy %>%
  filter(term=="b0")%>%
  select(estimate) 

rls.ls.int_ci <- rls.ls_tidy %>%
  filter(term=="b0")%>%
  select(std.error)

rls.ls.sl <- -1

rls.ls.bend <- rls.ls_tidy %>%
  filter(term=="b2")%>%
  select(estimate) 

rls.ls.bend_ci <- rls.ls_tidy %>%
  filter(term=="b2")%>%
  select(std.error)


# power bended
    
    # eDNA
edna.pb <- mle2(n_motus ~ dnbinom(mu=b0*(n^b1)*exp(-b2*n), size=exp(logdisp)),data=tab,control=list(maxit=1E5,trace=0),
                start=list(b0 = 10^(coef(linmod.motus)[1]), b1=coef(linmod.motus)[2],b2 = 0,logdisp=0))

confint(edna.pb,method="quad")

edna.pb_tidy <- (as.data.frame(tidy(edna.pb)))
# recuperer parametre
edna.pb.AIC <- AIC(edna.pb)

edna.pb.int <- edna.pb_tidy %>%
  filter(term=="b0")%>%
  select(estimate) 

edna.pb.int_ci <- edna.pb_tidy %>%
  filter(term=="b0")%>%
  select(std.error)

edna.pb.sl <- edna.pb_tidy %>%
  filter(term=="b1")%>%
  select(estimate)

edna.pb.sl_ci <- edna.pb_tidy %>%
  filter(term=="b1")%>%
  select(std.error) 

edna.pb.bend <- edna.pb_tidy %>%
  filter(term=="b2")%>%
  select(estimate) 

edna.pb.bend_ci <- edna.pb_tidy %>%
  filter(term=="b2")%>%
  select(std.error)

    # RLS

rls.pb <- mle2(Freq ~ dnbinom(mu=b0*(occ_RLS^b1)*exp(-b2*occ_RLS), size=exp(logdisp)),data=tab2,
               start=list(b0 = 10^(coef(linmod.rls)[1]), b1=coef(linmod.rls)[2],b2 = 0,logdisp=0))

confint(rls.pb,method="quad")
rls.pb_tidy <- (as.data.frame(tidy(rls.pb)))
# recuperer parametre
rls.pb.AIC <- AIC(rls.pb)

rls.pb.int <- rls.pb_tidy %>%
  filter(term=="b0")%>%
  select(estimate) 

rls.pb.int_ci <- rls.pb_tidy %>%
  filter(term=="b0")%>%
  select(std.error)

rls.pb.sl <- rls.pb_tidy %>%
  filter(term=="b1")%>%
  select(estimate)

rls.pb.sl_ci <- rls.pb_tidy %>%
  filter(term=="b1")%>%
  select(std.error) 

rls.pb.bend <- rls.pb_tidy %>%
  filter(term=="b2")%>%
  select(estimate) 

rls.pb.bend_ci <- rls.pb_tidy %>%
  filter(term=="b2")%>%
  select(std.error)



# Figures

# predict for each model
tab$pb <- predict(edna.pb)
tab$po <- predict(edna.po)
tab$ls <- predict(edna.ls)

tab2$pb <- predict(rls.pb)
tab2$po <- predict(rls.po)
tab2$ls <- predict(rls.ls)



# plot figure 4a edna --> automate parameters
edna_ls <- ggplot(tab, aes(x=logn, y=logmotus))+
  geom_point(colour="#d2981a", size=2, show.legend = TRUE)+
  geom_line(aes(x=logn, y=log10(ls)), linetype = "solid", size = 0.8)+
  xlim(0,2)+
  ylim(0,3)+
  annotate(geom="text", x=2, y=3, label="Log-series", hjust=1, size=3.5, fontface = "bold")+
  annotate(geom="text", x=2, y=2.5, label=paste("slope= -1\nBending=",edna.ls.bend, "\nCI=(",edna.ls.bend_ci,")\nAIC=",edna.ls.AIC), hjust=1, size=3.5)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        plot.title = element_text(size=12, face = "bold"))+
  labs(x="",y="log10(Nb of MOTUs)")

edna_pb <- ggplot(tab, aes(x=logn, y=logmotus))+
  geom_point(colour="#d2981a", size=2, show.legend = TRUE)+
  geom_line(aes(x=logn, y=log10(pb)), linetype = "solid", size = 0.8)+
  xlim(0,2)+
  ylim(0,3)+
  annotate(geom="text", x=2, y=3, label="Pareto-bended", hjust=1, size=3.5, fontface="bold")+
  annotate(geom="text", x=2, y=2.4, label=paste("slope=",edna.pb.sl,"\nCI=(",edna.pb.sl_ci,")\nBending=",edna.pb.bend,"\nCI=(",edna.pb.bend_ci,")\nAIC=",edna.pb.AIC), hjust=1, size=3.5)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        plot.title = element_text(size=12, face = "bold"))+
  labs(x="",y="")


edna_po <- ggplot(tab, aes(x=logn, y=logmotus))+
  geom_point(colour="#d2981a", size=2, show.legend = TRUE)+
  geom_line(aes(x=logn, y=log10(po)), linetype = "solid", size = 0.8)+
  xlim(0,2)+
  ylim(0,3)+
  annotate(geom="text", x=2, y=3, label="Pareto", hjust=1, size=3.5, fontface="bold")+
  annotate(geom="text", x=2, y=2.6, label=paste("slope=",edna.po.sl,"\nCI=(",edna.po.sl_ci,")\nAIC=",edna.po.AIC), hjust=1, size=3.5)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        plot.title = element_text(size=12, face = "bold"))+
  labs(x="",y="")

plot <- ggarrange(edna_ls, edna_po, edna_pb, ncol = 3, nrow=1)
x.grob <- textGrob("log10(Nb of stations)", 
                   gp=gpar(fontface="bold", col="black", fontsize=10), vjust = -0.5)


a <- grid.arrange(plot, bottom=x.grob)

# plot figure 4b rls
rls_ls <- ggplot(tab2, aes(x=logocc, y=logfreq))+
  geom_point(colour = "darkgrey", size=2, show.legend = TRUE)+
  geom_line(aes(x=logocc, y=log10(ls)), linetype = "solid", size = 0.8)+
  xlim(0,3)+
  ylim(0,3)+
  annotate(geom="text", x=3, y=3, label="Log-series", hjust=1, size=3.5, fontface="bold") +
  annotate(geom="text", x=3, y=2.5, label=paste("slope= -1\nBending=",rls.ls.bend, "\nCI=(",rls.ls.bend_ci,")\nAIC=",rls.ls.AIC), hjust=1, size=3.5) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        plot.title = element_text(size=12, face = "bold"))+
  labs(x="",y="log10(Number of species)")

rls_pb <- ggplot(tab2, aes(x=logocc, y=logfreq))+
  geom_point(colour = "darkgrey", size=2, show.legend = TRUE)+
  geom_line(aes(x=logocc, y=log10(pb)), linetype = "solid", size = 0.8)+
  xlim(0,3)+
  ylim(0,3)+
  annotate(geom="text", x=3, y=3, label="Pareto-bended", hjust=1, size=3.5, fontface="bold") +
  annotate(geom="text", x=3, y=2.4, label=paste("slope=",rls.pb.sl,"\nCI=(",rls.pb.sl_ci,")\nBending=",rls.pb.bend,"\nCI=(",rls.pb.bend_ci,")\nAIC=",rls.pb.AIC), hjust=1, size=3.5) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        plot.title = element_text(size=12, face = "bold"))+
  labs(x="",y="")

rls_po <- ggplot(tab2, aes(x=logocc, y=logfreq))+
  geom_point(colour = "darkgrey", size=2, show.legend = TRUE)+
  geom_line(aes(x=logocc, y=log10(po)), linetype = "solid", size = 0.8)+
  xlim(0,3)+
  ylim(0,3)+
  annotate(geom="text", x=3, y=3, label="Pareto", hjust=1, size=3.5, fontface="bold") +
  annotate(geom="text", x=3, y=2.6, label=paste("slope=",rls.po.sl,"\nCI=(",rls.po.sl_ci,")\nAIC=",rls.po.AIC), hjust=1, size=3.5) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        plot.title = element_text(size=12, face = "bold"))+
  labs(x="",y="")

plot2 <- ggarrange(rls_ls, rls_po, rls_pb, ncol = 3, nrow=1)
x.grob <- textGrob("log10(Nb of transects)", 
                   gp=gpar(fontface="bold", col="black", fontsize=10), vjust = -0.5)


b <- grid.arrange(plot2, bottom=x.grob)


# Plot all together

Fig4 <- ggarrange(a, b, nrow=2, labels = c("a", "b"), label.x =0, label.y=1)
Fig4
ggsave("outputs/00_Figures_for_paper/Figure4.png", width = 8, height = 8, unit = "in")

