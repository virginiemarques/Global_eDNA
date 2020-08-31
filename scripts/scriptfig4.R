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

# repertoire David
setwd("/Users/davidmouillot/Documents/articles/en cours/Global eDNA")

# fit log-log models

load("rarete_motu_station.rdata") # VM: cant open it here. I added a line with the proper path. 
load("Rdata/rarete_motu_station.rdata")

tab=as.data.frame(motu_station)

#logn=log10(tab$n)
#logmotus=log10(tab$n_motus)
#mlogn=-log10(tab$n)
#
#tab=cbind(tab,logn,logmotus,mlogn)

head(tab)

load("rarete_species_transects.rdata") # Same
load("Rdata/rarete_species_transects.rdata")

tab2 <- species_transects

# logocc=log10(tab2$occ_RLS)
# logfreq=log10(tab2$Freq)
# mlogocc=-log10(tab2$occ_RLS)
# 
# tab2=cbind(tab2,logocc,logfreq,mlogocc)

# -------- # Modif VM: add the logs sur tab et tab2 properly # 
tab$logn=log10(tab$n)
tab$logmotus=log10(tab$n_motus)
tab$mlogn=-log10(tab$n)

tab2$logocc <- log10(tab2$occ_RLS)
tab2$logfreq <- log10(tab2$Freq)
tab2$mlogocc <- -log10(tab2$occ_RLS)
# -------- #  End of modif VM

head(tab2)

# distribution

hist(tab$n_motus)

motus.poisson=fitdist(tab$n_motus, "pois")

motus.nb=fitdist(tab$n_motus, "nbinom")

gofstat(list(motus.poisson, motus.nb),fitnames = c("Poisson", "Negative Binomial"))



hist(tab2$Freq)

rls.poisson=fitdist(tab2$Freq, "pois")

rls.nb=fitdist(tab2$Freq, "nbinom")

gofstat(list(rls.poisson, rls.nb),fitnames = c("Poisson", "Negative Binomial"))


# fit non linear regression by maximum likelihood

# power law

# initial parameters

linmod.motus <- lm(log10(n_motus) ~ log10(n), data = tab)

coef(linmod.motus)

edna.po <- mle2(n_motus ~ dnbinom(mu=b0*n^b1, size=exp(logdisp)),data=tab,
                start=list(b0 = 10^(coef(linmod.motus)[1]), b1 = coef(linmod.motus)[2],logdisp=0))

confint(edna.po,method="quad")

plot(log10(n_motus) ~ log10(n), tab)
lines(log10(tab$n), log10(predict(edna.po)), col = 'blue')


linmod.rls <- lm(log10(Freq) ~ log10(occ_RLS), data = tab2)

coef(linmod.rls)

rls.po <- mle2(Freq ~ dnbinom(mu=b0*occ_RLS^b1, size=exp(logdisp)),data=tab2,
               start=list(b0 = 10^(coef(linmod.rls)[1]), b1 = coef(linmod.rls)[2],logdisp=0))

confint(rls.po,method="quad")

plot(log10(Freq) ~ log10(occ_RLS), tab2)
lines(log10(tab2$occ_RLS), log10(predict(rls.po)), col = 'blue')


# log series

edna.ls <- mle2(n_motus ~ dnbinom(mu=b0*(1/n)*exp(-b2*n), size=exp(logdisp)),data=tab,control=list(maxit=1E5,trace=0),
                start=list(b0 = 426, b2 = 0,logdisp=0))

plot(log10(n_motus) ~ log10(n), tab)
lines(log10(tab$n), log10(predict(edna.ls)), col = 'red')

confint(edna.ls,method="quad")

rls.ls <- mle2(Freq ~ dnbinom(mu=b0*(1/occ_RLS)*exp(-b2*occ_RLS), size=exp(logdisp)),data=tab2,
               start=list(b0 = 256, b2 = 0,logdisp=0))

plot(log10(Freq) ~ log10(occ_RLS), tab2)
lines(log10(tab2$occ_RLS), log10(predict(rls.ls )), col = 'red')

confint(rls.ls,method="quad")



# power bended

edna.pb <- mle2(n_motus ~ dnbinom(mu=b0*(n^b1)*exp(-b2*n), size=exp(logdisp)),data=tab,control=list(maxit=1E5,trace=0),
                start=list(b0 = 10^(coef(linmod.motus)[1]), b1=coef(linmod.motus)[2],b2 = 0,logdisp=0))


confint(edna.pb,method="quad")

plot(log10(n_motus) ~ log10(n), tab)
lines(log10(tab$n), log10(predict(edna.pb)), col = 'green')


rls.pb <- mle2(Freq ~ dnbinom(mu=b0*(occ_RLS^b1)*exp(-b2*occ_RLS), size=exp(logdisp)),data=tab2,
               start=list(b0 = 10^(coef(linmod.rls)[1]), b1=coef(linmod.rls)[2],b2 = 0,logdisp=0))

confint(rls.pb,method="quad")

plot(log10(Freq) ~ log10(occ_RLS), tab2)
lines(log10(tab2$occ_RLS), log10(predict(rls.pb )), col = 'green')


# model comparisons

AICtab(edna.po, edna.ls, edna.pb, weights=TRUE)

anova(edna.ls,edna.pb)
anova(edna.po,edna.pb)

AICtab(rls.po, rls.ls, rls.pb, weights=TRUE)

anova(rls.po,rls.pb)
anova(rls.ls,rls.pb)


# Figures

# predict for each model
tab$pb <- predict(edna.pb)
tab$po <- predict(edna.po)
tab$ls <- predict(edna.ls)

tab2$pb <- predict(rls.pb)
tab2$po <- predict(rls.po)
tab2$ls <- predict(rls.ls)



# plot figure 4a edna
edna_ls <- ggplot(tab, aes(x=logn, y=logmotus))+
  geom_point(colour="#d2981a", size=2, show.legend = TRUE)+
  geom_line(aes(x=logn, y=log10(ls)), linetype = "solid", size = 0.8)+
  xlim(0,2)+
  ylim(0,3)+
  annotate(geom="text", x=2, y=3, label="Log-series", hjust=1, size=3.5, fontface = "bold")+
  annotate(geom="text", x=2, y=2.5, label="slope= -1\nBending=0.04\nCI=[0.038:0.046]\nAIC=307", hjust=1, size=3.5)+
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
  annotate(geom="text", x=2, y=2.4, label="slope= -1.03\nCI=[-1.12:-0.94]\nBending=0.04\nCI=[0.029:0.053]\nAIC=309", hjust=1, size=3.5)+
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
  annotate(geom="text", x=2, y=2.6, label="slope= -1.64\nCI=[-1.69:-1.60]\nAIC=333", hjust=1, size=3.5)+
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
  annotate(geom="text", x=3, y=2.5, label="slope= -1\nBending=0.004\nCI=[-6e-04:3e-04]\nAIC=1078", hjust=1, size=3.5) +
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
  annotate(geom="text", x=3, y=2.4, label="slope= -0.96\nCI=[-9.9:-0.94]\nBending=0.0005\nCI=[-0.003:0.001]\nAIC=1076", hjust=1, size=3.5) +
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
  annotate(geom="text", x=3, y=2.6, label="slope= -0.97\nCI=[-0.99:-0.94]\nAIC=1074", hjust=1, size=3.5) +
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

ggarrange(a, b, nrow=2, labels = c("a", "b"), label.x =0, label.y=1)

conflict_prefer("ggsave", "ggplot2")

ggsave("outputs/Figures papier/Figure4.png")

