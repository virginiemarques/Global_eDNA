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

## Plot Figure 4 : dbRDA

# Panel a : eDNA

  # load rdata

load("Rdata/dbRDA_data_edna.rdata")

  # plot

a <- ggplot(site_scores_region, aes(x= CAP1, y = CAP2)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  geom_vline(xintercept = 0, lty = 2, col = "grey") +
  geom_encircle(data=site_scores_region[site_scores_region$province!="Southeast_Polynesia",], aes(group = province, linetype = province, fill= province), s_shape = 1, expand = 0,
                alpha = 0.4, show.legend = FALSE) + # hull area 
  geom_point(data = family_scores, aes(x= CAP1,y = CAP2), col = "grey", alpha = 0.5, cex = 0.5) +
  geom_point(data= family_scores_diff75, aes(x= CAP1, y=CAP2), col = "black", alpha = 1, cex = 0.5) +
  geom_point(aes(pch = province, fill = province), cex = 4, col = "black") +
  scale_fill_manual(values = c("#a6611a", "#E5A729","#b2182b", "#80cdc1", "#015462"),
                    name = "Region", labels = c("Southeast_Polynesia", "Tropical_Northwestern_Atlantic", "Tropical_Southwestern_Pacific", "Western_Coral_Triangle", "Western_Indian_Ocean")) +
  scale_shape_manual(values = c(25:21),
                     name = "Region", labels = c("Southeast_Polynesia", "Tropical_Northwestern_Atlantic", "Tropical_Southwestern_Pacific", "Western_Coral_Triangle", "Western_Indian_Ocean")) +
  labs(x = paste0("CAP1 (", CAP1, "%)"), y = paste0("CAP2 (", CAP2, "%)"),
       title = "") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.position = c(0,1),             # position in top left corner
        legend.justification = c(0, 1),        # correct legend justificaton
        legend.box.margin=margin(c(0.8,0.5,0.5,0.5)),  # add margin as to not overlap with axis box
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1)) 
a

# Panel b : rls

  #  load rdata

load("Rdata/dbRDA_data_rls.rdata")

  # plot


b <- ggplot(site_scores_region, aes(x= CAP1, y = CAP2)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  geom_vline(xintercept = 0, lty = 2, col = "grey") +
  geom_encircle(aes(group = province, linetype = province, fill= province), s_shape = 1, expand = 0,
                alpha = 0.4, show.legend = FALSE) + # hull area 
  geom_point(data = family_scores, aes(x= CAP1,y = CAP2), col = "grey", alpha = 0.5, cex = 0.5) +
  geom_point(data= family_scores_diff75, aes(x= CAP1, y=CAP2), col = "black", alpha = 1, cex = 0.5) +
  geom_point(aes(pch = province, fill = province), cex = 4, col = "black") +
  scale_fill_manual(values = c("#a6611a", "#E5A729","#b2182b", "#80cdc1", "#015462"),
                    name = "Region", labels = c("Southeast_Polynesia", "Tropical_Northwestern_Atlantic", "Tropical_Southwestern_Pacific", "Western_Coral_Triangle", "Western_Indian_Ocean")) +
  scale_shape_manual(values = c(25:21),
                     name = "Region", labels = c("Southeast_Polynesia", "Tropical_Northwestern_Atlantic", "Tropical_Southwestern_Pacific", "Western_Coral_Triangle", "Western_Indian_Ocean")) +
  labs(x = paste0("CAP1 (", CAP1, "%)"), y = paste0("CAP2 (", CAP2, "%)"),
       title = "") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.position = "none",             # position in top left corner
        legend.justification = c(1, 1),        # correct legend justificaton
        legend.box.margin=margin(c(1,1,1,1)),  # add margin as to not overlap with axis box
        legend.title = element_text(size=9),
        legend.text = element_text(size=9),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1)) 
b


# plot together

ggarrange(a, b, nrow=2, labels=c("a", "b"))
ggsave("outputs/00_Figures_for_paper/Figure4.png", width=6, height=10)
