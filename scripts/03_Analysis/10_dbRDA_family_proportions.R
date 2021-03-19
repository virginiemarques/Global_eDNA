library(dplyr)
library(forcats)
library(stringr)
# for (a) 
library(rsq)
library(margins)
# for (b)
library(betapart)
library(reshape)
library(tidyselect)
library(vegan)
library(ggplot2)
library(patchwork)
library(ggalt)
library(ggrepel)
library(grid)


load("Rdata/family_proportion_per_site.rdata")
load("Rdata/family_proportion_site_RLS.rdata")

met_eDNA <- df_all_site[,c("site", "province", "n_motus_total")] %>%
  distinct()

met_rls <- count_families_site_RLS[,c("site35", "province", "n_sp_total")]%>%
  distinct()

# prepare data : site x family with prop 
RLS_fam <- count_families_site_RLS[,c("site35", "family", "prop")]
RLS_fam <- as.data.frame(pivot_wider(RLS_fam, id_cols=site35, names_from=family, values_from = prop))
RLS_fam[is.na(RLS_fam)] <- 0
rownames(RLS_fam) <- RLS_fam[,1]
RLS_fam <- RLS_fam[,-1]

edna_fam <- df_all_site[,c("site", "family", "prop")]
edna_fam <- as.data.frame(pivot_wider(edna_fam, id_cols=site, names_from=family, values_from = prop))
edna_fam[is.na(edna_fam)] <- 0
rownames(edna_fam) <- edna_fam[,1]
edna_fam <- edna_fam[,-1]

#----------------------------------------------------------------------------------------------------------------
# dbRDA on eDNA
#----------------------------------------------------------------------------------------------------------------

# Total dbRDA
edna.tot.dbrda  <- capscale(edna_fam ~ province + n_motus_total, data = met_eDNA, distance = "bray")
RsquareAdj(edna.tot.dbrda)
anova(edna.tot.dbrda, by = "axis",   permutations = 9999)
anova(edna.tot.dbrda, by = "margin", permutations = 9999)

 # partial dbRDA for Province
edna.part.dbrda1  <- capscale(edna_fam ~ province + Condition(n_motus_total), data = met_eDNA, distance = "bray")
RsquareAdj(edna.part.dbrda1)
anova(edna.part.dbrda1,by="axis", permutations = 9999)


# partial dbRDA for site richness
edna.part.dbrda2  <- capscale(edna_fam ~ n_motus_total + Condition(province), data = met_eDNA, distance = "bray")
RsquareAdj(edna.part.dbrda2)
anova(edna.part.dbrda2, by="axis", permutations = 9999)


# get scores
site_scores <- scores(edna.part.dbrda1)$sites
family_scores <- scores(edna.part.dbrda1)$species %>% data.frame()

# get most differentiated species along first axis
quant75 <- quantile(abs(family_scores$CAP1), probs = c(0.75))
family_scores_diff75 <- family_scores[which(abs(family_scores$CAP1) > quant75["75%"]),]

# extract the percentage variability explained by axes
sumdbrda <- summary(edna.part.dbrda1)
CAP1 <- round(sumdbrda$cont$importance["Proportion Explained", "CAP1"]*100, 1)
CAP2 <- round(sumdbrda$cont$importance["Proportion Explained", "CAP2"]*100, 1)

# add metadata
identical(as.character(met_eDNA$site), rownames(site_scores)) # verify that data in same order
site_scores_region <- cbind(site_scores,select(met_eDNA, province))

# plot
grda_sites <- ggplot(site_scores_region, aes(x= CAP1, y = CAP2)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  geom_vline(xintercept = 0, lty = 2, col = "grey") +
  geom_encircle(aes(group = province, linetype = province, fill= province), s_shape = 1, expand = 0,
                alpha = 0.4, show.legend = FALSE) + # hull area 
  geom_point(data = family_scores, aes(x= CAP1,y = CAP2), col = "grey", alpha = 0.5, cex = 0.5) +
  geom_point(data= family_scores_diff75, aes(x= CAP1, y=CAP2), col = "black", alpha = 1, cex = 0.5) +
  geom_point(aes(pch = province, fill = province), cex = 4, col = "black") +
  scale_fill_manual(values = c("#a6611a", "#E5A729","#b2182b", "#80cdc1", "#015462"),
                    name = "Region", labels = c("Southeast_Pacific", "Tropical_Northwestern_Atlantic", "Tropical_Southwestern_Pacific", "Western_Coral_Triangle", "Western_Indian_Ocean")) +
  scale_shape_manual(values = c(25:21),
                     name = "Region", labels = c("Southeast_Pacific", "Tropical_Northwestern_Atlantic", "Tropical_Southwestern_Pacific", "Western_Coral_Triangle", "Western_Indian_Ocean")) +
  labs(x = paste0("CAP1 (", CAP1, "%)"), y = paste0("CAP2 (", CAP2, "%)"),
       title = "") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.position = c(0, 1),             # position in top left corner
        legend.justification = c(0, 1),        # correct legend justificaton
        legend.box.margin=margin(c(2,2,2,2)),  # add margin as to not overlap with axis box
        legend.title = element_text(size=11),
        legend.text = element_text(size=11),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1)) 
grda_sites

grda_family <- ggplot() + 
  geom_segment(data= family_scores, aes(x=0, xend=CAP1,y = 0, yend=CAP2), col = "grey",
               arrow=arrow(length=unit(0.01,"npc"))) + # all species
  geom_segment(data= family_scores_diff75, aes(x=0, xend=CAP1,y = 0, yend=CAP2), col = "black",
               arrow=arrow(length=unit(0.01,"npc"))) + # most differentiated species
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  geom_vline(xintercept = 0, lty = 2, col = "grey") +
  geom_label_repel(data= family_scores_diff75, 
                   aes(x= CAP1, y=CAP2, #hjust=0.5*(1-sign(CAP1)),vjust=0.5*(1-sign(CAP2)),
                   fontface=3), #size = 3,
                   label = rownames(family_scores_diff75),
                   show.legend = F) +
  labs(x = paste0("CAP1 (", CAP1, "%)"), y = paste0("CAP2 (", CAP2, "%)")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        legend.position = c(0, 1),              # position in top left corner
        legend.justification = c(0, 1),         # correct legend justificaton
        legend.box.margin=margin(c(2,2,2,2)),   # add margin as to not overlap with axis box
        legend.background = element_rect(fill =  alpha("white", 0.0)),
        legend.title = element_text(size=11),
        legend.text = element_text(size=11)) +
  # re-add legend and change text legend key by making invisible points and overriding its key shape
  geom_point(data= family_scores_diff75, 
             aes(x=CAP1, y=CAP2),
             size = 0, stroke = 0) + 
  guides(colour = guide_legend(override.aes = list(size = 5, 
                                                   shape = c(utf8ToInt("C"), utf8ToInt("B"), utf8ToInt("D"), utf8ToInt("P")))))
grda_family



#----------------------------------------------------------------------------------------------------------------
# dbRDA on RLS
#----------------------------------------------------------------------------------------------------------------

rls.tot.dbrda  <- capscale(RLS_fam ~ province + n_sp_total, data = met_rls, distance = "bray")
RsquareAdj(rls.tot.dbrda)
anova(rls.tot.dbrda, by = "axis",   permutations = 9999)
anova(rls.tot.dbrda, by = "margin", permutations = 9999)

# partial dbRDA for Province
rls.part.dbrda1  <- capscale(RLS_fam ~ province + Condition(n_sp_total), data = met_rls, distance = "bray")
RsquareAdj(rls.part.dbrda1)
anova(rls.part.dbrda1,by="axis", permutations = 99)


# partial dbRDA for site richness
rls.part.dbrda2  <- capscale(RLS_fam ~ n_sp_total + Condition(province), data = met_rls, distance = "bray")
RsquareAdj(rls.part.dbrda2)
anova(rls.part.dbrda2, by="axis", permutations = 99)


# get scores
site_scores <- scores(rls.part.dbrda1)$sites
family_scores <- scores(rls.part.dbrda1)$species %>% data.frame()

# get most differentiated species along first axis
quant75 <- quantile(abs(family_scores$CAP1), probs = c(0.75))
family_scores_diff75 <- family_scores[which(abs(family_scores$CAP1) > quant75["75%"]),]

# extract the percentage variability explained by axes
sumdbrda <- summary(rls.part.dbrda1)
CAP1 <- round(sumdbrda$cont$importance["Proportion Explained", "CAP1"]*100, 1)
CAP2 <- round(sumdbrda$cont$importance["Proportion Explained", "CAP2"]*100, 1)

# add metadata
identical(as.character(met_rls$site), rownames(site_scores)) # verify that data in same order
site_scores_region <- cbind(site_scores,select(met_rls, province))

# plot
grda_sites <- ggplot(site_scores_region, aes(x= CAP1, y = CAP2), col=province) +
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  geom_vline(xintercept = 0, lty = 2, col = "grey") +
  geom_encircle(aes(group = province, linetype = province, fill= province), s_shape = 1, expand = 0,
                alpha = 0.4, show.legend = FALSE) + # hull area 
  geom_point(data = family_scores, aes(x= CAP1,y = CAP2), col = "grey", alpha = 0.5, cex = 0.5) +
  geom_point(data= family_scores_diff75, aes(x= CAP1, y=CAP2), col = "black", alpha = 1, cex = 0.5) +
  geom_point(aes(pch = province, fill = province), cex = 4, col = "black") +
  #scale_fill_manual(palette=hue_pal(c(1,24)))+
                    #name = "Region", labels = province) +
  scale_shape_manual(values = c(25:1)) +
  labs(x = paste0("CAP1 (", CAP1, "%)"), y = paste0("CAP2 (", CAP2, "%)"),
       title = "") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.position = "none",             # position in top left corner
        legend.justification = c(0, 1),        # correct legend justificaton
        legend.box.margin=margin(c(2,2,2,2)),  # add margin as to not overlap with axis box
        legend.title = element_text(size=11),
        legend.text = element_text(size=11),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1)) 
grda_sites

grda_family <- ggplot() + 
  geom_segment(data= family_scores, aes(x=0, xend=CAP1,y = 0, yend=CAP2), col = "grey",
               arrow=arrow(length=unit(0.01,"npc"))) + # all species
  geom_segment(data= family_scores_diff75, aes(x=0, xend=CAP1,y = 0, yend=CAP2), col = "black",
               arrow=arrow(length=unit(0.01,"npc"))) + # most differentiated species
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  geom_vline(xintercept = 0, lty = 2, col = "grey") +
  geom_label_repel(data= family_scores_diff75, 
                   aes(x= CAP1, y=CAP2, #hjust=0.5*(1-sign(CAP1)),vjust=0.5*(1-sign(CAP2)),
                       fontface=3), #size = 3,
                   label = rownames(family_scores_diff75),
                   show.legend = F) +
  labs(x = paste0("CAP1 (", CAP1, "%)"), y = paste0("CAP2 (", CAP2, "%)")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        legend.position = c(0, 1),              # position in top left corner
        legend.justification = c(0, 1),         # correct legend justificaton
        legend.box.margin=margin(c(2,2,2,2)),   # add margin as to not overlap with axis box
        legend.background = element_rect(fill =  alpha("white", 0.0)),
        legend.title = element_text(size=11),
        legend.text = element_text(size=11)) +
  # re-add legend and change text legend key by making invisible points and overriding its key shape
  geom_point(data= family_scores_diff75, 
             aes(x=CAP1, y=CAP2),
             size = 0, stroke = 0) + 
  guides(colour = guide_legend(override.aes = list(size = 5, 
                                                   shape = c(utf8ToInt("C"), utf8ToInt("B"), utf8ToInt("D"), utf8ToInt("P")))))
grda_family


# filter same region as eDNA

site_scores_region <- site_scores_region %>%
  filter(province %in% c("Southeast_Polynesia", "Tropical_Northwestern_Atlantic", "Tropical_Southwestern_Pacific", "Western_Coral_Triangle", "Western_Indian_Ocean"))


grda_sites <- ggplot(site_scores_region, aes(x= CAP1, y = CAP2)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  geom_vline(xintercept = 0, lty = 2, col = "grey") +
  geom_encircle(aes(group = province, linetype = province, fill= province), s_shape = 1, expand = 0,
                alpha = 0.4, show.legend = FALSE) + # hull area 
  geom_point(data = family_scores, aes(x= CAP1,y = CAP2), col = "grey", alpha = 0.5, cex = 0.5) +
  geom_point(data= family_scores_diff75, aes(x= CAP1, y=CAP2), col = "black", alpha = 1, cex = 0.5) +
  geom_point(aes(pch = province, fill = province), cex = 4, col = "black") +
  scale_fill_manual(values = c("#a6611a", "#E5A729","#b2182b", "#80cdc1", "#015462"),
                    name = "Region", labels = c("Southeast_Pacific", "Tropical_Northwestern_Atlantic", "Tropical_Southwestern_Pacific", "Western_Coral_Triangle", "Western_Indian_Ocean")) +
  scale_shape_manual(values = c(25:21),
                     name = "Region", labels = c("Southeast_Pacific", "Tropical_Northwestern_Atlantic", "Tropical_Southwestern_Pacific", "Western_Coral_Triangle", "Western_Indian_Ocean")) +
  labs(x = paste0("CAP1 (", CAP1, "%)"), y = paste0("CAP2 (", CAP2, "%)"),
       title = "") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        legend.position = c(0, 1),             # position in top left corner
        legend.justification = c(0, 1),        # correct legend justificaton
        legend.box.margin=margin(c(2,2,2,2)),  # add margin as to not overlap with axis box
        legend.title = element_text(size=11),
        legend.text = element_text(size=11),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1)) 
grda_sites




