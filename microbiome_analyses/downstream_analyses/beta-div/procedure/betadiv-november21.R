##Betadiversity for Coral November 2021 only

setwd("~/Documents/OSUDocs/Projects/French_Polynesia/Tetiaroa/Island_Survey/TARP_motu_comparison/microbiome_analyses/downstream_analyses/alpha-div/procedure/")

library(phyloseq)
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyverse)
library(gridExtra)
library(vegan)
library(RColorBrewer)

#From Pre-processing file:
#Read in the rarefied november phyloseq object
physeq.nov.r <- readRDS("../output/downstream_analyses/november/nov_2021_rare.RDS")
data.nov.r <- as(sample_data(physeq.nov.r), "data.frame")

data.nov.r$dom.asv.name <- erich$dom.asv.name

##Let's just try an ordination to see
ord <- ordinate(physeq.nov.r, "NMDS", "bray", trymax = 500) #Run 20 stress 0.07221613
stressplot(ord)
scores <- as.data.frame(scores(ord))
scores <- cbind(scores, erich)

ord.wuf <- ordinate(physeq.nov.r, "NMDS", "wunifrac") #Run 20 stress 0.07221613
stressplot(ord.wuf)
scores.wuf <- as.data.frame(scores(ord.wuf))
scores.wuf <- cbind(scores.wuf, erich)

plot.wuf <- ggplot(scores.wuf, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = distance.along.transect, shape = motu), size = 4) +
  #geom_text(label = scores.outrm$sample.name) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = distance.along.transect), linetype = 2) +
  theme_bw()

plot.wuf+ facet_grid(cols = vars(motu))


##Jaccard
ord.jc <- ordinate(physeq.nov.r, "NMDS", "jaccard", trymax = 500) #Run 20 stress 0.07221613
stressplot(ord.jc)
scores.jc <- as.data.frame(scores(ord.jc))
scores.jc <- cbind(scores.jc, erich)

plot.jc <- ggplot(scores.jc, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = island.side, shape = dom.taxon), size = 4) +
  scale_shape_manual(values=11:18) +
  #geom_text(label = scores.outrm$sample.name) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = island.side), linetype = 2) +
  theme_bw() +
  facet_wrap(~motu)

disp.wuf <- betadisper(wuf, data.nov.r$motu, type = "centroid")
boxplot(disp.wuf)
boxplot(disp)

plot.br <- ggplot(scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = dom.taxon, shape = motu), size = 4) +
  #geom_text(label = scores.outrm$sample.name) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = dom.taxon), linetype = 2) +
  theme_bw()

ggsave("../output/downstream_analyses/november/figures/jaccard_nmds.pdf")

ggsave("../output/downstream_analyses/november/figures/coral_bray_nmds.pdf", plot = last_plot())

bc <- phyloseq::distance(physeq.nov.r, method = "bray")
disp <- betadisper(bc, data.nov.r$motu, type = "centroid")
#disp.motu <- betadisper(bc, data.r.coral$motu, type = "centroid")
boxplot(disp)


set.seed(47)
adonis2(bc ~ motu, data = data.nov.r)

#Yes they are significantly different by motu
#adonis2(formula = bc ~ motu, data = data.nov.r)
#Df SumOfSqs      R2      F Pr(>F)   
#motu      2   1.9319 0.09266 3.2681  0.002 **
#  Residual 64  18.9167 0.90734                 
#Total    66  20.8486 1.00000 
ggplot(scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = motu), size = 4) +
  #geom_text(label = scores.outrm$sample.name) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = motu), linetype = 2) +
  theme_bw()



##Since endo is driving the pattern, maybe we can map endo % onto a PCoA?
bray_curtis_pcoa <- ecodist::pco(bc)
bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], 
                                  pcoa2 = bray_curtis_pcoa$vectors[,2])


bray_curtis_plot <- ggplot(data = bray_curtis_pcoa_df, aes(x=pcoa1, y=pcoa2)) +
  geom_point() +
  labs(x = "PCo1",
       y = "PCo2", 
       title = "Bray-Curtis PCoA") +
  theme(title = element_text(size = 10)) # makes titles smaller


bray_curtis_pcoa_df <- cbind(bray_curtis_pcoa_df,
                             endo.relabund = erich$endo.relabund)
bray_curtis_pcoa_df <- cbind(bray_curtis_pcoa_df, motu = erich$motu)
bray_curtis_pcoa_df <- cbind(bray_curtis_pcoa_df, erich)

# Create a plot
bray_curtis_endo_plot <- ggplot(data = bray_curtis_pcoa_df, aes(x=pcoa1, y=pcoa2, color = endo.relabund)) + 
  geom_point(aes(shape = dom.taxon),size = 3) +
  scale_shape_manual(values=11:18) +
  scale_color_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = mean(erich$endo.relabund)) +
  labs(x = "PC1",
       y = "PC2",
       title = "PCoA of Microbial Communities vs  Endozoicomonas Relative Abundance") +
  theme_bw()

bray_curtis_endo_plot
ggsave("../output/downstream_analyses/november/figures/PCoA_Endo_domtaxa.pdf", plot = last_plot())

##Try another one with the algae stuff
bray_curtis_endo_plot <- ggplot(data = bray_curtis_pcoa_df, aes(x=pcoa1, y=pcoa2, color = algae.N.percent)) + 
  geom_point(size = 3) +
  #scale_shape_manual(values=11:18) +
  #scale_color_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = mean(erich$algae.N15)) +
  labs(x = "PC1",
       y = "PC2",
       title = "PCoA of Microbial Communities vs  Algae N15") +
  theme_bw()

bray_curtis_endo_plot
ggsave("../output/downstream_analyses/november/figures/PCoA_Endo_domtaxa.pdf", plot = last_plot())

##COuld we try a PCoA with wunifrac?
##Since endo is driving the pattern, maybe we can map endo % onto a PCoA?
wuf <- phyloseq::distance(physeq.nov.r, method = "wunifrac")
wunifrac_pcoa <- ecodist::pco(wuf)
wunifrac_pcoa_df <- data.frame(pcoa1 = wunifrac_pcoa$vectors[,1], 
                                  pcoa2 = wunifrac_pcoa$vectors[,2])


wunifrac_plot <- ggplot(data = wunifrac_pcoa_df, aes(x=pcoa1, y=pcoa2)) +
  geom_point() +
  labs(x = "PCo1",
       y = "PCo2", 
       title = "Weighted Unifrac PCoA") +
  theme(title = element_text(size = 10)) # makes titles smaller

wunifrac_pcoa_df <- cbind(wunifrac_pcoa_df, erich)

# Create a plot
wunifrac_endo_plot <- ggplot(data = wunifrac_pcoa_df, aes(x=pcoa1, y=pcoa2, color = endo.relabund)) + 
  geom_point(aes(shape = motu),size = 3) +
  #scale_shape_manual(values=11:18) +
  scale_color_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = mean(erich$endo.relabund)) +
  labs(x = "PC1",
       y = "PC2",
       title = "WUnifrac PCoA of Microbial Communities vs Endozoicomonas Relative Abundance") +
  theme_bw()

wunifrac_endo_plot
ggsave("../output/downstream_analyses/november/figures/wunifrac_PCoA_Endo_domtaxa.pdf", plot = last_plot())


##How about the barplots
percent.trial <- physeq.nov.r %>% 
  tax_glom(taxrank = "Phylum", NArm=TRUE) %>% 
  transform_sample_counts(function(x) {x/sum(x)} ) 
perc.melt <- psmelt(percent.trial)
sum <- ddply(perc.melt, c("Phylum", "motu", "island.side",),summarise,
             N = length(Abundance), 
             mean = mean(Abundance),
             sd = sd(Abundance), 
             se = sd/sqrt(N)
)

levels(sum$Phylum)
sum.as.factor <- sum
sum.as.factor$Phylum <- as.factor(sum.as.factor$Phylum)
levels(sum.as.factor$Phylum)

nb.cols <- 47
mycolors <- colorRampPalette(brewer.pal(11, "RdYlBu"))(nb.cols)
ggplot(sum, aes(x = island.side , y = mean, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=mycolors) +
  ylab("Relative Abundance") +
  xlab("Island Side") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~motu)

ggsave("../output/downstream_analyses/november/figures/phylum_barplot_motu.pdf", plot = last_plot())
