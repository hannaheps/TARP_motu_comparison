##Betadiversity for Coral November 2021 only

setwd("~/Documents/OSUDocs/Projects/French_Polynesia/Tetiaroa/Island_Survey/TARP_motu_comparison/microbiome_analyses/downstream_analyses/beta-div/procedure/")

library(phyloseq)
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyverse)
library(gridExtra)
library(vegan)
library(RColorBrewer)

#From Pre-processing file:
#Read in the november phyloseq object
physeq <- readRDS("../../pre-processing/output/physeq-water-nov21.RDS")
data <- as(sample_data(physeq), "data.frame")

#Check the rarefaction curve
mat <- otu_table(physeq)
class(mat) <- "matrix"
rarecurve((t(mat)), step = 1000)

#Based on the pre-processing script, we know we need to either cut-off or subsample to 17,959
#Unrarefied - cut out every sample with below 17959
physeq.nr <- prune_samples(sample_sums(physeq)>= 17959, physeq)
data.nr <- as(sample_data(physeq.nr), "data.frame")

#Rarefied - sub-sample to 17959
#Set seed for reproducibility pls
physeq.r <- rarefy_even_depth(physeq, sample.size = 17959, rngseed = 711) 
#4 samples removed = TW205, TW210
#3,404 ASVs were removed because they were no longer present in any sample after sub-sampling
data.r <- as(sample_data(physeq.r), "data.frame")


##Let's just try an ordination to see
ord <- ordinate(physeq.r, "NMDS", "bray", trymax = 500) #Run 20 stress 0.061
stressplot(ord)
scores <- (scores(ord))
scores <- (scores$sites)
scores <- cbind(scores, data.r)
scores$distance.along.transect <- as.factor(scores$distance.along.transect)

plot.bray <- ggplot(scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = island.side), size = 3, shape = 18) +
  scale_colour_manual(values = c("#63A022", "#2866AB")) +
  #geom_text(label = scores.outrm$sample.name) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = island.side), linetype = 6) +
  theme_bw()


plot.bray + facet_grid(cols = vars(motu))

ggsave("../output/plots/water_bray_ordination_islandside.pdf", plot = last_plot())








###Coral

physeq.r <- readRDS("../../alpha-div/output/physeq-coral-nov21-rarefied.RDS")
data.r <- as(sample_data(physeq.r), "data.frame")
physeq.r <- subset_samples(physeq.r, sample.type != "mock")

##Let's just try an ordination to see
ord <- ordinate(physeq.r, "NMDS", "bray", trymax = 500) #Run 20 stress 0.061
stressplot(ord)
scores <- (scores(ord))
scores <- (scores$sites)
scores <- cbind(scores, data.r)
scores$distance.along.transect <- as.factor(scores$distance.along.transect)

plot.bray <- ggplot(scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = island.side), size = 2) +
  scale_colour_manual(values = c("#63A022", "#2866AB")) +
  #geom_text(label = scores.outrm$sample.name) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = island.side), linetype = 6) +
  theme_bw()


plot.bray + facet_grid(cols = vars(motu))

ggsave("../output/plots/coral_bray_ordination_islandside.pdf", plot = last_plot())


