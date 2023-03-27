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


##We know from the barplots that the Aie communities are dominated
#by litoricola marina (heterotroph) and by Synecchococcus (cyanobacteria - bloomer)
##But the other two motus are mostly Prochlorococcus and SAR11 (oceanic bugs)

##Can we build a summarized stacked bar plot to show the top 10 most abundant genera??

percent.trial <- physeq.water.r %>% 
  tax_glom(taxrank = "Genus", NArm=TRUE) %>% 
  transform_sample_counts(function(x) {x/sum(x)} ) 
trial.top = prune_taxa(names(sort(taxa_sums(percent.trial), TRUE))[1:10], percent.trial)
perc.melt <- psmelt(trial.top)
sum <- ddply(perc.melt, c("Genus", "motu", "island.side"),summarise,
             N = length(Abundance), 
             mean = mean(Abundance),
             sd = sd(Abundance), 
             se = sd/sqrt(N)
)

sum$genus_order <- sum$Genus
sum$genus_order[sum$genus_order == "AEGEAN-169_marine_group" ] <- 'A-Rhodospirillales; AEGEAN-169'
sum$genus_order[sum$genus_order == "Clade_Ia"] <- 'B-SAR11; Clade Ia'
sum$genus_order[sum$genus_order == "Clade_Ib"] <- 'C-SAR11; Clade Ib'
sum$genus_order[sum$genus_order == "Clade_II"] <- 'D-SAR11; Clade II'
sum$genus_order[sum$genus_order == "Litoricola"] <- 'H-Litoricola'
sum$genus_order[sum$genus_order ==  "Prochlorococcus_MIT9313" ] <- 'G-Prochlorococcus'
sum$genus_order[sum$genus_order == "SAR116_clade"] <- 'E-SAR116 Clade'
sum$genus_order[sum$genus_order == "SAR86_clade"] <- "F-SAR86 Clade"
sum$genus_order[sum$genus_order == "Synechococcus_CC9902"] <- 'I-Synechococcus'
sum$genus_order[sum$genus_order == "uncultured"] <- 'J-Uncultured Alteromonadaceaea'


nb.cols <- 10
mycolors <- colorRampPalette(brewer.pal(10, "Set3"))(nb.cols)

ggplot(sum, aes(x = island.side , y = mean, fill = genus_order)) +
  geom_bar(stat = "identity") +
  #scale_fill_manual(values=mycolors) +
  scale_fill_manual(values = c('A-Rhodospirillales; AEGEAN-169' = "#C996CC", 'B-SAR11; Clade Ia' = "#7209B7", 'C-SAR11; Clade Ib' = "#170055", 'D-SAR11; Clade II' = "#14279B",
                               'E-SAR116 Clade' = "#185ADB", "F-SAR86 Clade" = "#77ACF1", 'G-Prochlorococcus' =  "#98BAE7",
                               'H-Litoricola' = "#FFB950", 
                               'I-Synechococcus' = "#FA5E1F", 
                                'J-Uncultured Alteromonadaceaea' = "#7A0103")) +
  #scale_fill_discrete() +
  ylab("Relative Abundance") +
  xlab("Island Side") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~motu)

ggsave("../output/plots/water_barplots_top10genera.pdf", plot = last_plot())
