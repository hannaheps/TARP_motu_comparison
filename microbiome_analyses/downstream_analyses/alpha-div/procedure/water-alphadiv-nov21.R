##Alpha Diversity metrics for both water and coral
##Combine the univariate data for downstream analyses

setwd("~/Documents/OSUDocs/Projects/French_Polynesia/Tetiaroa/Island_Survey/TARP_motu_comparison/microbiome_analyses/downstream_analyses/alpha-div/procedure/")


library(phyloseq)
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyverse)
library(gridExtra)
library(vegan)
library(RColorBrewer)
library(ggpubr)



##Water
#From Pre-processing file:
#Read in the november phyloseq object
physeq.water <- readRDS("../../pre-processing/output/physeq-water-nov21.RDS")
data.water <- as(sample_data(physeq.water), "data.frame")

#Check the rarefaction curve
mat.water <- otu_table(physeq.water)
class(mat.water) <- "matrix"
rarecurve((t(mat.water)), step = 1000)

#Based on the pre-processing script, we know we need to either cut-off or subsample to 17,959
#Rarefied - sub-sample to 17959
#Set seed for reproducibility pls
physeq.water.r <- rarefy_even_depth(physeq.water, sample.size = 17959, rngseed = 711) 
#2 samples removed = TW205, TW210
#3,404 ASVs were removed because they were no longer present in any sample after sub-sampling
data.water.r <- as(sample_data(physeq.water.r), "data.frame")

saveRDS(physeq.water.r, "../output/physeq-water-nov21-rarefied.RDS")


##Coral
physeq.coral <- readRDS("../../pre-processing/output/physeq-coral-nov21.RDS")
data.coral <- as(sample_data(physeq.coral), "data.frame")

#Check the rarefaction curve
mat.coral <- otu_table(physeq.coral)
class(mat.coral) <- "matrix"
rarecurve((t(mat.coral)), step = 50)

#Based on the pre-processing script, we know we need to either cut-off or subsample to 1,453
#Rarefied - sub-sample to 1,453
#Set seed for reproducibility pls
physeq.coral.r <- rarefy_even_depth(physeq.coral, sample.size = 1453, rngseed = 711) 
#3 samples removed = TC169, TC196, TC214
#8,350 ASVs were removed because they were no longer present in any sample after sub-sampling
data.coral.r <- as(sample_data(physeq.coral.r), "data.frame")

saveRDS(physeq.coral.r, "../output/physeq-coral-nov21-rarefied.RDS")

##Read these back in:
physeq.water.r <- read.RDS("../output/physeq-water-nov21-rarefied.RDS")
physeq.coral.r <- read.RDS("../output/physeq-coral-nov21-rarefied.RDS")


#Check that all the categorical and numerical values are correct in data frame
str(data.coral.r)
str(data.water.r)
#health.cat, distance.along.transect, site, collection year as character
data.coral.r$health.cat <- as.character(data.coral.r$health.cat)
data.coral.r$distance.along.transect <- as.character(data.coral.r$distance.along.transect)
data.coral.r$site <- as.character(data.coral.r$site)
data.coral.r$collection.year <- as.character(data.coral.r$collection.year)
data.water.r$distance.along.transect <- as.character(data.water.r$distance.along.transect)
data.water.r$site <- as.character(data.water.r$site)
data.water.r$collection.year <- as.character(data.water.r$collection.year)



##Estimating the alpha diversity metrics for each
#First need to input function estimate_richness_wPD (uses the estimate_richness() command in phyloseq
#but adds faith's phylogenetic distance metric to it as "FaithPD")
#Source code is located in procedure directory (aka your working directory)
source("estimate_richness_wPD.R")
erich.coral <- estimate_richness_wPD(physeq.coral.r, measures = c("Observed", "Shannon", "FaithPD"))
View(erich.coral)
erich.water <- estimate_richness_wPD(physeq.water.r, measure = c("Observed", "Shannon", "FaithPD"))
View(erich.water)
#add all the metadata
erich.coral <- cbind(erich.coral, data.coral.r)
erich.coral$sample.id <- rownames(erich.coral)
erich.water <- cbind(erich.water, data.water.r)
erich.water$sample.id <- rownames(erich.water)

erich.coral <- select(erich.coral, "sample.id", "Observed", "Shannon", "FaithPD", "collection.date", "motu", "island.side",
                  "station","site", "transect", "distance.along.transect", "lat", "long", "sample.type", "health.state",
                  "algae.N15", "algae.C13.not.acid", "algae.N.percent", "rat.status")
erich.water <- select(erich.water, "sample.id", "Observed", "Shannon", "FaithPD", "collection.date", "motu", "island.side",
                      "station","site", "transect", "distance.along.transect", "lat", "long", "sample_type", "mL_sterivex",
                      "algae.N15", "algae.C13.not.acid", "algae.N.percent", "rat.status")
                  
erich.water <- dplyr::rename(erich.water, "sample.type" = "sample_type")

erich.coral[setdiff(names(erich.water), names(erich.coral))] <- NA
erich.water[setdiff(names(erich.coral), names(erich.water))] <- NA

#Now ready to do an rbind because all column names are the same
erich <- rbind(erich.coral, erich.water)

#remove the mock community from the coral dataset
erich <- filter(erich, sample.type != "mock")

write.csv(erich, "../output/microbial_alphadiv.csv")

#Now we have a combined water and coral alpha diversity datasheet

ggplot(erich, aes(x = motu, y = Observed)) +
  geom_boxplot(aes(fill = motu)) +
  facet_wrap(~sample.type)


ggplot(erich, aes(x = motu, y =  Observed)) +
  geom_boxplot(aes(color = motu),alpha = 0.5) +
  geom_point(aes(colour = motu), size = 2, alpha = 3) +
  theme_bw() +
  facet_wrap(~sample.type)

ggsave("../output/plots/richness_motu.pdf", plot = last_plot())



