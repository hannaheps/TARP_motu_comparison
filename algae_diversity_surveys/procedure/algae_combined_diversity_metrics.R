algae.all <- read.csv("../../algae_diversity_surveys/div_abund_all_algae_nov21.csv")
#algae.all$sample.id <- c(1:221)

library(tidyverse)
library(plyr)
algae.sum <- ddply(algae.all, c("site.name", "transect.id.Maya", "algal.species"), summarise,
                      abundance = sum(abundance.cm2))

algae.matrix <- algae.sum %>% spread(algal.species,abundance)

#Set all NAs to O because they were observed 0 times
algae.matrix[is.na(algae.matrix)] <- 0

View(algae.matrix)
#remove spaces in the variable names in dataframe
algae.matrix <- algae.matrix %>% rename_all(make.names)
View(algae.matrix)
#Use Vegan to run alpha div metrics on the data matrix
library(vegan)

shannon <- diversity(algae.matrix[,3:33], index = "shannon")
algae.matrix$shannon <- shannon
richness <- specnumber(algae.matrix[,3:33])
algae.matrix$richness <- richness
evenness <- shannon/log(richness)
algae.matrix$evenness <- evenness

View(algae.matrix)

