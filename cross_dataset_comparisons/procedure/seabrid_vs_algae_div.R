##Combine seabirds and algae diversity
library(tidyverse)
library(plyr)
library(lme4)
library(car)
library(jtools)
library(emmeans)
library(multcompView)
library(multcomp)


#Bring in seabird data
seabirds.all <- read.csv("../../seabird_data/outputs/seabird_dens_bio_focal_data.csv")
View(seabirds.all)

#clean-up the dataframe by removing the X column (column 1) accidentally added in the csv file.
#At the same time, change the name to maintain consistency in downstream code
seabirds <- seabirds.all[,-1]
View(seabirds)
#add a site.name (Motu_Exposure) category separated by an underscore.
seabirds$site.name <- paste(seabirds$Motu, seabirds$Exposure, sep = "_")

#Bring in the algae diversity data
algae.div <- read.csv("../../algae_diversity_surveys/output/algae_div_summary_all.csv")
View(algae.div)

#Combine the algae and seabird data 
algae.div.seabirds <- merge(algae.div, seabirds, by = "site.name", all = TRUE, no.dups = TRUE)
View(algae.div.seabirds)
write.csv(algae.div.seabirds, "../output/seabird-algaediv/seabird_algaediv_noiti_totaldata.csv", row.names= F)

mod.algaediv <- lmer(richness ~ breeding_biomass_kgha_side + (1|site.name), 
            data = algae.div.seabirds)
Anova(mod.algaediv)


