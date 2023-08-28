##Combine the Seabird data vs. N15

##Bring in seabird data
seabirds.all <- read.csv("../../seabird_data/outputs/seabird_dens_bio_focal_data.csv")
View(seabirds.all)

#The following code under the #s will allow for the inclusion of iti birds, otherwise skip to line 25:
##########
#remove the density and biomass by side so that we can replace with the values that include iti
seabirds.all <- seabirds.all[, -c(1,6,10,14,18)] 
seabirds.all$site.name <- paste(seabirds.all$Motu, seabirds.all$Exposure, sep = "_") #create a site.name value for merging



#bring in the iti density data and combine with all other data
seabirds.iti <- read.csv("../../seabird_data/outputs/seabird_dens_bio_by_side_iti_rimatuu_exposed_combined.csv")
seabirds.iti$site.name <- paste(seabirds.iti$Motu, seabirds.iti$Exposure, sep = "_")

seabirds <- merge(seabirds.iti, seabirds.all, by = "site.name")

#reorder the column names
col_order <- c("site.name", "X", "Motu.x", "Exposure.x", "Motu.y", "Exposure.y", "breeding_biomass_kgha_100m", "breeding_biomass_kgha_200m",
               "breeding_biomass_kgha_side", "breeding_biomass_kgha_motu", "breeding_density_ha_100m", "breeding_density_ha_200m",
               "breeding_density_ha_side", "breeding_density_ha_motu", "adult_biomass_kgha_100m", "adult_biomass_kgha_200m",
               "adult_biomass_kgha_side", "adult_biomass_kgha_motu",  "adult_density_ha_100m", "adult_density_ha_200m",
               "adult_density_ha_side",  "adult_density_ha_motu")  
seabirds <- seabirds[, col_order]
colnames(seabirds)

seabirds <- seabirds[, -c(2:6)] #remove the non numerical columns aside from site.name needed for merging

write.csv(seabirds, "../output/seabirds.iti.combined.csv")

###########
##If NOT including iti run the following code to maintain consistent terminology for downstream analyses
seabirds <- seabirds.all
seabirds$site.name <- paste(seabirds$Motu, seabirds$Exposure, sep = "_")


##bring in algae data
algae <- read.csv("../../algae_isotopes/data/Tetiaroa_Turbinaria_Transects_November_2021_compiledMarch2023.csv")
View(algae)
algae$Distance_to_shore <- as.factor(algae$Distance_to_shore)

#manipulate algae data frame to get average N15 per motu, site and distance from shore
library(plyr)
sum.algae <- ddply(algae, c("Motu", "Site", "Distance_to_shore"),summarise,
                      mean.N15 = mean(N15)
)

sum.algae
#spread data to get values for each distance from shore
library(tidyverse)
sum.algae <- sum.algae %>% spread(Distance_to_shore, mean.N15)

sum.algae <- sum.algae %>% rename(c("N.15_at_10m" = "10", 
                                    "N.15_at_20m" = "20", "N.15_at_30m" = "30",
                                    "N.15_at_40m" = "40"))


#Rename Sites to match seabird data
sum.algae$Site <- as.factor(sum.algae$Site)
sum.algae <-sum.algae %>% mutate(Site = recode(Site, "1" = "Protected", "2" = 'Exposed'))

sum.algae$Exposure <- sum.algae$Site
sum.algae <- sum.algae[,-2]

#Create a site.name that combines motu and exposure in order to merge with seabrid data
sum.algae$site.name <-  paste(sum.algae$Motu, sum.algae$Exposure, sep = "_")


##Combine seabird and algae data
algae.seabirds <- merge(sum.algae, seabirds, by = "site.name", all = TRUE, no.dups = TRUE)
View(algae.seabirds)
algae.seabirds <- algae.seabirds[,-c(2, 7:10)] #This is just to remove duplicate columns
#algae.seabirds <- algae.seabirds[,-c(2,7)] #Use this one for Iti data, it has slightly different column #s 
View(algae.seabirds)

#save combined data
write.csv(algae.seabirds, "../output/n15_seabirds_combined_with_iti.csv") #change to "_no_iti" if removing iti 

##Create a matrix
data.matrix <- as.data.frame(algae.seabirds[,2:21])
View(data.matrix)


#Run a correlation test using the library corrplot
library(corrplot)

cor.mtest(data.matrix)
correlation.matrix <- cor(data.matrix, use = "pairwise.complete.obs")
write.csv(correlation.matrix, "../output/n15_seabirds_corrmatrix_with_iti.csv") #change to "_no_iti" if removing iti


pdf(file = "../output/seabird_v_N15_with_iti.pdf") #change to "_no_iti" if removing iti

corrplot(cor(data.matrix, use = "pairwise.complete.obs"), type = "upper", 
         addCoef.col = NULL, addCoefasPercent = FALSE, tl.col = "black", title = "seabirds vs. N15")

dev.off()
