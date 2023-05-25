##Combine the Seabird data vs. N15

##Bring in seabird data
seabirds.all <- read.csv("../../seabird_data/outputs/seabird_dens_bio_focal_data.csv")
View(seabirds.all)
#remove the density and biomass by side so that we can replace with the values that include iti
seabirds.all <- seabirds.all[, -c(1,6,10,14,18)]
seabirds.all$site.name <- paste(seabirds.all$Motu, seabirds.all$Exposure, sep = "_")



#bring in the iti density data and combine with all other data
seabirds.iti <- read.csv("../../seabird_data/outputs/seabird_dens_bio_by_side_iti_rimatuu_exposed_combined.csv")
seabirds.iti$site.name <- paste(seabirds.iti$Motu, seabirds.iti$Exposure, sep = "_")

seabirds <- merge(seabirds.iti, seabirds.all, by = "site.name")

seabirds <- seabirds[, -c(2:4,9:10)]

write.csv(seabirds, "../output/seabirds.iti.combined.csv")

##bring in algae data
algae <- read.csv("../../algae_isotopes/data/Tetiaroa_Turbinaria_Transects_November_2021_compiledMarch2023.csv")
View(algae)
algae$Distance_to_shore <- as.factor(algae$Distance_to_shore)

#manipulate algae data frame to get average N15 per mote, site and distance from shore
library(plyr)


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
algae.seabirds <- algae.seabirds[,-c(2, 7)]
View(algae.seabirds)

#save combined data
write.csv(algae.seabirds, "../output/n15_seabirds_combined.csv")

##Create a matrix
data.matrix <- as.data.frame(algae.seabirds[,2:21])


#Run a correlation test using the library corrplot
library(corrplot)

cor.mtest(data.matrix)
cor(data.matrix)

corrplot(cor(data.matrix), type = "upper", 
         addCoef.col = NULL, addCoefasPercent = FALSE, tl.col = "black", title = "seabirds vs. N15")

