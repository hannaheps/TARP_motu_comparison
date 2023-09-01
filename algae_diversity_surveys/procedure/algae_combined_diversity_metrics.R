algae.all <- read.csv("../input/div_abund_all_algae_nov21.csv")
head(algae.all)
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
algae.matrix <- algae.matrix %>% mutate(site.name = recode(site.name, "A1" = "Aie_Protected", "A2" = 'Aie_Exposed',
                                                   "Re1" = "Reiono_Protected", "Re2" = "Reiono_Exposed", 
                                                   "Rm1" = "Rimatuu_Protected", "Rm2" = "Rimatuu_Exposed"))

head(algae.matrix)
write.csv(algae.matrix, "../output/algae_div_summary_all.csv", row.names = F)

algae.by.site <- ddply(algae.matrix, c("site.name"), summarise,
                       mean.algae.richness = mean(richness),
                       mean.algae.evenness = mean(evenness), 
                       mean.algae.shannon = mean(shannon), 
                       mean.algae.turbinaria = mean(Turbinaria.ornata),
                       mean.algae.lobophora = mean(Lobophora.spp),
                       mean.algae.halimeda = mean(Halimeda.distorta),
                       mean.algae.chl.fastigita = mean(Chlorodesmis.fastigiata),
                       mean.algae.caul.serrulata = mean(Caulerpa.serrulata))
View(algae.by.site)

write.csv(algae.by.site, "../output/algae_div_summary_bysite.csv", row.names = F)

