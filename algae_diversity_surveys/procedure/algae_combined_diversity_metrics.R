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

#Stacked barchart
#1. Make the matrix into a long dataframe
colnames(algae.matrix[,3:33])
library(dplyr)
algae.long <- algae.matrix %>% gather(key=species,value=abundance, c(colnames(algae.matrix[,3:33])))

#2. Build a stacked bar by species relative abundance & save as pdf
pdf(file = "../output/algae_species_stackedbar.pdf")
ggplot(algae.long, aes(fill=species, y=abundance, x=site.name)) + 
  geom_bar(position="fill", stat="identity")
dev.off()

#3. differences among sites (PERMANOVAs)
vegdist.algae <- vegdist(algae.matrix[,3:33], method = "bray")
bc.algae <- as.matrix(vegdist.algae) 
perm.algae <- adonis2(bc.algae ~ site.name, data = algae.matrix)
#adonis2(formula = bc.algae ~ site.name, data = algae.matrix)
#           Df  SumOfSqs  R2     F    Pr(>F)    
#site.name  5   2.7515 0.32601 4.063  0.001 ***
#Residual  42   5.6885 0.67399                 
#Total     47   8.4399 1.00000 
disp.algae <- betadisper(vegdist.algae, algae.matrix$site.name, type = "centroid")
pdf(file = "../output/algae_dispersion.pdf")
boxplot(disp.algae)
dev.off()
disp.algae.test <- permutest(disp.algae)
#Response: Distances
#           Df  Sum Sq  Mean Sq   F     N.Perm Pr(>F)   
#Groups     5 0.31373 0.062746 4.5943    999  0.004 **
# Residuals 42 0.57361 0.013657 
