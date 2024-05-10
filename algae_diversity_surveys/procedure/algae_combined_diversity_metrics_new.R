#Fixing the issue with the relative abundance measurements from Maya

algae.all <- read.csv("../input/div_abund_all_algae_nov21.csv")
head(algae.all)
#algae.all$sample.id <- c(1:221)

#The abundance metric is the cm2 cover out of the total belt transect 
#This means we need to do a relative abundance count in some capacity...
algae.all$total.cm2 <- algae.all$transect.length
algae.all$total.cm2[algae.all$total.cm2=="25m"] <- 125000 #these are 2500cm x 50cm = 125k 
algae.all$total.cm2[algae.all$total.cm2 == "50m"] <- 250000 #these are 5000cm x 50cm = 250k
algae.all$total.cm2 <- as.integer(algae.all$total.cm2)

#Let's get a measure of relative percent cover that is abundance.cm2/total.cm2 

algae.all$rel.percent.cover <- (algae.all$abundance.cm2/algae.all$total.cm2)*100


library(tidyverse)
library(plyr)
algae.sum <- ddply(algae.all, c("site.name", "distance.along.transect", "transect.id.Maya", "algal.species"), summarise,
                   cover= sum(rel.percent.cover))

algae.matrix <- algae.sum %>% spread(algal.species,cover)

#Set all NAs to O because they were observed 0 times
algae.matrix[is.na(algae.matrix)] <- 0

View(algae.matrix)
#remove spaces in the variable names in dataframe
algae.matrix <- algae.matrix %>% rename_all(make.names)
View(algae.matrix)
#Add a non-algae cover category 
algae.matrix$non.algae.cover <- (100-(rowSums(algae.matrix[,4:34])))

#Use Vegan to run alpha div metrics on the data matrix
library(vegan)


shannon <- diversity(algae.matrix[,4:34], index = "shannon")
algae.matrix$shannon <- shannon
richness <- specnumber(algae.matrix[,4:34])
algae.matrix$richness <- richness
evenness <- shannon/log(richness)
algae.matrix$evenness <- evenness


View(algae.matrix)
algae.matrix <- algae.matrix %>% mutate(site.name = recode(site.name, "A1" = "Aie_Protected", "A2" = 'Aie_Exposed',
                                                           "Re1" = "Reiono_Protected", "Re2" = "Reiono_Exposed", 
                                                           "Rm1" = "Rimatuu_Protected", "Rm2" = "Rimatuu_Exposed"))
algae.matrix$distance.along.transect <- as.factor(algae.matrix$distance.along.transect)
write.csv(algae.matrix, "../output/algae_div_summary_all_new.csv", row.names = F)

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

write.csv(algae.by.site, "../output/algae_div_summary_bysite_new.csv", row.names = F)

#Stacked barchart
#1. Make the matrix into a long dataframe
colnames(algae.matrix[,4:35])
library(dplyr)
algae.long <- algae.matrix %>% gather(key=species,value=cover, c(colnames(algae.matrix[,4:35])))

#2. Build a stacked bar by species relative abundance & save as pdf
pdf(file = "../output/algae_species_stackedbar.pdf")
ggplot(algae.long, aes(fill=species, y=cover, x=site.name)) + 
  geom_bar(position="fill", stat = "identity") +
  ylab("Relative Percent Cover") +
  xlab("Site Name") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#3. differences among sites (PERMANOVAs)
vegdist.algae <- vegdist(algae.matrix[,4:34], method = "bray")
bc.algae <- as.matrix(vegdist.algae) 
perm.algae <- adonis2(bc.algae ~ site.name, data = algae.matrix)
#adonis2(formula = bc.algae ~ site.name, data = algae.matrix)
#Df SumOfSqs      R2      F Pr(>F)    
#site.name  5   2.5699 0.31118 3.7948  0.001 ***
#  Residual  42   5.6885 0.68882                  
#Total     47   8.2583 1.00000   
disp.algae <- betadisper(vegdist.algae, algae.matrix$site.name, type = "centroid")


pdf(file = "../output/algae_dispersion.pdf")
boxplot(disp.algae)
dev.off()
disp.algae.test <- permutest(disp.algae)
#Response: Distances
#Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)   
#Groups     5 0.31373 0.062746 4.5943    999  0.003 **
#  Residuals 42 0.57361 0.013657   
