###For TARP Data Integration Paper###
###Extraction of microbiome metrics (water and coral) for integration###

##What to include:
#1. extract endozoicomonas abundance in corals
#2. extract the bloomers/heteros abundances in water
#3. extract top 10 most abundant abundances in corals & water
#4. microbiome richness in corals & water
#5. microbiome eveness in corals & water
#6. microbiome phylogenetic diversity in corals & water
#7. microbiome Shannon diversity in corals & water
#8. microbiome PCoA or NMDS axis 1 in corals & water
#9. microbiome PCoA or NMDS axis 2 in corals & water
#10 microbiome dispersion in corals & water

##add metadata to this so that it can be combined with the other data!


##Set your working directory to integration/procedure folder
setwd("/Users/hannah/Documents/OSUDocs/Projects/French_Polynesia/Tetiaroa/Island_Survey/TARP_motu_comparison/microbiome_analyses/downstream_analyses/integration/procedure")

library(phyloseq)
library(dplyr)
library(tidyverse)

#pull in the CSV with the alpha diversity metrics
erich <- read.csv("../../alpha-div/output/microbial_alphadiv.csv", header = TRUE, strip.white = TRUE)
#erich now already includes all the alpha diversity that we want ASIDE from evenness

#Add an eveness column:
#evenness is just calculated by dividing Shannon div by the natural log of Richness
erich$Evenness <- erich$Shannon/(log(erich$Observed))


#pull in the two phyloseq objects that to grab the beta diversity metrics
physeq.coral.r <- readRDS("../../alpha-div/output/physeq-coral-nov21-rarefied.RDS")

physeq.water.r <- readRDS("../../alpha-div/output/physeq-water-nov21-rarefied.RDS")

#Water Top 10 Abundance (includes bloomers/heteros)
#######################
##Subset data to the top 10 most abundant taxa and add abundances to the erich dataframe
#first with water
percent.trial <- physeq.water.r %>% 
  tax_glom(taxrank = "Genus", NArm=TRUE) %>%
  transform_sample_counts(function(x) {x/sum(x)} ) 
trial.top = prune_taxa(names(sort(taxa_sums(percent.trial), TRUE))[1:10], percent.trial)
perc.melt <- psmelt(trial.top)
perc.melt$sample.id <- perc.melt$Sample

#summarize genera by sample id
sum.water.abund <- ddply(perc.melt, c("sample.id", "Genus"), summarise,
                         N = sum(Abundance)
)
#manipulate data frame so that there is an abundance per sample id for each genus
sum.water.abund <- sum.water.abund %>% spread(Genus, N)

#Rename the columns so that we know it's relative abundance
sum.water.abund <- sum.water.abund %>% dplyr::rename(RelAbund_Rhodospirillales_AEGEAN_169_marine_group = `AEGEAN-169_marine_group`)
sum.water.abund <- sum.water.abund %>% dplyr::rename(RelAbund_SAR11_CladeIa = Clade_Ia)
sum.water.abund <- sum.water.abund %>% dplyr::rename(RelAbund_SAR11_CladeIb = Clade_Ib)
sum.water.abund <- sum.water.abund %>% dplyr::rename(RelAbund_SAR11_CladeII = Clade_II)
sum.water.abund <- sum.water.abund %>% dplyr::rename(RelAbund_SAR116_Clade = SAR116_clade)
sum.water.abund <- sum.water.abund %>% dplyr::rename(RelAbund_SAR86_Clade = SAR86_clade)
sum.water.abund <- sum.water.abund %>% dplyr::rename(RelAbund_Prochlorococcus = Prochlorococcus_MIT9313)
sum.water.abund <- sum.water.abund %>% dplyr::rename(RelAbund_Litoricola = Litoricola)
sum.water.abund <- sum.water.abund %>% dplyr::rename(RelAbund_Synechococcus = Synechococcus_CC9902)
sum.water.abund <- sum.water.abund %>% dplyr::rename(RelAbund_Uncultured_Alteromonadaceaeae = uncultured)

View(sum.water.abund)

##Okay now we have water top abundance for each sample by sample.id



##Coral Top 10 Abundance (includes Endos)
#######################
##Now we need to repeat to get it for corals##
percent.cor <- physeq.coral.r %>% 
  tax_glom(taxrank = "Genus", NArm=TRUE) %>% 
  transform_sample_counts(function(x) {x/sum(x)} )
cor.top = prune_taxa(names(sort(taxa_sums(percent.cor), TRUE))[1:10], percent.cor)
perc.melt.cor <- psmelt(cor.top)
perc.melt.cor$sample.id <- perc.melt.cor$Sample

#summarize genera by sample id
sum.coral.abund <- ddply(perc.melt.cor, c("sample.id", "Genus"), summarise,
                         N = sum(Abundance)
)
#manipulate data frame so that there is an abundance per sample id for each genus
sum.coral.abund <- sum.coral.abund %>% spread(Genus, N)
print(colnames(sum.coral.abund))
#Rename the columns so that we know it's relative abundance
sum.coral.abund <- sum.coral.abund %>% dplyr::rename(RelAbund_Cellvibrionales_Aestuariicella = Aestuariicella)
sum.coral.abund <- sum.coral.abund %>% dplyr::rename(RelAbund_Alteromonadales_Agaribacter = Agaribacter)
sum.coral.abund <- sum.coral.abund %>% dplyr::rename(RelAbund_Cytophagales_Candidatus_Amoebophilus = Candidatus_Amoebophilus)
sum.coral.abund <- sum.coral.abund %>% dplyr::rename(RelAbund_Endozoicomonas = Endozoicomonas)
sum.coral.abund <- sum.coral.abund %>% dplyr::rename(RelAbund_Alteromonadales_Glaciecola = Glaciecola)
sum.coral.abund <- sum.coral.abund %>% dplyr::rename(RelAbund_Litoricola = Litoricola)
sum.coral.abund <- sum.coral.abund %>% dplyr::rename(RelAbund_Oceanospirillales_Neptuniibacter = Neptuniibacter)
sum.coral.abund <- sum.coral.abund %>% dplyr::rename(RelAbund_Cellvibrionales_Porticoccus = Porticoccus)
sum.coral.abund <- sum.coral.abund %>% dplyr::rename(RelAbund_Woesearchaeales_SCGC_AAA286E23 = `SCGC_AAA286-E23`)
sum.coral.abund <- sum.coral.abund %>% dplyr::rename(RelAbund_Uncultured_Alteromonadaceaeae = uncultured)

View(sum.coral.abund)

#merge with erich by sample.id

erich.abund <- merge(erich, sum.water.abund, by = "sample.id", all = TRUE)
erich.abund <- merge(erich.abund, sum.coral.abund, by = "sample.id", all = TRUE)
View(erich.abund)

#coalesce the two columns that have the same names & remove individual ones
erich.abund <- erich.abund %>%
  mutate(RelAbund_Litoricola = coalesce(RelAbund_Litoricola.x, RelAbund_Litoricola.y),
         RelAbund_Uncultured_Alteromonadaceaeae = coalesce(RelAbund_Uncultured_Alteromonadaceaeae.x, RelAbund_Uncultured_Alteromonadaceaeae.y)) %>%
  subset(select = -c(RelAbund_Litoricola.x, RelAbund_Litoricola.y, RelAbund_Uncultured_Alteromonadaceaeae.x, RelAbund_Uncultured_Alteromonadaceaeae.y))

View(erich.abund) #looks good so far! 


#Betadiversity Metrics
###########################
##Okay now for 8,9,& 10
##PCoA/NMDS axes and dispersion

ord.water <- ordinate(physeq.water.r, "NMDS", "bray", trymax = 500) #Run 20 stress 0.061
stressplot(ord.water)
str(ord.water)
scores.water <- (vegan::scores(ord.water))
scores.water <- as.data.frame(scores.water$sites)
scores.water$sample.id <- row.names(scores.water)
row.names(scores.water) <- NULL

ord.coral <- ordinate(physeq.coral.r, "NMDS", "bray", trymax = 1000) #Run 20 stress 0.061
stressplot(ord.coral)
scores.coral <- (vegan::scores(ord.coral))
scores.coral <- as.data.frame(scores.coral$sites)
scores.coral$sample.id <- row.names(scores.coral)
row.names(scores.coral) <- NULL

View(scores.coral)

#combine the two
scores <- rbind(scores.water, scores.coral)

##Merge the scores from coral and water into the erich abund dataframe
erich.abund.scores <- merge(erich.abund, scores, by = "sample.id", all = TRUE)

#just to view for the heck of it:
ggplot(erich.abund.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = island.side, shape = sample.type), size = 3) +
  scale_colour_manual(values = c("#63A022", "#2866AB")) +
  #geom_text(label = scores.outrm$sample.name) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = island.side)) +
  theme_bw() +
  facet_wrap(~motu)

##Merge the scores from coral and water into the 

##################################
###Can we get a measure of dispersion?
##This is a pain to extract from more than one group. AKA island side ON motu

#Since we used bray curtis above
bc.water <- phyloseq::distance(physeq.water.r, method = "bray")
data.water <- as(sample_data(physeq.water.r), "data.frame")
bc.coral <- phyloseq::distance(physeq.coral.r, method = "bray")
data.coral <- as(sample_data(physeq.coral.r), "data.frame")

#In order to do this, we need a variable that is a combined "motu_island.side"
data.water$motu_island.side <- (paste(data.water$motu, data.water$island.side, sep = "_"))
data.coral$motu_island.side <- paste(data.coral$motu, data.coral$island.side, sep = "_")

#Run the dispersion
disp.water <- betadisper(bc.water, data.water$motu_island.side, type = "centroid")
#view bc why not
boxplot(disp.water)

disp.coral <- betadisper(bc.coral, data.coral$motu_island.side, type = "centroid")
boxplot(disp.coral)
str(disp.coral)


#How to extract the values... which are the values? distances...
coral.dists <- as.data.frame(disp.coral$distances)
coral.dists$sample.id <- row.names(coral.dists)
row.names(coral.dists) <- NULL
coral.dists <- coral.dists %>% dplyr::rename(beta_dispersion_motu_islandside = `disp.coral$distances`)


water.dists <- as.data.frame(disp.water$distances)
water.dists$sample.id <- row.names(water.dists)
row.names(water.dists) <- NULL
water.dists <- water.dists %>% dplyr::rename(beta_dispersion_motu_islandside = `disp.water$distances`)


distances <- rbind(water.dists, coral.dists)


##Merge into the erich.abund.scores dataframe

erich.abund.scores.betadisp <- merge(erich.abund.scores, distances, by = "sample.id")
View(erich.abund.scores.betadisp)


##Add a site-name category properly to match the metadata
erich.abund.scores.betadisp$site.name <- paste(erich.abund.scores.betadisp$motu, erich.abund.scores.betadisp$site, sep = "")
erich.abund.scores.betadisp <- erich.abund.scores.betadisp %>%
  mutate(site.name = recode(site.name, aie1 = 'A1', aie2 = 'A2', reiono1 =  'Re1', reiono2 = "Re2", rimatuu1 ="Rm1", rimatuu2 = "Rm2"))


#####################
###VERY LAST STEP
write.csv(erich.abund.scores.betadisp, "../output/nov2021_microbiome_metrics.csv", row.names = FALSE)

