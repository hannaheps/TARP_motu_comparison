##This is the script for alpha diversity of November 2021 only
#This dataset happens to have the most algae nutrient data attached
#We will see...

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


##To start with provided RDS file, set your wd and begin here by uploading the following RDS file ##
##RDS file is pre-processed, using the pre-processing script (found in procedure)
physeq.prune <- readRDS("../output/downstream_analyses/island-survey-2021-pruned.RDS")

rarecurve(t(otu_table(physeq.prune)), step=50, cex=0.5)


#Based on the pre-processing script, we know the lowest valid sample had a read count of 3823
#Unrarefied - cut out every sample with below 1263 reads
physeq.nr <- prune_samples(sample_sums(physeq.prune)>= 1263, physeq.prune)
data.nr <- as(sample_data(physeq.nr), "data.frame")

#Rarefied - sub-sample to 1263
#Set seed for reproducibility pls
physeq.r <- rarefy_even_depth(physeq.prune, sample.size = 1263, rngseed = 711) 
#4 samples removed = TC109 and TC196
#3,404 ASVs were removed because they were no longer present in any sample after sub-sampling
data.r <- as(sample_data(physeq.r), "data.frame")

##Save the rarefied and unrarefied phyloseq objects
saveRDS(physeq.nr, "../output/downstream_analyses/isl_survey_2021_physeq_nonrare.RDS")
saveRDS(physeq.r, "../output/downstream_analyses/isl_survey_2021_physeq_rare.RDS")

#Okay great for comparing across time points but what about just the November 2021 data??

physeq.nov <- subset_samples(physeq.prune, collection.month == "november")
totalreads <- sample_sums(physeq.nov)
View(totalreads) #Lowest viable sample is 1934. We will lose two samples: TC196, TC169

#Rarefy - sub-sample to 1934
#Set seed for reproducibility pls
physeq.nov.r <- rarefy_even_depth(physeq.nov, sample.size = 1934, rngseed = 711) 
#4 samples removed = TC109 and TC196
#4,879 ASVs were removed because they were no longer present in any sample after sub-sampling
data.nov.r <- as(sample_data(physeq.nov.r), "data.frame")
saveRDS(physeq.nov.r, "../output/downstream_analyses/november/nov_2021_rare.RDS")



#Check that all the categorical and numerical values are correct in data frame
str(data.nov.r)
#health.cat, distance.along.transect, site, collection year as character
data.nov.r$health.cat <- as.character(data.nov.r$health.cat)
data.nov.r$distance.along.transect <- as.character(data.nov.r$distance.along.transect)
data.nov.r$site <- as.character(data.nov.r$site)
data.nov.r$collection.year <- as.character(data.nov.r$collection.year)

#Relative abundance
#physeq.ra <- transform_sample_counts(physeq.nr, function(x) x/ sum(x)) 
#data.ra <- as(sample_data(physeq.ra), "data.frame")

#Relative abundance with rarefied data
#physeq.r.ra <- transform_sample_counts(physeq.r, function(x) x/ sum(x)) 
#data.r.ra <- as(sample_data(physeq.r), "data.frame")

##Alpha Diversity using rarefied sample data for November only

?estimate_richness()

#First need to input function estimate_richness_wPD (uses the estimate_richness() command in phyloseq
#but adds faith's phylogenetic distance metric to it as "FaithPD")
#Source code is located in procedure directory (aka your working directory)
source("estimate_richness_wPD.R")
erich <- estimate_richness_wPD(physeq.nov.r, measures = c("Observed", "Shannon", "FaithPD"))
View(erich)
#add all the metadata
erich <- cbind(erich, data.nov.r)
erich$sample.id <- rownames(erich)
print(head(erich))

#check for normalization because microbe data rarely is! Lots of low diversity
hist(erich$Observed)
hist(erich$Shannon)
hist(erich$FaithPD)

hist(log10(erich$Observed))
shapiro.test(log10(erich$Observed))
#W = 0.94357, p-value = 0.004347***
hist(sqrt(erich$Shannon))
shapiro.test(sqrt(erich$Shannon))
#W = 0.97228, p-value = 0.1337
hist(sqrt(erich$FaithPD))
shapiro.test(sqrt(erich$FaithPD))
#W = 0.94773, p-value = 0.006458***

#Quick visualizations to see for patterns??

ggplot(erich, aes(x = algae.N15, y = (Observed))) +
  geom_boxplot() +
  facet_wrap(~motu)

ggplot(erich, aes(x = algae.N15, y = (Observed))) +
  geom_point(aes(color = distance.along.transect),size = 3) +
  geom_smooth(method=lm, aes(color = distance.along.transect)) +
  theme_bw()
  #facet_wrap(~ motu)

library(car)
scatterplot(endo.relabund~algae.N15, erich)
lm.arcsin <- lm(asin(sqrt(endo.relabund)) ~ algae.N15, data = erich)
plot(lm.arcsin)
library(effects)
plot(Effect("asin(sqrt(endo.relabund))",lm.arcsin, partial.residuals=TRUE))
summary(lm.arcsin)
mylogit <- glm(endo.relabund ~ algae.N15, data = erich, family = "binomial")
summary(mylogit)
AIC(mylogit, lm.log)
anova(lm)
summary(lm)
hist(erich$endo.relabund)

summary(lm(endo.relabund ~ motu*island.side, data = erich))
#Call:
#  lm(formula = endo.relabund ~ motu * island.side, data = erich)

#Residuals:
#  Min       1Q   Median       3Q      Max 
#-0.62164 -0.21185  0.06003  0.28491  0.66508 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                      0.56193    0.11554   4.863 8.46e-06 ***
#  motureiono                       0.02364    0.15285   0.155   0.8776    
#moturimatuu                      0.37339    0.15580   2.397   0.0196 *  
#  island.sidewindward             -0.34698    0.15580  -2.227   0.0296 *  
#  motureiono:island.sidewindward   0.45129    0.21047   2.144   0.0360 *  
#  moturimatuu:island.sidewindward  0.02830    0.21262   0.133   0.8946    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.3466 on 61 degrees of freedom
#Multiple R-squared:  0.2893,	Adjusted R-squared:  0.231 
#F-statistic: 4.965 on 5 and 61 DF,  p-value: 0.0007044


ggplot(erich, aes(x = motu, y = dom.relabund)) +
  geom_boxplot(aes(color = island.side), alpha = 0.5) +
  #geom_point(aes(color = island.side), size = 2, alpha = 0.8)
  theme_bw() #+
  #facet_wrap(~motu)
ggsave("../output/downstream_analyses/november/figures/endo_motu.pdf", plot = last_plot())

ggplot(erich, aes(x = motu, y =  endo.relabund)) +
  geom_boxplot(aes(color = motu),alpha = 0.5) +
  geom_point(aes(colour = motu, shape = island.side), size = 3, alpha = 1) +
  theme_bw() #+
  #facet_wrap(~motu)
ggsave("../output/downstream_analyses/november/figures/observed_motu_pooled.pdf", plot = last_plot())




##How about finding the abundances of the top most abundant taxon for each? By genus?
#Start with relative abundance
#Relative abundance with rarefied data
physeq.nov.ra <- transform_sample_counts(physeq.nov.r, function(x) x/ sum(x)) 
data.nov.ra <- as(sample_data(physeq.nov.ra), "data.frame")
str(data.nov.ra)
data.nov.ra$health.cat <- as.character(data.nov.ra$health.cat)
data.nov.ra$distance.along.transect <- as.character(data.nov.ra$distance.along.transect)
data.nov.ra$site <- as.character(data.nov.ra$site)
data.nov.ra$collection.year <- as.character(data.nov.ra$collection.year)

##Find most dominant taxon by iterating over each sample using the following for-loop code: 

df = data.frame()
names <- sample_names(physeq.nov.ra)

for (name in names){
  single_sample <- prune_samples(name, physeq.nov.ra) # This is your single sample
  percent <- single_sample %>% tax_glom(taxrank = "Genus")
  topOTU <- names(sort(taxa_sums(percent), TRUE)[1:1]) 
  newdat <- prune_taxa(topOTU, percent)
  newdatmelt <- psmelt(newdat)

  df <- rbind(df, data.frame(name, newdatmelt$OTU, paste("p_", newdatmelt$Phylum, "c_", newdatmelt$Class, "o_", newdatmelt$Order, "f_",
                                         newdatmelt$Family, "g_", newdatmelt$Genus), newdatmelt$Abundance))
  print(paste("Row added", name))

}

#You end with a df with 4 columns, sample id, the dominant taxon taxonomy string, the relative abundace and the ASV/feature id
View(df)
colnames(df)
df <- dplyr::rename(df, sample.id = name)
df <- dplyr::rename(df, dom.taxon = paste..p_...newdatmelt.Phylum...c_...newdatmelt.Class...o_...)
df <- dplyr::rename(df, dom.relabund = newdatmelt.Abundance)
df <- dplyr::rename(df, dom.asv.name = newdatmelt.OTU)
head(df)

#Add it to the erich data frame here:
erich <- cbind(erich, df)
head(erich)
#Note we have a duplicate column of"sample.id"
#To remove:
duplicated_names <- duplicated(colnames(erich))
erich <- erich[!duplicated_names]
head(erich)

##It is very clear that the most common dominant taxon is endozoicomonas,
#and probably it's this one single ASV: 79b1ad61f0860af468cb41e6ebe3f190
#And then add its relative abundance to our richness data frame

#First thing's first. Out of all taxa, the top most abundant taxon is:
table(erich['dom.taxon']) #52 of 67 individuals were dominated by Endozoicomonas
52/67*100 #That is 77.6% of all individuals were dominated by endo.

erich.endo <- subset(erich, dom.taxon == "p_ Proteobacteria c_ Gammaproteobacteria o_ Oceanospirillales f_ Endozoicomonadaceae g_ Endozoicomonas")
print(table(erich.endo['dom.asv.name'])) 
#4 ASVs dominate. 
#46ded9cc38b410dcfbf32b3d4222717f (12) 79b1ad61f0860af468cb41e6ebe3f190 (24)  908894451235b4ae5c7dcd54d8e5d3be (10)
#f6531c98c934029a9e9443af6c008e95  (6)


#Subset & extract endozoicomonas relative abundance
endo <- subset_taxa(physeq.nov.ra, Genus == "Endozoicomonas")
plot_bar(endo)
merged.endo <- tax_glom(endo, "Genus", NArm = TRUE)
endo.melt <- psmelt(endo)
sum.endo <- ddply(endo.melt,  "Sample", summarise,
             endo.relabund = sum(Abundance)
)
head(endo.melt)
#Insert into our alpha div dataframe
erich$endo.relabund <- sum.endo$endo.relabund
head(erich)

erich$algae.CN.ratio <- erich$algae.N.percent/erich$algae.C.percent

##Save erich with all the things here:
write.csv(erich, "../output/downstream_analyses/november/alphadiversity.csv")
erich <- read.csv("../output/downstream_analyses/november/alphadiversity.csv")
head(erich)



ggplot(erich, aes(x = algae.N15, y = FaithPD)) +
  geom_point(aes(color = motu), size = 3) +
  geom_smooth(method=lm) +
  theme_bw() #+ 
  #facet_wrap(~ motu)


library(car)
hist(log(erich$FaithPD))
scatterplot(FaithPD~algae.N15, erich)
lm <- lm(FaithPD ~ algae.N15, data = erich)
summary(lm)

lm.log <- lm(endo.relabund ~ algae.N15, data = erich)
summary(lm.log)


lm.arcsin <- lm(asin(sqrt(FaithPD)) ~ algae.N15, data = erich)
plot(lm.arcsin)
library(effects)
plot(Effect("asin(sqrt(endo.relabund))",lm.arcsin, partial.residuals=TRUE))
summary(lm.arcsin)
mylogit <- glm(endo.relabund ~ algae.N15, data = erich, family = "binomial")
summary(mylogit)

anova(mylogit)
AIC(mylogit, lm, lm.arcsin)
anova(lm)
summary(lm)
hist(erich$endo.relabund)

sum.n15 <- ddply(erich, c("motu", "island.side"), summarise,
                N = length(algae.N15), 
                mean = mean(algae.N15),
                sd = sd(algae.N15), 
                se = sd/sqrt(N)
)

ggplot(erich, aes(x = algae.N15, y = Observed)) +
  geom_point(size = 2.5) +
  geom_smooth(method=lm, aes(color = motu)) +
  theme_bw() #+ 
  #facet_wrap(~ motu)

ggsave("../output/downstream_analyses/november/figures/n15_vs_Observed.pdf")

library(MASS)
n15.nb <- glm.nb(Observed ~ algae.N15/island.side/motu, data = erich)
summary(n15.nb)
n15.poisson <- glm(Observed ~ algae.N15/island.side/motu, data = erich, family = "poisson")
summary(n15.poisson)
n15.loglm <- lm(log(Observed) ~ algae.N15/island.side/motu, data = erich)
summary(n15.loglm)
n15.lm <- lm(Observed ~ algae.N15/island.side/motu, data = erich)
summary(n15.lm)

AIC(n15.nb, n15.poisson, n15.loglm, n15.lm)
hist(log(erich$Observed))

library(vcd)
fit <- goodfit(erich$Observed, type='nbinomial')
summary(fit)
#Goodness-of-fit test for poisson distribution

#X^2 df    P(> X^2)
#Likelihood Ratio 25.33272 11 0.008146985
rootogram(fit)
Ord_plot(erich$Observed, tol=0.2)
?goodfit()

require("pscl")




n15.nb <- glm.nb()

erich$distance.along.transect <- as.character(erich$distance.along.transect)
sum.ob <- ddply(erich, c("motu", "island.side", "distance.along.transect"), summarise,
                 N = length(Observed), 
                 mean = mean(Observed),
                 sd = sd(Observed), 
                 se = sd/sqrt(Observed)
)


sum.ob$distance.along.transect <- as.character(sum.ob$distance.along.transect)
theme_set(theme_bw())
legend.title = "Island Side"
ggplot(data = sum.ob, aes(x = distance.along.transect, y = mean, group = island.side, color = island.side)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=0.1, show.legend = T) +
  geom_point(aes( color=Side), size=4, show.legend = T) +  #color points by symbiont
  geom_line(aes(color=Side),linetype = "dotted") + #dotted lilnes
  scale_color_manual(legend.title, values=c('#63A022', '#2866AB')) + #set colors manual
  geom_point( shape = 1,size = 4,colour = "black") +
  facet_wrap(~Motu)
ggsave("../output/downstream_analyses/november/figures/ob_motu_all.pdf", plot = last_plot())


