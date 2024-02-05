##Microbiome data vs. seabirds and N15
library(plyr)
library(tidyverse)

#Bring in Microbes and seabird data

#seabirds <- read.csv("../output/n15_seabirds_combined_with_iti.csv", strip.white = T, header = T)
seabirds <- read.csv("../output/n15_seabirds_combined_no_iti.csv", strip.white = T, header = T)
microbes <- read.csv("../../microbiome_analyses/downstream_analyses/integration/output/nov2021_microbiome_metrics.csv")

#recode the site.name category so we can combine with seabird data

microbes$site.name[microbes$site.name == "A1"] <- "Aie_Protected"
microbes$site.name[microbes$site.name == "A2"] <- "Aie_Exposed"
microbes$site.name[microbes$site.name == "Re1"] <- "Reiono_Protected"
microbes$site.name[microbes$site.name == "Re2"] <- "Reiono_Exposed"
microbes$site.name[microbes$site.name == "Rm1"] <- "Rimatuu_Protected"
microbes$site.name[microbes$site.name == "Rm2"] <- "Rimatuu_Exposed"

microbes$site.name <- as.factor(microbes$site.name)
levels(microbes$site.name)


#Combine with seabird data:
microbes.seabirds <- merge(microbes, seabirds, by = "site.name", all = TRUE, no.dups = TRUE)
head(microbes.seabirds)


#Slit into water and coral 
coral.microbes.seabirds <- microbes.seabirds %>% filter(sample.type == "coral")
water.microbes.seabirds <- microbes.seabirds %>% filter(sample.type ==  "water")

#This should allow for any downstream analyses 
library(lme4)
library(car)
library(jtools)
library(emmeans)

mod.micro.birds <- lmer(RelAbund_Endozoicomonas ~ algae.N.percent + (1|site.name), 
                     data = coral.microbes.seabirds)
mod.micro.birds <- lm(RelAbund_Endozoicomonas ~ breeding_biomass_kgha_side, 
                        data = coral.microbes.seabirds)
Anova(mod.micro.birds)
summary(mod.micro.birds)

#Run a loop now for coral and water : 

columns <- c( "Observed" , "Shannon", "FaithPD", "algae.N15", "Evenness", "algae.C13.not.acid", 
              "algae.N.percent", "RelAbund_Endozoicomonas",
              "NMDS1", "NMDS2","beta_dispersion_motu_islandside")
dat <- coral.microbes.seabirds
models <- list()
models2 <- list()
for (i in columns) {
  f <- formula(paste(i, "~ breeding_biomass_kgha_side + (1|site.name)"))
  models[[i]] <- lmer(f, data=dat)
  f2 <- formula(paste(i, "~ algae.N15 + (1|site.name)"))
  models2[[i]] <- lmer(f2, data = dat)
  print(Anova(models[[i]])) 
  print(Anova(models2[[i]]))
}


columns <- c( "Observed" , "Shannon", "FaithPD", "algae.N15", "Evenness", "algae.C13.not.acid", 
              "algae.N.percent","RelAbund_Synechococcus", "RelAbund_Prochlorococcus", "RelAbund_Litoricola",
              "NMDS1", "NMDS2","beta_dispersion_motu_islandside")
dat <- water.microbes.seabirds
models <- list()
models2 <- list()
for (i in columns) {
  f <- formula(paste(i, "~ breeding_biomass_kgha_side + (1|site.name)"))
  models[[i]] <- lmer(f, data=dat)
  f2 <- formula(paste(i, "~ algae.N15 + (1|site.name)"))
  models2[[i]] <- lmer(f2, data = dat)
  print(Anova(models[[i]])) 
  print(Anova(models2[[i]]))
}



##How about adding a seabird level to the dataset?

microbes.seabirds <-
  microbes.seabirds%>%
  mutate(seabird_level = case_when(breeding_biomass_kgha_side<10 ~"low",
                                   breeding_biomass_kgha_side>10&breeding_biomass_kgha_side <200 ~"mid",
                                   breeding_biomass_kgha_side>200 ~"high"))%>%
  mutate(seabird_level = as.factor(seabird_level))%>%
  mutate(seabird_level = fct_relevel(seabird_level, "low", "mid", "high"))

str(microbes.seabirds$distance.along.transect)
microbes.seabirds$Distance_to_shore <- as.factor(microbes.seabirds$distance.along.transect)

str(microbes.seabirds$seabird_level)

##plot----

#pdf(file = "../output/seabird-algaen15/N15vsBreedBiomass_levels_10_200_plus.pdf")


microbes.sum <- ddply(microbes.seabirds, c("sample.type", "Distance_to_shore", "seabird_level"), summarise,
                       mean_richness = mean(Observed),
                       n_richness = length(Observed),
                       se_richness = sd(Observed)/sqrt(n_richness)
)

microbes.sum%>%
  ggplot(aes(x = seabird_level, y = mean_richness, color = Distance_to_shore, fill = Distance_to_shore, group = Distance_to_shore)) +
  geom_point(alpha = 0.5)+
  geom_errorbar(aes(ymin= (mean_richness - se_richness), ymax = (mean_richness + se_richness)), alpha = 0.5, width = 0)+
  geom_line(alpha = 0.5) +
  theme_bw()+
  facet_wrap(~sample.type)



microbes.seabirds%>%
  group_by(Distance_to_shore, seabird_level)%>%
  summarize(mean_richness = mean(Observed),
            n_richness = length(Observed),
            se_richness = sd(Observed)/sqrt(n_richness))%>%
  ggplot(aes(x = seabird_level, y = mean_richness, color = Distance_to_shore, fill = Distance_to_shore, group = Distance_to_shore))+
  geom_point(alpha = .5)+
  geom_errorbar(aes(ymin = (mean_richness-se_richness), ymax = (mean_richness+se_richness)), alpha = .5, width = 0)+
  geom_line(alpha = .5)+
  theme_bw() 

dev.off()


###Test among seabird groups
#First break into coral vs. water again

coral.microbes.seabirds <- microbes.seabirds %>% filter(sample.type == "coral")
water.microbes.seabirds <- microbes.seabirds %>% filter(sample.type ==  "water")

mod.micro.birds <- lmer(RelAbund_Endozoicomonas ~ seabird_level*algae.N15 + (1|site.name), 
                        data = coral.microbes.seabirds)
mod.micro.birds <- lm(RelAbund_Endozoicomonas ~ seabird_level*algae.N15, 
                      data = coral.microbes.seabirds)
Anova(mod.micro.birds)
summary(mod.micro.birds)

#Loop it around the following columns & data

columns <- c( "Observed" , "Shannon", "FaithPD", "algae.N15", "Evenness", "algae.C13.not.acid", 
              "algae.N.percent", "RelAbund_Endozoicomonas","RelAbund_Cellvibrionales_Aestuariicella",
              "RelAbund_Alteromonadales_Agaribacter", "RelAbund_Cytophagales_Candidatus_Amoebophilus",
              "RelAbund_Alteromonadales_Glaciecola","RelAbund_Oceanospirillales_Neptuniibacter",
              "RelAbund_Cellvibrionales_Porticoccus","RelAbund_Woesearchaeales_SCGC_AAA286E23",
              "RelAbund_Litoricola","RelAbund_Uncultured_Alteromonadaceaeae",
              "NMDS1", "NMDS2","beta_dispersion_motu_islandside")
dat <- coral.microbes.seabirds
models <- list()
models2 <- list()

for (i in columns) {
  f <- formula(paste(i, "~ seabird_level + (1|site.name)"))
  models[[i]] <- lmer(f, data=dat)
  f2 <- formula(paste(i, "~ algae.N15 + (1|site.name)"))
  models2[[i]] <- lmer(f2, data = dat)
  print(Anova(models[[i]])) 
  print(Anova(models2[[i]]))
}


#Water
columns <- c( "Observed" , "Shannon", "FaithPD", "algae.N15", "Evenness", "algae.C13.not.acid", 
              "algae.N.percent", "RelAbund_Rhodospirillales_AEGEAN_169_marine_group","RelAbund_SAR11_CladeIa",
              "RelAbund_SAR11_CladeIb", "RelAbund_SAR11_CladeII", "RelAbund_Prochlorococcus", "RelAbund_SAR116_Clade",
              "RelAbund_SAR86_Clade", "RelAbund_Synechococcus",
              "RelAbund_Litoricola","NMDS1", "NMDS2","beta_dispersion_motu_islandside")
dat <- water.microbes.seabirds
models <- list()
models2 <- list()

for (i in columns) {
  f <- formula(paste(i, "~ seabird_level + (1|site.name)"))
  models[[i]] <- lmer(f, data=dat)
  f2 <- formula(paste(i, "~ algae.N15 + (1|site.name)"))
  models2[[i]] <- lmer(f2, data = dat)
  print(Anova(models[[i]])) 
  print(Anova(models2[[i]]))
}






#SAR11, SAR86, Synecchococcus all significant with N15 but not seabird level. 

#combine microbes into site averages
microbes.sum <- ddply(microbes, c("sample.type", "site.name"), summarise,
                                   mean.richness = mean(Observed),
                      mean.shannon = mean(Shannon),
                      mean.faithPD = mean(FaithPD),
                      mean.evenness = mean(Evenness),
                      mean.Rhodospirillales_AEGEAN_169 = mean(RelAbund_Rhodospirillales_AEGEAN_169_marine_group),
                      mean.SAR11.CladeIa = mean(RelAbund_SAR11_CladeIa),
                      mean.SAR11.CladeIb = mean(RelAbund_SAR11_CladeIb),
                      mean.SAR11.CladeII = mean(RelAbund_SAR11_CladeII),
                      mean.SAR116 = mean(RelAbund_SAR116_Clade),
                      mean.SAR86 = mean(RelAbund_SAR86_Clade),
                      mean.Cellvibrionales_Aestuariicella = mean(RelAbund_Cellvibrionales_Aestuariicella),
                      mean.Ca.Ameobophilus = mean(RelAbund_Cytophagales_Candidatus_Amoebophilus),
                      mean.Glaciecola = mean(RelAbund_Alteromonadales_Glaciecola),
                      mean.Meptuniibacter = mean(RelAbund_Alteromonadales_Glaciecola),
                      mean.Porticoccus = mean(RelAbund_Cellvibrionales_Porticoccus),
                      mean.Woesearchaeales = mean(RelAbund_Woesearchaeales_SCGC_AAA286E23),
                      mean.Unc.Alteromonadaceae = mean(RelAbund_Uncultured_Alteromonadaceaeae),
                      mean.Endozoicomonas = mean(RelAbund_Endozoicomonas),
                      mean.Litoricola = mean(RelAbund_Litoricola),
                      mean.Synechococcus = mean(RelAbund_Synechococcus),
                      mean.Prochlorococcus = mean(RelAbund_Prochlorococcus),
                      mean.betadisp = mean(beta_dispersion_motu_islandside)
                                   
)

microbes.sum.coral <- microbes.sum %>% filter(sample.type == "coral")
microbes.sum.coral <- microbes.sum.coral %>% rename(c("coral.mean.richness" = "mean.richness", "coral.mean.shannon" = "mean.shannon",
                                                      "coral.mean.faithPD" = "mean.faithPD", "coral.mean.evenness" = "mean.evenness",
                                                      "coral.mean.aestuariicella" = "mean.Cellvibrionales_Aestuariicella", "coral.mean.Ca.amoebophilus" = "mean.Ca.Ameobophilus",
                                                      "coral.mean.glaciecola" = "mean.Glaciecola", "coral.mean.neptuniibacter" = "mean.Meptuniibacter", 
                                                      "coral.mean.porticoccus" = "mean.Porticoccus", "coral.mean.woesearchaeales" = "mean.Woesearchaeales", 
                                                      "coral.mean.alteromonadaceae" = "mean.Unc.Alteromonadaceae",
                                                      "coral.mean.endozoicomonas" = "mean.Endozoicomonas", "coral.mean.litoricola" = "mean.Litoricola",
                                                      "coral.mean.betadisp" = "mean.betadisp"))
microbes.sum.coral <- microbes.sum.coral[, -c(1,7:12, 22, 23)]

microbes.sum.water <- microbes.sum %>% filter(sample.type == "water")
microbes.sum.water <- microbes.sum.water %>% rename(c("water.mean.richness" = "mean.richness", "water.mean.shannon" = "mean.shannon",
                                                      "water.mean.faithPD" = "mean.faithPD", "water.mean.evenness" = "mean.evenness",
                                                      "water.mean.rhodospirillales_aegean169" = "mean.Rhodospirillales_AEGEAN_169", 
                                                      "water.mean.SAR11.cladeIa" = "mean.SAR11.CladeIa",
                                                      "water.mean.SAR11.CladeIb" = "mean.SAR11.CladeIb",
                                                      "water.mean.SAR11.CladeII" = "mean.SAR11.CladeII", "water.mean.SAR116" = "mean.SAR116", "water.mean.SAR86" = "mean.SAR86",
                                                      "water.mean.alteromonadaceae" = "mean.Unc.Alteromonadaceae",
                                                      "water.mean.litoricola" = "mean.Litoricola",
                                                      "water.mean.Synechococcus" = "mean.Synechococcus", "water.mean.Prochlorococcus" = "mean.Prochlorococcus",
                                                      "water.mean.betadisp" = "mean.betadisp"))
microbes.sum.water <- microbes.sum.water[, -c(1,13:18,20)]

microbes.all <- merge(microbes.sum.coral, microbes.sum.water, by = "site.name", all = TRUE, no.dups = TRUE)
View(microbes.all)
write.csv(microbes.all, "../output/microbes_summary_bysite.csv", row.names = FALSE)

##Combine microbes with seabirds from the seabirds_vs_n15 script *** This is a change in the code - double check & come back

microbes.seabirds <- merge(microbes.all, seabirds, by = "site.name", all = TRUE, no.dups = TRUE)
colnames(microbes.seabirds) 

##Create a matrix
data.matrix <- as.data.frame(microbes.seabirds[,2:50]) #Use only the numerical values


#Run a correlation test using the library corrplot
library(corrplot)

cor.mtest(data.matrix)
correlation.matrix <- cor(data.matrix, use = "pairwise.complete.obs")
write.csv(correlation.matrix, "../output/seabirds_mmicrobes_corrmatrix_with_iti.csv") #change to "_no_iti" if removing iti

pdf(file = "../output/seabirds_v_microbes_with_iti.pdf")

corrplot(cor(data.matrix, use = "pairwise.complete.obs"), type = "upper",
         addCoef.col = NULL, addCoefasPercent = FALSE, tl.col = "black", tl.cex = 0.5, title = "seabirds vs. Microbes")

dev.off()

###### Linear Models ######

#Linear models, using breeding biomass as a predictor 
#this is going to be column 37 "breeding_biomass_kgha_side"

microbes <- read.csv("../../microbiome_analyses/downstream_analyses/integration/output/nov2021_microbiome_metrics.csv")

#Frankenstein the breeding biomass with full microbe data (so we can get confidence intervals)

microbes <- microbes %>% mutate(site.name = recode(site.name, "A1" = "Aie_Protected", "A2" = 'Aie_Exposed',
                                                   "Re1" = "Reiono_Protected", "Re2" = "Reiono_Exposed", 
                                                   "Rm1" = "Rimatuu_Protected", "Rm2" = "Rimatuu_Exposed"))

breeding.biomass <- seabirds[, c(1, 8)]

microbes.seabird.merge <- merge(microbes, breeding.biomass, by = "site.name", all = TRUE)

microbes.seabird.merge %>%
  ggplot(aes(x = site.name, y = alagae.N15, color = species, fill = species))+
  facet_wrap(~species, scales = "free")+
  geom_boxplot(alpha = .5)+
  geom_point(alpha = .7)

library(lme4)
library(car)
library(jtools)
library(emmeans)

microbesn15.mod <- lmer(RelAbund_Endozoicomonas ~ algae.N15 + (1|site.name), 
                            data = microbes.seabird.merge)

summary(microbesn15.mod)
vif(microbesn15.mod)
plot(microbesn15.mod)

plot_summs(microbesn15.mod)
summ(microbesn15.mod)
anova(microbesn15.mod)
Anova(microbesn15.mod)

microbesbirds.mod <- lmer(RelAbund_Endozoicomonas ~ breeding_biomass_kgha_side + (1|site.name), 
                        data = microbes.seabird.merge)

summary(microbesbirds.mod)
vif(microbesbirds.mod)
plot(microbesbirds.mod)

plot_summs(microbesbirds.mod)
summ(microbesbirds.mod)
anova(microbesbirds.mod)
Anova(microbesbirds.mod)


##Attempt a for loop to run through the columns I want
print(colnames(microbes.seabird.merge))
columns <- c( "Observed" , "Shannon", "FaithPD", "algae.N15", "Evenness", "algae.C13.not.acid", 
             "algae.N.percent", "RelAbund_Endozoicomonas","RelAbund_Synechococcus", "RelAbund_Prochlorococcus", "RelAbund_Litoricola",
             "NMDS1", "NMDS2","beta_dispersion_motu_islandside")
dat <- microbes.seabird.merge
models <- list()
models2 <- list()
for (i in columns) {
  f <- formula(paste(i, "~ breeding_biomass_kgha_side + (1|site.name)"))
  models[[i]] <- lmer(f, data=dat)
  f2 <- formula(paste(i, "~ algae.N15 + (1|site.name)"))
  models2[[i]] <- lmer(f2, data = dat)
  print(Anova(models[[i]])) 
  print(Anova(models2[[i]]))
}

##This is great but what about N15 @ 10m? Need to sort out data



