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

pdf(file = "../output/seabird-microbes/RichnessvsBreedBiomass_levels_10_200_plus.pdf")

microbes.sum <- ddply(microbes.seabirds, c("sample.type", "Distance_to_shore","seabird_level"), summarise,
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


dev.off()


###Test among seabird groups
#First break into coral vs. water again
names(microbes.seabirds)[names(microbes.seabirds) == "island.side"] <- "Exposure"
microbes.seabirds$Exposure <- as.factor(microbes.seabirds$Exposure)



coral.microbes.seabirds <- microbes.seabirds %>% filter(sample.type == "coral")
water.microbes.seabirds <- microbes.seabirds %>% filter(sample.type ==  "water")


hist(coral.microbes.seabirds$Observed)
plot(Observed ~ seabird_level, coral.microbes.seabirds)
#Not at all normally distributed

#Therefore, probably can't just loop over the same thing a bunch of times. 
#How to split this up? What do we want to look at? 

#load all the relevant packages
library(MASS) #for glmms
library(glmmTMB) #for glmms with a beta distribution
library(nlme) #for linear mixed effects models
library(lme4) #for linear mixed effects models
library(emmeans) #for multiple comparisons
library(car) #for qqplots
library(vcd) #for distplots of non gaussian distributions
 
####CORALS FIRST #######
#1. Observed Species Richness
hist(coral.microbes.seabirds$Observed) #not normal
qqPlot(coral.microbes.seabirds$Observed) #terrible fit

#Check poisson and nbinomial because these are count data
distplot(coral.microbes.seabirds$Observed, type="poisson") #looks pretty good!
distplot(coral.microbes.seabirds$Observed, typ = "nbinomial") #not good, let's go with Poisson
#using a Posison distribution

mod.glmm <- glmmPQL(Observed~seabird_level*Distance_to_shore + Exposure, random=~1|site.name, data=coral.microbes.seabirds, family="poisson")
summary(mod.glmm)
Anova(mod.glmm)
mod.glmm2 <- glmmPQL(Observed~algae.N15*Distance_to_shore + Exposure, random=~1|site.name, data=coral.microbes.seabirds, family="poisson")
summary(mod.glmm2)
Anova(mod.glmm2)

#Okay these are significant????!!!
#Multiple comparisons can't pull anything out. 
emmeans(mod.glmm, list(pairwise ~ seabird_level*Distance_to_shore), adjust = "fdr")
emmeans(mod.glmm2, list(pairwise~ algae.N15*Distance_to_shore), adjust = "fdr")


#Let's plot yeah

pdf(file = "../output/seabird-microbes/coral_richness_N15.pdf")
ggplot(coral.microbes.seabirds, aes(x = algae.N15, y = Observed)) +
  geom_point(aes(color = seabird_level), size = 4, alpha = 0.8) +
  geom_smooth(aes(group = seabird_level, fill = seabird_level, color = seabird_level), method = "lm", linewidth = 0.5, alpha = 0.1) +
  xlab(expression(italic(delta)^15*N))+
  ylab("Coral Microbiome Species Richness") +
  labs(color='Seabird Biomass \nLevel', fill = 'Seabird Biomass \nLevel') +
  theme_classic() +
  theme(legend.position = c(0.15, 0.8)) +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5)) +
  facet_wrap(~seabird_level)
dev.off()


library(ggeffects)
pdf(file = "../output/seabird-microbes/richness_by_n15_seabird.pdf")
p.predicted <- plot(ggpredict(mod.glmm2, terms = "algae.N15"))
pdf(file = "../output/seabird-microbes/richness_by_n15_seabird.pdf")
p.predicted +
  geom_point(data = coral.microbes.seabirds, aes(x = algae.N15, y = Observed, color = seabird_level), alpha = 0.5, size = 4)+
  scale_fill_manual(values=c("#CC0000", "#006600", "#669999")) +
  ylab("Coral Microbial Species Richness") +
  xlab(expression(italic(delta)^15*N)) +
  theme_classic() +
  ggtitle(NULL) +
  theme(legend.position = c(0.15, 0.82)) +
  labs(color='Seabird Biomass \nLevel') +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))
dev.off()
#Something is weird with the colour, but seems to have produced something okay. 

cor.obs.sum <- ddply(coral.microbes.seabirds, c("Distance_to_shore","seabird_level"), summarise,
                      mean_richness = mean(Observed),
                      n_richness = length(Observed),
                      se_richness = sd(Observed)/sqrt(n_richness)
)

pdf(file = "../output/seabird-microbes/richness_by_seabird_facet.pdf")
ggplot(cor.obs.sum, aes(x = seabird_level, y = mean_richness, color = Distance_to_shore, group = Distance_to_shore)) +
  geom_point(alpha = 0.5, size = 4)+
  geom_errorbar(aes(ymin= (mean_richness - se_richness), ymax = (mean_richness + se_richness)), alpha = 0.5, width = .1)+
  geom_line(alpha = 0.5) +
  ylab("Mean Microbial Species Richness") +
  xlab("Seabird Biomass Level") +
  theme_classic() +
  theme(legend.position = c(1.2, 0.8)) +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))  +
  facet_wrap (~Distance_to_shore)
dev.off()


#2. Shannon Div - Not significant
hist(coral.microbes.seabirds$Shannon)
plot(Shannon ~ seabird_level, coral.microbes.seabirds)
qqPlot(coral.microbes.seabirds$Shannon) #not good

##Can't use poisson though because it isn't count data
#can we use a log link?
y <- coral.microbes.seabirds$Shannon
log.gauss.glm <-glmmPQL(y~ seabird_level*Distance_to_shore + Exposure, random=~1|site.name, data=coral.microbes.seabirds, family=gaussian(link="log"))
summary(log.gauss.glm)
Anova(log.gauss.glm)
plot(log.gauss.glm)
library(sjPlot)
plot_grid(plot_model(log.gauss.glm, type = "diag")) #not bad

log.gauss.glm2<-glmmPQL(y ~ algae.N15*Distance_to_shore + Exposure, random=~1|site.name, data=coral.microbes.seabirds, family=gaussian(link="log"))
summary(log.gauss.glm2)
Anova(log.gauss.glm2)
plot(log.gauss.glm2)
library(sjPlot)
plot_grid(plot_model(log.gauss.glm2, type = "diag")) #not bad


#3. Faith's Phylogenetic Diversity
hist(coral.microbes.seabirds$FaithPD)
plot(FaithPD ~ seabird_level, coral.microbes.seabirds)
qqPlot(coral.microbes.seabirds$FaithPD )#not a good <- fit

#skewed toward 0 
#again can't use poisson because it is not count data. but perhaps a gaussian log link

y <- coral.microbes.seabirds$FaithPD
log.gauss.glm <-glmmPQL(y~ seabird_level*Distance_to_shore + Exposure, random=~1|site.name, data=coral.microbes.seabirds, family=gaussian(link="log"))
summary(log.gauss.glm)
Anova(log.gauss.glm)
plot(log.gauss.glm)
plot_grid(plot_model(log.gauss.glm, type = "diag")) #not bad

log.gauss.glm2<-glmmPQL(y ~ algae.N15*Distance_to_shore + Exposure, random=~1|site.name, data=coral.microbes.seabirds, family=gaussian(link="log"))
summary(log.gauss.glm2)
Anova(log.gauss.glm2)
plot(log.gauss.glm2)
plot_grid(plot_model(log.gauss.glm2, type = "diag")) #its okay

#significant effect of N15 on faith's PD
#residuals vs fitted look fine

#4. Evenness 
hist((coral.microbes.seabirds$Evenness))
plot(Eveness ~ seabird_level, coral.microbes.seabirds)
qqPlot(coral.microbes.seabirds$Evenness) #not bad

#use standard linear model
mod.lme<- lme(Evenness~seabird_level*Distance_to_shore + Exposure, random=~1|site.name, data=coral.microbes.seabirds, method = "REML")
summary(mod.lme)
Anova(mod.lme)
plot(mod.lme) #looks great

mod.lme2 <- lme(Evenness~algae.N15*Distance_to_shore + Exposure, random=~1|site.name, data=coral.microbes.seabirds, method = "REML", na.action = na.omit)
summary(mod.lme2)
Anova(mod.lme2)
plot(mod.lme2) #also looks great 

#Nothing significant

##Next is beta div: NMDS1, NMDS2
#5. NMDS1
hist(coral.microbes.seabirds$NMDS1)
plot(NMDS1 ~ seabird_level, coral.microbes.seabirds)
qqPlot(coral.microbes.seabirds$NMDS1) #horrible
#Problem with negative numbers since it's a coordinate????

mod.lme<- lme(NMDS1~seabird_level*Distance_to_shore + Exposure, random=~1|site.name, data=coral.microbes.seabirds, method = "REML")
summary(mod.lme)
Anova(mod.lme)
anova(mod.lme)
plot(mod.lme) #honestly the residuals look great

mod.lme2<- lme(NMDS1~algae.N15*Distance_to_shore + Exposure, random=~1|site.name, data=coral.microbes.seabirds, method = "REML", na.action = na.omit)
summary(mod.lme2)
Anova(mod.lme2)
anova(mod.lme2)
plot(mod.lme2)
library(sjPlot)
plot_grid(plot_model(mod.lme2, type = "diag")) #its okay

#6.NMDS2
hist(coral.microbes.seabirds$NMDS2) #looks normal
qqPlot(coral.microbes.seabirds$NMDS2) # not bad, not great though

mod.lme<- lme(NMDS2~seabird_level*Distance_to_shore + Exposure, random=~1|site.name, data=coral.microbes.seabirds, method = "REML")
summary(mod.lme)
Anova(mod.lme)
plot(mod.lme) #looks great
plot_grid(plot_model(mod.lme, type = "diag")) #nice

mod.lme2 <- lme(NMDS2~algae.N15*Distance_to_shore + Exposure, random=~1|site.name, data=coral.microbes.seabirds, method = "REML", na.action = na.omit)
summary(mod.lme2)
Anova(mod.lme2)
plot(mod.lme2) #also looks fine

#Nothing Significant

##Next are all of the individual taxa: 
#But this is super annoying because I can't loop. Try Endoz? 
hist(coral.microbes.seabirds$RelAbund_Endozoicomonas)
qqPlot(coral.microbes.seabirds$RelAbund_Endozoicomonas)#terrible
#Hmmm can't use log-link because it is skewed the opposite way. 
#What about proportional data?

##Since Endoz are skewed toward 1, what about using a beta distribution?
?glmmTMB()

glmmTMB(RelAbund_Endozoicomonas ~ seabird_level*Distance_to_shore + Exposure + (1|site.name), family = beta_family(), data = coral.microbes.seabirds, na.action = na.omit)
#says that the values need to be between zero and 1
coral.microbes.seabirds$RelAbund_Endozoicomonas #but they are....

#Maybe we can do beta but need a transformation first?
library(mgcv)
endo_trans <- mutate(coral.microbes.seabirds, RelAbund_Endozoicomonas = (RelAbund_Endozoicomonas*(n()-1) + 0.5)/n())

beta <- gam(RelAbund_Endozoicomonas ~ seabird_level*Distance_to_shore,data=endo_trans,
                    family=betar(link="logit"), na.action = na.omit)


plot.gam(beta)
anova(beta)
summary(beta)

beta1 <- glmmTMB(RelAbund_Endozoicomonas ~ seabird_level*Distance_to_shore + Exposure + (1|site.name), family = beta_family(), data = endo_trans, na.action = na.omit)
beta2 <- glmmTMB(RelAbund_Endozoicomonas ~ algae.N15*Distance_to_shore + Exposure + (1|site.name), family = beta_family(), data = endo_trans, na.action = na.omit)

summary(beta1)
Anova(beta1)

Anova(beta2)

##Maybe Ca. Amoebophilus
#Alteromonadaceae
#And Woesarchaeales, Litoricola

#This gives us 5 top abundant taxa 
#B) Ca. Amoebophilus
ameob_trans <- mutate(coral.microbes.seabirds, RelAbund_Cytophagales_Candidatus_Amoebophilus = (RelAbund_Cytophagales_Candidatus_Amoebophilus*(n()-1) + 0.5)/n())
amoeb.beta1 <- glmmTMB(RelAbund_Cytophagales_Candidatus_Amoebophilus ~ seabird_level*Distance_to_shore + Exposure + (1|site.name), family = beta_family(), data = ameob_trans, na.action = na.omit)
amoeb.beta2 <- glmmTMB(RelAbund_Cytophagales_Candidatus_Amoebophilus ~ algae.N15*Distance_to_shore + Exposure + (1|site.name), family = beta_family(), data = ameob_trans, na.action = na.omit)

Anova(amoeb.beta1)
Anova(amoeb.beta2)#does not explain


#C)Alteromonadaceae
hist(coral.microbes.seabirds$RelAbund_Uncultured_Alteromonadaceaeae)
alter_trans <- mutate(coral.microbes.seabirds, RelAbund_Uncultured_Alteromonadaceaeae = (RelAbund_Uncultured_Alteromonadaceaeae*(n()-1) + 0.5)/n())
alter.beta1 <- glmmTMB(RelAbund_Uncultured_Alteromonadaceaeae ~ seabird_level*Distance_to_shore + Exposure + (1|site.name), family = beta_family(), data = alter_trans, na.action = na.omit)
alter.beta2 <- glmmTMB(RelAbund_Uncultured_Alteromonadaceaeae ~ algae.N15*Distance_to_shore + Exposure + (1|site.name), family = beta_family(), data = alter_trans, na.action = na.omit)
Anova(alter.beta1)
Anova(alter.beta2)#woah! Does explain :) Let's plot 
library(ggeffects)
p <- plot(ggpredict(alter.beta2, terms = "algae.N15")) #this plots the predicted values

#add real data to this:
pdf(file = "../output/seabird-microbes/alteromon_true_betapredicted_data_by_n15.pdf")
p + geom_point(data = coral.microbes.seabirds, aes(x = algae.N15, y = RelAbund_Uncultured_Alteromonadaceaeae), alpha = .5, size = 4) +
  xlab(expression(italic(delta)^15*N))+
  ylab("Relative Abundance of Uncultured Alteromonadaceae \nin Coral Tissue") +
  ggtitle( label = NULL) +
  scale_y_continuous(n.breaks = 8) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))
dev.off()


#D)Woesarchaeales
coral.microbes.seabirds$RelAbund_Woesearchaeales_SCGC_AAA286E23
woes_trans <- mutate(coral.microbes.seabirds, RelAbund_Woesearchaeales_SCGC_AAA286E23 = (RelAbund_Woesearchaeales_SCGC_AAA286E23*(n()-1) + 0.5)/n())
woes.beta1 <- glmmTMB(RelAbund_Woesearchaeales_SCGC_AAA286E23 ~ seabird_level*Distance_to_shore + Exposure + (1|site.name), family = beta_family(), data = woes_trans, na.action = na.omit)
woes.beta2 <- glmmTMB(RelAbund_Woesearchaeales_SCGC_AAA286E23 ~ algae.N15*Distance_to_shore + Exposure + (1|site.name), family = beta_family(), data = woes_trans, na.action = na.omit)
Anova(woes.beta1)
Anova(woes.beta2)

#E)Litoricola
coral.microbes.seabirds$RelAbund_Litoricola
lit_trans <- mutate(coral.microbes.seabirds, RelAbund_Litoricola = (RelAbund_Litoricola*(n()-1) + 0.5)/n())
lit.beta1 <- glmmTMB(RelAbund_Litoricola ~ seabird_level*Distance_to_shore + Exposure + (1|site.name), family = beta_family(), data = lit_trans, na.action = na.omit)
lit.beta2 <- glmmTMB(RelAbund_Litoricola~ algae.N15*Distance_to_shore + Exposure + (1|site.name), family = beta_family(), data = lit_trans, na.action = na.omit)
Anova(lit.beta1)
Anova(lit.beta2)

####WATER NEXT#####
#Okay cool let's break and do this all for water 
#1. Observed Species Richness
hist(water.microbes.seabirds$Observed) #not normal
qqPlot(water.microbes.seabirds$Observed) #not bad

#Check poisson and nbinomial
distplot(water.microbes.seabirds$Observed, type="poisson") #looks literally perfect
#using a Poisson distribution

mod.glmm <- glmmPQL(Observed~seabird_level*Distance_to_shore + Exposure, random=~1|site.name, data=water.microbes.seabirds, family="poisson")
summary(mod.glmm)
Anova(mod.glmm)
mod.glmm2 <- glmmPQL(Observed~algae.N15*Distance_to_shore + Exposure, random=~1|site.name, data=water.microbes.seabirds, family="poisson")
summary(mod.glmm2)
Anova(mod.glmm2)

#not significant 


##2. Shannon 
hist(water.microbes.seabirds$Shannon) #not normal
qqPlot(water.microbes.seabirds$Shannon) #oof
#almost binomial - two peaks

y <- water.microbes.seabirds$Shannon
log.gauss.glm <-glmmPQL(y~ seabird_level*Distance_to_shore + Exposure, random=~1|site.name, data=water.microbes.seabirds, family=gaussian(link="log"))
summary(log.gauss.glm)
Anova(log.gauss.glm)
plot(log.gauss.glm)
plot_grid(plot_model(log.gauss.glm, type = "diag")) #not great but not bad

log.gauss.glm2<-glmmPQL(y ~ algae.N15*Distance_to_shore + Exposure, random=~1|site.name, data=water.microbes.seabirds, family=gaussian(link="log"))
summary(log.gauss.glm2)
Anova(log.gauss.glm2)
plot(log.gauss.glm2)
plot_grid(plot_model(log.gauss.glm2, type = "diag")) #also a little rough but okay



#

#3. Faith's Phylogenetic Diversity 
hist(water.microbes.seabirds$FaithPD)
qqPlot(water.microbes.seabirds$FaithPD) #also bad

y <- water.microbes.seabirds$FaithPD
log.gauss.glm <-glmmPQL(y~ seabird_level*Distance_to_shore + Exposure, random=~1|site.name, data=water.microbes.seabirds, family=gaussian(link="log"))
summary(log.gauss.glm)
Anova(log.gauss.glm)
plot(log.gauss.glm)
plot_grid(plot_model(log.gauss.glm, type = "diag")) #honestly not so bad

log.gauss.glm2<-glmmPQL(y ~ algae.N15*Distance_to_shore + Exposure, random=~1|site.name, data=water.microbes.seabirds, family=gaussian(link="log"))
summary(log.gauss.glm2)
Anova(log.gauss.glm2)
plot(log.gauss.glm2)
plot_grid(plot_model(log.gauss.glm2, type = "diag")) # okay



#4. Evenness
hist(water.microbes.seabirds$Evenness)
qqPlot(water.microbes.seabirds$Evenness) #NOOoooooo

y <- water.microbes.seabirds$Evenness
log.gauss.glm <-glmmPQL(y~ seabird_level*Distance_to_shore + Exposure, random=~1|site.name, data=water.microbes.seabirds, family=gaussian(link="log"))
summary(log.gauss.glm)
Anova(log.gauss.glm)
plot(log.gauss.glm)
plot_grid(plot_model(log.gauss.glm, type = "diag")) #Not good. 

log.gauss.glm2<-glmmPQL(y ~ algae.N15*Distance_to_shore + Exposure, random=~1|site.name, data=water.microbes.seabirds, family=gaussian(link="log"))
summary(log.gauss.glm2)
Anova(log.gauss.glm2)
plot(log.gauss.glm2)
plot_grid(plot_model(log.gauss.glm2, type = "diag")) # also not good. ugh


#5. NMDS1
hist(water.microbes.seabirds$NMDS1)
qqPlot(water.microbes.seabirds$NMDS1)

#let's just try a linear model and check residuals
mod.lme<- lme(NMDS1~seabird_level*Distance_to_shore + Exposure, random=~1|site.name, data=water.microbes.seabirds, method = "REML")
summary(mod.lme)
Anova(mod.lme)
plot(mod.lme) #looks great
plot_grid(plot_model(mod.lme, type = "diag")) #not great

mod.lme2 <- lme(NMDS1~algae.N15*Distance_to_shore + Exposure, random=~1|site.name, data=water.microbes.seabirds, method = "REML", na.action = na.omit)
summary(mod.lme2)
Anova(mod.lme2)
plot(mod.lme2) #also looks fine
plot_grid(plot_model(mod.lme2, type = "diag")) #its kinda weird


#6. NMDS2
hist(water.microbes.seabirds$NMDS2) #wow what a relief!
qqPlot(water.microbes.seabirds$NMDS2)

mod.lme<- lme(NMDS2~seabird_level*Distance_to_shore + Exposure, random=~1|site.name, data=water.microbes.seabirds, method = "REML")
summary(mod.lme)
Anova(mod.lme)
plot(mod.lme) #looks great
plot_grid(plot_model(mod.lme, type = "diag")) #nice nice

mod.lme2 <- lme(NMDS2~algae.N15*Distance_to_shore + Exposure, random=~1|site.name, data=water.microbes.seabirds, method = "REML", na.action = na.omit)
summary(mod.lme2)
Anova(mod.lme2)
plot(mod.lme2) #also looks fine
plot_grid(plot_model(mod.lme2, type = "diag")) #good good

#This was a good breather! 

#Okay done with the diversity metrics, how about let's look at some of the top taxa

##7. Synnechococcus 
hist(water.microbes.seabirds$RelAbund_Synechococcus)
qqPlot(water.microbes.seabirds$RelAbund_Synechococcus)
#use a beta I think

syn_trans <- mutate(water.microbes.seabirds, RelAbund_Synechococcus = (RelAbund_Synechococcus*(n()-1) + 0.5)/n())
syn.beta1 <- glmmTMB(RelAbund_Synechococcus ~ seabird_level*Distance_to_shore + Exposure + (1|site.name), family = beta_family(), data = syn_trans, na.action = na.omit)
syn.beta2 <- glmmTMB(RelAbund_Synechococcus ~ algae.N15*Distance_to_shore + Exposure + (1|site.name), family = beta_family(), data = syn_trans, na.action = na.omit)
Anova(syn.beta1)
Anova(syn.beta2)

#Can we remove protected side of Aie?
water.microbes.seabirds.noaiepro <- filter(water.microbes.seabirds, site.name != "Aie_Protected")

syn_trans <- mutate(water.microbes.seabirds.noaiepro, RelAbund_Synechococcus = (RelAbund_Synechococcus*(n()-1) + 0.5)/n())
syn.beta1 <- glmmTMB(RelAbund_Synechococcus ~ seabird_level*Distance_to_shore + Exposure + (1|site.name), family = beta_family(), data = syn_trans, na.action = na.omit)
syn.beta2 <- glmmTMB(RelAbund_Synechococcus ~ algae.N15*Distance_to_shore + Exposure + (1|site.name), family = beta_family(), data = syn_trans, na.action = na.omit)
Anova(syn.beta1)
Anova(syn.beta2)
# Synecchococcus significant with N15 but not seabird level. 

water.microbes.seabirds.noaiepro %>%
  group_by(seabird_level)%>%
  summarize(mean_syn = mean(RelAbund_Synechococcus),
            n_syn = length(RelAbund_Synechococcus),
            se_syn = sd(RelAbund_Synechococcus)/sqrt(n_syn))

ggplot(water.microbes.seabirds.noaiepro, aes(x = algae.N15, y = RelAbund_Synechococcus)) +
  geom_point(aes(color = seabird_level, shape = site.name))+
  theme_bw()


##What about Litoricola, Prochlorococcus, Uncultured Alteromonadaceae (same as corals) & SAR86

#8. Litoricola
hist(water.microbes.seabirds$RelAbund_Litoricola)
qqPlot(water.microbes.seabirds$RelAbund_Litoricola) #beta

syn_trans <- mutate(water.microbes.seabirds, RelAbund_Litoricola = (RelAbund_Litoricola*(n()-1) + 0.5)/n())
syn.beta1 <- glmmTMB(RelAbund_Litoricola ~ seabird_level*Distance_to_shore + Exposure + (1|site.name), family = beta_family(), data = syn_trans, na.action = na.omit)
syn.beta2 <- glmmTMB(RelAbund_Litoricola ~ algae.N15*Distance_to_shore + Exposure + (1|site.name), family = beta_family(), data = syn_trans, na.action = na.omit)
Anova(syn.beta1)
Anova(syn.beta2)



#9. Prochlorococcus

hist(water.microbes.seabirds$RelAbund_Prochlorococcus)
qqPlot(water.microbes.seabirds$RelAbund_Prochlorococcus) #beta

syn_trans <- mutate(water.microbes.seabirds, RelAbund_Prochlorococcus = (RelAbund_Prochlorococcus*(n()-1) + 0.5)/n())
syn.beta1 <- glmmTMB(RelAbund_Prochlorococcus ~ seabird_level*Distance_to_shore + Exposure + (1|site.name), family = beta_family(), data = syn_trans, na.action = na.omit)
syn.beta2 <- glmmTMB(RelAbund_Prochlorococcus ~ algae.N15*Distance_to_shore + Exposure + (1|site.name), family = beta_family(), data = syn_trans, na.action = na.omit)
Anova(syn.beta1)
Anova(syn.beta2)


#10. SAR86
hist(water.microbes.seabirds$RelAbund_SAR86_Clade)
qqPlot(water.microbes.seabirds$RelAbund_SAR86_Clade) #beta

syn_trans <- mutate(water.microbes.seabirds, RelAbund_SAR86_Clade = (RelAbund_SAR86_Clade*(n()-1) + 0.5)/n())
syn.beta1 <- glmmTMB(RelAbund_SAR86_Clade ~ seabird_level*Distance_to_shore + Exposure + (1|site.name), family = beta_family(), data = syn_trans, na.action = na.omit)
syn.beta2 <- glmmTMB(RelAbund_SAR86_Clade ~ algae.N15*Distance_to_shore + Exposure + (1|site.name), family = beta_family(), data = syn_trans, na.action = na.omit)
Anova(syn.beta1)
Anova(syn.beta2)

#Uncultured Alteromonadaceae

hist(water.microbes.seabirds$RelAbund_Uncultured_Alteromonadaceaeae)
qqPlot(water.microbes.seabirds$RelAbund_Uncultured_Alteromonadaceaeae) #beta

syn_trans <- mutate(water.microbes.seabirds, RelAbund_Uncultured_Alteromonadaceaeae = (RelAbund_Uncultured_Alteromonadaceaeae*(n()-1) + 0.5)/n())
syn.beta1 <- glmmTMB(RelAbund_Uncultured_Alteromonadaceaeae ~ seabird_level*Distance_to_shore + Exposure + (1|site.name), family = beta_family(), data = syn_trans, na.action = na.omit)
syn.beta2 <- glmmTMB(RelAbund_Uncultured_Alteromonadaceaeae ~ algae.N15*Distance_to_shore + Exposure + (1|site.name), family = beta_family(), data = syn_trans, na.action = na.omit)
Anova(syn.beta1)
Anova(syn.beta2)



#Plot endo - no differences according to distance from shore but definitely trending toward higher at low seabird_levels (confounded by site though)

pdf(file = "../output/seabird-microbes/endo_true_data_by_seabirds.pdf")
endo.sum <- ddply(coral.microbes.seabirds, c("seabird_level"), summarise,
                  mean_endo = mean(RelAbund_Endozoicomonas),
                  n_endo = length(RelAbund_Endozoicomonas),
                  se_endo = sd(RelAbund_Endozoicomonas)/sqrt(n_endo) )
pdf(file = "../output/seabird-microbes/endo_truedata_seabirdlevel.pdf")
ggplot(endo.sum, aes(x = seabird_level, y = mean_endo))+
  geom_point(alpha = .5, size = 4, aes(color = seabird_level))+
  geom_errorbar(aes(ymin = (mean_endo-se_endo), ymax = (mean_endo+se_endo), color = seabird_level), alpha = .5, width = 0.2)+
  xlab("Seabird Biomass Level")+
  ylab("Mean Relative Abundance of Endozoicomonas \nin Coral Tissue") +
  labs(color='Seabird Biomass \nLevel') +
  theme_classic() +
  theme(legend.position = c(0.85, 0.82)) +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))

dev.off()

#What about a plot that is like a barplot or boxplot or something
ggplot(endo.sum, aes(x = seabird_level, y = mean_endo, fill = Distance_to_shore))+
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = (mean_endo-se_endo), ymax = (mean_endo+se_endo)), position = "dodge", alpha = .5, width = 0.2)+
  #geom_line(alpha = .5)+
  theme_bw() 
#omg totally not interesting 







######Deprecated###########
#Can we plot the model
##Just be aware this model has bad fit!!!
plot(mod.micro.birds)

emm.summary.sb.endo<-
  mod.micro.birds %>% 
  emmeans(~ seabird_level) 

emm_plot<-emmip(mod.micro.birds,  ~ seabird_level, CIs = TRUE, level = c(0.75),
                CIarg = list(lwd = .9, alpha = .25), lwd = 2, linearg = list(lwd = 1.2, alpha = .8), dotarg = list(lwd = 1.2, alpha = .8))


emm_plot_sb_endo<-
  emm_plot+
  ylab("Endozoicomonas Rel Abundance")+
  xlab("seabid biomass")+
  #labs(colour = "fish species")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), # remove gridlines
        panel.grid.minor = element_blank(), #remove gridlines
        strip.background = element_blank(), 
        # text = element_text(size = 14, family = "Arial"),
        rect = element_rect(fill = "transparent"),  
        plot.background = element_rect(fill = "transparent", color = NA)
  )

pdf(file = "../output/seabird-microbes/model_endo_by_seabirds.pdf")
emm_plot_sb_endo
dev.off()





########DEPRECATED MAYBE - INCLUDES FOR LOOP######
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
  f2 <- formula(paste(i, "~ algae.N15 + Distance_to_shore + island.side  + (1|site.name)"))
  models2[[i]] <- lmer(f2, data = dat)
  print(Anova(models[[i]])) 
  print(Anova(models2[[i]]))
}

Anova(lmer(NMDS1 ~ algae.N15 + Distance_to_shore + island.side + (1|site.name), data = coral.microbes.seabirds))

#Beta dispersion_motu_islandside is significant but this is maybe irrelevant 
#need to make a new NMDS and colour by seabird level??
ggplot(coral.microbes.seabirds, aes (x = NMDS1, y = NMDS2))+
  geom_point(aes(color = seabird_level, shape = Exposure), size = 4, alpha = 0.8)+
  stat_ellipse(aes(color = seabird_level)) +
  #facet_wrap(~motu) +
  theme_bw()

ggplot(coral.microbes.seabirds, aes(x = NMDS1, y = algae.N15)) +
  geom_point(aes(color = seabird_level, shape =Exposure)) +
  geom_smooth(method = lm) +
  theme_bw()

pdf(file = "../output/seabird-microbes/coral_richness_N15.pdf")
ggplot(coral.microbes.seabirds, aes(x = algae.N15, y = Observed)) +
  geom_point(aes(color = seabird_level), size = 4, alpha = 0.8) +
  geom_smooth(aes(group = seabird_level, fill = seabird_level, color = seabird_level), method = "lm", linewidth = 0.5, alpha = 0.1) +
  xlab(expression(italic(delta)^15*N))+
  ylab("Coral Microbiome Species Richness") +
  labs(color='Seabird Biomass \nLevel', fill = 'Seabird Biomass \nLevel') +
  theme_classic() +
  theme(legend.position = c(0.15, 0.8)) +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5)) 
dev.off()



emm_plot_sb_endo<-
  emm_plot+
  ylab("Endozoicomonas Rel Abundance")+
  xlab("seabid biomass")+
  #labs(colour = "fish species")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), # remove gridlines
        panel.grid.minor = element_blank(), #remove gridlines
        strip.background = element_blank(), 
        # text = element_text(size = 14, family = "Arial"),
        rect = element_rect(fill = "transparent"),  
        plot.background = element_rect(fill = "transparent", color = NA)
  )

pdf(file = "../output/seabird-microbes/model_endo_by_seabirds.pdf")
emm_plot_sb_endo
dev.off()


ggplot(water.microbes.seabirds, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = seabird_level, shape = seabird_level)) +
  stat_ellipse(aes(color = seabird_level)) +
  facet_wrap(~motu) +
  theme_bw()



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
  f2 <- formula(paste(i, "~ algae.N15 + Distance_to_shore + island.side + (1|site.name)"))
  models2[[i]] <- lmer(f2, data = dat)
  print(Anova(models[[i]])) 
  print(Anova(models2[[i]]))
}


anova(lmer(RelAbund_Synechococcus ~ algae.N15 + Distance_to_shore + island.side + (1|site.name), data = water.microbes.seabirds))
# Synecchococcus significant with N15 but not seabird level. 


water.microbes.seabirds %>%
  group_by(seabird_level)%>%
  summarize(mean_syn = mean(RelAbund_Synechococcus),
            n_syn = length(RelAbund_Synechococcus),
            se_syn = sd(RelAbund_Synechococcus)/sqrt(n_syn))
pdf(file = "../output/synechococcus_N15.pdf")
ggplot(water.microbes.seabirds, aes(x = algae.N15, y = RelAbund_Synechococcus)) +
  geom_point(aes(color = seabird_level, shape = site.name), size = 4)+
  ylab("Mean Relative Abundance Synechococcus") +
  xlab(expression(italic(delta)^15*N))+
  labs(color='Seabird Biomass \nLevel', fill = 'Seabird Biomass \nLevel') +
  theme_classic() +
  #theme(legend.position = c(0.15, 0.8)) +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5)) 

dev.off()
##What about removing Aie Protected and trying again if there is a strong relationship?
water.microbes.seabirds.noaiepro <- filter(water.microbes.seabirds, site.name != "Aie_Protected")

anova(lmer(RelAbund_Synechococcus ~ algae.N15 + Distance_to_shore + island.side + (1|site.name), data = water.microbes.seabirds.noaiepro))
# Synecchococcus significant with N15 but not seabird level. 


water.microbes.seabirds.noaiepro %>%
  group_by(seabird_level)%>%
  summarize(mean_syn = mean(RelAbund_Synechococcus),
            n_syn = length(RelAbund_Synechococcus),
            se_syn = sd(RelAbund_Synechococcus)/sqrt(n_syn))

ggplot(water.microbes.seabirds.noaiepro, aes(x = algae.N15, y = RelAbund_Synechococcus)) +
  geom_point(aes(color = seabird_level, shape = site.name))+
  theme_bw()


#Can we do some NMDS at least for coral? 


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


