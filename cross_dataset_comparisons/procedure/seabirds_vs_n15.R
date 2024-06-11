##Combine the Seabird data vs. N15

##Bring in seabird data
seabirds <- read.csv("../../seabird_data/outputs/seabird_dens_bio_focal_data.csv")

seabirds$site.name <- paste(seabirds$Motu, seabirds$Exposure, sep = "_")


##bring in algae data
algae <- read.csv("../../algae_isotopes/data/Tetiaroa_Turbinaria_Transects_November_2021_compiledMarch2023.csv")
View(algae)
algae$Distance_to_shore <- as.factor(algae$Distance_to_shore)
algae$Side[algae$Side=="leeward"] <- "Protected"
algae$Side[algae$Side == "windward"] <- "Exposed"
names(algae)[names(algae) == 'Side'] <- 'Exposure'

#manipulate algae data frame to get average N15 per motu, site and distance from shore

library(plyr)
sum.algae <- ddply(algae, c("Motu", "Exposure", "Distance_to_shore"),summarise,
                      mean.N15 = mean(N15)
)

sum.algae
#spread data to get values for each distance from shore
library(tidyverse)
sum.algae <- sum.algae %>% spread(Distance_to_shore, mean.N15)

sum.algae <- sum.algae %>% rename(c("10" = "N.15_at_10m", 
                                    "20" = "N.15_at_20m", "30" = "N.15_at_30m",
                                    "40" = "N.15_at_40m"))


#Create a site.name that combines motu and exposure in order to merge with seabrid data
sum.algae$site.name <-  paste(sum.algae$Motu, sum.algae$Exposure, sep = "_")


##Combine seabird and algae data
algae.seabirds <- right_join(seabirds, sum.algae)
View(algae.seabirds)
algae.seabirds <- algae.seabirds[,-1] #This is just to remove duplicate columns
#algae.seabirds <- algae.seabirds[,-c(2,7)] #Use this one for Iti data, it has slightly different column #s 
View(algae.seabirds)

#save combined data
write.csv(algae.seabirds, "../output/n15_seabirds_combined_no_iti.csv") #change to "_no_iti" if removing iti 

##Create a matrix
#data.matrix <- as.data.frame(algae.seabirds[ , 2:21])
data.matrix <- as.data.frame(algae.seabirds[,c(3:18,20:23)])
View(data.matrix)


#Run a correlation test using the library corrplot
library(corrplot)

cor.mtest(data.matrix)
correlation.matrix.psn <- cor(data.matrix, use = "pairwise.complete.obs",) 
write.csv(correlation.matrix.psn, "../output/n15_seabirds_corrmatrix_no_iti_pearson.csv") #change to "_no_iti" if removing iti
correlation.matrix.sp <- cor(data.matrix, use = "pairwise.complete.obs", method  = "spearman") 
write.csv(correlation.matrix.sp, "../output/n15_seabirds_corrmatrix_no_iti_spearman.csv")

pdf(file = "../output/correlation_seabird_N15_spearman.pdf") #change to "_no_iti" if removing iti

corrplot(cor(data.matrix, use = "pairwise.complete.obs", method = "spearman"), type = "upper", 
         addCoef.col = NULL, addCoefasPercent = FALSE, tl.col = "black", title = "seabirds vs. N15")

dev.off()



##Below are the linear models
#Re-organisation of the dataframe was necessary to include not only the summary
#data for the algae
################################
##Combine the Seabird data vs. N15 again
##Analyses using full data (no algal summary)
#start over in your environment: 
rm(list = ls())

##Bring in seabird data again
seabirds.all <- read.csv("../../seabird_data/outputs/seabird_dens_bio_focal_data.csv")
View(seabirds.all)

#clean-up the dataframe by removing the X column (column 1) accidentally added in the csv file.
#At the same time, change the name to maintain consistency in downstream code
seabirds <- seabirds.all[,-1]
View(seabirds)
#add a site.name (Motu_Exposure) category separated by an underscore.
seabirds$site.name <- paste(seabirds$Motu, seabirds$Exposure, sep = "_")

##Bring in algae data
algae <- read.csv("../../algae_isotopes/data/Tetiaroa_Turbinaria_Transects_November_2021_compiledMarch2023.csv")
View(algae)

#Clean-up the dataframe by making distance_to_shore a factor 
algae$Distance_to_shore <- as.factor(algae$Distance_to_shore)
algae$Site[algae$Site=="1"] <- "Protected"
algae$Site[algae$Site == "2"] <- "Exposed"
#add a site.name (Motu_Site) category separated by an underscore to allow merging with seabird data
algae$site.name <- paste(algae$Motu, algae$Site, sep = "_")

#Combine the algae and seabird data 
algae.seabirds <- merge(algae, seabirds, by = "site.name", all = TRUE, no.dups = TRUE)
View(algae.seabirds)

##Save the combined data, unsummarized for future analyses (or to come back to)
write.csv(algae.seabirds, "../output/seabird-algaen15/seabird_algaen15_noiti_totaldata.csv")
algae.seabirds <- read.csv("../output/seabird-algaen15/seabird_algaen15_noiti_totaldata.csv", header = TRUE, strip.white = TRUE)
algae.seabirds$Distance_to_shore <- as.factor(algae.seabirds$Distance_to_shore)

##Supplementary correlation plots:

#we want data from all the seabirds and all the N15 at each distance from shore
#so we need to spread the data according to distance to shore & 15

library(plyr)
sum.n15<- ddply(algae.seabirds, c("site.name", "Transect", "Distance_to_shore"),summarise,
                   mean.N15 = mean(N15),
                  breeding_biomass_kgha_side = breeding_biomass_kgha_side
)

sum.n15 <- sum.n15 %>% spread(Distance_to_shore, mean.N15)

sum.n15 <- sum.n15 %>% rename(c("10" = "N.15_at_10m", 
                                    "20" = "N.15_at_20m", "30" = "N.15_at_30m" ,
                                    "40" = "N.15_at_40m"))
data.matrix <- as.data.frame(sum.n15[,3:7])
View(data.matrix)
library(corrplot)

cor.mtest(data.matrix)
correlation.matrix <- cor(data.matrix, use = "pairwise.complete.obs")
write.csv(correlation.matrix, "../output/n15_seabirds_corrmatrix_with_iti.csv") #change to "_no_iti" if removing iti


pdf(file = "../output/seabird-algaen15/corrplot_breedingbiomass_n15_no_iti.pdf") #change to "_no_iti" if removing iti

corrplot(cor(data.matrix, use = "pairwise.complete.obs"), type = "upper", 
         addCoef.col = NULL, addCoefasPercent = FALSE, tl.col = "black", title = "seabirds vs. N15")

dev.off()


##Run lmers 
library(lme4)
library(car)
library(jtools)
library(emmeans)

mod <- lmer(N15 ~ breeding_biomass_kgha_side * Distance_to_shore + (1|site.name), 
            data = algae.seabirds)
Anova(mod)

#Response: N15
#Chisq Df Pr(>Chisq)    
#breeding_biomass_kgha_side                    3.5491  1    0.05958 .  
#Distance_to_shore                            68.0129  3  1.137e-14 ***
#breeding_biomass_kgha_side:Distance_to_shore 24.2504  3  2.215e-05 ***


#Yes there is a significant interaction between breeding biomass and distance to shore. 

summary(mod)
vif(mod)
plot(mod)

plot_summs(mod)
summ(mod)
anova(mod)
Anova(mod)

#Can we tease apart what's going on?
library(multcompView)
library(emmeans)
marginal <- lsmeans(mod, ~Distance_to_shore | site.name)
cld(marginal, alpha = 0.05, Letters = letters, adjust = "sidak")

#Distance_to_shore lsmean    SE   df lower.CL upper.CL .group
#40                  5.40 0.518 4.82     3.40     7.40  a    
#30                  5.75 0.512 4.60     3.73     7.77  ab   
#20                  6.23 0.512 4.60     4.20     8.25   b   
#10                  7.25 0.513 4.66     5.23     9.26    c 


#Plot

pdf(file = "../output/seabird-algaen15/N15vsBreedBiomass.pdf")

algae.seabirds%>%
  group_by(Distance_to_shore, breeding_biomass_kgha_side)%>%
  summarize(mean_n15 = mean(N15),
            n_n15 = length(N15),
            se_n15 = sd(N15)/sqrt(n_n15))%>%
  ggplot(aes(x = breeding_biomass_kgha_side, y = mean_n15, color = Distance_to_shore, fill = Distance_to_shore))+
  geom_point( alpha = .5)+
  geom_errorbar(aes(ymin = (mean_n15-se_n15), ymax = (mean_n15+se_n15)), alpha = .5)+
  geom_line(alpha = .5)+
  theme_bw() 

dev.off()

algae.seabirds %>%
  #mutate(site.name = fct_relevel(site.name, "Rimatuu_Protected", "Reiono_Exposed", "Aie_Exposed", "Aie_Protected"))%>%
  ggplot(aes(x = breeding_biomass_kgha_side, y = N15, colour = site.name, fill = site.name, shape = site.name)) + 
  geom_point(size = 3, alpha= .6)+
  stat_ellipse(geom = "polygon", alpha = .05, level = .75)+
  #facet_wrap(~Distance_to_shore, scales = 'free')+
  theme_bw()+
  theme(panel.grid.major = element_blank(), # remove gridlines
        panel.grid.minor = element_blank(), #remove gridlines
        strip.background = element_blank()
  )


##Because of variability, combine into low medium high (following the fish vs seabird data)
algaen15.seabirds <-
  algae.seabirds%>%
  mutate(seabird_level = case_when(breeding_biomass_kgha_side<10 ~"low",
                                   breeding_biomass_kgha_side>10&breeding_biomass_kgha_side <200 ~"mid",
                                   breeding_biomass_kgha_side>200 ~"high"))%>%
  mutate(seabird_level = as.factor(seabird_level))%>%
  mutate(seabird_level = fct_relevel(seabird_level, "low", "mid", "high"))


##re-plot----

sum <- ddply(algaen15.seabirds, c("Distance_to_shore", "seabird_level"), summarise,
             mean_n15 = mean(N15), 
             n_n15 = length(N15), 
             se_n15 = sd(N15)/sqrt(n_n15))

pdf(file = "../output/seabird-algaen15/N15vsBreedBiomass_levels_10_200_plus.pdf")

ggplot(sum, aes(x = seabird_level, y = mean_n15, color = Distance_to_shore, fill = Distance_to_shore, group = Distance_to_shore))+
  geom_point(alpha = .5, size = 4)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  geom_errorbar(aes(ymin = (mean_n15-se_n15), ymax = (mean_n15+se_n15)), alpha = .5, width = 0.1)+
  geom_line(alpha = .5)+
  ylab(expression(italic(delta)^15*N))+
  xlab("\n Seabird Biomass Level") +
  labs(color='Distance \nFrom Shore (m)', fill = 'Distance \nFrom Shore (m)') +
  theme_classic() +
  theme(legend.position = c(0.15, 0.8)) +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5)) 

dev.off()

#What about stats

hist(algaen15.seabirds$N15)
plot(N15 ~ seabird_level, algaen15.seabirds)
with(algaen15.seabirds, lines(lowess(N15 ~ seabird_level)))

#Trial the standard mixed effects model with seabird level/distance interaction effect and site as a random factor (REML)
library(nlme)
mod.lme <- lme(N15 ~ seabird_level*Distance_to_shore + Exposure, random = ~1|site.name, algaen15.seabirds, method = "REML")
#Trial an optimized version (more complex)
mod.lme1 <- lme(N15 ~ seabird_level*Distance_to_shore + Exposure, random = ~ seabird_level | site.name, algaen15.seabirds, 
                         method = "REML",control = lmeControl(opt = "optim", method = "BFGS"))
#Check AIC values
anova(mod.lme, mod.lme1) 

library(sjPlot)
plot_grid(plot_model(mod.lme, type = "diag"))
plot(mod.lme) #model fit = quite good

anova(mod.lme)
intervals(mod.lme)
#Look at the pairwise comparisons using glht
emmeans(mod.lme, list(pairwise ~ seabird_level*Distance_to_shore), adjust = "fdr")
