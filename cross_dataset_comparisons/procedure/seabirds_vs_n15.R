##Combine the Seabird data vs. N15

##Bring in seabird data
seabirds.all <- read.csv("../../seabird_data/outputs/seabird_dens_bio_focal_data.csv")
View(seabirds.all)

#The following code under the #s will allow for the inclusion of iti birds, otherwise skip to line 32:
##########
#remove the density and biomass by side so that we can replace with the values that include iti
#seabirds.all <- seabirds.all[, -c(1,6,10,14,18)] 
#seabirds.all$site.name <- paste(seabirds.all$Motu, seabirds.all$Exposure, sep = "_") #create a site.name value for merging

#bring in the iti density data and combine with all other data
#seabirds.iti <- read.csv("../../seabird_data/outputs/seabird_dens_bio_by_side_iti_rimatuu_exposed_combined.csv")
#seabirds.iti$site.name <- paste(seabirds.iti$Motu, seabirds.iti$Exposure, sep = "_")

#seabirds <- merge(seabirds.iti, seabirds.all, by = "site.name")

#reorder the column names
#col_order <- c("site.name", "X", "Motu.x", "Exposure.x", "Motu.y", "Exposure.y", "breeding_biomass_kgha_100m", "breeding_biomass_kgha_200m",
#               "breeding_biomass_kgha_side", "breeding_biomass_kgha_motu", "breeding_density_ha_100m", "breeding_density_ha_200m",
#               "breeding_density_ha_side", "breeding_density_ha_motu", "adult_biomass_kgha_100m", "adult_biomass_kgha_200m",
#               "adult_biomass_kgha_side", "adult_biomass_kgha_motu",  "adult_density_ha_100m", "adult_density_ha_200m",
 #              "adult_density_ha_side",  "adult_density_ha_motu")  
#seabirds <- seabirds[, col_order]
#colnames(seabirds)

#seabirds <- seabirds[, -c(2:6)] #remove the non numerical columns aside from site.name needed for merging

#write.csv(seabirds, "../output/seabirds.combined.csv")

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
#data.matrix <- as.data.frame(algae.seabirds[ , 2:21])
data.matrix <- as.data.frame(algae.seabirds[,2:16])
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

sum.n15 <- sum.n15 %>% rename(c("N.15_at_10m" = "10", 
                                    "N.15_at_20m" = "20", "N.15_at_30m" = "30",
                                    "N.15_at_40m" = "40"))
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

pdf(file = "../output/seabird-algaen15/N15vsBreedBiomass_levels_10_200_plus.pdf")

algaen15.seabirds%>%
  group_by(Distance_to_shore, seabird_level)%>%
  summarize(mean_n15 = mean(N15),
            n_n15 = length(N15),
            se_n15 = sd(N15)/sqrt(n_n15))%>%
  ggplot(aes(x = seabird_level, y = mean_n15, color = Distance_to_shore, fill = Distance_to_shore, group = Distance_to_shore))+
  geom_point(alpha = .5)+
  geom_errorbar(aes(ymin = (mean_n15-se_n15), ymax = (mean_n15+se_n15)), alpha = .5, width = 0)+
  geom_line(alpha = .5)+
  theme_bw() 

dev.off()

#What about stats

mod.sb.level <- lmer(N15 ~ seabird_level*Distance_to_shore + (1|site.name), 
            data = algaen15.seabirds)
Anova(mod.sb.level)
summary(mod.sb.level)
vif(mod.sb.level)
plot(mod.sb.level)
plot_summs(mod.sb.level)
summ(mod.sb.level)
anova(mod.sb.level)

#At 10m only: 
algaen15.seabirds.10m <- filter(algaen15.seabirds, Distance_to_shore == "10")
mod.10m.level <- lmer(N15 ~ seabird_level + (1|site.name), data = algaen15.seabirds.10m)
Anova(mod.10m.level)
summary(mod.10m.level)
anova(mod.10m.level)


algaen15.seabirds.20m <- filter(algaen15.seabirds, Distance_to_shore == "20")
mod.20m.level <- lmer(N15 ~ seabird_level + (1|site.name), data = algaen15.seabirds.20m)
Anova(mod.20m.level)
summary(mod.20m.level)
anova(mod.20m.level)

algaen15.seabirds.30m <- filter(algaen15.seabirds, Distance_to_shore == "30")
mod.30m.level <- lmer(N15 ~ seabird_level + (1|site.name), data = algaen15.seabirds.30m)
Anova(mod.30m.level)

algaen15.seabirds.40m <- filter(algaen15.seabirds, Distance_to_shore == "40")
mod.40m.level <- lmer(N15 ~ seabird_level + (1|site.name), data = algaen15.seabirds.40m)
Anova(mod.40m.level)

#Can we tease apart what's going on?
library(multcompView)
library(multcomp)
library(emmeans)
marginal.level <- lsmeans(mod.sb.level, ~Distance_to_shore*seabird_level)
cld(marginal.level, alpha = 0.05, Letters = letters, adjust = "sidak")




