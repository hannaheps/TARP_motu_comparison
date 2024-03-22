##Combine seabirds and algae diversity
library(tidyverse)
library(plyr)
library(lme4)
library(car)
library(jtools)
library(emmeans)
library(multcompView)
library(multcomp)

#Bring in seabird data
seabirds.all <- read.csv("../../seabird_data/outputs/seabird_dens_bio_focal_data.csv")
View(seabirds.all)

#clean-up the dataframe by removing the X column (column 1) accidentally added in the csv file.
#At the same time, change the name to maintain consistency in downstream code
seabirds <- seabirds.all[,-1]
View(seabirds)
#add a site.name (Motu_Exposure) category separated by an underscore.
seabirds$site.name <- paste(seabirds$Motu, seabirds$Exposure, sep = "_")

#Bring in the algae diversity data
algae.div <- read.csv("../../algae_diversity_surveys/output/algae_div_summary_all.csv")
algae.div$distance.along.transect <- as.factor(algae.div$distance.along.transect)

#Combine the algae and seabird data 
algae.div.seabirds <- merge(algae.div, seabirds, by = "site.name", all = TRUE, no.dups = TRUE)
View(algae.div.seabirds)
write.csv(algae.div.seabirds, "../output/seabird-algaediv/seabird_algaediv_noiti_totaldata.csv", row.names= F)


##trial model
mod.algaediv <- lmer(richness ~ breeding_biomass_kgha_side*distance.along.transect + (1|site.name), 
            data = algae.div.seabirds)
Anova(mod.algaediv)



##Because of variability, combine into low medium high (following the fish vs seabird data)
algae.div.seabirds.lvl <-
  algae.div.seabirds%>%
  mutate(seabird_level = case_when(breeding_biomass_kgha_side<10 ~"low",
                                   breeding_biomass_kgha_side>10&breeding_biomass_kgha_side <200 ~"mid",
                                   breeding_biomass_kgha_side>200 ~"high"))%>%
  mutate(seabird_level = as.factor(seabird_level))%>%
  mutate(seabird_level = fct_relevel(seabird_level, "low", "mid", "high"))

algae.div.seabirds.lvl$richness <- as.numeric(algae.div.seabirds.lvl$richness)
##re-plot----


algae.sum <- ddply(algae.div.seabirds.lvl, c("distance.along.transect", "seabird_level"),summarise, 
                    mean_richness = mean(richness), 
                    n_richness = length(richness),
                    se_richness = sd(richness)/sqrt(n_richness))
             
pdf(file = "../output/seabird-algaediv/RichnessvsSeabirdLevel.pdf")
ggplot(algae.sum, aes(x = seabird_level, y = mean_richness, color = distance.along.transect, 
              group = distance.along.transect))+
  geom_point(alpha = 0.5, size  = 4)+
  geom_errorbar(aes(ymin = (mean_richness-se_richness), ymax = (mean_richness+se_richness)), alpha = .5, width = 0.1)+
  geom_line(alpha = 0.5)+
  ylab("Mean Species Richness \n") +
  xlab("\n Seabird Biomass") +
  labs(color='Distance \nFrom Shore (m)') +
  scale_color_manual(labels = c("10", "40")
                     ,values = c("orangered", "purple3")) +
  theme_classic() +
  theme(legend.position = c(0.8, 0.82)) +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5)) 

dev.off()

##Let's model this

library(car)
library(nlme)

#View the data
ggplot(algae.div.seabirds.lvl, aes(y = richness, x = seabird_level, fill = distance.along.transect, color = distance.along.transect)) +
  geom_smooth(method = "lm") + 
  geom_point() + 
  theme_classic()

#Look at the distribution of data
hist(algae.div.seabirds.lvl$richness)
plot(richness ~ seabird_level, algae.div.seabirds.lvl)
with(algae.div.seabirds.lvl, lines(lowess(richness ~ seabird_level)))

#Trial the standard mixed effects model with seabird level/distance interaction effect and site as a random factor (REML)
mod.algaediv.lme <- lme(richness ~ seabird_level*distance.along.transect + Exposure, random = ~1|site.name, algae.div.seabirds.lvl, method = "REML")
#Trial an optimized version (more complex)
mod.algaediv.lme1 <- lme(richness ~ seabird_level*distance.along.transect + Exposure, random = ~ seabird_level | site.name, algae.div.seabirds.lvl, 
                         method = "REML",control = lmeControl(opt = "optim", method = "BFGS"))
#Check AIC values
anova(mod.algaediv.lme, mod.algaediv.lme1) #Looks like the simpler model is better
library(sjPlot)
plot_grid(plot_model(mod.algaediv.lme, type = "diag"))
plot(mod.algaediv.lme)

anova(mod.algaediv.lme)
intervals(mod.algaediv.lme)
#Look at the pairwise comparisons using glht
emmeans(mod.algaediv.lme, list(pairwise ~ seabird_level*distance.along.transect), adjust = "fdr")

##Can we bring in the n15 data so that we can also look at difference in N15
#Is it possible to bring in 10 and 40m algae?

n15.seabirds <- read.csv("../output/seabird-algaen15/seabird_algaen15_noiti_totaldata.csv", header = TRUE, strip.white = TRUE)
n15.seabirds <- n15.seabirds[,-1]

n15.seabirds$Distance_to_shore <- as.factor(n15.seabirds$Distance_to_shore)
#filter to distance from shore
n15.seabirds <- n15.seabirds %>% filter(Distance_to_shore != "20")
n15.seabirds <- n15.seabirds[n15.seabirds$Distance_to_shore != '20',]
str(n15.seabirds$Distance_to_shore)
#change name of distance_to_shore to the same as in algal div
names(n15.seabirds)[names(n15.seabirds) == "Distance_to_shore"] <- "distance.along.transect"
View(n15.seabirds)

n15.only <- n15.seabirds[,c(1,10,16)]

#Merge the two data frames by two columns
algae.div.seabirds.lvl.n15 <- merge(n15.only, algae.div.seabirds.lvl, by = c("site.name", "distance.along.transect"), all = TRUE, no.dups = TRUE)
View(algae.div.seabirds.lvl.n15)


str(algae.div.seabirds.lvl.n15$distance.along.transect)
algae.div.seabirds.lvl.n15 <- filter(algae.div.seabirds.lvl.n15, distance.along.transect == c("10", "40"))

#remove rows with NAs
algae.div.seabirds.lvl.n15 <- algae.div.seabirds.lvl.n15[!is.na(algae.div.seabirds.lvl.n15$N15), ]

#Now we can do the same by n15
hist(algae.div.seabirds.lvl.n15$richness)
plot(richness ~ N15, algae.div.seabirds.lvl.n15)
with(algae.div.seabirds.lvl.n15, lines(lowess(richness ~ N15)))

#Trial the standard mixed effects model with seabird level/distance interaction effect and site as a random factor (REML)
mod.algaediv.lme <- lme(richness ~ N15*distance.along.transect + Exposure, random = ~1|site.name, algae.div.seabirds.lvl.n15, method = "REML")
#Trial an optimized version (more complex)
mod.algaediv.lme1 <- lme(richness ~ N15*distance.along.transect + Exposure, random = ~ seabird_level | site.name, algae.div.seabirds.lvl.n15, 
                         method = "REML",control = lmeControl(opt = "optim", method = "BFGS"))
#Check AIC values
anova(mod.algaediv.lme, mod.algaediv.lme1) #Looks like the simpler model is better
library(sjPlot)
plot_grid(plot_model(mod.algaediv.lme, type = "diag"))
plot(mod.algaediv.lme)

anova(mod.algaediv.lme)
intervals(mod.algaediv.lme)
#Look at the pairwise comparisons using emmeans
emmeans(mod.algaediv.lme, list(pairwise ~ distance.along.transect), adjust = "fdr")


pdf(file = "../output/seabird-algaediv/RichnessvsN15.pdf")

ggplot(algae.div.seabirds.lvl.n15, aes(x = N15, y = richness, color = distance.along.transect, 
                      group = distance.along.transect))+
  geom_point(stat = "summary", fun = "mean", alpha = 0.5, size  = 3) +
  geom_smooth(method = lm, alpha = 0.2, aes(group = distance.along.transect))+
  ylab("Mean Species Richness \n") +
  xlab(expression(italic(delta)^15*N))+
  labs(color='Distance \nFrom Shore (m)') +
  scale_color_manual(labels = c("10", "40")
                     ,values = c("orangered", "purple3")) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.82)) +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5)) 


dev.off()
##Not perfect because this removes the points at 40m from Aie exposed where there were no algal collections



##How do we look at whether one taxon is driving the relationship? 
#I think I need to make the data long version:
#colnames(algae.div.seabirds.lvl.n15)
#algae.div.long <- algae.div.seabirds.lvl.n15 %>% gather(key=species.names,value=abundance, colnames(algae.div.seabirds[,c(4:34)]))
#head(algae.div.long)

#Okay that works BUT need to set up the data for relative abundances

#I need relative abundances first? Or at least rowsums
relabund <- algae.div.seabirds.lvl.n15
relabund$sums <- rowSums(relabund[,5:35])

relabund.long <- relabund %>% gather(key = species.names, value = abundance, colnames(relabund[, 5:35]))

#Then make a relative abundance column
relabund.long$rel.abund <- relabund.long$abundance/relabund.long$sums

#Now we have a relative abundance column and we can split into species and genus names!

#Can I split the species into species names and genus names?
#1. make a new duplicate column 

relabund.long$species.name.full <- relabund.long$species.names
relabund.long <- separate(data = relabund.long, col = species.name.full, into = c("genus", "species"), sep = "\\.")
head(relabund.long)

##Can I plot against seabird_level and color by genus?

algae.sum <- ddply(relabund.long, c("seabird_level", "genus"),summarise, 
                   mean_abund = mean(rel.abund),
                   n_abund = length(rel.abund),
                   se_abund = sd(rel.abund)/sqrt(n_abund))

#Remove any genera that have 0 mean
algae.sum <- algae.sum[algae.sum$mean_abund != 0, ]

ggplot(algae.sum, aes(x = seabird_level, y = mean_abund, color = genus))+
  geom_point( alpha = .5)+
  geom_errorbar(aes(ymin = (mean_abund-se_abund), ymax = (mean_abund+se_abund)), alpha = .5)+
  geom_line(alpha = .5) +
  ylab("Mean Abundance \n") +
  xlab("\n Seabird Biomass") +
  labs(color='Distance \nFrom Shore (m)') +
  theme_classic() 

pdf(file = "../output/seabird-algaediv/algal_genus_abundance.pdf")
ggplot(algae.sum, aes(x = genus, y = mean_abund)) +
  geom_point(position = position_dodge(width = 1), alpha = 0.5, aes(color = algae.sum$seabird_level)) +
  geom_errorbar(position = position_dodge(width = 1), alpha = 0.5, aes(ymin = (mean_abund-se_abund), ymax = (mean_abund + se_abund), color = algae.sum$seabird_level)) +
  ylab("Mean Relative Abundance \n") +
  xlab("\n Macroalgal Genus") +
  labs(color = 'Seabird Biomass', fill = 'Seabird Biomass') +
  #scale_color_manual(labels = c("low", "mid", "high")
                     #,values = c("blue3", "gold", "red3")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()


#Can we do any tests? 


##What about NMDS like what Casey has for fish

library(vegan)
nmds.trial <-metaMDS(relabund[, 5:35], distance = "bray", trymax=200, k=2, autotransform = FALSE)
stressplot(nmds.trial)
plot(nmds.trial)
scores(nmds.trial, display="species")
nmds.scores = as.data.frame(scores(nmds.trial)$sites)

nmds.scores$seabird_level <- relabund$seabird_level
nmds.scores$N15 <- relabund$N15
nmds.scores$site.name <- relabund$site.name
nmds.scores$distance.along.transect <- relabund$distance.along.transect

hull_nmds <- nmds.scores %>%
  group_by(seabird_level) %>%
  slice(chull(NMDS1, NMDS2))

pdf(file = "../output/seabird-algaediv/nmds_seabirdlevel.pdf")
ggplot(nmds.scores, aes(x = NMDS1, y= NMDS2)) + 
  geom_point(aes(colour = nmds.scores$seabird_level, shape = nmds.scores$distance.along.transect), size = 3) +
  geom_polygon(data = hull_nmds, aes(x=NMDS1,y=NMDS2,fill=seabird_level),alpha=0.30) +
  labs(color = 'Seabird Biomass', fill = 'Seabird Biomass', shape = 'Distance from Shore') +
  theme_classic() +
  theme(legend.position = c(0.8, 0.75)) +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5)) 
dev.off()

#how about a permanova?
relabund.nozero <- relabund[, colSums(relabund[,5:35] !=0) >0]
colSums(relabund.nozero[,5:32]) #doesn't remove all the 0s, no idea why
#index: 31,28, 16
relabund.nozero <- relabund.nozero[, -c(16,28,31)]
colSums(relabund.nozero[,5:29]) #SO ANNOYING WTF. but it works now

bc <- vegdist(relabund.nozero[,5:29], method = "bray")
bc.mat <- as.matrix(bc)
adonis2(bc.mat ~ seabird_level*distance.along.transect + Exposure, data = relabund.nozero, strata = relabund.nozero$site.name)
##Weird output, all p values are 0.001
#Can't actually use the strata = function because it doesn't treat it correctly

adonis2(bc.mat ~ seabird_level *distance.along.transect, data = relabund.nozero)
adonis2(bc.mat ~ N15*distance.along.transect, data = relabund.nozero)

#PERMDISP##
disp <- betadisper(bc, relabund.nozero$site.name, type = "centroid")
par(pty = 's')
boxplot(disp) 
disp.test <- permutest(disp)
disp.test

##Add N15 for visualisation??
envfit <- envfit(nmds.trial, relabund$N15, strata = relabund$site.name)
n15.arrow <-  as.data.frame(scores(envfit, display = "vectors"))
n15.arrow <- cbind(n15.arrow, Species = rownames(n15.arrow))

pdf(file = "../output/seabird-algaediv/nmds_seabirdlevel_envfit.pdf")
ggplot(nmds.scores, aes(x = NMDS1, y= NMDS2)) + 
  geom_point(aes(colour = nmds.scores$seabird_level, shape = nmds.scores$distance.along.transect), size = 3) +
  geom_polygon(data = hull_nmds, aes(x=NMDS1,y=NMDS2,fill=seabird_level),alpha=0.30) +
  geom_segment(data = n15.arrow, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), color = "black")+
  geom_text(data = n15.arrow, aes(x = NMDS1, y = NMDS2 - 0.1, label = "N15"), size = 4) +
  labs(color = 'Seabird Biomass', fill = 'Seabird Biomass', shape = 'Distance from Shore') +
  theme_classic() +
  #theme(legend.position = c(0.8, 0.75)) +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5)) 

ggsave("../output/seabird-algaediv/nmds_seabirdlevel_envfit.pdf", plot = last_plot())
dev.off()


#Very cool
#Can we plot a small stacked barplot that's averaged by seabird_level and genus?
library(viridis)
library(RColorBrewer)
pdf(file = "../output/seabird-algaediv/barplot_by_seabirdlevel.pdf")
ggplot(relabund.long, aes(fill=genus, y=rel.abund, x=seabird_level)) + 
  geom_bar(position="fill", stat="identity") +
  #scale_fill_viridis(discrete = T) +
  ylab("Relative Abundance\n") +
  xlab("\nSeabird Biomass Level") +
  labs(fill = "Genus") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


#Finally, what are the relative abundances of turbs, lobophoras,and chlorodesmis?

filter(algae.sum, genus == "Turbinaria")
filter(algae.sum, genus == "Lobophora")
filter(algae.sum, genus == "Chlorodesmis")
