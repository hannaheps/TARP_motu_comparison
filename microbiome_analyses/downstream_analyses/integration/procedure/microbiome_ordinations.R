##Microbiome Ordinations
library(plyr)
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(vegan)

#pull in the two phyloseq objects that to grab the beta diversity metrics
physeq.coral.r <- readRDS("../../alpha-div/output/physeq-coral-nov21-rarefied.RDS")
physeq.water.r <- readRDS("../../alpha-div/output/physeq-water-nov21-rarefied.RDS")

#Export the data
coral.micro.data <- as(sample_data(physeq.coral.r), "data.frame")
water.micro.data <- as(sample_data(physeq.water.r), "data.frame")

coral.micro.data$site <- as.factor(coral.micro.data$site)
coral.micro.data$site.name <- paste(coral.micro.data$motu, coral.micro.data$site)

#recode the site.name category so we can combine with seabird data

coral.micro.data$site.name[coral.micro.data$site.name == "aie 1"] <- "Aie_Protected"
coral.micro.data$site.name[coral.micro.data$site.name == "aie 2"] <- "Aie_Exposed"
coral.micro.data$site.name[coral.micro.data$site.name == "reiono 1"] <- "Reiono_Protected"
coral.micro.data$site.name[coral.micro.data$site.name == "reiono 2"] <- "Reiono_Exposed"
coral.micro.data$site.name[coral.micro.data$site.name == "rimatuu 1"] <- "Rimatuu_Protected"
coral.micro.data$site.name[coral.micro.data$site.name == "rimatuu 2"] <- "Rimatuu_Exposed"

coral.micro.data$site.name <- as.factor(coral.micro.data$site.name)


water.micro.data$site <- as.factor(water.micro.data$site)
water.micro.data$site.name <- paste(water.micro.data$motu, water.micro.data$site)

water.micro.data$site.name[water.micro.data$site.name == "aie 1"] <- "Aie_Protected"
water.micro.data$site.name[water.micro.data$site.name == "aie 2"] <- "Aie_Exposed"
water.micro.data$site.name[water.micro.data$site.name == "reiono 1"] <- "Reiono_Protected"
water.micro.data$site.name[water.micro.data$site.name == "reiono 2"] <- "Reiono_Exposed"
water.micro.data$site.name[water.micro.data$site.name == "rimatuu 1"] <- "Rimatuu_Protected"
water.micro.data$site.name[water.micro.data$site.name == "rimatuu 2"] <- "Rimatuu_Exposed"

water.micro.data$site.name <- as.factor(water.micro.data$site.name)

#distance along transect needs to be made into a factor
coral.micro.data$distance.along.transect <- as.factor(coral.micro.data$distance.along.transect)
water.micro.data$distance.along.transect <- as.factor(water.micro.data$distance.along.transect)
#algae.N15 is there !


#need to add seabird level. 
seabirds <- read.csv("../../../../cross_dataset_comparisons/output/n15_seabirds_combined_no_iti.csv", strip.white = T, header = T)
coral.microbes.seabirds <- merge(coral.micro.data, seabirds, by = "site.name", all = TRUE, no.dups = TRUE)
head(coral.microbes.seabirds)
water.microbes.seabirds <- merge(water.micro.data, seabirds, by = "site.name", all = TRUE, no.dups = TRUE)
head(water.microbes.seabirds)


coral.microbes.seabirds <-
  coral.microbes.seabirds%>%
  mutate(seabird_level = case_when(breeding_biomass_kgha_side<10 ~"low",
                                   breeding_biomass_kgha_side>10&breeding_biomass_kgha_side <200 ~"mid",
                                   breeding_biomass_kgha_side>200 ~"high"))%>%
  mutate(seabird_level = as.factor(seabird_level))%>%
  mutate(seabird_level = fct_relevel(seabird_level, "low", "mid", "high"))
str(coral.microbes.seabirds$seabird_level)

water.microbes.seabirds <-
  water.microbes.seabirds%>%
  mutate(seabird_level = case_when(breeding_biomass_kgha_side<10 ~"low",
                                   breeding_biomass_kgha_side>10&breeding_biomass_kgha_side <200 ~"mid",
                                   breeding_biomass_kgha_side>200 ~"high"))%>%
  mutate(seabird_level = as.factor(seabird_level))%>%
  mutate(seabird_level = fct_relevel(seabird_level, "low", "mid", "high"))
str(water.microbes.seabirds$seabird_level)






#okay we have the right sample data and we have the right phyloseq object.
#Let's pull the ordination stuff from our other script. 

nmds.coral <- ordinate(physeq.coral.r, "NMDS", "bray", trymax = 500, k = 2, autotransform = FALSE)
nmds.coral
#2 axes stress = 0.1239, a little weak but still below 0.2

plot(nmds.coral)
scores(nmds.coral, display="sites")

#look at plot and centroids with ordfit
plot(nmds.coral)

ord.fit.level<-envfit(nmds.coral~seabird_level + island.side, data = coral.microbes.seabirds, na.rm=TRUE)
ord.fit.level 
plot(ord.fit.level)

plot(nmds.coral)
ord.fit.algae<-envfit(nmds.coral~algae.N15 + island.side, data = coral.microbes.seabirds, na.rm=TRUE)
ord.fit.algae 
plot(ord.fit.algae)

#rotate by algae n15 so can more easily interpret----
nmds.rotate<-MDSrotate(nmds.coral, coral.microbes.seabirds$algae.N15, na.rm = TRUE)
nmds.rotate
plot(nmds.rotate)

ord.fit<-envfit(nmds.rotate~algae.N15, data = coral.microbes.seabirds, na.rm=TRUE)
ord.fit #p and r^2 are the same, just rotated so all with NMDS1 now.


#Let's make it pretty
#extract scores (non-rotated)
nmds.coral.scores <- as.data.frame(scores(nmds.coral)$sites)

#combine with other relevant data
nmds.coral.scores$seabird_level <- coral.microbes.seabirds$seabird_level
nmds.coral.scores$N15 <- coral.microbes.seabirds$algae.N15
nmds.coral.scores$site.name <- coral.microbes.seabirds$site.name
nmds.coral.scores$distance.along.transect <- coral.microbes.seabirds$distance.along.transect

#NMDS plots - rotated:----
#get hull data - seabird level:
sb.low.ro <- nmds.coral.scores[nmds.coral.scores$seabird_level == "low", ][chull(nmds.coral.scores[nmds.coral.scores$seabird_level == 
                                                                                      "low", c("NMDS1", "NMDS2")]), ]

sb.mid.ro <- nmds.coral.scores[nmds.coral.scores$seabird_level == "mid", ][chull(nmds.coral.scores[nmds.coral.scores$seabird_level == 
                                                                                      "mid", c("NMDS1", "NMDS2")]), ]

sb.high.ro <- nmds.coral.scores[nmds.coral.scores$seabird_level == "high", ][chull(nmds.coral.scores[nmds.coral.scores$seabird_level == 
                                                                                        "high", c("NMDS1", "NMDS2")]), ]

hull.data.ro <- rbind(sb.low.ro, sb.mid.ro, sb.high.ro)  #combine
hull.data.ro


##get arrows:
ord.fit.sb.ro<-envfit(nmds.coral ~ algae.N15, data = coral.microbes.seabirds, na.rm=TRUE)
ord.fit.sb.ro

sb.fit.scores.ro <- as.data.frame(scores(ord.fit.sb.ro, display = "vectors"))
sb.fit.scores.ro <- cbind(sb.fit.scores.ro, variable = rownames(sb.fit.scores.ro))
sb.fit.scores.ro<-
  sb.fit.scores.ro %>%
  mutate(variable2 = c("algae n15"))


coral.nmds.plot<-
  ggplot() + 
  #geom_text(data=nmds.scores,aes(x=NMDS1,y=NMDS2), alpha = .7, size = 3, hjust = .1) + 
  geom_point(data=nmds.coral.scores,aes(x=NMDS1,y=NMDS2,colour=seabird_level, fill = seabird_level, pch = seabird_level), stat="identity", size=5, alpha = .9) +
  geom_polygon(data=hull.data.ro,aes(x=NMDS1,y=NMDS2,fill=seabird_level,group=seabird_level),alpha=0.30) + # add the convex hulls
  #  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = sb.fit.scores.ro,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text(data = sb.fit.scores.ro, aes(x = NMDS1, y = NMDS2, label = sb.fit.scores.ro$variable2),
            size = 3, vjust = "outward")+
  scale_fill_manual(values = c( "low" = "#CD1913", "mid" = "#F2BB05", "high" ="#2F9D3E"))+ #low-mid-high = red-yellow-green
  scale_colour_manual(values = c( "low" = "#CD1913", "mid" = "#F2BB05", "high" ="#2F9D3E"))+#low-mid-high = red-yellow-green
  labs(colour = "seabird biomass", fill = "seabird biomass", shape = "seabird biomass")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), # remove gridlines
        panel.grid.minor = element_blank(), #remove gridlines
        strip.background = element_blank(),
        legend.position = c(.89, .85)) 

coral.nmds.plot

#save plot
ggsave(coral.nmds.plot, file = "../output/coral_micro_nmds_plot_norotate.jpg",
       width = 7, height = 5)





##Alright water now
#okay we have the right sample data and we have the right phyloseq object.
#Let's pull the ordination stuff from our other script. 

nmds.water <- ordinate(physeq.water.r, "NMDS", "bray", trymax = 500, k = 2, autotransform = FALSE)
nmds.water
#2 axes stress = 0.034, looks good

plot(nmds.water)
scores(nmds.coral, display="sites")

#look at plot and centroids with ordfit
plot(nmds.water)

ord.fit.level<-envfit(nmds.water~seabird_level + island.side, data = water.microbes.seabirds, na.rm=TRUE)
ord.fit.level #sig
plot(ord.fit.level)

plot(nmds.water)
ord.fit.water<-envfit(nmds.water~algae.N15 + island.side, data = water.microbes.seabirds, na.rm=TRUE)
ord.fit.water
plot(ord.fit.water)

#rotate by algae n15 so can more easily interpret----
nmds.rotate<-MDSrotate(nmds.water, water.microbes.seabirds$algae.N15, na.rm = TRUE)
nmds.rotate
plot(nmds.rotate)

ord.fit<-envfit(nmds.rotate~algae.N15, data = water.microbes.seabirds, na.rm=TRUE)
ord.fit


#Let's make it pretty
#extract scores (non-rotated)
nmds.water.scores <- as.data.frame(scores(nmds.water)$sites)

#combine with other relevant data
nmds.water.scores$seabird_level <- water.microbes.seabirds$seabird_level
nmds.water.scores$N15 <- water.microbes.seabirds$algae.N15
nmds.water.scores$site.name <- water.microbes.seabirds$site.name
nmds.water.scores$distance.along.transect <- water.microbes.seabirds$distance.along.transect

#NMDS plots - rotated:----
#get hull data - seabird level:
sb.low.ro <- nmds.water.scores[nmds.water.scores$seabird_level == "low", ][chull(nmds.water.scores[nmds.water.scores$seabird_level == 
                                                                                                     "low", c("NMDS1", "NMDS2")]), ]

sb.mid.ro <- nmds.water.scores[nmds.water.scores$seabird_level == "mid", ][chull(nmds.water.scores[nmds.water.scores$seabird_level == 
                                                                                                     "mid", c("NMDS1", "NMDS2")]), ]

sb.high.ro <- nmds.water.scores[nmds.water.scores$seabird_level == "high", ][chull(nmds.water.scores[nmds.water.scores$seabird_level == 
                                                                                                       "high", c("NMDS1", "NMDS2")]), ]

hull.data.ro <- rbind(sb.low.ro, sb.mid.ro, sb.high.ro)  #combine
hull.data.ro


##get arrows:
ord.fit.sb.ro<-envfit(nmds.rotate ~ algae.N15, data = water.microbes.seabirds, na.rm=TRUE)
ord.fit.sb.ro

sb.fit.scores.ro <- as.data.frame(scores(ord.fit.sb.ro, display = "vectors"))
sb.fit.scores.ro <- cbind(sb.fit.scores.ro, variable = rownames(sb.fit.scores.ro))
sb.fit.scores.ro<-
  sb.fit.scores.ro %>%
  mutate(variable2 = c("algae n15"))


water.nmds.plot<-
  ggplot() + 
  #geom_text(data=nmds.scores,aes(x=NMDS1,y=NMDS2), alpha = .7, size = 3, hjust = .1) + 
  geom_point(data=nmds.water.scores,aes(x=NMDS1,y=NMDS2,colour=seabird_level, fill = seabird_level, pch = seabird_level), shape = water.microbes.seabirds$site.name, stat="identity", size=5, alpha = .9) +
  geom_polygon(data=hull.data.ro,aes(x=NMDS1,y=NMDS2,fill=seabird_level,group=seabird_level),alpha=0.30) + # add the convex hulls
  #  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = sb.fit.scores.ro,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text(data = sb.fit.scores.ro, aes(x = NMDS1, y = NMDS2, label = sb.fit.scores.ro$variable2),
            size = 3, vjust = "outward")+
  scale_fill_manual(values = c( "low" = "#CD1913", "mid" = "#F2BB05", "high" ="#2F9D3E"))+ #low-mid-high = red-yellow-green
  scale_colour_manual(values = c( "low" = "#CD1913", "mid" = "#F2BB05", "high" ="#2F9D3E"))+#low-mid-high = red-yellow-green
  labs(colour = "seabird biomass", fill = "seabird biomass", shape = "seabird biomass")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), # remove gridlines
        panel.grid.minor = element_blank(), #remove gridlines
        strip.background = element_blank(),
        legend.position = c(.89, .85)) 

water.nmds.plot

#save plot
ggsave(water.nmds.plot, file = "../output/water_micro_nmds_plot_rotate.jpg",
       width = 7, height = 5)

#Weirdddddd - probably is by site. 




