## Code for wind data
library(ggplot2)
library(plyr)

#Read in data from Meteo France (open access online: https://meteo.data.gouv.fr/form)
tet_weather <- read.csv("../input/Tet_Weather_2020_to_2022.csv")

#Subset to november 2021 only
tet_weather_nov <- tet_weather %>% filter(AAAAMM == "2021_11")

#Explore data
sum.FF <- ddply(tet_weather_nov, "JJ", summarise,
                 N = length(FF),
                 mean = mean(FF),
                 sd = sd(FF), 
                 se = sd/sqrt(N)
)

pdf(file = "../output/windspeed_tet_november2021.pdf")
ggplot(sum.FF, aes(x = JJ, y = mean)) +
  geom_bar(stat = "identity") +
  ggtitle("FF in November at Tetiaroa") +
  #guides(fill=NULL) +
  xlab("Day (in November)") +
  ylab("FF") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))
dev.off()


sum.DD <- ddply(tet_weather_nov, "JJ", summarise,
                N = length(DD),
                mean = mean(DD),
                sd = sd(DD), 
                se = sd/sqrt(N)
)

pdf(file = "../output/winddirection_tet_november2021.pdf")
ggplot(sum.DD, aes(x = JJ, y = mean)) +
  geom_bar(stat = "identity") +
  ggtitle("DD in November at Tetiaroa") +
  #guides(fill=NULL) +
  xlab("Day (in November)") +
  ylab("DD") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))
dev.off()

##Use the windrose code
source(windrose.R)

plot.windrose(spd = tet_weather_nov$FF,
                    dir = tet_weather_nov$DD)

#Remove NAs in the speed (FF) and direction (DD) columns
tet_weather_nov_narm <- tet_weather_nov[!is.na(tet_weather_nov$FF) | !is.na(tet_weather_nov_narm$DD), ]

#build the windrose, using customised binning for speed
pdf(file = "../output/windrose_tet_november2021_custbin.pdf")
plot.windrose(data = tet_weather_nov_narm, 
              spd = "FF",
              dir = "DD",
              spdseq = c(0.0,3,6,9,12))
dev.off()


#Can make it ggplot compatible: 
p.wind <- plot.windrose(data = tet_weather_nov_narm, 
                        spd = "FF",
                        dir = "DD",
                        spdseq = c(0.0,3,6,9,12))

p.wind + theme(axis.title.x = element_blank(), axis.text.x = element_blank())

#How about try to facet wrap for the different time periods of each metric: 
#Algae isotopes - 7-11 November 2021
#Algae surveys - 7-11 November 2021
#Fish collections  (Aie - 17, 24, 29 30 July & 20 Jun, Reiono - 12 Nov, Rimatuu - 27 Nov, 10 & 13 Dec)  )
#Microbes - 7-11 November 2021
#Benthic - 6,8,10,14 November 2021
#Fish Surveys - 6,8,10,14 November 2021


tet_weather_transects<- tet_weather %>% filter(AAAAMM == "2021_11")
tet_weather_transects<- tet_weather_transects %>% filter(JJ == c("7", "8", "9", "10", "11"))

tet_weather_surveys <- tet_weather %>% filter(AAAAMM == "2021_11" )
tet_weather_surveys <- tet_weather_surveys %>% filter(JJ == c("6", "8", "10", "14"))

tet_weather_fishcoll <- tet_weather %>% filter(AAAAMM == c("2021_7", "2021_6", "2021_11", "2021_12"))
tet_weather_fishcoll$MMJJ <- as.factor(paste(tet_weather_fishcoll$MM, tet_weather_fishcoll$JJ, sep = "_"))
str(tet_weather_fishcoll$MMJJ)
tet_weather_fishcoll <- tet_weather_fishcoll %>% filter(tet_weather_fishcoll$MMJJ %in% c("6_20", "7_17", "7_24", "7_29", "7_30", "11_12", "11_27", "12_10", "12_13"))
tet_weather_fishcoll$MMJJ


#okay! Now we have our separate datasets, can we build some roses
#transects, which includes algae isotopes, algae surveys, microbes
p.wind.trans <- plot.windrose(data = tet_weather_transects, 
                        spd = "FF",
                        dir = "DD",
                        spdseq = c(0.0,3,6,9,12))
pdf(file = "../output/windrose_tet_november2021_transects_byday.pdf")
p.wind.trans + theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +facet_wrap(~JJ)
dev.off()

#benthic and fish surveys
p.wind.surv <- plot.windrose(data = tet_weather_surveys, 
                              spd = "FF",
                              dir = "DD",
                              spdseq = c(0.0,3,6,9,12))
pdf(file = "../output/windrose_tet_november2021_survey_byday.pdf")
p.wind.surv + theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +facet_wrap(~JJ)
dev.off()

#fish collection
p.wind.fishcoll <- plot.windrose(data = tet_weather_fishcoll, 
                              spd = "FF",
                              dir = "DD",
                              spdseq = c(0.0,3,6,9,12))

pdf(file = "../output/windrose_tet_november2021_fishcollection_bymonthday.pdf")
p.wind.fishcoll + theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + facet_wrap(~ MMJJ)
dev.off()

##Can we facet wrap the three things?
#First combine the dataframes and then facet_wrap by the dataset
tet_weather_transects<- tet_weather_transects %>% filter(JJ == c("7", "8", "9", "10", "11"))
tet_weather_transects$MMJJ <- as.factor(paste(tet_weather_transects$MM, tet_weather_transects$JJ, sep = "_"))
tet_weather_transects$dataset <- "transects"

tet_weather_surveys <- tet_weather_surveys %>% filter(JJ == c("6", "8", "10", "14"))
tet_weather_surveys$MMJJ <- as.factor(paste(tet_weather_surveys$MM, tet_weather_surveys$JJ, sep = "_"))
tet_weather_surveys$dataset <- "benthic_fish_surveys"

tet_weather_fishcoll$dataset <- "fish_collection"

#combine
tet_weather_datsets <- rbind(tet_weather_transects, tet_weather_surveys)
tet_weather_datsets <- rbind(tet_weather_datsets, tet_weather_fishcoll)

p.wind.all <- plot.windrose(data = tet_weather_datsets, 
                                 spd = "FF",
                                 dir = "DD",
                                 spdseq = c(0.0,3,6,9,12))

pdf(file = "../output/windrose_tet_all_datasets.pdf")
p.wind.all + theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + facet_wrap(~ dataset)
dev.off()

#How about a table with the summarised data
sum.speed.dir <- ddply(tet_weather_datsets, c("dataset", "MMJJ"), summarise,
                Speed_N = length(FF),
                Speed_mean = mean(FF),
                Speed_sd = sd(FF), 
                Speed_se = Speed_sd/sqrt(Speed_N),
                Direction_N = length(DD),
                Direction_mean = mean(DD),
                Direction_sd = sd(DD), 
                Direction_se = Direction_sd/sqrt(Direction_N)
)

sum.speed.dir

write.csv(sum.speed.dir, "../output/data_summary_speed_direction.csv", row.names = F)

