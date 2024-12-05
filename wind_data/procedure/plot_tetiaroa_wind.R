## Code for wind data
library(ggplot2)
library(plyr)
library(dplyr)

##Pull in the windrose code
source("windrose.R")

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



#What about adding in different time scales of wind roses

#1. immediate Nov 6-14 
tet_weather_imm <- tet_weather %>% filter(AAAAMM == "2021_11")
tet_weather_imm <- tet_weather_imm %>% filter(JJ == c("6", "7", "8", "9", "10", "11", "12", "13"))
#tet_weather_imm$dataset <- "immediate"
immo.extract <- tet_weather_imm[,c(1,16,18)]
#immo.extract <- na.omit(immo.extract)
View(immo.extract)
p.imm <- plot.windrose(data = immo.extract, 
              spd = "FF",
              dir = "DD",
              spdseq = c(0.0,3,6,9,12))
p.imm <- p.imm + theme(axis.title.x = element_blank(), axis.text.x = element_blank())

pdf(file = "../output/1_windrose_tet_immediate.pdf")
p.imm
dev.off()
#2. 3 months Sep-Nov

tet_weather_3mo <- tet_weather %>% filter(AAAAMM == c("2021_9", "2021_10", "2021_11"))
#tet_weather_3mo$dataset <- "3_month"
threemo.extract <- tet_weather_3mo[,c(1,16,18)]
threemo.extract <- na.omit(threemo.extract)
View(threemo.extract)
#length(threemo.extract$FF)
range(threemo.extract$FF)

p.3mo <- plot.windrose(data = threemo.extract, 
                       spd = "FF",
                       dir = "DD",
                       spdseq = c(0,3,6,9,12))

p.3mo <- p.3mo + theme(axis.title.x = element_blank(), axis.text.x = element_blank())

pdf(file = "../output/2_windrose_tet_3month.pdf")
p.3mo
dev.off()
#3. 6 months June-Nov

#tet_weather$AAAAMM <- as.factor(tet_weather$AAAAMM)
#levels(tet_weather$AAAAMM)
tet_weather_6mo <- tet_weather %>% filter(AAAAMM %in% c("2021_6", "2021_7", "2021_8", 
                                                      "2021_9", "2021_10", "2021_11"))

View(tet_weather_6mo)
#tet_weather_6mo$dataset <- "6_month"
sixmo.extract <- tet_weather_6mo[,c(1,11,16,18)]
sixmo.extract <- na.omit(sixmo.extract)
View(sixmo.extract)
length(sixmo.extract$FF)

p.6mo <- plot.windrose(data = sixmo.extract, 
                       spd = "FF",
                       dir = "DD",
                       spdseq = c(0.0,3,6,9,12))

p.6mo <- p.6mo + theme(axis.title.x = element_blank(), axis.text.x = element_blank())

pdf(file = "../output/3_windrose_tet_6month.pdf")
p.6mo
dev.off()

#4. annual -2021

tet_weather_yr <- tet_weather %>% filter(AAAA == "2021")
View(tet_weather_yr)
#tet_weather_yr$dataset <- "annual"

mean(na.omit(tet_weather_yr$FF)) #5.528576
sd(na.omit(tet_weather_yr$FF))/sqrt(length(na.omit(tet_weather_yr$FF))) #0.028131
mean(na.omit(tet_weather_yr$DD)) #118.9599
sd(na.omit(tet_weather_yr$DD))/sqrt(length(na.omit(tet_weather_yr$DD))) #0.8329114
range(na.omit(tet_weather_yr$DD))
hist(na.omit(tet_weather_yr$DD))

tet_weather_yr_0.170 <- tet_weather_yr %>% filter(tet_weather_yr$DD < 171)
tet_weather_yr_171.359 <- tet_weather_yr %>% filter(tet_weather_yr$DD >170)

length(na.omit(tet_weather_yr_0.170$DD)) #7528
length(na.omit(tet_weather_yr_171.359$DD)) #1231 
median(na.omit(tet_weather_yr$DD))

#to get rid of NAs maybe need to extract just the FF and the DD
yr.extract <- tet_weather_yr[,c(1,16,18)]
#yr.extract <- na.omit(yr.extract)

p.yr <- plot.windrose(data = yr.extract, 
                       spd = "FF",
                       dir = "DD",
                       spdseq = c(0.0,3,6,9,12))
p.yr <- p.yr + 
  labs(y = "Number of Observations (n)") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

pdf(file = "../output/4_windrose_tet_annual.pdf")
p.yr
dev.off()


#Something is going on with this one - should only have ~150-200 observations but it has >700: 
tet_weather_nov <- tet_weather %>% filter(AAAAMM == "2021_11")
tet_weather_nov$AAAA
nov.extract <- tet_weather_nov[,c(1,16,18)]
View(nov.extract)
p.nov <- plot.windrose(data = nov.extract, 
                       spd = "FF",
                       dir = "DD",
                       spdseq = c(0.0,3,6,9,12))
p.nov <- p.nov + theme(axis.title.x = element_blank(), axis.text.x = element_blank())

pdf(file = "../output/5_windrose_tet_november.pdf")
p.nov
dev.off()


##How about pulling out the data we want for the table 
immo.extract$dataset <- "immediate"

sum.imm <- ddply(immo.extract, "dataset", summarise,
                       Speed_N = length(FF),
                       Speed_mean = mean(FF),
                       Speed_se = sd(FF)/sqrt(Speed_N),
                       Direction_N = length(DD),
                       Direction_mean = mean(DD),
                       Direction_median= median(DD), 
                       Direction_se = sd(DD)/sqrt(Direction_N)
)

sum.imm
easterly.immo <- immo.extract %>% filter(immo.extract$DD < 171)
westerly.immo <- immo.extract %>% filter(immo.extract$DD >170)
length(na.omit(easterly.immo$DD))
length(na.omit(westerly.immo$DD))
(17/(17+7))*100


threemo.extract$dataset <- "3month"

sum.three <- ddply(threemo.extract, "dataset", summarise,
                 Speed_N = length(FF),
                 Speed_mean = mean(FF),
                 Speed_se = sd(FF)/sqrt(Speed_N),
                 Direction_N = length(DD),
                 Direction_mean = mean(DD),
                 Direction_median= median(DD), 
                 Direction_se = sd(DD)/sqrt(Direction_N)
)

sum.three

easterly.three <- threemo.extract %>% filter(threemo.extract$DD < 171)
westerly.three <- threemo.extract %>% filter(threemo.extract$DD >170)
length(na.omit(easterly.three$DD))
length(na.omit(westerly.three$DD))
(623/(623+105))*100

sixmo.extract$dataset <- "6month"
sum.six <- ddply(sixmo.extract, "dataset", summarise,
                   Speed_N = length(FF),
                   Speed_mean = mean(FF),
                   Speed_se = sd(FF)/sqrt(Speed_N),
                   Direction_N = length(DD),
                   Direction_mean = mean(DD),
                   Direction_median= median(DD), 
                   Direction_se = sd(DD)/sqrt(Direction_N)
)

sum.six
easterly.six <- sixmo.extract %>% filter(sixmo.extract$DD < 171)
westerly.six <- sixmo.extract %>% filter(sixmo.extract$DD >170)
length(na.omit(easterly.six$DD))
length(na.omit(westerly.six$DD))
(3788/(3788+604))*100

yr.extract$dataset <- "annual"
yr.extract <- na.omit(yr.extract)
sum.yr <- ddply(yr.extract, "dataset", summarise,
                 Speed_N = length(FF),
                 Speed_mean = mean(FF),
                 Speed_se = sd(FF)/sqrt(Speed_N),
                 Direction_N = length(DD),
                 Direction_mean = mean(DD),
                 Direction_median= median(DD),
                 Direction_se = sd(DD)/sqrt(Direction_N)
)

sum.yr
easterly.yr <- yr.extract %>% filter(yr.extract$DD < 171)
westerly.yr <- yr.extract %>% filter(yr.extract$DD >170)
length(na.omit(easterly.yr$DD))
length(na.omit(westerly.yr$DD))
(7528/(7528+1231))*100



nov.extract$dataset <- "one-month"

sum.nov <- ddply(nov.extract, "dataset", summarise,
                 Speed_N = length(FF),
                 Speed_mean = mean(FF),
                 Speed_se = sd(FF)/sqrt(Speed_N),
                 Direction_N = length(DD),
                 Direction_mean = mean(DD),
                 Direction_median= median(DD), 
                 Direction_se = sd(DD)/sqrt(Direction_N)
)

sum.nov
easterly.nov <- nov.extract %>% filter(nov.extract$DD < 171)
westerly.nov <- nov.extract %>% filter(nov.extract$DD >170)
length(na.omit(easterly.nov$DD))
lengt(na.omit(westerly.nov$DD))
(610/(610+110))*100


#Get the modes rather than the medians

Mode <- function(x) {
  a <- table(x)
  as.numeric(names(a)[a == max(a)])
}

Mode(yr.extract$DD)
Mode(immo.extract$DD)
Mode(threemo.extract$DD)
Mode(sixmo.extract$DD)

write.csv(sum.speed.dir, "../output/data_summary_speed_direction.csv", row.names = F)

