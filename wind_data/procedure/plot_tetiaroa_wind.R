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


