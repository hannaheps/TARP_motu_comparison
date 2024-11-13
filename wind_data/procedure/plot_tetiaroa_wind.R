## Code for wind data

tet_weather <- read.csv("../input/Tet_Weather_2020_to_2022.csv")


tet_weather_nov <- tet_weather %>% filter(AAAAMM == "2021_11")

sum.FF <- ddply(tet_weather_nov, "JJ", summarise,
                 N = length(FF),
                 mean = mean(FF),
                 sd = sd(FF), 
                 se = sd/sqrt(N)
)

ggplot(sum.FF, aes(x = JJ, y = mean)) +
  geom_bar(stat = "identity") +
  ggtitle("FF in November at Tetiaroa") +
  #guides(fill=NULL) +
  xlab("Day (in November)") +
  ylab("FF") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))



sum.DD <- ddply(tet_weather_nov, "JJ", summarise,
                N = length(DD),
                mean = mean(DD),
                sd = sd(DD), 
                se = sd/sqrt(N)
)

ggplot(sum.DD, aes(x = JJ, y = mean)) +
  geom_bar(stat = "identity") +
  ggtitle("DD in November at Tetiaroa") +
  #guides(fill=NULL) +
  xlab("Day (in November)") +
  ylab("DD") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linetype = 1, linewidth = 0.5))
library(ggplot2)

plot.windrose()

plot.windrose(spd = tet_weather_nov$FF,
                    dir = tet_weather_nov$DD)
