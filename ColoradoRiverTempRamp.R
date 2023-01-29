###########################################
## Code to create a plot for a ramping temp increase at Lees Ferry
###########################################
library(purrr)
library(tidyverse)
library(lubridate)
library(plyr)
library(dplyr)
library(ggplot2)
# data retrieval tool from USGS
library(dataRetrieval)
#Lees Ferry Temps
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2021-09-30")

# historical average for each 2 weeks timestep concatenated into a yearlong timeseries
temps <- average.yearly.temp(temp, "X_00010_00003", "Date")

# n is number of years I want the final time series to have
n = 13
# qr is the temp ramps I want to increase the average Lees Ferry temp by 
qr <- c(0, 1, 2.5, 5, 7.5)
# how many years I want each temp ramp to last
r <- c(5, 2, 2, 2, 2)

# manipulating temps to make a temp ramp
temps <- rep.avg.year(temps, n, change.in.temp = qr, years.at.temp = r)

# some code for plotting
arrows <- tibble(
  x1 = c("2005-01-07", "2007-01-07", "2009-01-07", "2011-01-07"),
  x2 = c("2005-01-07", "2007-01-07", "2009-01-07", "2011-01-07"),
  y1 = c(20, 20, 20, 20), 
  y2 = c(7, 7, 7, 7)
)
arrows$x1 <- as.Date(arrows$x1)
arrows$x2 <- as.Date(arrows$x2)

ggplot(data = temp_seq, mapping = aes(x = dts, y = Temperature ))+
  geom_line()+
  ylab("Water Temperature (°C)")+
  xlab(" ")+
  theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13))+
  scale_x_date(date_labels="%B", date_breaks  ="6 months")+
  annotate("segment", x = arrows$x1, y = arrows$y1, xend = arrows$x2, yend = arrows$y2, color = "red" )+
  annotate("text", x = arrows$x1[1]+7, y = 21, label = "+1°C", size = 5)+
  annotate("text", x = arrows$x1[2]+7, y = 21, label = "+2.5°C", size = 5)+
  annotate("text", x = arrows$x1[3]+7, y = 21, label = "+5°C", size = 5)+
  annotate("text", x = arrows$x1[4]+7, y = 21, label = "+7.5°C", size = 5 )





