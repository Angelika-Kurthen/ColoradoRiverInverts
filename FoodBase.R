library(purrr)
library(tidyverse)
library(lubridate)
library(plyr)
library(dplyr)
library(ggplot2)
# data retrieval tool from USGS
library(dataRetrieval)
#install.packages("devtools")
library(devtools)
#install_github(repo = "jmuehlbauer-usgs/R-packages", subdir = "foodbase")
library(foodbase)

# we will say that Lees Ferry = -2 km to 2 km (4k reach)
source("NZMS_1sp_Model.R")

discharge <- readNWISdv("09380000", "00060", "2007-10-01", "2023-05-16")
flow.magnitude <- TimestepDischarge(discharge, 85000)
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-16")
temps <- TimestepTemperature(temp, "Colorado River")
out <- NZMSmodel(flow.data = flow.magnitude$Discharge, temp.data = temps, disturbanceK = 10000, baselineK = 1000, Qmin = 0.2, extinct = 0, iteration = 10, peaklist = 0.27, peakeach = length(temps$Temperature))

means.list.NZMS <- mean.data.frame(out,burnin = 1, iteration= 1)
means.list.NZMS <- cbind(means.list.NZMS, temps$dts[1:length(means.list.NZMS$timesteps)])
means.list.NZMS$`temps$dts` <- as.Date(means.list.NZMS$`temps$dts`)
# plot abundance over time

abund.trends.NZMS <- ggplot(data = means.list.NZMS, aes(x = `temps$dts`,
                                                        y = mean.abund, group = 1)) +
  geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
                   ymax = mean.abund + 1.96 * se.abund),
               colour = 'transparent',
               alpha = .5,
               show.legend = FALSE) +
  geom_line(show.legend = FALSE, linewidth = 0.7) +
  geom_point()+
  geom_line(data = NZMS.samp, aes(x = Group.1, y = x), color = "red")+
  geom_line(data = flow.magnitude, aes(x = as.Date(dts), y = X_00060_00003), color = "blue") +
  geom_line(data = temps, aes(x = as.Date(dts), y = Temperature*1000), color = "green")+
  coord_cartesian(ylim = c(0,70000)) +
  ylab('New Zealand Mudsnail Abundance/Recruitment Limit') +
  xlab("")+
  theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13))+
  scale_x_date(date_labels="%B", date_breaks  ="6 months")



drift.data.total <- readDB(gear = "Drift", type = "Sample", updater = TRUE)

drift.LF <- drift.data.total[which(drift.data.total$RiverMile >= -3 & drift.data.total$RiverMile <= 3),]

#Baet.samp.LF <- sampspec(samp = drift.LF, species = "BAET", stats = T)

#View(Baet.samp.LF$Samples)
#m <- merge(Baet.samp.LF$Statistics, Baet.samp.LF$Samples, by = "BarcodeID", all = T)

NZMS.samp.LF <- sampspec(samp = drift.LF, species = "NZMS", stats = T)
NZMS.samp <- merge(NZMS.samp.LF$Statistics, NZMS.samp.LF$Samples, by = "BarcodeID", all = T)
NZMS.samp <- aggregate(NZMS.samp$CountTotal, list(NZMS.samp$Date), FUN = sum)
