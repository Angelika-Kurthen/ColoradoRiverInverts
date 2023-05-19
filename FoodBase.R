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

discharge <- readNWISdv("09380000", "00060", "2007-10-01", "2023-05-01")
flow.magnitude <- TimestepDischarge(discharge, 85000)
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
temps <- TimestepTemperature(temp)
out <- NZMSmodel(flow.data = flow.magnitude$Discharge, temp.data = temps, disturbanceK = 500, baselineK = 1.25, Qmin = 0.19, extinct = 0, iteration = 100, peaklist = 0.27, peakeach = length(temps$Temperature))


adults<-as.data.frame(cbind(as.Date(temps$dts), out[1:length(temps$dts),2:3,1]))
colnames(adults) <- c("Time","Adult")
adults$Time <- as.Date(adults$Time, origin = "1970-01-01")

means.list.NZMS <- mean.data.frame(out,burnin = 1, iteration= 1)
means.list.NZMS <- cbind(means.list.NZMS, temps$dts[1:length(means.list.NZMS$timesteps)])
means.list.NZMS$`temps$dts` <- as.Date(means.list.NZMS$`temps$dts`)
# plot abundance over time

abund.trends.NZMS <- ggplot(data = means.list.NZMS, aes(x = `temps$dts`,
                                                        y = mean.abund, group = 1, color = "Model")) +
  geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
                   ymax = mean.abund + 1.96 * se.abund),
               colour = 'transparent',
               alpha = .5,
               show.legend = T) +
  geom_line(show.legend = T, linewidth = 0.7) +
  geom_line(data = NZMS.samp.sum, aes(x = as.Date(V1, origin = "1970-01-01"), y = sums, color = "Empirical"), show.legend = T)+
  #geom_line(data = flow.magnitude, aes(x = as.Date(dts), y = X_00060_00003), color = "blue") +
  #geom_line(data = temps, aes(x = as.Date(dts), y = Temperature*1000), color = "green")+
  coord_cartesian(ylim = c(0,600)) +
  ylab('New Zealand Mudsnail Abundance Density (m2)') +
  xlab("")+
  labs(colour=" ")
  theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), )+
  scale_x_date(date_labels="%B", date_breaks  ="6 months")

ggplot(data = adults, aes(x = Time,
                                   y = Adult, group = 1)) +
  # #geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
  #                 ymax = mean.abund + 1.96 * se.abund),
  #             colour = 'transparent',
  #             alpha = .5,
  #             show.legend = FALSE) +
  geom_line(show.legend = FALSE, linewidth = 0.7) +
  geom_point()+
  geom_line(data = NZMS.samp, aes(x = V1, y = exp(x)), color = "red")+
  geom_line(data = flow.magnitude, aes(x = as.Date(dts), y = X_00060_00003), color = "blue") +
  geom_line(data = temps, aes(x = as.Date(dts), y = Temperature*1000), color = "green")+
  coord_cartesian(ylim = c(0,70000)) +
  ylab('New Zealand Mudsnail Abundance') +
  xlab("")+
  theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13))+
  scale_x_date(date_labels="%B", date_breaks  ="6 months")


drift.data.total <- readDB(gear = "Drift", type = "Sample", updater = TRUE)

drift.LF <- drift.data.total[which(drift.data.total$RiverMile >= -3 & drift.data.total$RiverMile <= 3),]

#Baet.samp.LF <- sampspec(samp = drift.LF, species = "BAET", stats = T)

#View(Baet.samp.LF$Samples)
#m <- merge(Baet.samp.LF$Statistics, Baet.samp.LF$Samples, by = "BarcodeID", all = T)
discharge
NZMS.samp.LF <- sampspec(samp = drift.LF, species = "NZMS", stats = T)
NZMS.samp <- merge(NZMS.samp.LF$Statistics, NZMS.samp.LF$Samples, by = "BarcodeID", all = T)
NZMS.samp <- NZMS.samp[which(NZMS.samp$GearID == 4),]
NZMS.samp$Density <- NZMS.samp$CountTotal/NZMS.samp$Volume
#NZMS.samp <- merge(NZMS.samp, discharge[, 3:4], by = "Date")
#NZMS.samp$Density <- NZMS.samp$Density/((NZMS.samp$X_00060_00003 * 0.028316832)^4.1)
NZMS.samp <- aggregate(NZMS.samp$Density, list(NZMS.samp$Date), FUN = sum)

calculate_sum_within_interval <- function(start_time, end_time, dates, values) {
  sum(values[dates >= start_time & dates <= end_time])
}

# Calculate the sum within each time interval
interval_sums <- mapply(calculate_sum_within_interval, start_time = temps$dts[-length(temps$dts)], 
                        end_time = temps$dts[-1], dates = NZMS.samp$Group.1, values = NZMS.samp$x)


calculate_sum_within_interval(start_time = temps$dts[31], end_time = temps$dts[32], dates = NZMS.samp$Group.1, values = NZMS.samp$x)


sums <- vector()
for (i in 1:length(temps$dts)){
  d <- NZMS.samp[which(NZMS.samp$Group.1 %within% interval(temps$dts[i], temps$dts[i+1]) == T),]
  if (sum(d$x) == 0) {
    s = NA
  } else {
  s<- sum(d$x)}
  sums <- append(sums, s)
}
NZMS.samp.sum <- na.omit(as.data.frame(cbind(as.Date(temps$dts), sums)))
NZMS.samp.sum$V1 <- as.Date(NZMS.samp.sum$V1, origin = "1970-01-01")

cor.df <- left_join(NZMS.samp.sum, means.list.NZMS, by=c('V1'='temps$dts'))
cor.test(cor.df$sums, cor.df$mean.abund, method = "pearson")
