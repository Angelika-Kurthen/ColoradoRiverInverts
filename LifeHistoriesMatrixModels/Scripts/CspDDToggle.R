#############################
# C sp degree day toggle
############################
library(lubridate)
#Code for HPC
#library(lubridate, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")



source("C_1sp_Model.R")
source("1spFunctions.R")
source("NegExpSurv.R")

Time <- c(1:36500)
Date <- rep(1:365, times = 100)
Day <- seq(as.Date("2022-01-01"), as.Date("2121-12-31"), by="days")
# find and remove leap days
find_leap = function(x){
  day(x) == 29 & month(x) == 2 
}
Day <- Day[which(find_leap(Day) == F)]


Temperature <-  -7.374528  * (cos(((2*pi)/365)*Date))  +  (-1.649263* sin(2*pi/(365)*Date)) + 13.956243

temp <- as.data.frame(cbind(Time, Day, Temperature))
temp$Day <- as.Date(temp$Day, origin= "1970-01-01")
colnames(temp) <- c("Time", "Date", "Temperature")
temp <- TimestepTemperature(temp)
temp <- temp[c(1,3)]

Year <- year(temp$dts)
uYear <- unique(Year)
Month <- month(temp$dts)
discharge <- rep(0.1, time = length(temp$dts))

dd_seq <- seq(900*0.9, 900*1.1, length.out = 21)


# Run C sp (itermediate response to different scenarios) with no HFE, no hydropeaking (so no disturbance), temperature of North American temperate stream between 3 and 18 C, and fecundities from between 20 and 2000 eggs per brood. 
dd_means <- vector()
for (dd in 1:length(dd_seq)){
  print(dd)
  out <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), dds = dd_seq[dd])
  means.list.C <- mean.data.frame(out, burnin = 250, iteration = 2) 
  dd_means[dd] <- mean(means.list.C$mean.abund)
}

cdd_df <- as.data.frame(cbind(dd_seq, dd_means, rep("C", length(dd_means))))
cdd_df$dd_seq <- as.numeric(cdd_df$dd_seq)
cdd_df$dd_means <- as.numeric(cdd_df$dd_means)
cdd_lm <- lm((dd_means/10000) ~ dd_seq, data = cdd_df)
# cdd <- ggplot(data = dd_df, mapping = aes(x = dd_seq, y = dd_means/10000))+
#   geom_point(size = 1, col = "#EE6677")+
#   xlab("Degree Day Requirement")+
#   ylab("C sp Abundance Relative to K")+
#   theme_bw()
