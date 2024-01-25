###########################
## D sp Degree Day toggle
##########################
library(lubridate)


source("D_1sp_Model.R")
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

dd_seq <- seq(1500*0.9, 1500*1.1, length.out = 21)


# Run D sp (itermediate response to different scenarios) with no HFE, no hydropeaking (so no disturbance), temperature of North American temperate stream between 3 and 18 C, and fecundities from between 20 and 2000 eggs per brood. 
dd_means <- vector()
for (dd in 1:length(dd_seq)){
  print(dd)
  out <- Dmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), dds = dd_seq[dd])
  means.list.D <- mean.data.frame(out, burnin = 250, iteration = 2) 
  dd_means[dd] <- mean(means.list.D$mean.abund)
}

dd_df <- as.data.frame(cbind(dd_seq, dd_means))
ddd_lm <- lm((dd_means/10000) ~ dd_seq, data = dd_df)
ddd <- ggplot(data = dd_df, mapping = aes(x = dd_seq, y = dd_means/10000))+
  geom_point(size = 1, col = "#AA3377")+
  xlab("Degree Day Requirement")+
  ylab("D sp Abundance Relative to K")+
  theme_bw()
