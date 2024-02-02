##################################
## Sp C Press Magnitude
##################################
library(lubridate)
source("C_1sp_Model.R")
source("1spFunctions.R")
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

discharge <- rep(0.1, length(temp$dts))

magnitudes <- seq(0, 1, by = 0.05)
mag_response <- vector()
short_response <- vector() 
for (i in 1:length(magnitudes)){
  out <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 10, peaklist = magnitudes[i], peakeach = length(temp$Temperature))
  m <- mean.data.frame(out, burnin = 250, iteration = 10)
  mag_response[i] <- mean(m$mean.abund)
  }
# calculate immediate response to the different magnitudes

c_magnitude_df <- as.data.frame(cbind(magnitudes, mag_response, rep("C", times = length(mag_response))))
#short_df <- as.data.frame(cbind(magnitudes, short_response))
