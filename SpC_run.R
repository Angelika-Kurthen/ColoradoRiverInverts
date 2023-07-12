

# Time <- c(1:365)
# Date <- rep(c(1:365), times = 100)
# Day <- seq(as.Date("2022-01-01"), as.Date("2222-12-31"), by="days")
# Day <- Day[-which(Day == "2024-02-29" & Day = "2028-02-29" & Day = "20")]
# 
Temperature <-  -7.374528  * (cos(((2*pi)/365)*Date))  +  (-1.649263* sin(2*pi/(365)*Date)) + 10.956243

temp <- as.data.frame(cbind(Time, Day, Temperature))
temp$Day <- as.Date(temp$Day, origin= "1970-01-01")
colnames(temp) <- c("Time", "Date", "Temperature")
temp <- TimestepTemperature(temp)
temp <- temp[c(1,3)]
discharge <- rep(0, times = 131)
flow.magnitude <- as.data.frame(cbind(temp$dts, discharge))

source("C_1sp_Model.R")


out <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 9, peaklist = 0.13, peakeach = 
                length(temps$Temperature))
means.list.C <- mean.data.frame(out,burnin = 50, iteration= 1)
means.list.C <- cbind(means.list.C, temp$dts[1:length(means.list.C$timesteps)])
means.list.C$`temp$dts` <- as.Date(means.list.C$`temp$dts`)

plot(means.list.C$`temp$dts`, means.list.C$mean.abund, type = "l", xlab = "Time", ylab = "Type C species")
