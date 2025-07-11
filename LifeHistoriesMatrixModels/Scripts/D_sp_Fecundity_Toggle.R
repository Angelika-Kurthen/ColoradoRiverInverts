#############################
# B sp fecundity toggle
############################
library(lubridate)
#Code for HPC
#library(lubridate, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")

source("LifeHistoriesMatrixModels/Scripts/D_1sp_Model.R")
source("LifeHistoriesMatrixModels/Scripts/1spFunctions.R")
source("LifeHistoriesMatrixModels/Scripts/NegExpSurv.R")

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

fec_seq <- seq(300*0.9, 300*1.1, length.out = 21)


# Run C sp (itermediate response to different scenarios) with no HFE, no hydropeaking (so no disturbance), temperature of North American temperate stream between 3 and 18 C, and fecundities from between 20 and 2000 eggs per brood. 
fec_means <- vector()
for (fec in 1:length(fec_seq)){
  print(fec)
  out <- Dmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), fecundity = fec_seq[fec])
  means.list.D <- mean.data.frame(out, burnin = 250, iteration = 2) 
  fec_means[fec] <- mean(means.list.D$mean.abund)
}

d_fec_df <- as.data.frame(cbind(fec_seq, fec_means, rep("D", times = length(fec_means))))
d_fec_df$fec_seq <- as.numeric(d_fec_df$fec_seq)
d_fec_df$fec_means <- as.numeric(d_fec_df$fec_means)

dfec <- dfec_lm <- lm((fec_means/10000) ~ fec_seq, data = d_fec_df)
summary(dfec_lm)
# dfec <- ggplot(data = fec_df, mapping = aes(x = fec_seq, y = fec_means/10000))+
#   geom_point(size = 1, col = "#AA3377")+
#   xlab("Egg Mass Size")+
#   ylab("D sp Abundance Relative to K")+
#   theme_bw()


