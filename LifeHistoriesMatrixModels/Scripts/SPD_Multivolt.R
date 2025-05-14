######################
# Multivoltinism Sp D
######################

source("LifeHistoriesMatrixModels/Scripts/D_1sp_Model.R")
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


discharge <- rep(0.1, time = length(temp$dts))

# run model under different temp regimes

runs <- c(-2, 6)

# for adults

for (i in 1:length(runs)){
  temp$Temperature <- temp$Temperature + runs[i]
  out <- Dmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), stage_output = "3")
  m <- rowMeans(out)
  m <- cbind.data.frame(temp$dts, m[-1], rep(i, length(temp$dts)))
  assign(paste0("m",i), m)
  temp$Temperature <- temp$Temperature - runs[i]
}

mlist <- rbind(m1, m2)
colnames(mlist) <- c("Date", "Abund", "MeanTemp")
D.oneyear <- mlist[which(mlist$Date >= "2035-01-01" & mlist$Date <= "2035-12-31"), ]
peaks <- D.oneyear[c(16, 40, 46, 61, 66, 71, 80, 86, 91, 95, 100), ]
arrows <- tibble(
  x1 = peaks$Date,
  x2 = peaks$Date,
  y1 = c(8000, 16900, 14000, 6100, 14800, 11000, 4200, 8000 , 6300, 8500, 20200), 
  y2 = c(6000, 14900, 12000, 4100, 12800, 9000, 2200, 6000 , 4300, 6500, 18200)
)
arrowcols <- c("#4477AA", "#EE6677","#EE6677", "#228833","#228833","#228833","#CCBB44","#CCBB44","#CCBB44", "#CCBB44","#CCBB44")


arrows$x1 <- as.POSIXct(arrows$x1, format = "%Y-%m-%d")

arrows$x1 <- as.Date(arrows$x1)
arrows$x2 <- as.POSIXct(arrows$x2, fomrat = "%Y-%m_%d")
arrows$x2 <- as.Date(arrows$x2)
D.oneyear$Date <- as.Date(D.oneyear$Date)
D.oneyear$Strategy <- rep("Slow", times = length(D.oneyear$Date))