############################
# Code to run Species C code
############################

# Time <- c(1:365)
# Date <- rep(c(1:365), times = 100)
# Day <- seq(as.Date("2022-01-01"), as.Date("2222-12-31"), by="days")
# Day <- Day[-which(Day == "2024-02-29" & Day = "2028-02-29" & Day = "20")]
# 
# Temperature <-  -7.374528  * (cos(((2*pi)/365)*Date))  +  (-1.649263* sin(2*pi/(365)*Date)) + 10.956243
# 
# temp <- as.data.frame(cbind(Time, Day, Temperature))
# temp$Day <- as.Date(temp$Day, origin= "1970-01-01")
# colnames(temp) <- c("Time", "Date", "Temperature")
# temp <- TimestepTemperature(temp)
# temp <- temp[c(1,3)]
# discharge <- rep(0, times = 131)
# flow.magnitude <- as.data.frame(cbind(temp$dts, discharge))
# 
library(lubridate)
library(parallel)



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


Temperature <-  -7.374528  * (cos(((2*pi)/365)*Date))  +  (-1.649263* sin(2*pi/(365)*Date)) + 10.956243

temp <- as.data.frame(cbind(Time, Day, Temperature))
temp$Day <- as.Date(temp$Day, origin= "1970-01-01")
colnames(temp) <- c("Time", "Date", "Temperature")
temp <- TimestepTemperature(temp)
temp <- temp[c(1,3)]

# pull one random day in March for each 
Year <- year(temp$dts)
uYear <- unique(Year)
Month <- month(temp$dts)

season <- as.numeric(unlist(c(1, 4, 7, 10)))

for (s in season){
  results <- lapply(uYear, function(x){
    sample_dates <- temp$dts[Month == s & Year == x]
    return(sample(sample_dates, size = 1))
  })
  
  discharge <- rep(0.1, time = length(temp$dts))
  results <- do.call("c", results)
  discharge[match(results, temp$dts)] <- 0.3
  flow.magnitude <- as.data.frame(cbind(temp$dts, discharge))

system.time( 
  
out <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 99, peaklist = 0, peakeach = length(temp$Temperature))

)
m <- mean.data.frame(out, burnin = 250, iteration = 99)
assign(paste0("means.list.", s), m)
}



means.list <- as.data.frame(cbind(means.list.1$mean.abund, means.list.4$mean.abund, means.list.7$mean.abund, means.list.10$mean.abund, temp$dts[250:length(means.list.1$timesteps)]))
means.list[,5] <- as.Date(as.POSIXct(means.list[,5], origin = "1970-01-01"))
colnames(means.list) <- c("January", "April", "July", "October", "Time")
x11()
means.list.abridged <- means.list[106:158,]
ggplot(data = means.list.abridged, aes(x = Time,
                              y = January/10000, group = 1, color = "January")) +
  #geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
  #                ymax = mean.abund + 1.96 * se.abund),
  #            colour = 'transparent',
  #            alpha = .15,
  #            show.legend = T) +
  geom_line(show.legend = T) +
  geom_line(data = means.list.abridged, aes(x = Time, y = April/10000, color = "April"), show.legend = T)+
  geom_line(data = means.list.abridged, aes(x = Time, y = July/10000, color = "July"), show.legend = T)+
  geom_line(data = means.list.abridged, aes(x = Time, y = October/10000, color = "October"), show.legend = T)+
  coord_cartesian(ylim = c(0,2.5)) +
  ylab('Caddisfly Abundance Relative to Baseline Recuritment Limit') +
  xlab("")+
  labs(colour=" ")+
  theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), )+
  scale_x_date(date_labels="%B", date_breaks  ="4 months")


means.list.C <- mean.data.frame(out,burnin = 50, iteration= 1)
means.list.C <- cbind(means.list.C, temp$dts[1:length(means.list.C$timesteps)])
means.list.C$`temp$dts` <- as.Date(means.list.C$`temp$dts`)

plot(means.list.C$`temp$dts`, means.list.C$mean.abund, type = "l", xlab = "Time", ylab = "Type C species")
