#############################
# A sp temperature toggle
############################


#Code for HPC
#library(lubridate, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")

library(lubridate)


source("A_1sp_Model.R")
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

temp_regime <- vector()
temp_means <- vector()
temp_seq <- seq(-10, 10, by = 1)
for (te in 1:length(temp_seq)){
  temp$Temperature <- temp$Temperature + temp_seq[te]
  temp_regime[te] <- mean(temp$Temperature)
  out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
  temp$Temperature <- temp$Temperature - temp_seq[te]
  means.list.A <- mean.data.frame(out, burnin = 250, iteration = 2) 
  temp_means[te] <- mean(means.list.A$mean.abund)
}


a_temp_adjust_df <- as.data.frame(cbind(temp_regime, temp_means, rep("A", times = length(temp_means))))
a_temp_adjust_df$temp_regime <- as.numeric(a_temp_adjust_df$temp_regime)
a_temp_adjust_df$temp_means <- as.numeric(a_temp_adjust_df$temp_means)


size_means <- vector()
for (te in 1:length(temp_seq)){
  temp$Temperature <- temp$Temperature + temp_seq[te]
  temp_regime[te] <- mean(temp$Temperature)
  out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature), stage_output = "size")
  temp$Temperature <- temp$Temperature - temp_seq[te]
  size_means[te] <- mean(out)
}

size_means <- 0.0094*(size_means)^2.754  # multiply relative size (which is also biologically plausible) by Benke et al 1999 Table 2 a and b params (M(mg) = aL^b) 

a_size_df <- as.data.frame(cbind(temp_regime, size_means, rep("A", times = length(temp_means))))
a_size_df$temp_regime <- as.numeric(a_size_df$temp_regime)
a_size_df$size_means <- as.numeric(a_size_df$size_means)

# atemp <- ggplot(data = temp_adjust_df, mapping = aes(x = temp_seq, y = temp_means/10000))+
#   geom_line(size= 1, col = "#66CCEE")+
#   xlab("Degree C Change")+
#   ylab("Sp A Abundance Relative to K")+
#   theme_bw()
# 


# Disturbance by temp

temp_regime <- vector()
temp_means <- vector()
temp_seq <- seq(-10, 10, by = 1)
short <- vector()
discharge[259] <- 1
for (te in 1:length(temp_seq)){
  temp$Temperature <- temp$Temperature + temp_seq[te]
  temp_regime[te] <- mean(temp$Temperature)
  out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
  temp$Temperature <- temp$Temperature - temp_seq[te]
  means.list.A<- mean.data.frame(out, burnin = 250, iteration = 2)
  short[te] <- mean(means.list.A$mean.abund[10:16])
}

winter <- as.data.frame(cbind(temp_regime, log(short)))

size_means <- vector()
for (te in 1:length(temp_seq)){
  temp$Temperature <- temp$Temperature + temp_seq[te]
  temp_regime[te] <- mean(temp$Temperature)
  out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature), stage_output = "size")
  temp$Temperature <- temp$Temperature - temp_seq[te]
  size_means[te] <- mean(out)
}

size_means <- 0.0094*(size_means)^2.754  # multiply relative size (which is also biologically plausible) by Benke et al 1999 Table 2 a and b params (M(mg) = aL^b) 
winter_size_means <- as.data.frame(cbind(temp_regime, size_means))


temp_regime <- vector()
temp_means <- vector()
temp_seq <- seq(-10, 10, by = 1)
short <- vector()
discharge[272] <- 1
for (te in 1:length(temp_seq)){
  temp$Temperature <- temp$Temperature + temp_seq[te]
  temp_regime[te] <- mean(temp$Temperature)
  out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
  temp$Temperature <- temp$Temperature - temp_seq[te]
  means.list.A<- mean.data.frame(out, burnin = 250, iteration = 2)
  short[te] <- mean(means.list.A$mean.abund[23:29])
}
summer <- as.data.frame(cbind(temp_regime,log(short)))

size_means <- vector()
for (te in 1:length(temp_seq)){
  temp$Temperature <- temp$Temperature + temp_seq[te]
  temp_regime[te] <- mean(temp$Temperature)
  out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature), stage_output = "size")
  temp$Temperature <- temp$Temperature - temp_seq[te]
  size_means[te] <- mean(out)
}
size_means <- 0.0094*(size_means)^2.754  # multiply relative size (which is also biologically plausible) by Benke et al 1999 Table 2 a and b params (M(mg) = aL^b) 
summer_size_means <- as.data.frame(cbind(temp_regime, size_means))

# bind together, 1 = winter 2 = summer
# bind together, 1 = winter 2 = summer
# bind together, 1 = winter 2 = summer
temp_dist_a <- bind_rows(winter, summer, .id = "season")
sizes <- rbind(winter_size_means, summer_size_means)
deltatemp_a <- as.data.frame(cbind(rep(3, times = length(temp_regime)),temp_regime, summer[,2]-winter[,2]))
temp_size_a <- bind_rows(winter_size_means, summer_size_means, .id = "season")
deltasize_a <- as.data.frame(cbind(rep(3, times = length(temp_regime)), temp_regime, (summer_size_means[,2]*summer[,2])-(winter_size_means[,2]*winter[2])))
temp_size_a <- mutate(.data = temp_size_a, size_means = temp_dist_a$short * sizes$size_means )

deltatemp_a <- setNames(deltatemp_a, names(temp_dist_a))
temp_dist_a <- rbind(temp_dist_a, deltatemp_a)
deltasize_a <- setNames(deltasize_a, names(temp_size_a))
temp_size_a <- rbind(temp_size_a, deltasize_a)