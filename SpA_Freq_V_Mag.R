##################################
## Sp A Flood Pulse Magnitude
##################################
#library(lubridate)
library(lubridate, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
library(tidyr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
source("A_1sp_Model.R")
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

Year <- year(temp$dts)
uYear <- unique(Year)
Month <- month(temp$dts)

discharge <- rep(0.1, length(temp$dts))

magnitudes <- seq(0.11, 1, by = 0.04 )

iterations <- 99
# Filter temp$dts to only include dates in the year 2035
dates_2035 <- temp$dts[format(temp$dts, "%Y") == "2035"]

# create empty arrays to store results
imm_array <- array(data = NA, dim = c(27, iterations))
short_array <- array(data = NA, dim = c(27,iterations))
long_array <-array(data = NA, dim = c(27,iterations))
immediate  <- array(data = NA, dim = c(27, length(magnitudes)))
short <- array(data = NA, dim = c(27, length(magnitudes)))
long <- array(data = NA, dim = c(27, length(magnitudes)))



# Initializes the progress bar
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = length(magnitudes) * iterations * 26, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar
steps <- 0

for (mts in 1:length(magnitudes)){ # iterate over each magnitude
  
  for (j in 1:iterations){ # will want to repeat quite a few times
    # this is how often the disturbance occurs from once to every other week
    for (i in 1:27){# Randomly select dates without replacement
      random_dates <- sample(dates_2035, i)
      # create a list of non-disturbance discharges
      if (i == 1) { # unless 0 disturbance regime
        # Randomly select dates without replacement
        discharge <- rep(0.1, time = length(temp$dts))
        random_date <- sample(dates_2035, i)
        discharge[match(random_date, temp$dts)] <- 0.1
      }else{
        # Randomly select dates without replacement
        discharge <- rep(0.1, time = length(temp$dts))
        random_date <- sample(dates_2035, i-1)
        discharge[match(random_dates, temp$dts)] <- magnitudes[mts]
      }
      # run model
      out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
      # create summary dataframe 
      m <- mean.data.frame(data = out, burnin = 250, iteration = 2)
      m <- cbind.data.frame(temp$dts[250:last(m$timesteps)], m)  
      colnames(m) <- c("time", 'timestep', "abund", "sd", "se")
      imm_array[i,j] <- m$abund[which(m$time == last(random_date))+1]
      short_array[i,j] <- max(m$abund[(which(m$time == last(random_date))):(which(m$time == last(random_date)) + 6)])
      
      steps<-steps+1
      setTxtProgressBar(pb, steps)
    }
  }
  immediate[,mts] <- rowMeans(imm_array)
  short[,mts] <- rowMeans(short_array)
  long[,mts] <- rowMeans(long_array)
}

immediate_df <- cbind.data.frame(immediate, seq(0,26))
colnames(immediate_df) <- c(magnitudes, "frequency")
a_immediate_df <- pivot_longer(immediate_df, cols = 1:length(magnitudes), names_to = "magnitude", values_to = "abundance")
a_immediate_df$magnitude <- as.numeric(a_immediate_df$magnitude)
write.csv(a_immediate_df, file = "SpA_FreqVMag_immediate.csv")
# 
# ggplot(data = immediate_df, aes(x = frequency, y = magnitude))+
#   geom_raster(aes(fill = abundance/10000), interpolate = T)+
#   scale_fill_viridis_c(option = "magma") +
#   xlab("Pulse Frequency")+
#   ylab("Pulse Magnitude")+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#   axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
#   guides(fill=guide_legend(title="Sp A Relativized Abundance"))
# 
# 
  
short_df <- cbind.data.frame(short, seq(0,26))
colnames(short_df) <- c(magnitudes, "frequency")
a_short_df <- pivot_longer(short_df, cols = 1:length(magnitudes), names_to = "magnitude", values_to = "abundance")
a_short_df$magnitude <- as.numeric(a_short_df$magnitude)
write.csv(a_short_df, file = "SpA_FreqVMag_short.csv")

# png("SpA_freq_v_mag_short.png")
# plot(a_short)
# dev.off()

# 
# long_df <- cbind.data.frame(long, seq(1:26))
# colnames(long_df) <- c(magnitudes, "frequency")
# a_long_df <- pivot_longer(long_df, cols = 1:length(magnitudes), names_to = "magnitude", values_to = "abundance")
# a_long_df$magnitude <- as.numeric(a_long_df$magnitude)
# write.csv(a_long_df, file = "SpA_FreqVMag_long.csv")

# a_long <- ggplot(data = long_df, aes(x = frequency, y = magnitude))+
#   geom_raster(aes(fill = abundance/10000), interpolate = F)+
#   scale_fill_viridis_c() +
#   xlab("Pulse Frequency")+
#   ylab("Pulse Magnitude")+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#   axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
#   guides(fill=guide_legend(title="Sp A Relativized Abundance"))
# 
# png("SpA_freq_v_mag_long.png")
# plot(a_long)
# dev.off()
