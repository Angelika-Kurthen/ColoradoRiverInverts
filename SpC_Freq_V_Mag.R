##################################
## Sp C Flood Pulse Magnitude
##################################
library(lubridate)
# library(lubridate, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")

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

Year <- year(temp$dts)
uYear <- unique(Year)
Month <- month(temp$dts)

discharge <- rep(0.1, length(temp$dts))

magnitudes <- seq(0.11, 1, by = 0.04 )

iterations <- 99
# Filter temp$dts to only include dates in the year 2035
dates_2035 <- temp$dts[format(temp$dts, "%Y") == "2035"]

# create empty arrays to store results
imm_array <- array(data = NA, dim = c(26, iterations))
short_array <- array(data = NA, dim = c(26,iterations))
long_array <-array(data = NA, dim = c(26,iterations))
immediate  <- array(data = NA, dim = c(26, length(magnitudes)))
short <- array(data = NA, dim = c(26, length(magnitudes)))
long <- array(data = NA, dim = c(26, length(magnitudes)))



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
    for (i in 1:26){# Randomly select dates without replacement
      random_dates <- sample(dates_2035, i)
      # create a list of non-disturbance discharges
      discharge <- rep(0.1, time = length(temp$dts))
      # from that list of dates from above, assign a disturbance discharge date
      discharge[match(random_dates, temp$dts)] <- magnitudes[mts]
      # run model
      out <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature))
      # create summary dataframe 
      sums <- rowSums(out)
      m <- cbind(sums[250:402], discharge[250:402])
      last <- max(which(m[,2] == magnitudes[mts]))
      immediate_response <- m[last+1,1]
      short_term <- mean(m[(last + 2):(last + 6), 1])
      long_term <- mean(m[(last + 26):(last + 36),1])
      imm_array[i,j] <- immediate_response
      short_array[i,j] <- short_term
      long_array[i,j] <- long_term
      
      steps<-steps+1
      setTxtProgressBar(pb, steps)
    }
  }
  immediate[,mts] <- rowMeans(imm_array)
  short[,mts] <- rowMeans(short_array)
  long[,mts] <- rowMeans(long_array)
}

immediate_df <- cbind.data.frame(immediate, seq(1:26))
colnames(immediate_df) <- c(magnitudes, "frequency")
c_immediate_df <- pivot_longer(immediate_df, cols = 1:length(magnitudes), names_to = "magnitude", values_to = "abundance")
c_immediate_df$magnitude <- as.numeric(c_immediate_df$magnitude)
write.csv(c_immediate_df, file = "SpC_FreqVMag_immediate.csv")

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
  
short_df <- cbind.data.frame(short, seq(1:26))
colnames(short_df) <- c(magnitudes, "frequency")
c_short_df <- pivot_longer(short_df, cols = 1:length(magnitudes), names_to = "magnitude", values_to = "abundance")
c_short_df$magnitude <- as.numeric(c_short_df$magnitude)
write.csv(c_short_df, file = "SpC_FreqVMag_short.csv")

# png("SpA_freq_v_mag_short.png")
# plot(a_short)
# dev.off()


long_df <- cbind.data.frame(long, seq(1:26))
colnames(long_df) <- c(magnitudes, "frequency")
c_long_df <- pivot_longer(long_df, cols = 1:length(magnitudes), names_to = "magnitude", values_to = "abundance")
c_long_df$magnitude <- as.numeric(c_long_df$magnitude)
write.csv(c_long_df, file = "SpC_FreqVMag_long.csv")
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
