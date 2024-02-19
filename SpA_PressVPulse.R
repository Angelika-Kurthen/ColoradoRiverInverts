####################################
## Sp A Press vs Pulse Magnitude
####################################

library(lubridate)
library(tidyr)

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

discharge <- rep(0.1, length(temp$dts))
selected_date <- temp$dts[temp$dts >= as.Date("2035-05-01") & temp$dts <= as.Date("2035-05-15")]


press_magnitudes <- seq(0, 1, by = 0.1)
pulse_magnitudes <- seq(0.1, 1, by = 0.05 )

immediate_response <- array(data= NA, dim = c(length(press_magnitudes), length(pulse_magnitudes)))
short_response <-  array(data= NA, dim = c(length(press_magnitudes), length(pulse_magnitudes)))
max_response <- array(data= NA, dim = c(length(press_magnitudes), length(pulse_magnitudes)))

for (j in 1:length(pulse_magnitudes)){
for (i in 1:length(press_magnitudes)){
  discharge[which(temp$dts == selected_date)] <- pulse_magnitudes[j]
  out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = press_magnitudes[i], peakeach = length(temp$Temperature))
  m <- mean.data.frame(out, burnin = 250, iteration = 2)
  m <- cbind.data.frame(temp$dts[250:last(m$timesteps)], m)
  colnames(m) <- c("time", 'timestep', "abund", "sd", "se")
  immediate_response[i, j] <- m$abund[which(m$time == selected_date)+1]
  short_response[i, j] <- mean(m$abund[(which(m$time == selected_date) + 2):(which(m$time == selected_date) + 6)])
  max_response[i,j] <- max(m$abund[(which(m$time == selected_date)):(which(m$time == selected_date) + 6)])
}
}


#short_df <- as.data.frame(cbind(magnitudes, short_response))
immediate_df <- as.data.frame(cbind(immediate_response, press_magnitudes))
colnames(immediate_df) <- c(pulse_magnitudes, "Press_mag")
a_immediate_df <- pivot_longer(immediate_df, cols = 1:length(pulse_magnitudes), names_to = "Pulse_mag", values_to = "abundance")
a_immediate_df <- cbind.data.frame(a_immediate_df, rep("A", times = length(a_immediate_df$abundance)))
colnames(a_immediate_df) <- c("Press_mag", "Pulse_mag", "abundance", "Taxa")
# imm <- ggplot(data = immediate_df, aes(x = Press_mag, y = Pulse_mag))+
#   geom_raster(aes(fill = abundance/10000), interpolate = TRUE)+
#   scale_fill_viridis_c() +
#   xlab("Press Magnitude")+
#   ylab("Pulse Magnitude")+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
#   guides(fill=guide_legend(title="Sp A Relativized Abundance"))
short_df <- as.data.frame(cbind(short_response, press_magnitudes))
colnames(short_df) <- c(pulse_magnitudes, "Press_mag")
a_short_df <- pivot_longer(short_df, cols = 1:length(pulse_magnitudes), names_to = "Pulse_mag",  values_to = "abundance")
a_short_df <- cbind.data.frame(a_short_df, rep("A", times = length(a_short_df$abundance)))
colnames(a_short_df) <- c("Press_mag", "Pulse_mag", "abundance", "Taxa")


max_df <- cbind.data.frame(max_response, press_magnitudes)
colnames(max_df) <- c(pulse_magnitudes, "Press_mag")
a_max_df <- pivot_longer(max_df, cols = 1:length(pulse_magnitudes), names_to = "Pulse_mag", values_to = "abundance")
a_max_df <- cbind.data.frame(a_max_df, rep("A", times = length(a_max_df$abundance)))
colnames(a_max_df) <- c("Press_mag", "Pulse_mag", "abundance", "Taxa")


ggplot(data = a_immediate_df, aes(x = Press_mag, y = Pulse_mag))+
  geom_raster(aes(fill = abundance/10000), interpolate = TRUE)+
  scale_fill_viridis_c() +
  xlab("Press Magnitude")+
  ylab("Pulse Magnitude")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5),
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
  guides(fill=guide_legend(title="Sp A Relativized Abundance"))
