##################################
## Sp B Flood Pulse Magnitude
##################################
library(lubridate)
source("B_1sp_Model.R")
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

selected_date <- temp$dts[temp$dts >= as.Date("2035-05-01") & temp$dts <= as.Date("2035-05-15")]
selected_date <- sample(selected_date, 1)

discharge <- rep(0.1, length(temp$dts))

magnitudes <- seq(0.1, 2, by = 0.05)
immediate_response <- vector()
short_response <- vector() 

for (i in 1:length(magnitudes)){
    discharge[which(temp$dts == selected_date)] <- magnitudes[i]
    # calculate the response to the different magnitudes
    out <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
    m <- mean.data.frame(data = out, burnin = 250, iteration = 2)
    m <- cbind.data.frame(temp$dts[250:last(m$timesteps)], m)   
    colnames(m) <- c("time", 'timestep', "mean.abund", "sd", "se")
    immediate_response[i] <- m$mean.abund[which(m$time == selected_date)+1]
    short_response[i] <- max(m$mean.abund[(which(m$time == selected_date)):(which(m$time == selected_date) + 8)])
  }
# calculate immediate response to the different magnitudes

b_magnitude_df <- as.data.frame(cbind(magnitudes, immediate_response, rep("B", times = length(immediate_response))))
b_short_df <- as.data.frame(cbind(magnitudes, short_response, rep("B", times = length(short_response))))
b_short_df$magnitudes <- as.numeric (b_short_df$magnitudes) 
b_short_df$short_response <- as.numeric(b_short_df$short_response)
bmag <- ggplot(data = b_short_df, aes(x = magnitudes, y = short_response/10000))+
  geom_line(size = 1, col = "#228833")+
  xlab("Discharge Magnitude (Proportion Bankfull)")+
  ylab("Sp B post pulse abundance relative to K")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5),
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

# # 
# alph <- seq(0.01,1, by = 0.05)
# tim <- seq(343, 349, by = 1)
# 
# for(ti in tim){
# for (al in alph){
#   selected_date <- temp$dts[ti]
#   discharge <- rep(0.1, length(temp$dts))
#   magnitudes <- seq(0.1, 1, by = 0.05)
#   immediate_response <- vector()
#   short_response <- vector()
#   for (i in 1:length(magnitudes)){
#     discharge[which(temp$dts == selected_date)] <- magnitudes[i]
#     # calculate the response to the different magnitudes
#     out <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), alph = al)
#     m <- mean.data.frame(data = out, burnin = 250, iteration = 2)
#     m <- cbind.data.frame(temp$dts[250:last(m$timesteps)], m)
#     colnames(m) <- c("time", 'timestep', "mean.abund", "sd", "se")
#     #immediate_response[i] <- m$mean.abund[which(m$time == selected_date)+1]
#     short_response[i] <- max(m$mean.abund[(which(m$time == selected_date)):(which(m$time == selected_date) + 7)])
#   }
#   # calculate immediate response to the different magnitudes
# 
#   #b_magnitude_df <- as.data.frame(cbind(magnitudes, immediate_response, rep("B", times = length(immediate_response))))
#   b_short_df <- as.data.frame(cbind(magnitudes, short_response, rep("B", times = length(short_response))))
#   b_short_df$magnitudes <- as.numeric (b_short_df$magnitudes)
#   b_short_df$short_response <- as.numeric(b_short_df$short_response)
#   bmag <- ggplot(data = b_short_df, aes(x = magnitudes, y = short_response/10000))+
#     geom_line(size = 1, col = "#228833")+
#     xlab("Discharge Magnitude (Proportion Bankfull)")+
#     ylab("Sp B post pulse abundance relative to K")+
#     theme_bw()+
#     theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5),
#           axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
#   png(filename = paste0(al,"Bmag_himort", selected_date, ".png"))
#   plot(bmag)
#   dev.off()
# }}
# #-----------------------
# selected_date <- temp$dts[temp$dts >= as.Date("2035-05-01") & temp$dts <= as.Date("2035-05-15")]
# selected_date <- sample(selected_date, 1)
# 
# discharge <- rep(0.1, length(temp$dts))
# 
# magnitudes <- seq(0.1, 1, by = 0.05)
# immediate_response <- vector()
# short_response <- vector()
# for (i in 1:length(magnitudes)){
#   discharge[which(temp$dts == selected_date)] <- magnitudes[i]
#   # calculate the response to the different magnitudes
#   out <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature))
#   m <- cbind.data.frame(temp$dts[250:length(out)], out[-c(1:249)])
#   colnames(m) <- c("time", "abund")
#   immediate_response[i] <- m$abund[which(m$time == selected_date)+1]
#   short_response[i] <- max(m$abund[(which(m$time == selected_date)):(which(m$time == selected_date) + 26)])
# }
# # calculate immediate response to the different magnitudes
# 
# b_S3magnitude_df <- as.data.frame(cbind(magnitudes, immediate_response, rep("B", times = length(immediate_response))))
# b_S3short_df <- as.data.frame(cbind(magnitudes, short_response, rep("B", times = length(immediate_response))))
# 
# b_S3short_df$short_response <- as.numeric(b_S3short_df$short_response)
# 
# bmag <- ggplot(data = b_S3short_df, aes(x = magnitudes, y = short_response, group = 1))+
#     geom_line(size = 1, col = "#228833")+
#     xlab("Discharge Magnitude (Proportion Bankfull)")+
#     ylab("Sp B post pulse abundance relative to K")+
#     theme_bw()+
#     theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5),
#           axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# 
