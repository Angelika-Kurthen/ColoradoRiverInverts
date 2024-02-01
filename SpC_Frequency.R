#############################
## Sp C Pulse Frequency 
############################
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

Year <- year(temp$dts)
uYear <- unique(Year)
Month <- month(temp$dts)

iterations <- 100
# Filter temp$dts to only include dates in the year 2035
dates_2035 <- temp$dts[format(temp$dts, "%Y") == "2035"]
m_array <- array(data = NA, dim = c(26, iterations))
short_array <- array(data = NA, dim = c(26,iterations))
long_array <-array(data = NA, dim = c(26,iterations))
# this is how often the disturbance occurs from once to every other week
for (j in 1:iterations){ # will want to repeat quite a few times
for (i in 1:26){# Randomly select dates without replacement
  random_dates <- sample(dates_2035, i)
  # create a list of non-disturbance discharges
  discharge <- rep(0.1, time = length(temp$dts))
  # from that list of dates from above, assign a disturbance discharge date
  discharge[match(random_dates, temp$dts)] <- 0.3
  # run model
  out <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature))
  # create summary dataframe 
  sums <- rowSums(out)
  m <- cbind(sums[250:402], discharge[250:402])
  last <- max(which(m[,2] == 0.3))
  immediate_response <- m[last+1,1]
  short_term <- mean(m[(last + 2):(last + 6), 1])
  long_term <- mean(m[(last + 26):(last + 36),1])
  m_array[i,j] <- immediate_response
  short_array[i,j] <- short_term
  long_array[i,j] <- long_term
}
}

immediate <- rowMeans(m_array)
c_immediate_df <- as.data.frame(cbind(immediate, seq(1:26), rep("C", length(immediate))))
short <- rowMeans(short_array)
c_short_df <- as.data.frame(cbind(short, seq(1:26), rep("C", length(short))))
long <- rowMeans(long_array)
c_long_df <- as.data.frame(cbind(long, seq(1:26), rep("C", length(long))))

# c_imm <- ggplot(data = immediate_df, aes(x = V2, y = immediate/10000))+
#   geom_line(size=1 ,col = "#EE6677")+
#   xlab("Frequency of pulse disturbance per year")+
#   ylab("Sp C abundance relative to K")+
#   theme_bw()+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# 
# 
# c_short <- ggplot(data = short_df, aes(x = V2, y = short/10000))+
#   geom_line(size = 1, col = "#EE6677")+
#   xlab("Frequency of pulse disturbance per year")+
#   ylab("Sp C abundance relative to K")+
#   theme_bw()+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# 
# 
# c_long <- ggplot(data = long_df, aes(x = V2, y = long/10000))+
#   geom_line(size = 1, col = "#EE6677")+
#   xlab("Frequency of pulse disturbance per year")+
#   ylab("Sp C abundance relative to K")+
#   theme_bw()+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# 
