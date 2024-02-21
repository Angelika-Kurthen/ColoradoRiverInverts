#############################
## Sp D Pulse Frequency 
############################
library(lubridate)
source("D_1sp_Model.R")
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

iterations <- 99
# Filter temp$dts to only include dates in the year 2035
dates_2035 <- temp$dts[format(temp$dts, "%Y") == "2035"]
m_array <- array(data = NA, dim = c(27, iterations))
short_array <- array(data = NA, dim = c(27,iterations))
immediate_response <- vector()
short_response <- vector()
# this is how often the disturbance occurs from once to every other week
for (j in 1:iterations){ # will want to repeat quite a few times
for (i in seq(1,27)){
  # create a list of non-disturbance discharges
  discharge <- rep(0.1, time = length(temp$dts))
  # from that list of dates from above, assign a disturbance discharge date
  if (i == 1) { # unless 0 disturbance regime
    # Randomly select dates without replacement
    random_date <- sample(dates_2035, i)
    discharge[match(random_date, temp$dts)] <- 0.1
  }else{
    # Randomly select dates without replacement
    random_date <- sample(dates_2035, i-1)
    discharge[match(random_date, temp$dts)] <- 0.3
    }
  # run model
  out <- Dmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
  # create summary dataframe 
  m <- mean.data.frame(data = out, burnin = 250, iteration = 2)
  m <- cbind.data.frame(temp$dts[250:last(m$timesteps)], m)   
  colnames(m) <- c("time", 'timestep', "abund", "sd", "se")
  m_array[i,j] <- m$abund[which(m$time == last(random_date))+1]
  short_array[i,j] <- mean(m$abund[(which(m$time == last(random_date))):(which(m$time == last(random_date)) + 6)])
}
}

immediate <- rowMeans(m_array)
d_immediate_df <- as.data.frame(cbind(immediate, seq(0,26), rep("D", length(immediate))))
short <- rowMeans(short_array)
d_short_df <- as.data.frame(cbind(short, seq(0,26), rep("D", length(short))))

d_immediate_df$immediate <- as.numeric(d_immediate_df$immediate)
# a_imm <- ggplot(data = a_immediate_df, aes(x = V2, y = immediate/10000, group = 1))+
#   geom_line(size=1 ,col = "#EE6677")+
#   xlab("Frequency of pulse disturbance per year")+
#   ylab("Sp C abundance relative to K")+
#   theme_bw()+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5),
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))



d_short_df$short <- as.numeric(d_short_df$short)
# a_short <- ggplot(data = a_short_df, aes(x = V2, y = short/10000, group = 1))+
#   geom_line(col = "#EE6677")+
#   xlab("Frequency of pulse disturbance per year")+
#   ylab("Sp C abundance relative to K")+
#   theme_bw()+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5),
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# 

# c_long <- ggplot(data = long_df, aes(x = V2, y = long/10000))+
#   geom_line(size = 1, col = "#EE6677")+
#   xlab("Frequency of pulse disturbance per year")+
#   ylab("Sp C abundance relative to K")+
#   theme_bw()+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# 
