#############################
## Sp A Pulse Frequency 
############################
library(lubridate)
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
  out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature))
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
a_immediate_df <- as.data.frame(cbind(immediate, seq(1:26), rep("A", length(immediate))))
short <- rowMeans(short_array)
a_short_df <- as.data.frame(cbind(short, seq(1:26), rep("A", length(short))))
long <- rowMeans(long_array)
a_long_df <- as.data.frame(cbind(long, seq(1:26), rep("A", length(long))))

# a_imm <- ggplot(data = immediate_df, aes(x = V2, y = immediate/10000))+
#   geom_line(size=1 ,col = "#66CCEE")+
#   xlab("Frequency of pulse disturbance per year")+
#   ylab("Sp A abundance relative to K")+
#   theme_bw()+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# 
# 
# a_short <- ggplot(data = short_df, aes(x = V2, y = short/10000))+
#   geom_line(size = 1, col = "#66CCEE")+
#   xlab("Frequency of pulse disturbance per year")+
#   ylab("Sp A abundance relative to K")+
#   theme_bw()+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# 
# 
# a_long <- ggplot(data = long_df, aes(x = V2, y = long/10000))+
#   geom_line(size = 1, col = "#66CCEE")+
#   xlab("Frequency of pulse disturbance per year")+
#   ylab("Sp A abundance relative to K")+
#   theme_bw()+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# 
# 
# 
# means[i] <- list(m)
#     short_abund <- means[[i]][which(discharge[250:402] == "0.3"), 2]
#     resil_abund <- means[[i]][(which(discharge[250:402] == "0.3") + 1):(which(discharge[250:402] == "0.3") + 6), 2]
#     long_abund <- means[[i]][(which(discharge[250:402] == "0.3") + 26):(which(discharge[250:402] == "0.3") + 36), 2]
#     #post <- seq(0, 10, by = 1)
#     #dat <- rep(d, times = length(post))
#     #jday_data <- as.data.frame(rbind(jday_data, cbind(dat, post, abund)))
#     short_mean[i] <- short_abund
#     resil_mean[i] <- mean(resil_abund)
#     long_mean[i] <- mean(long_abund)
# 
# 
# which(discharge == 0.3)
# 
# 
# means <- list()
# #jday_data <- data.frame(date = NULL, post_dist = NULL, abund = NULL)
# jday_max <- vector()
# 
# short_mean <- vector()
# resil_mean <- vector()
# long_mean <- vector()
# 
# for (d in 1:length(all.dates)){ # 30 reps takes 60 mins
#   sample_dates <- temp$dts[which(format(temp$dts, "%m-%d")== all.dates[d])]
#   samp <- which(temp$dts == sample(sample_dates[which(sample_dates > temp$dts[300] & sample_dates < temp$dts[2508])], size = 1))
#   dates <- temp[(samp-300):(samp+100),]
#   discharge <- rep(0.1, time = length(dates$dts)) # create a list of non-disturbance discharges
#   discharge[match(temp$dts[samp], dates$dts)] <- 0.3 # from that list of dates from above, assign a disturbance discharge to that date
#   
#   # run model
#   out <- Amodel(discharge, dates, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
#   # create summary dataframe 
#   m <- mean.data.frame(out, burnin = 250, iteration = 2)
#   m <- cbind(m, discharge[250:402])
#   means[d] <- list(m)
#   short_abund <- means[[d]][which(discharge[250:402] == "0.3"), 2]
#   resil_abund <- means[[d]][(which(discharge[250:402] == "0.3") + 1):(which(discharge[250:402] == "0.3") + 6), 2]
#   long_abund <- means[[d]][(which(discharge[250:402] == "0.3") + 26):(which(discharge[250:402] == "0.3") + 36), 2]
#   #post <- seq(0, 10, by = 1)
#   #dat <- rep(d, times = length(post))
#   #jday_data <- as.data.frame(rbind(jday_data, cbind(dat, post, abund)))
#   short_mean[d] <- short_abund
#   resil_mean[d] <- mean(resil_abund)
#   long_mean[d] <- mean(long_abund)
# }
