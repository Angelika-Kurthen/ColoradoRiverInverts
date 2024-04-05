###################################
## Sp A Julian Date Pulse Response
###################################

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
#Calculate Julian Date Effect
all.dates <- unique(format(temp$dts, "%m-%d"))[order(unique(format(temp$dts, "%m-%d")))]
# loop to select a date from a Week-Month combo from each unique year
means <- list()
#jday_data <- data.frame(date = NULL, post_dist = NULL, abund = NULL)
jday_max <- vector()

short_mean <- vector()
resil_mean <- vector()
long_mean <- vector()

for (d in 1:length(all.dates)){ # 30 reps takes 60 mins
  sample_dates <- temp$dts[which(format(temp$dts, "%m-%d")== all.dates[d])]
  samp <- which(temp$dts == sample(sample_dates[which(sample_dates > temp$dts[300] & sample_dates < temp$dts[2508])], size = 1))
  dates <- temp[(samp-300):(samp+100),]
  discharge <- rep(0.1, time = length(dates$dts)) # create a list of non-disturbance discharges
  discharge[match(temp$dts[samp], dates$dts)] <- 1 # from that list of dates from above, assign a disturbance discharge to that date
  
  # run model
  out <- Amodel(discharge, dates, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
  # create summary dataframe 
  m <- mean.data.frame(out, burnin = 250, iteration = 2)
  m <- cbind(m, discharge[250:402])
  means[d] <- list(m)
  short_abund <- means[[d]][which(discharge[250:402] == 1.0), 2]
  resil_abund <- means[[d]][(which(discharge[250:402] == 1.0) + 1):(which(discharge[250:402] == 1) + 6), 2]
  long_abund <- means[[d]][(which(discharge[250:402] == 1.0) + 26):(which(discharge[250:402] == 1) + 36), 2]
  #post <- seq(0, 10, by = 1)
  #dat <- rep(d, times = length(post))
  #jday_data <- as.data.frame(rbind(jday_data, cbind(dat, post, abund)))
  short_mean[d] <- short_abund
  resil_mean[d] <- mean(resil_abund)
  long_mean[d] <- mean(long_abund)
}

ashort_df <- as.data.frame(cbind(short_mean, all.dates, rep("A", length(all.dates))))
aresil_df <- as.data.frame(cbind(resil_mean, all.dates, rep("A", length(all.dates))))
along_df <- as.data.frame(cbind(long_mean, all.dates, rep("A", length(all.dates))))

ashort_df$short_mean <- as.numeric(ashort_df$short_mean)
aresil_df$resil_mean <- as.numeric(aresil_df$resil_mean)
along_df$long_mean <- as.numeric(along_df$long_mean)

ashort_df$all.dates <- as.Date(ashort_df$all.dates, format = "%m-%d")
aresil_df$all.dates <- as.Date(aresil_df$all.dates, format = "%m-%d")
along_df$all.dates <- as.Date(along_df$all.dates, format = "%m-%d")

# ashort <- ggplot(data = short_df, aes(all.dates, short_mean/10000, group = 1))+
#   geom_line(size = 1, col = "#66CCEE")+
#   theme_bw()+
#   ylab("Sp A abundance relative to K")+
#   xlab("Date of one-time Pulse")+
#   scale_x_date(date_labels="%B", date_breaks  ="1 month")+
#   theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# 
# 
# aresil <- ggplot(data = resil_df, aes(all.dates, resil_mean/10000, group =1))+
#   geom_line(size = 1,  col = "#66CCEE")+
#   theme_bw()+
#   ylab("Sp A abundance relative to K")+
#   xlab("Date of one-time Pulse")+
#   scale_x_date(date_labels="%B", date_breaks  ="1 month")+
#   theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# 
# 
# 
# along <- ggplot(data = long_df, aes(all.dates, long_mean/10000, group = 1))+
#   geom_line(size = 1, col = "#66CCEE")+
#   theme_bw()+
#   ylab("Sp A abundance relative to K")+
#   xlab("Date of one-time Pulse")+
#   scale_x_date(date_labels="%B", date_breaks  ="1 month")+
#   theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# 

# can also do rolling mean
#l <- rollmean(long_df$long_mean, 30)
#plot(as.Date(all.dates[29:364], format = "%m-%d"), l, type = "l")
