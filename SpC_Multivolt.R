######################
# Multivoltinism Sp C
######################

source("C_1sp_Model.R")
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


discharge <- rep(0.1, time = length(temp$dts))

# run model under different temp regimes

runs <- c(-5, -2.5, 0, 5)

# for adults

for (i in 1:length(runs)){
  temp$Temperature <- temp$Temperature + runs[i]
  out <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
  m <- rowMeans(out)
  m <- cbind.data.frame(temp$dts, m[-1], rep(i, length(temp$dts)))
  assign(paste0("m",i), m)
  temp$Temperature <- temp$Temperature - runs[i]
}

mlist <- rbind(m1, m2, m3, m4)
colnames(mlist) <- c("Date", "Abund", "MeanTemp")
oneyear <- mlist[which(mlist$Date >= "2035-01-01" & mlist$Date <= "2035-12-31"), ]
peaks <- oneyear[c(19, 41, 47, 64, 69, 75, 82, 89, 94, 100), ]
arrows <- tibble(
  x1 = peaks$Date,
  x2 = peaks$Date,
  y1 = c(5200, 6300, 7700, 4300, 7200, 6200, 4800, 9800, 7500, 26100), 
  y2 = c(3200, 4300, 5700, 2300, 5200, 4200,2800, 7800, 5500, 24100)
)
arrowcols <- c("#4477AA", "#EE6677","#EE6677", "#228833","#228833","#228833","#CCBB44","#CCBB44","#CCBB44", "#CCBB44")
arrows$x1 <- as.POSIXct(arrows$x1, format = "%Y-%m-%d")
arrows$x1 <- as.Date(arrows$x1)
arrows$x2 <- as.POSIXct(arrows$x2, fomrat = "%Y-%m_%d")
arrows$x2 <- as.Date(arrows$x2)
oneyear$Date <- as.Date(oneyear$Date)
#create list of peak points
ggplot(data = oneyear, aes(x = Date, y = Abund, group = as.factor(MeanTemp), color = as.factor(MeanTemp)))+
        geom_line(size = 1 )+
        scale_color_manual(name = "Mean Temperature (C)", labels = c("9", "11.5", "14", "19"), values = c("#4477AA", "#EE6677", "#228833", "#CCBB44"))+
        annotate("segment", x = arrows$x1, y = arrows$y1, xend = arrows$x2, yend = arrows$y2, color = arrowcols,
            arrow = arrow(type = "closed", length = unit(0.02, "npc")))+
        ylab("Adult Abundance") + 
        ylim(c(0, 26100))+
        theme_bw()+
        scale_x_date(date_labels="%B", date_breaks  ="1 month")+
        theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle=45, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

# for all 
# for (i in 1:length(runs)){
#   temp$Temperature <- temp$Temperature + runs[i]
#   out <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
#   m <- mean.data.frame(out, burnin = 250, iteration = 2)
#   m <- cbind.data.frame(temp$dts[250:last(m$timesteps)], m, rep(i, length(m$mean.abund)))
#   assign(paste0("m",i), m)
#   temp$Temperature <- temp$Temperature - runs[i]
# }
# 
# mlist <- rbind(m1, m2, m3, m4)
# colnames(mlist) <- c("Date", "Timesteps", "Abund", "SD", "SE", "MeanTemp")
# oneyear <- mlist[which(mlist$Date >= "2035-01-01" & mlist$Date <= "2035-12-31"), ]
# index <- seq(1:length(oneyear$Date))
# oneyear <- cbind(oneyear, index)
# peaks <- oneyear[c(18, 23, 27, 41, 48, 55, 62, 67, 73, 81, 91, 94, 99), ]
# arrows <- tibble(
#   x1 = peaks$Date,
#   x2 = peaks$Date,
#   y1 = c(15, 29, 15, 12, 31, 21, 11, 15, 37, 18, 23, 21, 34), 
#   y2 = c(13, 27, 13, 10, 29, 19, 9, 13, 35, 16, 21, 19, 32)
# )
# arrowcols <- c("#4477AA","#4477AA", "#EE6677","#EE6677","#EE6677", "#228833","#228833","#228833","#228833","#CCBB44","#CCBB44","#CCBB44", "#CCBB44")
# arrows$x1 <- as.POSIXct(arrows$x1, format = "%Y-%m-%d")
# arrows$x1 <- as.Date(arrows$x1)
# arrows$x2 <- as.POSIXct(arrows$x2, fomrat = "%Y-%m_%d")
# arrows$x2 <- as.Date(arrows$x2)
# oneyear$Date <- as.Date(oneyear$Date)
# #create list of peak points
# ggplot(data = oneyear, aes(x = Date,  y = Abund/10000, group = as.factor(MeanTemp), color = as.factor(MeanTemp)))+
#   geom_line(size = 1 )+
#   scale_color_manual(name = "Mean Temperature (C)", labels = c("9", "14", "19", "21"), values = c("#4477AA", "#EE6677", "#228833", "#CCBB44"))+
#   annotate("segment", x = arrows$x1, y = arrows$y1, xend = arrows$x2, yend = arrows$y2, color = arrowcols,
#            arrow = arrow(type = "closed", length = unit(0.02, "npc")))+
#   ylab(" Abundance") +   
#   theme_bw()+
#   scale_x_date(date_labels="%B", date_breaks  ="1 month")+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle=45, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# 
source("C_1sp_Model.R")
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


discharge <- rep(0.1, time = length(temp$dts))

# run model under different temp regimes

runs <- c(-5,  5)

# for adults

for (i in 1:length(runs)){
  temp$Temperature <- temp$Temperature + runs[i]
  out <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
  m <- rowMeans(out)
  m <- cbind.data.frame(temp$dts, m[-1], rep(i, length(temp$dts)))
  assign(paste0("m",i), m)
  temp$Temperature <- temp$Temperature - runs[i]
}

mlist <- rbind(m1, m2)
colnames(mlist) <- c("Date", "Abund", "MeanTemp")


oneyear <- mlist[which(mlist$Date >= "2035-01-01" & mlist$Date <= "2035-12-31"), ]
peaks <- oneyear[c(17, 22, 28, 35, 40, 46), ]
arrows <- tibble(
  x1 = peaks$Date,
  x2 = peaks$Date,
  y1 = c(4.2, 5, 5.3000, 2.9000, 3.5000, 10.700), 
  y2 = c(3.8, 4.6000, 4.9000, 2.5000, 3.1000, 10.3000)
)
arrowcols <- c("#4477AA", "#4477AA", "#EE6677","#EE6677", "#EE6677","#EE6677")
arrows$x1 <- as.POSIXct(arrows$x1, format = "%Y-%m-%d")
arrows$x1 <- as.Date(arrows$x1)
arrows$x2 <- as.POSIXct(arrows$x2, fomrat = "%Y-%m_%d")
arrows$x2 <- as.Date(arrows$x2)
oneyear$Date <- as.Date(oneyear$Date)
#create list of peak points
ggplot(data = oneyear, aes(x = Date, y = Abund/10000, group = as.factor(MeanTemp), color = as.factor(MeanTemp)))+
  geom_line(size = 1 )+
  scale_color_manual(name = "Mean Temperature (C)", labels = c("9", "19"), values = c("#4477AA", "#EE6677"))+
  annotate("segment", x = arrows$x1, y = arrows$y1, xend = arrows$x2, yend = arrows$y2, color = arrowcols,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")))+
  ylab("Relativized Abundance") + 
  theme_bw()+
  scale_x_date(date_labels="%B", date_breaks  ="1 month")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle=45, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
