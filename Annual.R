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


discharge <- rep(0.1, time = length(temp$dts))

A_out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature))
A_annual <- cbind.data.frame(as.Date(temp$dts[340:365]), rowSums(A_out[341:366, ]), rep("A", length(temp$dts[340:365])))
B_out <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature))
B_annual <- cbind.data.frame(as.Date(temp$dts[340:365]), rowSums(B_out[341:366, ]), rep("B", length(temp$dts[340:365])))
C_out <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature))
C_annual <- cbind.data.frame(as.Date(temp$dts[340:365]), rowSums(C_out[341:366, ]), rep("C", length(temp$dts[340:365])))
D_out <- Dmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature))
D_annual <- cbind.data.frame(as.Date(temp$dts[340:365]), rowSums(D_out[341:366, ]), rep("D", length(temp$dts[340:365])))
colnames(A_annual) <- c("Date", "Abundance", "Taxa")
colnames(B_annual) <- c("Date", "Abundance", "Taxa")
colnames(C_annual) <- c("Date", "Abundance", "Taxa")
colnames(D_annual) <- c("Date", "Abundance", "Taxa")


annual <- rbind(A_annual, B_annual, C_annual, D_annual)

ggplot(data = annual, aes(x = Date, y  =Abundance/10000, color = Taxa))+
  geom_point(size = 1, alpha = 0.5)+
  stat_smooth(size= 1, span = 0.4, se =F)+
  xlab("Month")+
  ylab("Relativized Abundance")+
  scale_color_manual(name = "Taxa", values=c("#66CCEE", "#228833", "#EE6677", "#AA3377"))+
  scale_x_date(date_labels="%B", date_breaks  ="1 month")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle = 45, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))


  
selected_date <- temp$dts[temp$dts >= as.Date("2035-04-15") & temp$dts <= as.Date("2035-04-30")]
discharge[match(selected_date, temp$dts)] <- 0.3
A_out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature))
A_pulse <- cbind.data.frame(as.Date(temp$dts[340:365]), rowSums(A_out[341:366, ]), rep("A", length(temp$dts[340:365])))
B_out <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature))
B_pulse <- cbind.data.frame(as.Date(temp$dts[340:365]), rowSums(B_out[341:366, ]), rep("B", length(temp$dts[340:365])))
C_out <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature))
C_pulse <- cbind.data.frame(as.Date(temp$dts[340:365]), rowSums(C_out[341:366, ]), rep("C", length(temp$dts[340:365])))
D_out <- Dmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature))
D_pulse <- cbind.data.frame(as.Date(temp$dts[340:365]), rowSums(D_out[341:366, ]), rep("D", length(temp$dts[340:365])))
colnames(A_pulse) <- c("Date", "Abundance", "Taxa")
colnames(B_pulse) <- c("Date", "Abundance", "Taxa")
colnames(C_pulse) <- c("Date", "Abundance", "Taxa")
colnames(D_pulse) <- c("Date", "Abundance", "Taxa")

pulse <- rbind(A_pulse, B_pulse, C_pulse, D_pulse)

ggplot(data = pulse, aes(x = Date, y  =Abundance/10000, color = Taxa))+
  geom_point(size = 1, alpha = 0.5)+
  #geom_line(size = 1)+
  stat_smooth(size= 1, span = 0.4, se =F)+
  xlab("Month")+
  ylab("Relativized Abundance")+
  geom_vline(xintercept = as.numeric(as.Date("2035-04-24")), 
             color = "black", 
             lwd = 1,
             linetype = "dotted") +
  scale_color_manual(name = "Taxa", values=c("#66CCEE", "#228833", "#EE6677", "#AA3377"))+
  scale_x_date(date_labels="%B", date_breaks  ="1 month")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle = 45, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))



