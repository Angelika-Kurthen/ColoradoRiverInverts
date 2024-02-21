#############################
## Sp A Pulse Frequency 
############################
library(lubridate)
source("A_1sp_Model.R")
source("B_1sp_Model.R")
source("C_1sp_Model.R")
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


discharge <- rep(0.1, time = length(temp$dts))

A_out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
A_out <- mean.data.frame(A_out, burnin = 250, iteration = 2)
A_annual <- cbind.data.frame(temp$dts[250:last(A_out$timesteps)], A_out, rep("A", length(A_out$mean.abund)))
B_out <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
B_out <- mean.data.frame(B_out, burnin = 250, iteration = 2)
B_annual <- cbind.data.frame(temp$dts[250:last(B_out$timesteps)], B_out, rep("B", length(B_out$mean.abund)))
C_out <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
C_out <- mean.data.frame(C_out, burnin = 250, iteration = 2)
C_annual <- cbind.data.frame(temp$dts[250:last(C_out$timesteps)], C_out, rep("C", length(C_out$mean.abund)))
D_out <- Dmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
D_out <- mean.data.frame(D_out, burnin = 250, iteration = 2)
D_annual <- cbind.data.frame(temp$dts[250:last(C_out$timesteps)], D_out, rep("D", length(B_out$mean.abund)))


colnames(A_annual) <- c("Date", "timesteps", "Abundance", "sd", "se", "Taxa")
colnames(B_annual) <- c("Date", "timesteps", "Abundance", "sd", "se", "Taxa")
colnames(C_annual) <- c("Date", "timesteps", "Abundance", "sd", "se", "Taxa")
colnames(D_annual) <- c("Date", "timesteps", "Abundance", "sd", "se", "Taxa")


annual <- rbind(A_annual, B_annual, C_annual, D_annual)
annual$Date <- as.Date(annual$Date)
annual <- subset(annual, Date >= "2035-01-01" & Date <= "2035-12-31")
ggplot(data = annual, aes(x = Date, y  =Abundance/10000, color = Taxa))+
  geom_point(size = 1, alpha = 0.5)+
  geom_line(size = 1)+
  #stat_smooth(size= 1, span = 0.4, se =F)+
  xlab("Month")+
  ylab("Relativized Abundance")+
  scale_color_manual(name = "Taxa", values=c("#66CCEE", "#228833", "#EE6677", "#AA3377"))+
  scale_x_date(date_labels="%B", date_breaks  ="1 month")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle = 45, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))


  
selected_date <- temp$dts[temp$dts >= as.Date("2035-05-01") & temp$dts <= as.Date("2035-05-15")]
discharge[match(selected_date, temp$dts)] <- 0.3

A_pout <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
A_pout <- mean.data.frame(A_pout, burnin = 250, iteration = 2)
A_pulse <-  cbind.data.frame(temp$dts[250:last(A_pout$timesteps)], A_pout, rep("A", length(A_pout$mean.abund)))
B_pout <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
B_pout <- mean.data.frame(B_pout, burnin = 250, iteration = 2)
B_pulse <-  cbind.data.frame(temp$dts[250:last(B_pout$timesteps)], B_pout, rep("B", length(B_pout$mean.abund)))
C_pout <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
C_pout <- mean.data.frame(C_pout, burnin = 250, iteration = 2)
C_pulse <-  cbind.data.frame(temp$dts[250:last(C_pout$timesteps)], C_pout, rep("C", length(C_pout$mean.abund)))
D_pout <- Dmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
D_pout <- mean.data.frame(D_pout, burnin = 250, iteration = 2)
D_pulse <-  cbind.data.frame(temp$dts[250:last(D_pout$timesteps)], D_pout, rep("D", length(D_pout$mean.abund)))
colnames(A_pulse) <- c("Date", "timesteps", "Abundance", "sd", "se", "Taxa")
colnames(B_pulse) <- c("Date", "timesteps", "Abundance", "sd", "se", "Taxa")
colnames(C_pulse) <- c("Date", "timesteps", "Abundance", "sd", "se", "Taxa")
colnames(D_pulse) <- c("Date", "timesteps", "Abundance", "sd", "se", "Taxa")

pulse <- rbind(A_pulse, B_pulse, C_pulse, D_pulse)
pulse$Date <- as.Date(pulse$Date)
pulse <- subset(pulse, Date >= "2035-01-01" & Date <= "2035-12-31")
# 
# max(B_pulse$Abundance[(which(B_pulse$Date == selected_date)):(which(B_pulse$Date == selected_date) + 6)])
# 
# ggplot(data = pulse, aes(x = Date, y  =Abundance/10000, color = Taxa))+
#   geom_point(size = 1, alpha = 0.5)+
#   geom_line(size = 1)+
#   ylim(c(0,60))+
#   #stat_smooth(size= 1, span = 0.4, se =F)+
#   xlab("Month")+
#   ylab("Relativized Abundance")+
#   geom_vline(xintercept = as.numeric(as.Date("2035-03-13")), 
#              color = "black", 
#              lwd = 1,
#              linetype = "dotted") +
#   scale_color_manual(name = "Taxa", values=c("#66CCEE", "#228833", "#EE6677", "#AA3377"))+
#   scale_x_date(date_labels="%B", date_breaks  ="1 month")+
#   theme_bw()+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle = 45, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# 
# 
# 
# 
# annual1 <- subset(annual, Date >= "2035-01-01" & Date <= "2035-12-31" & Taxa == "A")
# annual2 <- subset(annual, Date >= "2036-01-01" & Date <= "2036-12-31" & Taxa == "A")
# 
# annual3 <- subset(annual, Date >= "2035-01-01" & Date <= "2035-12-31")
# annual4 <- subset(annual, Date >= "2035-01-01" & Date <= "2035-12-31")
# annual5 <- subset(annual, Date >= "2035-01-01" & Date <= "2035-12-31")
# annual6 <- subset(annual, Date >= "2035-01-01" & Date <= "2035-12-31")
# annual7 <- subset(annual, Date >= "2035-01-01" & Date <= "2035-12-31")
# annual8 <- subset(annual, Date >= "2035-01-01" & Date <= "2035-12-31")
# annual9 <- subset(annual, Date >= "2035-01-01" & Date <= "2035-12-31")
# annual10 <- subset(annual, Date >= "2035-01-01" & Date <= "2035-12-31")
# annual11<- subset(annual, Date >= "2035-01-01" & Date <= "2035-12-31")
# annual12 <- subset(annual, Date >= "2035-01-01" & Date <= "2035-12-31")
# annual13 <- subset(annual, Date >= "2035-01-01" & Date <= "2035-12-31")
# annual14 <- subset(annual, Date >= "2049-01-01" & Date <= "2049-12-31" & Taxa == "A")
# annual <- subset(annual, Date >= "2035-01-01" & Date <= "2035-12-31")
# annual <- subset(annual, Date >= "2035-01-01" & Date <= "2035-12-31")
# 
discharge <- rep(0.1, time = length(temp$dts))

A_out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature))
A_annual <- cbind.data.frame(temp$dts[250:length(A_out)], A_out[250:length(A_out)], rep("A", length(A_out)-249))
B_out <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature))
B_annual <- cbind.data.frame(temp$dts[250:length(B_out)], B_out[250:length(B_out)], rep("B", length(B_out)-249))
C_out <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature))
C_annual <- cbind.data.frame(temp$dts[250:length(C_out)], C_out[250:length(C_out)], rep("C", length(C_out)-249))
D_out <- Dmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature))
D_annual <- cbind.data.frame(temp$dts[250:length(D_out)], D_out[250:length(D_out)], rep("D", length(D_out)-249))


colnames(A_annual) <- c("Date", "Adults", "Taxa")
colnames(B_annual) <- c("Date", "Adults", "Taxa")
colnames(C_annual) <- c("Date", "Adults", "Taxa")
colnames(D_annual) <- c("Date", "Adults", "Taxa")


annual <- rbind(A_annual, B_annual, C_annual, D_annual)
annual$Date <- as.Date(annual$Date)
annual <- subset(annual, Date >= "2035-01-01" & Date <= "2035-12-31")
ggplot(data = annual, aes(x = Date, y  =Adults, color = Taxa))+
  geom_point(size = 1, alpha = 0.5)+
  geom_line(size = 1)+
  ylim(c(0,20000))+
  #stat_smooth(size= 1, span = 0.4, se =F)+
  xlab("Month")+
  ylab("Relativized Abundance")+
  scale_color_manual(name = "Taxa", values=c("#66CCEE", "#228833", "#EE6677", "#AA3377"))+
  scale_x_date(date_labels="%B", date_breaks  ="1 month")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle = 45, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))



selected_date <- temp$dts[temp$dts >= as.Date("2035-05-01") & temp$dts <= as.Date("2035-05-15")]
discharge[match(selected_date, temp$dts)] <- 0.3

A_out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature))
A_annual <- cbind.data.frame(temp$dts[250:length(A_out)], A_out[250:length(A_out)], rep("A", length(A_out)-249))
B_out <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature))
B_annual <- cbind.data.frame(temp$dts[250:length(B_out)], B_out[250:length(B_out)], rep("B", length(B_out)-249))
C_out <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature))
C_annual <- cbind.data.frame(temp$dts[250:length(C_out)], C_out[250:length(C_out)], rep("C", length(C_out)-249))
D_out <- Dmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature))
D_annual <- cbind.data.frame(temp$dts[250:length(D_out)], D_out[250:length(D_out)], rep("D", length(D_out)-249))
colnames(A_annual) <- c("Date", "Adults", "Taxa")
colnames(B_annual) <- c("Date", "Adults", "Taxa")
colnames(C_annual) <- c("Date", "Adults", "Taxa")
colnames(D_annual) <- c("Date", "Adults", "Taxa")
pulse <- rbind(A_annual, B_annual, C_annual, D_annual)
pulse$Date <- as.Date(pulse$Date)
pulse <- subset(pulse, Date >= "2035-01-01" & Date <= "2035-12-31")


ggplot(data = pulse, aes(x = Date, y  =Adults, color = Taxa))+
  geom_point(size = 1, alpha = 0.5)+
  geom_line(size = 1)+
  ylim(c(0,20000))+
  #stat_smooth(size= 1, span = 0.4, se =F)+
  xlab("Month")+
  ylab("Relativized Abundance")+
  geom_vline(xintercept = as.numeric(as.Date("2035-05-08")),
             color = "black",
             lwd = 1,
             linetype = "dotted") +
  scale_color_manual(name = "Taxa", values=c("#66CCEE", "#228833", "#EE6677", "#AA3377"))+
  scale_x_date(date_labels="%B", date_breaks  ="1 month")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle = 45, size = 12.5),
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
