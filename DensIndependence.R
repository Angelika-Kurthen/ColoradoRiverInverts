####################
# Code for density independent runs
####################

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


Temperature <-  rep(14, times = length(Time))

temp <- as.data.frame(cbind(Time, Day, Temperature))
temp$Day <- as.Date(temp$Day, origin= "1970-01-01")
colnames(temp) <- c("Time", "Date", "Temperature")
temp <- TimestepTemperature(temp)
temp <- temp[c(1,3)]


discharge <- rep(0.1, time = length(temp$dts))

A_out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), dens.dep = F)
A_out <- mean.data.frame(A_out, burnin = 250, iteration = 2)
A_dens.ind <- cbind.data.frame(temp$dts[250:last(A_out$timesteps)], A_out, rep("A", length(A_out$mean.abund)))
B_out <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), dens.dep = F)
B_out <- mean.data.frame(B_out, burnin = 250, iteration = 2)
B_dens.ind <- cbind.data.frame(temp$dts[250:last(B_out$timesteps)], B_out, rep("B", length(B_out$mean.abund)))
C_out <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), dens.dep = F)
C_out <- mean.data.frame(C_out, burnin = 250, iteration = 2)
C_dens.ind <- cbind.data.frame(temp$dts[250:last(C_out$timesteps)], C_out, rep("C", length(C_out$mean.abund)))
D_out <- Dmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), dens.dep = F)
D_out <- mean.data.frame(D_out, burnin = 250, iteration = 2)
D_dens.ind <- cbind.data.frame(temp$dts[250:last(C_out$timesteps)], D_out, rep("D", length(B_out$mean.abund)))


colnames(A_dens.ind)<- c("Date", "timesteps", "Abundance", "sd", "se", "Taxa")
colnames(B_dens.ind) <- c("Date", "timesteps", "Abundance", "sd", "se", "Taxa")
colnames(C_dens.ind) <- c("Date", "timesteps", "Abundance", "sd", "se", "Taxa")
colnames(D_dens.ind) <- c("Date", "timesteps", "Abundance", "sd", "se", "Taxa")


dens.ind <- rbind(A_dens.ind, B_dens.ind, C_dens.ind, D_dens.ind)
dens.ind$Date <- as.Date(dens.ind$Date)
dens.ind <- subset(dens.ind, Date >= "2030-01-01" & Date <= "2033-12-31")
# ggplot(data = dens.ind, aes(x = Date, y  =log(Abundance), color = Taxa))+
#   geom_point(size = 1, alpha = 0.5)+
#   geom_line(size = 1)+
#   #stat_smooth(size= 1, span = 0.4, se =F)+
#   xlab("Month")+
#   ylab("Log Abundance")+
#   scale_color_manual(name = "Taxa", values=c("#66CCEE", "#228833", "#EE6677", "#AA3377"))+
#   scale_x_date(date_labels="%B", date_breaks  ="4 month")+
#   theme_bw()+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle = 45, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
