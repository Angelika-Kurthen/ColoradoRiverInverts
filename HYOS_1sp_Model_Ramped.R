##########################
# HYOS 1 sp model
###########################


library(purrr)
library(tidyverse)
library(lubridate)
library(plyr)
library(dplyr)
library(ggplot2)
# data retrieval tool from USGS
library(dataRetrieval)


source("HYOS_1sp.R")

source("1spFunctions.R")
source("HYOS_1sp.R")


# load Water Temperature data from above Diamond Creek Confluence (RM226)
temp <- read.delim("CRaboveDC_Temp.tsv", header=T)
colnames(temp) <- c("Date", "Temperature")
temp$Date <- as.Date(temp$Date, format = "%Y-%m-%d")
temp <- aggregate(temp$Temperature, by = list(temp$Date), FUN = mean)
colnames(temp) <- c("Date", "Temperature")
temp$Temperature[which(temp$Temperature < 0)] <- NA
temp <- na.omit(temp)
#want last 10 years of data (data ends in 2023)
temps <- temp[which(temp$Date >= "2013-01-01"),]
tempslist <- average.yearly.temp(temps, "Temperature", "Date")

years <- seq(23, 30, by = 1)
temps <- average.yearly.temp(temps, "Temperature", "Date")
# add average yearly temps to 
# add new years so we can progress forwards in time, then concatonate
for(year in years){
  year(tempslist$dts) <- (2000 + year)
  temps <- rbind(temps, tempslist)
}
temps <- temps[-(1:26),-c(1,4, 5, 6)]

# but also want those average temps to line up seamlessly with out emprical data
# before 2000, lots of missing data
temp <- temp[which(temp$Date >= "2000-01-05"),]
temp <- TimestepTemperature(temp)
lastval <- last(temp)

colnames(temp) <- c("dts", "Temperature")
temperature <- rbind(temp, temps)
plot(temperature$dts, temperature$Temperature)



# want to check what water temp was when lees ferry temp over 15.5 and 13
lftemp <- readNWISdv("09380000", "00010", "2014-05-01", "2024-05-01")
lftemp <- lftemp[, c(3,4)]
colnames(lftemp) <- c("Date", "Temperature")
lftemps <- (lftemp$Date[which(lftemp$Temperature == 15.5)]) # in 2021, it was from 09-03 and 09-30; in 2022 it was 06-09, 06-16, and 11-03; and in 2023 it was 06-26, 06-27, and 11-1
lftempscold <- lftemp$Date[which(lftemp$Temperature == 13)]
# match the temperatures and calculate mean when LF is at that temp, what is RM 226?
# these values will be what we put into the temperature model, while flow will be the same

value15.5 <- mean(temps$Temperature[temps$Date %in% as.Date(lftemps)])
value13 <- mean(temps$Temperature[temps$Date %in% as.Date(lftempscold)])

# load discharge data from above Diamond Creek Confluence (RM226)
discharge <- readNWISdv("09404200", "00060", "2000-01-07", "2024-01-01")
discharge <- TimestepDischarge(discharge, 85000)

flow10yr <- readNWISdv("09404200", "00060", "2014-01-01", "2024-01-01")
flows10yr <- average.yearly.flows(flowdata = flow10yr, "X_00060_00003", "Date")

# do the same with out yearly flows
flow.df <- as.data.frame(cbind(temps$dts[which(year(temps$dts) > 2023)], rep(flows$Discharge/85000, times = 101)))
colnames(flow.df) <- c("dts", "flow.magnitude")
flow.df$dts <- as.Date(flow.df$dts)
rbind(discharge, flow.df)

# Increase in summer temps
# means <- vector()
# increase <- seq(0, 5, by = 0.5)
# for (inc in 1:length(increase)){
#   temps$Temperature[which(month(temps$dts) == 6 | month(temps$dts) == 7 | month(temps$dts) == 8)] <- temps$Temperature[which(month(temps$dts) == 6 | month(temps$dts) == 7 | month(temps$dts) == 8)] + increase[inc]
#   out <- HYOSmodel(flow.data = flow.df$flow.magnitude ,temp.data = temps, disturbanceK = 40000, baselineK = 5000, Qmin = 0.2, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
#   HYOS.df <- mean.data.frame(out, 250, 9)
#   means[inc] <- mean(HYOS.df$mean.abund)
#   temps$Temperature[which(month(temps$dts) == 6 | month(temps$dts) == 7 | month(temps$dts) == 8)] <- temps$Temperature[which(month(temps$dts) == 6 | month(temps$dts) == 7 | month(temps$dts) == 8)] - increase[inc]
# }
# 
# 
# summer.HYOS <- as.data.frame(cbind(increase, means))
# 
# 
# ggplot(data = summer.HYOS, aes(x = increase, y = means))+
#   geom_point()+
#   geom_line()+
#   xlab("")

# sensitivity to hydropeaking
# 
# hydropeak <- seq(0, 0.5, by = 0.025)
# means <- vector()
# for (hydr in 1:length(hydropeak)){
#   out <- HYOSmodel(flow.data = flow.df$flow.magnitude, temp.data = temps, disturbanceK = 40000, baselineK = 5000, Qmin = 0.2, extinct = 50, iteration = 9, peaklist = hydropeak[hydr], peakeach = length(temps$Temperature))
#   HYOS.df <- mean.data.frame(out, 250, 9)
#   means[hydr] <- mean(HYOS.df$mean.abund)
# }
# 
# 
# hydropeak <- as.data.frame(cbind(hydropeak, means))
# 
# ggplot(data = hydropeak, aes(x = hydropeak, y = means))+
#   geom_point()+
#   geom_line()+
#   theme_bw()+
#   labs(x = "Hydropeaking Index", y = "Average Annual Abundance")+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"), plot.margin = margin(5,5,5,20))



#BAU scenario
# regular temps and regular flows - we will look at year 2024 thru year 
bauflow <- flow.df
bauflow[which(month(bauflow$dts) == 11 & day(bauflow$dts) == 9),] <- 0.18
#bauflow$flow.magnitude[which(bauflow$dts == "2024-11-09" | bauflow$dts == "2025-11-09" | bauflow$dts == "2026-11-09" | bauflow$dts == "2027-11-09")] <- 0.28
out <- HYOSmodel(flow.data = bauflow$flow.magnitude, temp.data = temperature, disturbanceK = 40000, baselineK = 5000, Qmin = 0.2, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temperature$Temperature))
HYOS.df <- mean.data.frame(out, 250, 9)
HYOS.df <- as.data.frame(cbind(temperature$dts[249:length(temperature$dts)], HYOS.df$mean.abund))
colnames(HYOS.df) <- c("Date", "Abundance")
HYOS.df$Date <- as.Date(as.POSIXct(HYOS.df$Date))
# 2024 - 2027

# HFE of around 37000 cfs allowed Oct 1 through Nov 30. Will split the difference with 11-09 HFE day
# assume 7 days at high cfs and 7 days at baseline, so magnitude = 0.28 

bau <- HYOS.df[which(HYOS.df$Date >= "2024-01-01" & HYOS.df$Date < "2028-01-01"),]
plot(bau$Date, bau$Abundance, type = "both")
lines(temps$dts, temps$Temperature, col = "red")
lines(bauflow$dts, bauflow$flow.magnitude*100)
# cool mix scenario
coolmix <- temps
# all temperatures are supposed to be 15.5 or below - this won't affect Lees Ferry because we don't really have super warm average temps (BUT we did in 2023)
coolmix[coolmix$Temperature > value15.5,] <- value15.5

# will use bauflow, since HFE still allowed
out <- HYOSmodel(flow.data = bauflow$flow.magnitude, temp.data = coolmix, disturbanceK = 40000, baselineK = 5000, Qmin = 0.2, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
HYOS.df <- mean.data.frame(out, 250, 9)
HYOS.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], HYOS.df$mean.abund))
colnames(HYOS.df) <- c("Date", "Abundance")
HYOS.df$Date <- as.Date(HYOS.df$Date)

cool <- HYOS.df[which(HYOS.df$Date >= "2024-01-01" & HYOS.df$Date < "2028-01-01"),]
cool2 <- HYOS.df[which(HYOS.df$Date >= "2010-01-01" & HYOS.df$Date < "2014-01-01"),]

# cool mix + flow spike (suggested is 32000) in Late May Early June
# average daily flow is around 11000 cfs per fortnight, based on proposed hydrograph, with a 3 day increase to 32000 (this average flow would be 12500)... not such a large jump. Could 
# so we have tw0 options
# a) one time increase using mean (0.1823529) on 5-25-2024, 2 in 2025, 2 in 2025, 0 in 2027 (pg 3-10)
# b) one time increase using max (0.3764706) on 5-25-2025,2 in 2025, 2 in 2025, 0 in 2027 (pg 3-10)


# will do scenario based on suggested # of flow spikes, HFEs, and proposed hydrograph
coolflow <- bauflow
# add in flow spikes
coolflow$flow.magnitude[which(coolflow$dts == "2024-05-25")] <- 0.18
coolflow$flow.magnitude[which(coolflow$dts == "2025-05-25" | coolflow$dts == "2025-06-08")] <- 0.18
coolflow$flow.magnitude[which(coolflow$dts == "2026-05-25" | coolflow$dts == "2026-06-08")] <- 0.18


out <- HYOSmodel(flow.data = coolflow$flow.magnitude, temp.data = coolmix, disturbanceK = 40000, baselineK = 5000, Qmin = 0.2, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
HYOS.df <- mean.data.frame(out, 250, 9)
HYOS.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], HYOS.df$mean.abund))
colnames(HYOS.df) <- c("Date", "Abundance")
HYOS.df$Date <- as.Date(HYOS.df$Date)

coolspike <- HYOS.df[which(HYOS.df$Date >= "2024-01-01" & HYOS.df$Date < "2028-01-01"),]
plot(coolspike$Date, coolspike$Abundance, type = "l" )

# cold schock alternative
# keep water below 13 C
coldtemps <- temps
coldtemps$Temperature[which(coldtemps$Temperature > value13)] <-value13
# use bauflows

out <- HYOSmodel(flow.data = bauflow$flow.magnitude, temp.data = coldtemps, disturbanceK = 40000, baselineK = 5000, Qmin = 0.2, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
HYOS.df <- mean.data.frame(out, 250, 9)
HYOS.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], HYOS.df$mean.abund))
colnames(HYOS.df) <- c("Date", "Abundance")
HYOS.df$Date <- as.Date(HYOS.df$Date)

coldshock <- HYOS.df[which(HYOS.df$Date >= "2024-01-01" & HYOS.df$Date < "2028-01-01"),]
plot(coldshock$Date, coldshock$Abundance, type = "l" )


# scenario 2
# coolflow$flow.magnitude[which(coolflow$dts == "2024-05-25")] <- 0.376
# coolflow$flow.magnitude[which(coolflow$dts == "2025-05-25" | coolflow$dts == "2025-06-08")] <- 0.376
# coolflow$flow.magnitude[which(coolflow$dts == "2026-05-25" | coolflow$dts == "2026-06-08")] <- 0.376
# 
# 
# out <- HYOSmodel(flow.data = coolflow$flow.magnitude, temp.data = coolmix, disturbanceK = 40000, baselineK = 5000, Qmin = 0.1, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
# HYOS.df <- mean.data.frame(out, 250, 9)
# HYOS.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], HYOS.df$mean.abund))
# colnames(HYOS.df) <- c("Date", "Abundance")
# HYOS.df$Date <- as.Date(HYOS.df$Date)
# 
# coolspike <- HYOS.df[which(HYOS.df$Date >= "2024-01-01" & HYOS.df$Date < "2028-01-01"),]
# plot(coolspike$Date, coolspike$Abundance, type = "l" )

# cold schock alternative
# keep water below 13 C
coldtemps <- temps
coldtemps$Temperature[which(coldtemps$Temperature > 13)] <- 13
# use bauflows

out <- HYOSmodel(flow.data = bauflow$flow.magnitude, temp.data = coldtemps, disturbanceK = 40000, baselineK = 5000, Qmin = 0.2, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
HYOS.df <- mean.data.frame(out, 250, 9)
HYOS.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], HYOS.df$mean.abund))
colnames(HYOS.df) <- c("Date", "Abundance")
HYOS.df$Date <- as.Date(HYOS.df$Date)

coldshock <- HYOS.df[which(HYOS.df$Date >= "2024-01-01" & HYOS.df$Date < "2028-01-01"),]
plot(coldshock$Date, coldshock$Abundance, type = "l" )

# cold schock with flow spike
# we will add flow spikes on same dates as for cool spike alt
# magnitude will be either a) mean given hydrograph or b) max
# mean = (9*11000)+(3*14000)+(2*32000)/14/85000 = 0.1722689
# max = 0.376
coldflow <- bauflow
# add in flow spikes
coldflow$flow.magnitude[which(coldflow$dts == "2024-05-25")] <- 0.17
coldflow$flow.magnitude[which(coldflow$dts == "2025-05-25" | coldflow$dts == "2025-06-08")] <- 0.17
coldflow$flow.magnitude[which(coldflow$dts == "2026-05-25" | coldflow$dts == "2026-06-08")] <- 0.17

out <- HYOSmodel(flow.data = coldflow$flow.magnitude, temp.data = coldtemps, disturbanceK = 40000, baselineK = 5000, Qmin = 0.2, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
HYOS.df <- mean.data.frame(out, 250, 9)
HYOS.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], HYOS.df$mean.abund))
colnames(HYOS.df) <- c("Date", "Abundance")
HYOS.df$Date <- as.Date(HYOS.df$Date)

coldspike <- HYOS.df[which(HYOS.df$Date >= "2024-01-01" & HYOS.df$Date < "2028-01-01"),]
plot(coldspike$Date, coldspike$Abundance, type = "l" )

# scenario 2
# coldflow$flow.magnitude[which(coldflow$dts == "2024-05-25")] <- 0.376
# coldflow$flow.magnitude[which(coldflow$dts == "2025-05-25" | coldflow$dts == "2025-06-08")] <- 0.376
# coldflow$flow.magnitude[which(coldflow$dts == "2026-05-25" | coldflow$dts == "2026-06-08")] <- 0.376
# 
# out <- HYOSmodel(flow.data = coldflow$flow.magnitude, temp.data = coldtemps, disturbanceK = 40000, baselineK = 5000, Qmin = 0.1, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
# HYOS.df <- mean.data.frame(out, 250, 9)
# HYOS.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], HYOS.df$mean.abund))
# colnames(HYOS.df) <- c("Date", "Abundance")
# HYOS.df$Date <- as.Date(HYOS.df$Date)
# 
# coldshock <- HYOS.df[which(HYOS.df$Date >= "2024-01-01" & HYOS.df$Date < "2028-01-01"),]
# plot(coldshock$Date, coldshock$Abundance, type = "l" )

# nonbypass drop
# this one is a little tricky
# mean would be 0.1416807
# but we can divide up that is split so that one fortnight has 3 (2000s) and next has 3 (27300)
# so one timestep 0.12
# and the next is 0.16
# other option is to just do max 0.32
# seems like June would be the time for this


nonbypass <- bauflow
nonbypass$flow.magnitude[which(nonbypass$dts == "2024-06-08" )] <- 0.12
nonbypass$flow.magnitude[which(nonbypass$dts == "2025-06-08" )] <- 0.12
nonbypass$flow.magnitude[which(nonbypass$dts == "2026-06-08" )] <- 0.12
nonbypass$flow.magnitude[which(nonbypass$dts == "2027-06-08" )] <- 0.12
nonbypass$flow.magnitude[which(nonbypass$dts == "2024-06-22" )] <- 0.16
nonbypass$flow.magnitude[which(nonbypass$dts == "2025-06-22" )] <- 0.16
nonbypass$flow.magnitude[which(nonbypass$dts == "2026-06-22" )] <- 0.16
nonbypass$flow.magnitude[which(nonbypass$dts == "2027-06-22" )] <- 0.16

out <- HYOSmodel(flow.data = nonbypass$flow.magnitude, temp.data = temps, disturbanceK = 40000, baselineK = 5000, Qmin = 0.2, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
HYOS.df <- mean.data.frame(out, 250, 9)
HYOS.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], HYOS.df$mean.abund))
colnames(HYOS.df) <- c("Date", "Abundance")
HYOS.df$Date <- as.Date(HYOS.df$Date)

nbypass <- HYOS.df[which(HYOS.df$Date >= "2024-01-01" & HYOS.df$Date < "2028-01-01"),]
plot(nbypass$Date, nbypass$Abundance, type = "l" )

HYOS.scenarios <- as.data.frame(cbind(bau, cool$Abundance, coolspike$Abundance, coldshock$Abundance, coldspike$Abundance, nbypass$Abundance))


colors <- c("No Action" = "#FF7F00", "Cool Mix" = "#A6CEE3", "Cool Mix with Flow Spike" = "#1F78B4")
ggplot(data = HYOS.scenarios, aes(x = Date, y = Abundance))+
  geom_line(aes(color = "No Action" ), alpha = 0.8, size = 1)+  
  geom_line(aes(y = cool$Abundance,  color = "Cool Mix"),linetype = "dashed", alpha = 0.8, size = 1)+
  geom_line(aes(y = coolspike$Abundance, color = "Cool Mix with Flow Spike"), size = 1, linetype = "dashed", alpha = 0.8)+
  scale_color_manual(values = colors)+
  theme_bw()+
  labs(x = "Year", y = "Abundance", color = "Scenario")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"), plot.margin = margin(5,5,5,20))

colors <- c("No Action" = "#FF7F00", "Cold Shock" = "#CAB2D6", "Cold Shock with Flow Spike" = "#6A3D9A")
ggplot(data = HYOS.scenarios, aes(x = Date, y = Abundance))+
  geom_line(aes(color = "No Action" ), alpha = 0.8, size = 1)+  
  geom_line(aes(y = coldshock$Abundance, color = "Cold Shock"), size = 1, linetype = "dotdash",alpha = 0.8)+
  geom_line(aes(y = coldspike$Abundance, color = "Cold Shock with Flow Spike"), size = 1, linetype = "longdash", alpha = 0.8)+
  scale_color_manual(values = colors)+
  theme_bw()+
  labs(x = "Year", y = "Abundance", color = "Scenario")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"), plot.margin = margin(5,5,5,20))

colors <- c("No Action" = "#FF7F00", "NonBypass" = "#33A02C")
ggplot(data = HYOS.scenarios, aes(x = Date, y = Abundance))+
  geom_line(aes(color = "No Action" ), alpha = 0.8, size = 1)+ 
  geom_line(aes(y = nbypass$Abundance, color = "NonBypass"), linetype = "dashed", size = 1, alpha = 0.8)+
  scale_color_manual(values = colors)+
  theme_bw()+
  labs(x = "Year", y = "Abundance", color = "Scenario")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"), plot.margin = margin(5,5,5,20))



# what about down by Diamond Creek



# #
# degreedays <- as.data.frame(cbind(temp_seq$dts, temps$Temperature * 14))
# colnames(degreedays) <- c("dts", "DegreeDay")
# 
# # specify iterations
# iterations <- 50
# 
# # baseline K in the absence of disturbance
# Kb <- 10000
# # max K after a big disturbance
# Kd <- 40000
# 
# 
# # specify baseline transition probabilities for each species
# # 3 stages - we have egg - larval instar V, pupae, and adult
# 
# G1_HYOS = 0.1 # move onto stage 2
# G2_HYOS = 0.445 # move onto stage 3
# P1_HYOS = 0.7 # remain in stage 1
# P2_HYOS = 0.0 # remain in stage 2
# 
# # want to run this for one year, in 14 day timesteps 
# #timestep <- seq(2, (length(flow.magnitude$Discharge) + 1), by = 1) # OR
# timestep <- seq(2, (length(temps$Temperature) + 1), by = 1)
# #timestep <- seq(2, (length(out_sample) + 1), by = 1)
# 
# # create an array to put our output into
# #output.N.array <- array(0, dim = c(length(timestep) + 1, length(species)))
# output.N.array <- array(0, dim = c(length(timestep) + 1))
# 
# output.N.list <- list(output.N.array)
# 
# 
# ## Assigning names to each array from sppnames vector
# #names(output.N.list) <- species
# 
# # create array to put the total N of all species into
# Total.N <- array(0,
#                  dim  <-c((length(timestep) +1 ), iterations),
#                  dimnames <- list(1:(length(timestep) + 1), 1:iterations))
# 
# # create list of arrays w/ abundance data for each spp
# reparray <- array(0,
#                   
#                   dim = c(length(timestep) + 1, 3, iterations),
#                   dimnames = list(1:(length(timestep)+1), c("S1", "S2", "S3"), 1:iterations)
# )
# 
# output.N.list <- reparray
# ## Repeating the array 7 times 
# #replist <- rep(list(reparray), 3)
# #names(replist) <- species
# 
# 
# # Q is equal to average discharge over 14 days
# #Q <- flow.magnitude$Discharge
# Q <- rep(0.1, length(temps$Temperature))
# Qmin <- 0.25
# a <- 0.1
# g <- 0.1
# h <- surv.fit.HYOS$m$getPars()[2]   
# k <- surv.fit.HYOS$m$getPars()[1] 
# extinction <- 500
# #-------------------------
# # Outer Loop of Iterations
# #--------------------------
# 
# 
# for (iter in c(1:iterations)) {
#   K = 10000 # need to reset K for each iteration
#   # we can also create a random flow scenario by sampleing flows
#   #out_sample <- sample(out$Discharge,length(out$Discharge), replace=TRUE)
#   #Q <- out_sample
#   
#   # another option is to keep flow the same each timestep (to test temp effects)
#   #out_same <- rep(10000, length(out$Discharge))
#   #Q <- out_same
#   
#   
#   # need to assign starting value
#   # we can start with 10 S1 individuals for each species
#   #for (sp in species){
#   #  output.N.list[[sp]][1,1] <- 10
#   #}
#   
#   # or we can pull randomw values from a uniform distribution 
#   output.N.list[1,1:3, iter]<- runif(3, min = 1, max = (0.5*K))
#   
#   # we often want to look at different parameter values after we run code, so we create some lists
#   
#   # list to input Ks
#   Klist <- vector()
#   Klist[1] <- 10000
#   
#   # list to imput flow morts
#   flowmortlist <- vector()
#   
#   Flist <- vector()
#   
#   emergetime <- vector()
#   
#   sizelist <- vector()
#   #-------------------------
#   # Inner Loop of Timesteps
#   #-------------------------
#   
#   for (t in timestep) {
#     
#     #----------------------------------------------------------
#     # Calculate how many timesteps emerging adults have matured
#     # Calculate how many timesteps emerging adults have matured
#     
#     emergetime <- append(emergetime, back.count.degreedays(t, 1680))
#     #---------------------------------------------------------
#     # Calculate fecundity per adult
#     
#       F_HYOS = rnorm(1, mean = 235.6, sd = 11.05102 ) * 0.5  #from Willis Jr & Hendricks, sd calculated from 95% CI = 21.66 = 1.96*sd
#     # * 0.5 assuming 50% female
#       
#       # we can scale fecundity based on the 95% CI of 21.66 (min = 213.94, max = 257.26) 
#       if (t > 15) {
#         size <- emergetime[t-1]
#         sizelist <- append(sizelist, size)
#         F_HYOS <- ((8.664 * size) + 127.3) * 0.5
#       }
#       
#     #---------------------------------------------------
#     # Calculate the disturbance magnitude-K relationship 
#     # Sets to 0 if below the Qmin
#     # Calculate the disturbance magnitude-K relationship 
#     # Sets to 0 if below the Qmin
#     Qf <- Qf.Function(Q[t-1], Qmin, a)
#   
#     #-------------------------------------------------------------------
#     # Calculate K arrying capacity immediately following the disturbance
#     # Calculate K arrying capacity immediately following the disturbance
#     K0 <- K + ((Kd-K)*Qf)
#     
#     # Calculate final K for timestep, including relationship between K and time since disturbance
#     K <- post.dist.K(K0, Kb, g)
#     Klist <- append(Klist, K)
# 
# 
#     #---------------------------------------------
#     # Calculate effect of density dependence on fecundity
#     
#     # Logistic via Rogosch et al. Fish Model
#     F_HYOS <- Logistic.Dens.Dependence(F_HYOS, K, Total.N[t-1, iter]) * 0.5
#     Flist <- append(Flist, F_HYOS)
#     #-----------------------------------------------
#     # Create Lefkovitch Matrix
#     
#     HYOS1 <- c(P1_HYOS, 0, F_HYOS)
#     HYOS2 <- c(G1_HYOS, P2_HYOS, 0)
#     HYOS3 <- c(0, G2_HYOS, 0) 
#     
#     AHYOS <- rbind( HYOS1, HYOS2, HYOS3)
#     
#     #-----------------------------------------------
#     # Calculate new transition probabilities based on temperature
#     # This is the growth v development tradeoff 
#     # don't know if this exists for HYOS - they can emerge under a wide temp gradient (<5 - 25+ C) but relationship between growth and temp 
#     
#     # development measures (basically, if below 10C, no development, if between 10 and 12, follows a function, if above 12, prob of transition to next stage is 0.6395)
#     # if (10 > temps$Temperature[t-1]) AHYOS[3,2] <- 0.001
#     # if (temps$Temperature[t-1] > 13) AHYOS[3,2] <- 0.55
#     # if (10 <= temps$Temperature[t-1] & temps$Temperature[t-1] <= 13) AHYOS[3,2] <- (0.183 * temps$Temperature[t-1]) -1.829 #(0.2 * temps$Temperature[t-1]) -2
#     # AHYOS[2,1] <- AHYOS[3,2] 
#     # 
#     # # growth (if below 10C, no growth can occur - everything basically freezes, if between 10 and 11, prob of remaining in same stage = 0.6395, if above 13, prob of transition to next stage is 0 )
#     # if (10 >= temps$Temperature[t-1]) AHYOS[2,2] <- 0.55
#     # if (temps$Temperature[t-1] > 13) AHYOS[2,2] <- 0.001
#     # if (10 < temps$Temperature[t-1] & temps$Temperature[t-1] <=  13) AHYOS[2,2] <- (-0.183 * temps$Temperature[t-1]) + 2.38 #(-0.1 temps$Temperature[t-1]) - 2.6
#     # 
#     # AHYOS[1,1] <- AHYOS[2,2] 
#     
#     # 
#     #Glist <-append(Glist, AHYOS[3,2])
#     #Plist <- append(Plist, AHYOS[2,2])
#     
#     #--------------------------------------
#     # Calculate abundances for each stage
#     
#     output.N.list[t, 1:3, iter] <- AHYOS %*% output.N.list[t-1, 1:3, iter] 
#     
#     #------------------------------------------
#     #Calculate immediate mortality due to flows
#     # mortality due to flooding follows N0 = Nz*e^-hQ
#     # but what about using a sigmoidal logistic function so we can control the threshold point and rate of growth
#     # following m = 1/1+e^-h*(x-xf)
#     # where h is is shape value
#     # x is Q, and xf is threshold point (100% of pop dies)
#     #plot(x = Q, y = 1/(1+exp(-0.02*(Q-100000))))
#     
#     #s1
#     output.N.list[t, 1, iter] <- flood.mortality(output.N.list[t, 1, iter], k, h, Q[t-1], Qmin)
#     #s2
#     output.N.list[t,2,iter] <- flood.mortality(output.N.list[t,2,iter], k, h, Q[t-1], Qmin)
#     #3
#     output.N.list[t,3,iter] <- flood.mortality(output.N.list[t,3,iter], k, h, Q[t-1], Qmin)
#     
#     flowmortlist <- append(flowmortlist, flood.mortality(1, k, h, Q[t-1], Qmin))
#     #replist[[1]][,,1] <- output.N.list[[1]]
#     Total.N[,iter] <- apply(output.N.list[,,iter],1,sum)
#     # check extinction threshold
#     if (Total.N[t, iter] < extinction){
#       output.N.list[t,,iter] <- 0
#       Total.N[t, iter] <- 0  
#       } #-------------------------
# }# End Inner Loop  
#   #------------------------- 
# } #----------------------
# 
# # End Outer Loop
# #----------------------

#------------------
# Analyzing Results
#-------------------
# summarizing iterations
means.list.HYOS <- mean.data.frame(out, burnin = 550, iteration = 9)
means.list.HYOS <- cbind(means.list.HYOS, temps$dts[549:1300])
means.list.HYOS$`temps$dts` <- as.Date(means.list.HYOS$`temps$dts`)
# plot abundance over time

arrows <- tibble(
  x1 = c("2030-01-07", "2035-01-07", "2040-01-07", "2045-01-07"),
  x2 = c("2030-01-07", "2035-01-07", "2040-01-07", "2045-01-07"),
  y1 = c(0.23, 0.23, 0.23, 0.23), 
  y2 = c(0.1, 0.1, 0.1, 0.1)
)



arrows$x1 <- as.Date(arrows$x1)
arrows$x2 <- as.Date(arrows$x2)

abund.trends.HYOS <- ggplot(data = means.list.HYOS, aes(x =  `temps$dts`,
                                              y = mean.abund/10000, group = 1)) +
  geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
                  ymax = mean.abund + 1.96 * se.abund),
              colour = 'transparent',
              alpha = .5,
              show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  coord_cartesian(ylim = c(0,0.25)) +
  ylab('Hydrospyche spp. Relative Abundance') +
  xlab(" ")+
  theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13))+
  scale_x_date(date_labels="%Y", date_breaks  ="1 year")+
  annotate("segment", x = arrows$x1, y = arrows$y1, xend = arrows$x2, yend = arrows$y2,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")), color = "red")+
  annotate("text", x = arrows$x1[1], y = 0.25, label = "+1°C", size = 5)+
  annotate("text", x = arrows$x1[2], y = 0.25, label = "+2.5°C", size = 5)+
  annotate("text", x = arrows$x1[3], y = 0.25, label = "+5°C", size = 5)+
  annotate("text", x = arrows$x1[4], y = 0.25, label = "+7.5°C", size = 5 )

a<- mean(means.list.HYOS$mean.abund[which(means.list.HYOS$`temps$dts` < "2030-01-07")])
b <- mean(means.list.HYOS$mean.abund[which(means.list.HYOS$`temps$dts` >= "2030-01-07" & means.list.HYOS$`temps$dts` < "2035-01-07")])
c <- mean(means.list.HYOS$mean.abund[which(means.list.HYOS$`temps$dts` >= "2035-01-07" & means.list.HYOS$`temps$dts` < "2040-01-07")])
d <- mean(means.list.HYOS$mean.abund[which(means.list.HYOS$`temps$dts` >= "2040-01-07" & means.list.HYOS$`temps$dts` < "2045-01-07")])
e <- mean(means.list.HYOS$mean.abund[which(means.list.HYOS$`temps$dts` >= "2045-01-07")])


plot(qr, c(a, b, c, d, e))
# take a look at results
# 
# par(mfrow = c(1,1))
plot(timestep[9:(length(timestep)+1)], output.N.list[9:(length(timestep)+1), 3, 1], type = "l", ylab = "Hydropsyche spp. Adults", xlab = "Timestep (1 fortnight)")
plot(timestep[9:length(timestep)], Total.N[10:(length(timestep)+1)], type= "l", ylab = "Hydropsyche spp. Total N", xlab = "Timestep (1 fortnight)")
# 
# 
# 
# #creating plots to analyze how temp relationship is working
#create dataframe with timestemp, s1, s2, and s3 abundances, and tempterature
 data <- as.data.frame(cbind(timestep, output.N.list[2:(length(timestep)+1) ,1, 1], output.N.list[2:(length(timestep)+1) ,2, 1], output.N.list[2:(length(timestep)+1) ,3, 1], temps$Temperature))
 colnames(data) <- c("timestep", "Stage1", "Stage2", "Stage3", "Temperature")

 data <- data[210:260, ]
ggplot(data = data, aes(x = timestep, y = Stage1, color = "Stage1"))+
   geom_path()+
   geom_path(aes(x = timestep, y = Stage2, color = "Stage2"))+
   geom_path(aes(x = timestep, y = Stage3, color = "Stage3"))+
   geom_path(aes(x = timestep, y = Temperature*200, color = "Temperature"))+
   scale_y_continuous(

     # Features of the first axis
     name = "Abundance",

     # Add a second axis and specify its features
     sec.axis = sec_axis( ~.*0.005, name="Temperature C")
   )
# 
#  # plot to show relationship between temp and fecundity
#  plot(temps$Temperature, Flist, ylab = "Fecundity per individual", xlab = "Temperature (C)", pch = 19)
# 
#  plot(timestep[10:60], output.N.list[10:60, 3,1] * 1, type = "l", ylim  = c(0, 500), ylab = " ")
#  lines(timestep[10:60], ((Flist[10:60]*output.N.list[10:60, 3,1])/10), col = "blue", ylim = c(0, 450))
#  lines(timestep[10:60], temps$Temperature[10:60]*15, col = "red")
# lines(timestep[10:60], emergetime[4:54]*18, col = "green")
#  lines(timestep[10:60], sizelist[4:54]*100, col = "magenta")
#  legend(25, 430, legend = c("Adult Abundance", "Realized Fecundity per female * 0.1", "Temperature (C) * 15", "Dry Weight (mg) * 100"),
#         col = c("black", "blue", "red", "magenta"), lty = 1, cex = 0.8)
# 
#  plot(timestep, temps$Temperature, type = "l", col = "blue")
# 
#  par(mfrow = c(1,2))
#  plot(timestep[10:37], Glist[10:37], type = "l", col = "blue")
#  lines(timestep[10:37], Plist[10:37], type = "l", col = "black")
#  legend(18, 0.1, legend=c("Growth (remain)", "Development (transition)"),
#         col=c("Black", "Blue"), lty=1, cex=0.8)
#  plot(timestep[10:37], temps$Temperature[10:37], type = "l", col = "red")
# 
#  plot(timestep[200:210], Total.N[201:211], type= "l", ylab = "HYOSis spp. Total N", xlab = "Timestep (1 fortnight")
#  par(new=TRUE)
#  lines(timestep[200:210],temps$Temperature[201:211],col="green")
# 
#  Total.N
# 
#  r <-Total.N[2:(length(timestep)+1)]/Total.N[1:length(timestep)]
#  plot(timestep, r, type = "l")
# 
#  plot(Q, Klist[1:940])
# 
#  
#  
