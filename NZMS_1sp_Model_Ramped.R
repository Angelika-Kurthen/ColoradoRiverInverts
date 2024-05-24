##########################
# NZMS 1 sp model
###########################
library(purrr)
library(tidyverse)
library(lubridate)
library(plyr)
library(dplyr)
library(ggplot2)
# data retrieval tool from USGS
library(dataRetrieval)

# Code for HPC - tidyverse has some issues on our HPC because one of the packages is deprecated
# We have to manually load all tidyverse packages
# library(purrr, lib.loc = "/home/ib/kurthena/R_libs/4Phosphorus-mediated changes in life history traits of the invasive New Zealand mudsnail (Potamopyrgus antipodarum) .2.1")
# library(tibble, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(tidyr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(readr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(stringr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(forcats, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(lubridate, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(plyr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(dplyr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(ggplot2, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(dataRetrieval, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")

# source functions
source("NZMSSurvivorship.R")
source("1spFunctions.R")
source("NZMS_1sp_Model.R")

#read in flow data from USGS gauge at Lees Ferry, AZ between 1985 to the end of the last water year
flow <- readNWISdv("09380000", "00060", "2014-01-01", "2024-09-30")
flows <- average.yearly.flows(flowdata = flow, "X_00060_00003", "Date")

# read in temperature data from USGS gauge at Lees Ferry, AZ between _ to the end of last water year
temp <- readNWISdv("09380000", "00010", "2014-05-01", "2024-05-01")
temp <- temp[, c(3,4)]
colnames(temp) <- c("Date", "Temperature")
tempslist <- average.yearly.temp(temp, "Temperature", "Date")
years <- seq(0, 100, by = 1)
temps <- average.yearly.temp(temp, "Temperature", "Date")

# add new years so we can progress forwards in time, then concatonate
for(year in years){
  year(tempslist$dts) <- (2000 + year)
  temps <- rbind(temps, tempslist)
}
temps <- temps[-(1:26),-c(1,4, 5, 6)]

# do the same with out yearly flows
flow.df <- as.data.frame(cbind(temps$dts, rep(flows$Discharge/85000, times = 101)))
colnames(flow.df) <- c("dts", "flow.magnitude")
flow.df$dts <- as.Date(flow.df$dts)

# Increase in summer temps
means <- vector()
increase <- seq(0, 5, by = 0.5)
for (inc in 1:length(increase)){
  temps$Temperature[which(month(temps$dts) == 6 | month(temps$dts) == 7 | month(temps$dts) == 8)] <- temps$Temperature[which(month(temps$dts) == 6 | month(temps$dts) == 7 | month(temps$dts) == 8)] + increase[inc]
  out <- NZMSmodel(flow.data = flow.df$flow.magnitude ,temp.data = temps, disturbanceK = 1000000, baselineK = 5000, Qmin = 0.25, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
  NZMS.df <- mean.data.frame(out, 250, 9)
  means[inc] <- mean(NZMS.df$mean.abund)
  temps$Temperature[which(month(temps$dts) == 6 | month(temps$dts) == 7 | month(temps$dts) == 8)] <- temps$Temperature[which(month(temps$dts) == 6 | month(temps$dts) == 7 | month(temps$dts) == 8)] - increase[inc]
}


summer.NZMS <- as.data.frame(cbind(increase, means))


ggplot(data = summer.NZMS, aes(x = increase, y = means))+
  geom_point()+
  geom_line()+
  xlab("")

# sensitivity to hydropeaking

hydropeak <- seq(0, 0.5, by = 0.025)
means <- vector()
for (hydr in 1:length(hydropeak)){
  out <- NZMSmodel(flow.data = flow.df$flow.magnitude, temp.data = temps, disturbanceK = 1000000, baselineK = 5000, Qmin = 0.25, extinct = 50, iteration = 9, peaklist = hydropeak[hydr], peakeach = length(temps$Temperature))
  NZMS.df <- mean.data.frame(out, 250, 9)
  means[hydr] <- mean(NZMS.df$mean.abund)
}


hydropeak <- as.data.frame(cbind(hydropeak, means))

ggplot(data = hydropeak, aes(x = hydropeak, y = means))+
  geom_point()+
  geom_line()+
  theme_bw()+
  labs(x = "Hydropeaking Index", y = "Average Annual Abundance")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"), plot.margin = margin(5,5,5,20))



#BAU scenario
# regular temps and regular flows - we will look at year 2024 thru year 
bauflow <- flow.df
bauflow[which(month(bauflow$dts) == 11 & day(bauflow$dts) == 9),] <- 0.28
out <- NZMSmodel(flow.data = bauflow$flow.magnitude, temp.data = temps, disturbanceK = 1000000, baselineK = 5000, Qmin = 0.25, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
NZMS.df <- mean.data.frame(out, 250, 9)
NZMS.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], NZMS.df$mean.abund))
colnames(NZMS.df) <- c("Date", "Abundance")
NZMS.df$Date <- as.Date(NZMS.df$Date)
# 2024 - 2027

# HFE of around 37000 cfs allowed Oct 1 through Nov 30. Will split the difference with 11-09 HFE day
# assume 7 days at high cfs and 7 days at baseline, so magnitude = 0.28 

bau <- NZMS.df[which(NZMS.df$Date >= "2024-01-01" & NZMS.df$Date < "2028-01-01"),]
plot(bau$Date, bau$Abundance, type = "l")

# cool mix scenario
coolmix <- temps
# all temperatures are supposed to be 15.5 or below - this won't affect Lees Ferry because we don't really have super warm average temps (BUT we did in 2023)
coolmix[coolmix$Temperature > 15.5,] <- 15.5

# will use bauflow, since HFE still allowed
out <- NZMSmodel(flow.data = bauflow$flow.magnitude, temp.data = coolmix, disturbanceK = 1000000, baselineK = 5000, Qmin = 0.25, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
NZMS.df <- mean.data.frame(out, 250, 9)
NZMS.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], NZMS.df$mean.abund))
colnames(NZMS.df) <- c("Date", "Abundance")
NZMS.df$Date <- as.Date(NZMS.df$Date)

cool <- NZMS.df[which(NZMS.df$Date >= "2024-01-01" & NZMS.df$Date < "2028-01-01"),]
cool2 <- NZMS.df[which(NZMS.df$Date >= "2010-01-01" & NZMS.df$Date < "2014-01-01"),]

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


out <- NZMSmodel(flow.data = coolflow$flow.magnitude, temp.data = coolmix, disturbanceK = 1000000, baselineK = 5000, Qmin = 0.25, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
NZMS.df <- mean.data.frame(out, 250, 9)
NZMS.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], NZMS.df$mean.abund))
colnames(NZMS.df) <- c("Date", "Abundance")
NZMS.df$Date <- as.Date(NZMS.df$Date)

coolspike <- NZMS.df[which(NZMS.df$Date >= "2024-01-01" & NZMS.df$Date < "2028-01-01"),]
plot(coolspike$Date, coolspike$Abundance, type = "l" )

# cold schock alternative
# keep water below 13 C
coldtemps <- temps
coldtemps$Temperature[which(coldtemps$Temperature > 13)] <- 13
# use bauflows

out <- NZMSmodel(flow.data = bauflow$flow.magnitude, temp.data = coldtemps, disturbanceK = 1000000, baselineK = 5000, Qmin = 0.25, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
NZMS.df <- mean.data.frame(out, 250, 9)
NZMS.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], NZMS.df$mean.abund))
colnames(NZMS.df) <- c("Date", "Abundance")
NZMS.df$Date <- as.Date(NZMS.df$Date)

coldshock <- NZMS.df[which(NZMS.df$Date >= "2024-01-01" & NZMS.df$Date < "2028-01-01"),]
plot(coldshock$Date, coldshock$Abundance, type = "l" )


# scenario 2
# coolflow$flow.magnitude[which(coolflow$dts == "2024-05-25")] <- 0.376
# coolflow$flow.magnitude[which(coolflow$dts == "2025-05-25" | coolflow$dts == "2025-06-08")] <- 0.376
# coolflow$flow.magnitude[which(coolflow$dts == "2026-05-25" | coolflow$dts == "2026-06-08")] <- 0.376
# 
# 
# out <- NZMSmodel(flow.data = coolflow$flow.magnitude, temp.data = coolmix, disturbanceK = 1000000, baselineK = 5000, Qmin = 0.1, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
# NZMS.df <- mean.data.frame(out, 250, 9)
# NZMS.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], NZMS.df$mean.abund))
# colnames(NZMS.df) <- c("Date", "Abundance")
# NZMS.df$Date <- as.Date(NZMS.df$Date)
# 
# coolspike <- NZMS.df[which(NZMS.df$Date >= "2024-01-01" & NZMS.df$Date < "2028-01-01"),]
# plot(coolspike$Date, coolspike$Abundance, type = "l" )

# cold schock alternative
# keep water below 13 C
coldtemps <- temps
coldtemps$Temperature[which(coldtemps$Temperature > 13)] <- 13
# use bauflows

out <- NZMSmodel(flow.data = bauflow$flow.magnitude, temp.data = coldtemps, disturbanceK = 1000000, baselineK = 5000, Qmin = 0.25, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
NZMS.df <- mean.data.frame(out, 250, 9)
NZMS.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], NZMS.df$mean.abund))
colnames(NZMS.df) <- c("Date", "Abundance")
NZMS.df$Date <- as.Date(NZMS.df$Date)

coldshock <- NZMS.df[which(NZMS.df$Date >= "2024-01-01" & NZMS.df$Date < "2028-01-01"),]
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

out <- NZMSmodel(flow.data = coldflow$flow.magnitude, temp.data = coldtemps, disturbanceK = 1000000, baselineK = 5000, Qmin = 0.25, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
NZMS.df <- mean.data.frame(out, 250, 9)
NZMS.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], NZMS.df$mean.abund))
colnames(NZMS.df) <- c("Date", "Abundance")
NZMS.df$Date <- as.Date(NZMS.df$Date)

coldspike <- NZMS.df[which(NZMS.df$Date >= "2024-01-01" & NZMS.df$Date < "2028-01-01"),]
plot(coldspike$Date, coldspike$Abundance, type = "l" )

# scenario 2
# coldflow$flow.magnitude[which(coldflow$dts == "2024-05-25")] <- 0.376
# coldflow$flow.magnitude[which(coldflow$dts == "2025-05-25" | coldflow$dts == "2025-06-08")] <- 0.376
# coldflow$flow.magnitude[which(coldflow$dts == "2026-05-25" | coldflow$dts == "2026-06-08")] <- 0.376
# 
# out <- NZMSmodel(flow.data = coldflow$flow.magnitude, temp.data = coldtemps, disturbanceK = 1000000, baselineK = 5000, Qmin = 0.1, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
# NZMS.df <- mean.data.frame(out, 250, 9)
# NZMS.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], NZMS.df$mean.abund))
# colnames(NZMS.df) <- c("Date", "Abundance")
# NZMS.df$Date <- as.Date(NZMS.df$Date)
# 
# coldshock <- NZMS.df[which(NZMS.df$Date >= "2024-01-01" & NZMS.df$Date < "2028-01-01"),]
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

out <- NZMSmodel(flow.data = nonbypass$flow.magnitude, temp.data = temps, disturbanceK = 1000000, baselineK = 5000, Qmin = 0.25, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
NZMS.df <- mean.data.frame(out, 250, 9)
NZMS.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], NZMS.df$mean.abund))
colnames(NZMS.df) <- c("Date", "Abundance")
NZMS.df$Date <- as.Date(NZMS.df$Date)

nbypass <- NZMS.df[which(NZMS.df$Date >= "2024-01-01" & NZMS.df$Date < "2028-01-01"),]
plot(nbypass$Date, nbypass$Abundance, type = "l" )

NZMS.scenarios <- as.data.frame(cbind(bau, cool$Abundance, coolspike$Abundance, coldshock$Abundance, coldspike$Abundance, nbypass$Abundance))


colors <- c("No Action" = "#FF7F00", "Cool Mix" = "#A6CEE3", "Cool Mix with Flow Spike" = "#1F78B4")
ggplot(data = NZMS.scenarios, aes(x = Date, y = Abundance))+
  geom_line(aes(color = "No Action" ), alpha = 0.8, size = 1)+  
  geom_line(aes(y = cool$Abundance,  color = "Cool Mix"),linetype = "dashed", alpha = 0.8, size = 1)+
  geom_line(aes(y = coolspike$Abundance, color = "Cool Mix with Flow Spike"), size = 1, linetype = "dashed", alpha = 0.8)+
  scale_color_manual(values = colors)+
  theme_bw()+
  labs(x = "Year", y = "Abundance", color = "Scenario")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"), plot.margin = margin(5,5,5,20))

colors <- c("No Action" = "#FF7F00", "Cold Shock" = "#CAB2D6", "Cold Shock with Flow Spike" = "#6A3D9A")
ggplot(data = NZMS.scenarios, aes(x = Date, y = Abundance))+
  geom_line(aes(color = "No Action" ), alpha = 0.8, size = 1)+  
  geom_line(aes(y = coldshock$Abundance, color = "Cold Shock"), size = 1, linetype = "dotdash",alpha = 0.8)+
  geom_line(aes(y = coldspike$Abundance, color = "Cold Shock with Flow Spike"), size = 1, linetype = "longdash", alpha = 0.8)+
  scale_color_manual(values = colors)+
  theme_bw()+
  labs(x = "Year", y = "Abundance", color = "Scenario")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"), plot.margin = margin(5,5,5,20))

colors <- c("No Action" = "#FF7F00", "NonBypass" = "#33A02C")
ggplot(data = NZMS.scenarios, aes(x = Date, y = Abundance))+
  geom_line(aes(color = "No Action" ), alpha = 0.8, size = 1)+ 
  geom_line(aes(y = nbypass$Abundance, color = "NonBypass"), linetype = "dashed", size = 1, alpha = 0.8)+
  scale_color_manual(values = colors)+
  theme_bw()+
  labs(x = "Year", y = "Abundance", color = "Scenario")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"), plot.margin = margin(5,5,5,20))

# # specify iterations
# iterations <- 5
# 
# # baseline K in the absence of disturbance
# Kb <- 10000
# # max K after a big disturbance
# Kd <- 1000000
# 
# # specify baseline transition probabilities
# # its speculated (Cross et al 2010) that survivorship is between 80 - 100% for NZMS in Grand Canyon - will say 90% survive, 10% baseline mortality
# # from timestep to timestep, we expect 90% to survive so for stage 1 (which lasts aproximately 14 timesteps, survival should be approx 09^14 = 0.2287679
# # from that, only 1/14th will transition out, 13/14 remain in stage
# 
# # stage 1 G1 (prob transition to stage 2) = 0.95^14/14
# # stage 1 P1 (prob remaining in stage 2) = 0.95^14 -0.95^14/14
# # stage 2 G2 (prob transition to stage 3) = 0.95^7/7
# # stage 2 P2 (prob remaining in stage 2) = 0.95^7 - 0.95^7/7
# # stage 3 P3 (prob remaining in stage 3) = 0.95 ^5 - 0.95^5/5
# 
# G1_NZMS = 0.1
# G2_NZMS = 0.1
# P1_NZMS = 0.4
# P2_NZMS = 0.5
# P3_NZMS = 0.8
# 
# # want to run this for one year, in 14 day timesteps 
# #timestep <- seq(2, (length(flow.magnitude$Discharge) + 1), by = 1) # OR
# timestep <- seq(2, (length(temps$Temperature) + 1), by = 1)
# # create an array to put our output into
# output.N.array <- array(0, dim = c(length(timestep) + 1))
# 
# output.N.list <- list(output.N.array)
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
# 
# # Q is equal to average discharge over 14 days
# #Q <- flow.magnitude$Discharge
# 
# Q <- rep(0.1, length(temps$Temperature))
# 
# Qmin <- 0.25
# a <- 0.1
# g <- 0.1
# h <- surv.fit.NZMS$m$getPars()[2]   
# k <- surv.fit.NZMS$m$getPars()[1]
# e = 2.71828
# extinction <- 1
# 
# #-------------------------
# # Outer Loop of Iterations
# #--------------------------
# 
# 
# for (iter in c(1:iterations)) {
#   
#   K = 10000 # need to reset K for each iteration
#   
#   # we can also create a random flow scenario by sampling flows
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
#   # or we can pull random values from a uniform distribution 
#   output.N.list[1,1:3, iter]<- runif(3, min = 1, max = (0.5*K))
#   
#   # we often want to look at different parameter values after we run code, so we create some lists
#   # list to input Ks
#   Klist <- vector()
#   Klist[1] <- 10000
#   Flist <- vector()
#   # list to imput flow morts
#   flowmortlist <- vector()
#   Flist <- vector()
#   emergetime <- vector()
#   sizelist <- vector()
#   
#   #-------------------------
#   # Inner Loop of Timesteps
#   #-------------------------
#   
#   for (t in timestep) {
#     
#     #---------------------------------------------------------
#     # Calculate starting fecundity per adult
#     
#     # fecundities estimated from McKenzie et al. 2013; 
#     if (temps$Temperature[t-1] <= 10) { 
#       F2_NZMS <- 2
#       F3_NZMS <- 2  
#     } else {
#       F2_NZMS <- 8.87473
#       F3_NZMS <- 27.89665
#     }
#     
#     #---------------------------------------------------
#     # Calculate the disturbance magnitude-K relationship 
# 
#     # Calculate the disturbance magnitude-K relationship
#     # Sets to 0 if below the Qmin
#     Qf <- Qf.Function(Q[t-1], Qmin, a)
#     
#     #-------------------------------------------------------------------
#     # Calculate K arrying capacity immediately following the disturbance
#     K0 <- K + ((Kd-K)*Qf)
#     
#     # Calculate final K for timestep, including relationship between K and time since disturbance
#     K <- post.dist.K(K0, Kb, g)
#     Klist <- append(Klist, K)
#     #---------------------------------------------
#     # Calculate effect of density dependnce on fecundity 
#     
#     # Logistic Density Dependence on Fecundity via Rogosch et al. Fish Model \
#     # assume 97% Female poplation, and 30% instant mortality for eggs laid
#     F2_NZMS <- Logistic.Dens.Dependence(F2_NZMS, K, Total.N[t-1, iter]) * 0.97 * 0.7
#     F3_NZMS <- Logistic.Dens.Dependence(F3_NZMS, K, Total.N[t-1, iter]) * 0.97 * 0.7
#     
#     Flist <- append(Flist, (F2_NZMS + F3_NZMS))
#     #-----------------------------------------------
#     # Create Lefkovitch Matrix
#     
#     NZMS1 <- c(P1_NZMS, F2_NZMS, F3_NZMS)
#     NZMS2 <- c(G1_NZMS, P2_NZMS, 0)
#     NZMS3 <- c(0, G2_NZMS, P3_NZMS) 
#     
#     ANZMS <- rbind( NZMS1, NZMS2, NZMS3)
#     
#     #-----------------------------------------------
#     # Calculate new transition probabilities based on temperature
#     # This is the growth v development tradeoff
#     #--------------------------------------
#     # Calculate abundances for each stage
#     
#     output.N.list[t, 1:3, iter] <- ANZMS %*% output.N.list[t-1, 1:3, iter] 
#     
#     #------------------------------------------
#     #Calculate immediate mortality due to flows
# 
#     #s1
#     output.N.list[t, 1, iter] <- flood.mortality(output.N.list[t,1,iter], k, h, Q[t-1], Qmin)
#     #s2
#     output.N.list[t,2,iter] <- flood.mortality(output.N.list[t,2,iter], k, h, Q[t-1], Qmin)
#     #3
#     output.N.list[t,3,iter] <- flood.mortality(output.N.list[t,3,iter], k, h, Q[t-1], Qmin)
#     
#     # in case we want to look back on what the flow mortality rates were
#     flowmortlist <- append(flowmortlist, flood.mortality(1, k, h, Q[t-1], Qmin))
# 
#     #-------------------------------------------------
#     # Calculate sum of all stages (total population)
#     Total.N[,iter] <- apply(output.N.list[,,iter],1,sum)
#     # check extinction threshold
#     if (Total.N[t, iter] < extinction){
#       output.N.list[t,,iter] <- 0
#       Total.N[t, iter] <- 0
#   
#     } #-------------------------
#   # End Inner Loop  
#   #------------------------- 
# } #----------------------
# # End Outer Loop
# #----------------------

#------------------
# Analyzing Results
#-------------------
# summarizing iterations



## turning replist into a d
means.list.NZMS <- mean.data.frame(out, burnin = 550, iteration  = 9)
means.list.NZMS <- cbind(means.list.NZMS, temps$dts[549:1300])
means.list.NZMS$`temps$dts` <- as.Date(means.list.NZMS$`temps$dts`)
# plot abundance over time

arrows <- tibble(
  x1 = c("2030-01-07", "2035-01-07", "2040-01-07", "2045-01-07"),
  x2 = c("2030-01-07", "2035-01-07", "2040-01-07", "2045-01-07"),
  y1 = c(0.9, 0.9, 0.9, 0.9), 
  y2 = c(0.7, 0.7, 0.7, 0.7)
)

arrows$x1 <- as.Date(arrows$x1)
arrows$x2 <- as.Date(arrows$x2)
abund.trends.NZMS <- ggplot(data = means.list.NZMS, aes(x = `temps$dts`,
                                              y = mean.abund/10000, group = 1)) +
  geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
                  ymax = mean.abund + 1.96 * se.abund),
              colour = 'transparent',
              alpha = .5,
              show.legend = FALSE) +
  geom_line(show.legend = FALSE, linewidth = 0.7) +
  coord_cartesian(ylim = c(0,1)) +
  ylab('New Zealand Mudsnail Relative Abundance') +
  xlab(" ")+
  theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13))+
  scale_x_date(date_labels="%Y", date_breaks  ="1 year")+
  annotate("segment", x = arrows$x1, y = arrows$y1, xend = arrows$x2, yend = arrows$y2,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")), color = "red")+
  annotate("text", x = arrows$x1[1], y = 1, label = "+1°C", size = 5)+
  annotate("text", x = arrows$x1[2], y = 1, label = "+2.5°C", size = 5)+
  annotate("text", x = arrows$x1[3], y = 1, label = "+5°C", size = 5)+
  annotate("text", x = arrows$x1[4], y = 1, label = "+7.5°C", size = 5 )


a <- mean(means.list.NZMS$mean.abund[which(means.list.NZMS$`temps$dts` < "2030-01-07")])
b <- mean(means.list.NZMS$mean.abund[which(means.list.NZMS$`temps$dts` >= "2030-01-07" & means.list.NZMS$`temps$dts` < "2035-01-07")])
c <- mean(means.list.NZMS$mean.abund[which(means.list.NZMS$`temps$dts` >= "2035-01-07" & means.list.NZMS$`temps$dts` < "2040-01-07")])
d <- mean(means.list.NZMS$mean.abund[which(means.list.NZMS$`temps$dts` >= "2040-01-07" & means.list.NZMS$`temps$dts` < "2045-01-07")])
e <- mean(means.list.NZMS$mean.abund[which(means.list.NZMS$`temps$dts` >= "2045-01-07")])


plot(qr, c(a, b, c, d, e))


