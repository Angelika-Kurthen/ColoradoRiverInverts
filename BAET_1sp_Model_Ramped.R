#--------------------------------
# Baet ramped
#--------------------------------
library(purrr)
library(tidyverse)
library(lubridate)
library(plyr)
library(dplyr)
library(ggplot2)
# data retrieval tool from USGS
library(dataRetrieval)


source("BAETSurvivorship.R")
source("1spFunctions_CR.R")
source("BAET_1sp_Model.R")

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
  out <- BAETmodel(flow.data = flow.df$flow.magnitude ,temp.data = temps, disturbanceK = 40000, baselineK = 5000, Qmin = 0.15, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
  BAET.df <- mean.data.frame(out, 250, 9)
  means[inc] <- mean(BAET.df$mean.abund)
  temps$Temperature[which(month(temps$dts) == 6 | month(temps$dts) == 7 | month(temps$dts) == 8)] <- temps$Temperature[which(month(temps$dts) == 6 | month(temps$dts) == 7 | month(temps$dts) == 8)] - increase[inc]
}


summer.BAET <- as.data.frame(cbind(increase, means))


ggplot(data = summer.BAET, aes(x = increase, y = means))+
  geom_point()+
  geom_line()+
  xlab("")

# sensitivity to hydropeaking

hydropeak <- seq(0, 0.5, by = 0.025)
means <- vector()
for (hydr in 1:length(hydropeak)){
  out <- BAETmodel(flow.data = flow.df$flow.magnitude, temp.data = temps, disturbanceK = 40000, baselineK = 5000, Qmin = 0.1, extinct = 50, iteration = 9, peaklist = hydropeak[hydr], peakeach = length(temps$Temperature))
  BAET.df <- mean.data.frame(out, 250, 9)
  means[hydr] <- mean(BAET.df$mean.abund)
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
out <- BAETmodel(flow.data = bauflow$flow.magnitude, temp.data = temps, disturbanceK = 40000, baselineK = 5000, Qmin = 0.1, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
BAET.df <- mean.data.frame(out, 250, 9)
BAET.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], BAET.df$mean.abund))
colnames(BAET.df) <- c("Date", "Abundance")
BAET.df$Date <- as.Date(BAET.df$Date)
# 2024 - 2027

# HFE of around 37000 cfs allowed Oct 1 through Nov 30. Will split the difference with 11-09 HFE day
# assume 7 days at high cfs and 7 days at baseline, so magnitude = 0.28 

bau <- BAET.df[which(BAET.df$Date >= "2024-01-01" & BAET.df$Date < "2028-01-01"),]
plot(bau$Date, bau$Abundance, type = "l")

# cool mix scenario
coolmix <- temps
# all temperatures are supposed to be 15.5 or below - this won't affect Lees Ferry because we don't really have super warm average temps (BUT we did in 2023)
coolmix[coolmix$Temperature > 15.5,] <- 15.5

# will use bauflow, since HFE still allowed
out <- BAETmodel(flow.data = bauflow$flow.magnitude, temp.data = coolmix, disturbanceK = 40000, baselineK = 5000, Qmin = 0.1, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
BAET.df <- mean.data.frame(out, 250, 9)
BAET.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], BAET.df$mean.abund))
colnames(BAET.df) <- c("Date", "Abundance")
BAET.df$Date <- as.Date(BAET.df$Date)

cool <- BAET.df[which(BAET.df$Date >= "2024-01-01" & BAET.df$Date < "2028-01-01"),]


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


out <- BAETmodel(flow.data = coolflow$flow.magnitude, temp.data = coolmix, disturbanceK = 40000, baselineK = 5000, Qmin = 0.1, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
BAET.df <- mean.data.frame(out, 250, 9)
BAET.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], BAET.df$mean.abund))
colnames(BAET.df) <- c("Date", "Abundance")
BAET.df$Date <- as.Date(BAET.df$Date)

coolspike <- BAET.df[which(BAET.df$Date >= "2024-01-01" & BAET.df$Date < "2028-01-01"),]
plot(coolspike$Date, coolspike$Abundance, type = "l" )

# cold schock alternative
# keep water below 13 C
coldtemps <- temps
coldtemps$Temperature[which(coldtemps$Temperature > 13)] <- 13
# use bauflows

out <- BAETmodel(flow.data = bauflow$flow.magnitude, temp.data = coldtemps, disturbanceK = 40000, baselineK = 5000, Qmin = 0.1, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
BAET.df <- mean.data.frame(out, 250, 9)
BAET.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], BAET.df$mean.abund))
colnames(BAET.df) <- c("Date", "Abundance")
BAET.df$Date <- as.Date(BAET.df$Date)

coldshock <- BAET.df[which(BAET.df$Date >= "2024-01-01" & BAET.df$Date < "2028-01-01"),]
plot(coldshock$Date, coldshock$Abundance, type = "l" )


# scenario 2
# coolflow$flow.magnitude[which(coolflow$dts == "2024-05-25")] <- 0.376
# coolflow$flow.magnitude[which(coolflow$dts == "2025-05-25" | coolflow$dts == "2025-06-08")] <- 0.376
# coolflow$flow.magnitude[which(coolflow$dts == "2026-05-25" | coolflow$dts == "2026-06-08")] <- 0.376
# 
# 
# out <- BAETmodel(flow.data = coolflow$flow.magnitude, temp.data = coolmix, disturbanceK = 40000, baselineK = 5000, Qmin = 0.1, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
# BAET.df <- mean.data.frame(out, 250, 9)
# BAET.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], BAET.df$mean.abund))
# colnames(BAET.df) <- c("Date", "Abundance")
# BAET.df$Date <- as.Date(BAET.df$Date)
# 
# coolspike <- BAET.df[which(BAET.df$Date >= "2024-01-01" & BAET.df$Date < "2028-01-01"),]
# plot(coolspike$Date, coolspike$Abundance, type = "l" )

# cold schock alternative
# keep water below 13 C
coldtemps <- temps
coldtemps$Temperature[which(coldtemps$Temperature > 13)] <- 13
# use bauflows

out <- BAETmodel(flow.data = bauflow$flow.magnitude, temp.data = coldtemps, disturbanceK = 40000, baselineK = 5000, Qmin = 0.1, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
BAET.df <- mean.data.frame(out, 250, 9)
BAET.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], BAET.df$mean.abund))
colnames(BAET.df) <- c("Date", "Abundance")
BAET.df$Date <- as.Date(BAET.df$Date)

coldshock <- BAET.df[which(BAET.df$Date >= "2024-01-01" & BAET.df$Date < "2028-01-01"),]
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

out <- BAETmodel(flow.data = coldflow$flow.magnitude, temp.data = coldtemps, disturbanceK = 40000, baselineK = 5000, Qmin = 0.1, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
BAET.df <- mean.data.frame(out, 250, 9)
BAET.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], BAET.df$mean.abund))
colnames(BAET.df) <- c("Date", "Abundance")
BAET.df$Date <- as.Date(BAET.df$Date)

coldspike <- BAET.df[which(BAET.df$Date >= "2024-01-01" & BAET.df$Date < "2028-01-01"),]
plot(coldspike$Date, coldspike$Abundance, type = "l" )

# scenario 2
# coldflow$flow.magnitude[which(coldflow$dts == "2024-05-25")] <- 0.376
# coldflow$flow.magnitude[which(coldflow$dts == "2025-05-25" | coldflow$dts == "2025-06-08")] <- 0.376
# coldflow$flow.magnitude[which(coldflow$dts == "2026-05-25" | coldflow$dts == "2026-06-08")] <- 0.376
# 
# out <- BAETmodel(flow.data = coldflow$flow.magnitude, temp.data = coldtemps, disturbanceK = 40000, baselineK = 5000, Qmin = 0.1, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
# BAET.df <- mean.data.frame(out, 250, 9)
# BAET.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], BAET.df$mean.abund))
# colnames(BAET.df) <- c("Date", "Abundance")
# BAET.df$Date <- as.Date(BAET.df$Date)
# 
# coldshock <- BAET.df[which(BAET.df$Date >= "2024-01-01" & BAET.df$Date < "2028-01-01"),]
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

out <- BAETmodel(flow.data = nonbypass$flow.magnitude, temp.data = temps, disturbanceK = 40000, baselineK = 5000, Qmin = 0.1, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
BAET.df <- mean.data.frame(out, 250, 9)
BAET.df <- as.data.frame(cbind(temps$dts[249:length(temps$dts)], BAET.df$mean.abund))
colnames(BAET.df) <- c("Date", "Abundance")
BAET.df$Date <- as.Date(BAET.df$Date)

nbypass <- BAET.df[which(BAET.df$Date >= "2024-01-01" & BAET.df$Date < "2028-01-01"),]
plot(nbypass$Date, nbypass$Abundance, type = "l" )

BAET.scenarios <- as.data.frame(cbind(bau, cool$Abundance, coolspike$Abundance, coldshock$Abundance, coldspike$Abundance, nbypass$Abundance))


colors <- c("No Action" = "#FF7F00", "Cool Mix" = "#A6CEE3", "Cool Mix with Flow Spike" = "#1F78B4")
ggplot(data = BAET.scenarios, aes(x = Date, y = Abundance))+
  geom_line(aes(color = "No Action" ), alpha = 0.8, size = 1)+  
  geom_line(aes(y = cool$Abundance,  color = "Cool Mix"),linetype = "dashed", alpha = 0.8, size = 1)+
  geom_line(aes(y = coolspike$Abundance, color = "Cool Mix with Flow Spike"), size = 1, linetype = "dashed", alpha = 0.8)+
  scale_color_manual(values = colors)+
  theme_bw()+
  labs(x = "Year", y = "Abundance", color = "Scenario")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"), plot.margin = margin(5,5,5,20))

colors <- c("No Action" = "#FF7F00", "Cold Shock" = "#CAB2D6", "Cold Shock with Flow Spike" = "#6A3D9A")
ggplot(data = BAET.scenarios, aes(x = Date, y = Abundance))+
  geom_line(aes(color = "No Action" ), alpha = 0.8, size = 1)+  
  geom_line(aes(y = coldshock$Abundance, color = "Cold Shock"), size = 1, linetype = "dotdash",alpha = 0.8)+
  geom_line(aes(y = coldspike$Abundance, color = "Cold Shock with Flow Spike"), size = 1, linetype = "longdash", alpha = 0.8)+
  scale_color_manual(values = colors)+
  theme_bw()+
  labs(x = "Year", y = "Abundance", color = "Scenario")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"), plot.margin = margin(5,5,5,20))

colors <- c("No Action" = "#FF7F00", "NonBypass" = "#33A02C")
ggplot(data = BAET.scenarios, aes(x = Date, y = Abundance))+
  geom_line(aes(color = "No Action" ), alpha = 0.8, size = 1)+ 
  geom_line(aes(y = nbypass$Abundance, color = "NonBypass"), linetype = "dashed", size = 1, alpha = 0.8)+
  scale_color_manual(values = colors)+
  theme_bw()+
  labs(x = "Year", y = "Abundance", color = "Scenario")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"), plot.margin = margin(5,5,5,20))



 
# #read in flow data from USGS gauge at Lees Ferry, AZ between 1985 to the end of the last water year
# flow <- readNWISdv("09380000", "00060", "1985-10-01", "2021-09-30")
# 
# # read in temperature data from USGS gauge at Lees Ferry, AZ between _ to the end of last water year
# temp <- readNWISdv("09380000", "00010", "2007-10-01", "2021-09-30")
# temps <- average.yearly.temp(temp, "X_00010_00003", "Date")
# n <- 50
# # qr is the temp ramps I want to increase the average Lees Ferry temp by 
# qr <- c(0,1, 2.5, 5, 7.5)
# # how many years I want each temp ramp to last
# r <- c(30, 5, 5, 5, 5)
# temps <- rep.avg.year(temps, n, change.in.temp = qr, years.at.temp = r)
# 
# temps <- temps[2:3]
# flow.magnitude <- rep(mean(flow$X_00060_00003)/85000, times = length(temps$Temperature))
# 
# out <- BAETmodel(flow.data = flow.magnitude,temp.data = temps, disturbanceK = 40000, baselineK = 5000, Qmin = 0.17, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))
# 






# 
# # if desired, can pre-count how many degreedays each cohort should experience
# #forward.count.degreedays(559)
# 
# #----------------------------------------------
# # set model parameters
# #----------------------------------------------
# # specify iterations
# iterations <- 1
# 
# # baseline K in the absence of disturbance
# Kb <- 10000
# # max K after a big disturbance
# Kd <- 40000
# 
# # specify baseline transition probabilities for each species
# G1_BAET = 0.25 #move to Stage2 (subimago)
# G2_BAET = 0.25 #move to Stage3 (adult)
# P1_BAET = 0.3 #stay in Stage1 (larvae)
# P2_BAET = 0.3 #stay in Stage2 (subimago)
# 
# # want this to run as long as our temperature timeseries
# timestep <- seq(2, (length(temps$Temperature) + 1), by = 1) 
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
# # Q is equal to average discharge over 14 days/bankful discharge for the system
# # in the Colorado River, bankful discharge = 85000 cfs (personal communication with Theodore Kennedy)
# #Q <- flow.magnitude$Discharge 
# 
# # We want to create a model with no flow mortality, so set Q less than Qmin
# Q <- rep(0.1, length(temps$Temperature))
# 
# Qmin <- 0.25 # Q min is the minimum flow required to impact mortality and carryin capactity (K)
# a <- 0.1 #shape param for flow magnitude and K
# g <- 0.1 #shape param for relationship between K and time since disturbance
# h <- surv.fit.BAET$m$getPars()[2] # shape param for flood mortality
# k <- surv.fit.BAET$m$getPars()[1] # shape param for flood mortality  
# extinction <- 50 # extinction threshold - if Total abundance below this, population goes extinct
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
#   # pull random values from a uniform distribution for starting pop
#   output.N.list[1,1:3, iter]<- runif(3, min = 1, max = (0.3*K))
#   
#   # we often want to look at different parameter values after we run code, so we create some lists
#   # list to input Ks
#   Klist <- vector()
#   Klist[1] <- 10000 #set first K
#   
#   # list to input flow mortality
#   flowmortlist <- vector()
#   
#   # list to input fecundities
#   Flist <- vector()
#   
#   # list to input back-looking emergence times
#   emergetime <- vector()
#   
#   # list to input size
#   sizelist <- vector()
#   
#   # list to input probs to remaining in same stage
#   Glist <- vector()
#   
#   # list to input eigenvalue
#   eigenlist <- vector()
#   #-------------------------
#   # Inner Loop of Timesteps
#   #-------------------------
#   
#   for (t in timestep) {
#     
#     #----------------------------------------------------------
#     # Calculate how many timesteps emerging adults have matured
#     
#     emergetime <- append(emergetime, back.count.degreedays(t, 559)) # value from Perry and Kennedy, 2016 
#     #---------------------------------------------------------
#     # Calculate fecundity per adult
#     
#     # we start by pulling fecundities from normal distribution
#     F_BAET = rnorm(1, mean = 1104.5, sd = 42.75) * 0.5 * 0.5  #Baetidae egg minima and maxima from Degrange, 1960, assuming 1:1 sex ratio and 50% egg mortality
#     
#     # we can also relate fecundities to body mass.
#     # Sweeney and Vannote 1980 have recorded dry body weight between 0.9 and 2.0 mg. 
#     # That weight is related to fecundity Y = 614X - 300
#     # we can "convert" emergetime to mg by multiplying to get dry weights between 0.9 - 2 mg, and then convert to fecunity
#     # Issue: this data is for Ephemerella spp, not Baetidae spp
#     
#     if (t > 15) {
#       size <- (emergetime[t-1] * 0.55)-0.75
#       sizelist <- append(sizelist, size)
#       F_BAET <- ((614 * size) - 300)* 0.5 * 0.5 
#     }
#     #--------------------------------------------------
#     # Calculate the disturbance magnitude-K relationship
#     # Sets to 0 if below the Qmin
#     Qf <- Qf.Function(Q[t-1], Qmin, a)
#     
#     #-------------------------------------------------------------------
#     # Calculate K carrying capacity immediately following the disturbance
#     K0 <- K + ((Kd-K)*Qf)
#     
#     # Calculate final K for timestep, including relationship between K and time since disturbance
#     K <- post.dist.K(K0, Kb, g)
#     
#     Klist <- append(Klist, K)
#     #---------------------------------------------
#     # Calculate effect of density dependence on fecundity 
#     
#     # Logistic via Rogosch et al. Fish Model
#     # no immediate egg mortality incorporated
#     F_BAET <- Logistic.Dens.Dependence(F_BAET, K, Total.N[t-1, iter])
#     # add F_BAET to list
#     Flist <- append(Flist, F_BAET)
#     #-----------------------------------------------
#     # Calculate new transition probabilities based on temperature
#     # This is the growth v development tradeoff
#     
#     # development measures
#     # in this function, we assume that if below the min temp threshold (9) no maturation occurs (slow maturation, large growth)
#     # if above the max temp threshold (15), no one remains more than 1 timestep in each stage (fast maturation, small growth)
#     
#     # Probabilities of remaining in stages (when temps low, high prob of remaining)
#     P1_BAET <- growth.development.tradeoff(temps$Temperature[t-1],  9, 13, 0.43, 0.0)
#     P2_BAET <- growth.development.tradeoff(temps$Temperature[t-1], 9, 13, 0.43, 0)
#     
#     # growth (if below 10C, no growth can occur - everything basically freezes, if between 10 and 11, prob of remaining in same stage = 0.6395, if above 13, prob of transition to next stage is 0 )
#     G1_BAET <- 0.43 - P1_BAET
#     G2_BAET <- 0.43 - P2_BAET
#     #-----------------------------------------------
#     # Create Lefkovitch Matrix
#     
#     BAET1 <- c(P1_BAET, 0, F_BAET)
#     BAET2 <- c(G1_BAET, P2_BAET, 0)
#     BAET3 <- c(0, G2_BAET, 0) 
#     
#     ABAET <- rbind( BAET1, BAET2, BAET3)
#     eigenlist <- append(eigenlist, eigen(ABAET)$values[1])
#     
#     #--------------------------------------
#     # Calculate abundances for each stage
#     
#     output.N.list[t, 1:3, iter] <- ABAET %*% output.N.list[t-1, 1:3, iter] 
#     
#     #------------------------------------------
#     #Calculate immediate mortality due to flows
#     # mortality due to flooding follows N0 = Nz*e^-hQ
#     
#     #s1
#     output.N.list[t, 1, iter] <- flood.mortality(output.N.list[t, 1, iter], k, h, Q[t-1], Qmin)
#     #s2Q
#     output.N.list[t,2,iter] <- flood.mortality(output.N.list[t,2,iter], k, h, Q[t-1], Qmin)
#     
#     output.N.list[t,3,iter] <- flood.mortality(output.N.list[t,3,iter], k, h, Q[t-1], Qmin)
#     
#     flowmortlist <- append(flowmortlist, flood.mortality(1, k, h, Q[t-1], Qmin))
#     
#     #------------------------------------------------------
#     # check extinction threshold and if below set to 0
#     Total.N[,iter] <- apply(output.N.list[,,iter],1,sum)
#     if (Total.N[t, iter] < extinction){
#       output.N.list[t,,iter] <- 0
#       Total.N[t, iter] <- 0}
#     
#     
#   } #-------------------------
#   # End Inner Loop  
#   #------------------------- 
# } #----------------------
# # End Outer Loop
# #----------------------
# #------------------
# # Analyzing Results
# #-------------------
# summarizing iterations

## turning replist into a df
means.list.BAET <- mean.data.frame(out, burnin = 550, iteration = 9)
means.list.BAET <- cbind(means.list.BAET, temps$dts[549:1300])
means.list.BAET$`temps$dts` <- as.Date(means.list.BAET$`temps$dts`)

arrows <- tibble(
  x1 = c("2030-01-07", "2035-01-07", "2040-01-07", "2045-01-07"),
  x2 = c("2030-01-07", "2035-01-07", "2040-01-07", "2045-01-07"),
  y1 = c(0.23, 0.23, .23, .23), 
  y2 = c(0.1,.1,.1,.1)
)

arrows$x1 <- as.Date(arrows$x1)
arrows$x2 <- as.Date(arrows$x2)

# note how boom and bust this model is - K is set to be 10,000 not 100,000
abund.trends.BAET <- ggplot(data = means.list.BAET, aes(x = `temps$dts`,
                                                        y = mean.abund/40000, group = 1)) +
  geom_line(show.legend = FALSE) +
  coord_cartesian(ylim = c(0,0.25)) +
  ylab('Baetidae Relative Abundance') +
  xlab(' ')+
  theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13))+
  scale_x_date(date_labels="%Y", date_breaks  ="1 year")+
  annotate("segment", x = arrows$x1, y = arrows$y1, xend = arrows$x2, yend = arrows$y2,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")), color = "red")+
  annotate("text", x = arrows$x1[1], y = 0.25, label = "+1°C", size = 5)+
  annotate("text", x = arrows$x1[2], y = 0.25, label = "+2.5°C", size = 5)+
  annotate("text", x = arrows$x1[3], y = 0.25, label = "+5°C", size = 5)+
  annotate("text", x = arrows$x1[4], y = 0.25, label = "+7.5°C", size = 5 )


