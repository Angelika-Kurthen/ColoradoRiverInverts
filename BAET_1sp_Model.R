##########################
# Baetis 1 sp model
###########################


library(purrr)
library(tidyverse)
library(lubridate)
library(dplyr)
library(demogR)

# data retrieval tool from USGS
install.packages("dataRetrieval")
library(dataRetrieval)

source("EggMortalityIntegration.R")

#read in flow data from USGS gauge at Lees Ferry, AZ between 1985 to the end of the last water year
flow <- readNWISdv("09380000", "00060", "1985-10-01", "2021-09-30")

# Make an index to be used for aggregating
ID <- as.numeric(as.factor(flow$Date))-1
# want it to be every 14 days, hence the 14
ID <- ID %/% 14
# aggregate over ID and TYPE for all numeric data.
out <- aggregate(flow[sapply(flow,is.numeric)],
                 by=list(ID,flow$X_00060_00003),
                 FUN=mean)
# format output
names(out)[1:2] <-c("dts","Discharge")
# add the correct dates as the beginning of every period
out$dts <- as.POSIXct(flow$Date[(out$dts*14)+1])
# order by date in chronological order
out <- out[order(out$dts),]
# get mean Discharge data for every 14 days
out <- aggregate(out, by = list(out$dts), FUN = mean)


# read in temperature data from USGS gauge at Lees Ferry, AZ between _ to the end of last water year
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2021-09-30")

# Make an index to be used for aggregating
ID <- as.numeric(as.factor(temp$Date))-1
# want it to be every 14 days, hence the 14
ID <- ID %/% 14
# aggregate over ID and TYPE for all numeric data.
outs <- aggregate(temp[sapply(temp,is.numeric)],
                  by=list(ID,temp$X_00010_00003),
                  FUN=mean)
# format output
names(outs)[1:2] <-c("dts","Temperature")
# add the correct dates as the beginning of every period
outs$dts <- as.POSIXct(temp$Date[(outs$dts*14)+1])
# order by date in chronological order
outs <- outs[order(outs$dts),]
# get mean Discharge data for every 14 days
outs <- aggregate(outs, by = list(outs$dts), FUN = mean)
mean_temp <- mean(outs$Temperature)

# there are less temperatures than discharge readings, so we will just repeat the same temp sequence for this exercise
temps <- rbind(outs, outs, outs)
temps <- temps[1:length(out$Discharge), ]

#make matrix with fecundities

# specify iterations
iterations <- 1
#species <- c("BAET", "SIMU", "CHIRO")


# set carrying capacity
K = 10000

# specify baseline transition probabilities for each species
G1_BAET = 0.032
G2_BAET = 0.032
P1_BAET = 0.032
P2_BAET = 0.032
# # transition probabilites when there is lowered flow (Q<8000)
# DG1_BAET = 0.75
# DG2_BAET = 0.7
# DP1_BAET = 0.15
# DP2_BAET = 0.15
# # transition probabilities when there is a higher flow (15000 > Q >  10000)
# SG1_BAET = 0.45
# SG2_BAET = 0.4
# SP1_BAET = 0.15
# SP2_BAET = 0.15

# transition probabilities when thre is a high flow (Q>15000)
# LG1_BAET = 0.29
# LG2_BAET = 0.2
# LP1_BAET = 0.15
# LP2_BAET = 0.15


# want to run this for one year, in 14 day timesteps 
timestep <- seq(2, (length(out$Discharge) + 1), by = 1) # OR
timestep <- seq(2, (length(out_sample) + 1), by = 1)

# create an array to put our output into
#output.N.array <- array(0, dim = c(length(timestep) + 1, length(species)))
output.N.array <- array(0, dim = c(length(timestep) + 1))

output.N.list <- list(output.N.array)


## Assigning names to each array from sppnames vector
#names(output.N.list) <- species

# create array to put the total N of all species into
Total.N <- array(0,
                 dim  <-c((length(timestep) +1 ), iterations),
                 dimnames <- list(1:(length(timestep) + 1), 1:iterations))

# create list of arrays w/ abundance data for each spp
reparray <- array(0,
                  
                  dim = c(length(timestep) + 1, 3, iterations),
                  dimnames = list(1:(length(timestep)+1), c("S1", "S2", "S3"), 1:iterations)
)


## Repeating the array 7 times 
#replist <- rep(list(reparray), 3)
#names(replist) <- species


# Q is equal to average discharge over 14 days
Q <- out$Discharge #OR

Qmin <- 20000
a <- 100
g <- 0.01

# beverton holt is Nt+1 = rNt/1-Nt(r-1)/K
# it is supposed to be depensatory, so as t -> inf, Nt+1 -> K, BUT 
# the discrete nature of this causes it overshoot by a lot, 
# meaning it isn't any better or worse than traditional logistric growth models

# ricker model Nt+1= Nt* e^r(1-Nt/K)
# overcompensatory, but could represent pops better? used in fisheries a lot
# r is instrisice rate of increase when N is small 
# r can be species specific, or the same for all species
# in this case, we will use the r that McMullen et al 2017 used for Beatis
# r = 1.23
e = 2.71828
b = 0.005

for (iter in iterations){
  
  # we can also create a random flow scenario by sampleing flows
  out_sample <- sample(out$Discharge,length(out$Discharge), replace=TRUE)
  Q <- out_sample
  
  # need to assign starting value
  # in the future, we can pull these #s from a randomly selected date in the Colorado River data
  # for now, will start with 10 S1 individuals for each species
  #for (sp in species){
  #  output.N.list[[sp]][1,1] <- 10
  #}
  
  output.N.list <- reparray
  output.N.list[1,1:3,]<- runif(3, min = 1, max = (0.5*K))
  
  for (t in timestep){
    F_BAET = rnorm(1, mean = 1104.5, sd = 42.75) # * H_BAET #Baetidae egg minima and maxima from Degrange, 1960
    
    
    # Calculate the disturbance magnitude-K relationship. Sets to 0 if below the Qmin
    if (Q[t-1] < Qmin) {
      Qf <- 0
    } else {
      Qf <- (Q[t-1] - Qmin)/(a + Q[t-1]- Qmin)
    }
    
    
    
    # Calculate K arrying capacity immediately following the disturbance
    K <- 10000 + ((40000-10000)*Qf)
    
    # Function to calc. K as a function of time post-disturbance at a particular disturbance intensity
    K <- 10000 + ((K - 10000)*exp(-g*14))
    
    # Ricker model - pro = doesn't go negative
    #F_BAET <- F_BAET*exp(1.23*(1-(Total.N[t-1]/K)))
    
    #Ricker model from Mathmatica
    F_BAET <- F_BAET*exp(-b * Total.N[t-1])
    
    # Beverton Holt from Mathmatica - can't go negative
    #if (Total.N[t-1] < K){
    #  F_BAET <- F_BAET*((K - Total.N[t-1])/K)
    #} else{
    #  F_BAET <- F_BAET*((K - (K-1))/K)
    #}
    # Beverton Holt - issue, can go negative
    #if (Total.N[t-1] < K){
    #F_BAET <- F_BAET*((r*Total.N[t-1])/(1 - (Total.N[t-1]*(r-1)/K)))
    #} else{
    #F_BAET <- F_BAET*((K - (K-1))/K)
    #}
    
    # regular old logistic equation - issue, can go negative
    #if (Total.N[t-1] < K){# we can also just imagine 'normal' logistic growth of the K - N/K
    #F_BAET <- F_BAET*(1 - Total.N[t-1]/K)
    #} else {
    # F_BAET <- F_BAET*((1 - (K-1))/K)
    #}
    
    
    #Baetidae
    
    # if drought
    # if (Q[t-1] < 8000){
    #   BAET1 <- c(DP1_BAET, 0, F_BAET)
    #   BAET2 <- c(DG1_BAET, DP2_BAET, 0)
    #   BAET3 <- c(0,DG2_BAET, 0)
    # }
    # 
    
    BAET1 <- c(P1_BAET, 0, F_BAET)
    BAET2 <- c(G1_BAET, P2_BAET, 0)
    BAET3 <- c(0, G2_BAET, 0) 
    # 
    # if (Q[t-1] >= 10000 & Q[t-1] < 15000){
    #   BAET1 <- c(SP1_BAET, 0, F_BAET)
    #   BAET2 <- c(SG1_BAET, SP2_BAET, 0 )
    #   BAET3 <- c(0, SG2_BAET, 0)
    # }
    # 
    # if (Q[t-1] >= 15000){
    #   BAET1 <- c(LP1_BAET, 0, F_BAET)
    #   BAET2 <- c(LG1_BAET, LP2_BAET , 0)
    #   BAET3 <- c(0, LG2_BAET, 0)
    # }
    
    ABAET <- rbind( BAET1, BAET2, BAET3)
    
    # if the water temp is warmer than the average water temp, then development is favored over growth
    # if (temps$Temperature[t-1] > mean_temp){
    #   ABAET[3,2] = ((temps$Temperature[t-1] - mean_temp)/mean_temp * ABAET[3,2]) + ABAET[3,2]
    #   ABAET[2,1] = ABAET[2,1] - ((temps$Temperature[t-1] - mean_temp)/mean_temp * ABAET[2,1]) 
    # }
    # 
    # # if the water temp is cooler than the average water temp, then growth is favored over development
    # if (temps$Temperature[t-1] < mean_temp){
    #   ABAET[2,1] = ((mean_temp - temps$Temperature[t-1])/mean_temp * ABAET[2,1]) + ABAET[2,1]
    #   ABAET[3,2] = ABAET[3,2] - ((temps$Temperature[t-1] - mean_temp)/mean_temp * ABAET[3,2])
    # }
    # 
    output.N.list[t, 1:3, 1] <- output.N.list[t-1, 1:3,1] %*% ABAET
    
    # immediate mortality due to flows
    # assume that mortality 
    # mortality due to flooding follows N0 = Nz*e^-hQ
    # but what about using a sigmoidal logistic function so we can control the threshold point and rate of growth
    # following m = 1/1+e^-h*(x-xf)
    # where h is is shape value
    # x is Q, and xf is threshold point (100% of pop dies)
    #plot(x = Q, y = 1/(1+exp(-0.02*(Q-100000))))
    
    
    #s1
    output.N.list[t, 1, 1] <- output.N.list[t, 1, 1] - (Q[t-1] * 1/(1+exp(-0.02*(Q[t-1]-100000))))
    #s2
    output.N.list[t,2,1] <- output.N.list[t,2,1] - (Q[t-1] * 1/(1+exp(-0.02*(Q[t-1]-100000))))
    #3
    output.N.list[t,3,1] <- output.N.list[t,3,1] - (Q[t-1] * 1/(1+exp(-0.02*(Q[t-1]-100000))))
    #replist[[1]][,,1] <- output.N.list[[1]]
    Total.N[,iter] <- apply(output.N.list,1,sum)
  }
}

# take a look at results

plot(timestep, Total.N[2:(length(timestep)+1)], type= "l", ylab = "Baetis spp. Total N", xlab = "Timestep (1 fortnight")

plot(timestep[200:length(timestep)], Total.N[201:(length(timestep)+1)], type= "l", ylab = "Baetis spp. Total N", xlab = "Timestep (1 fortnight")

Total.N
