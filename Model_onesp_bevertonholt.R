#######
# 
#######
library(purrr)
library(tidyverse)
library(lubridate)
library(dplyr)

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

# set carrying capacity
K = 100

# specify baseline transition probabilities
G1 = 0.9
G2 = 0.89

# transition probabilites when there is lowered flow (Q<8000)
DG1 = 0.85
DG2 = 0.8

# transition probabilities when there is a higher flow (15000 > Q >  10000)
SG1 = 0.55
SG2 = 0.5

# transition probabilities when thre is a high flow (Q>15000)
LG1 = 0.29
LG2 = 0.2



# want to run this for one year, in 14 day timesteps 
timestep <- seq(2, (length(out$Discharge) + 1), by = 1)

# create an array to put our output into
output.N.array <- array(0, dim = c(length(timestep) + 1))

# create array to put the total N of all species into
Total.N <- array(0,
                 dim  <-c((length(timestep) +1 ), iterations),
                 dimnames <- list(1:(length(timestep) + 1), 1:iterations)
)

# create list of arrays w/ abundance data for each spp
reparray <- array(0,
                  
                  dim = c(length(timestep) + 1, 3, iterations),
                  dimnames = list(1:(length(timestep)+1), c("S1", "S2", "S3"), 1:iterations)
)
output.N.list <- reparray

#insert N0 abundance
output.N.list[1,1,] <- 10

# insert TotalNO abundace
Total.N[1] <- 10


# Q is equal to average discharge over 14 days
Q <- out$Discharge
Qmin <- 8000
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
  for (t in timestep){
    F_BAET = rnorm(1, mean = 1104.5, sd = 42.75) * H_BAET #Baetidae egg minima and maxima from Degrange, 1960
    
    
    # Calculate the disturbance magnitude-K relationship. Sets to 0 if below the Qmin
    if (Q[t-1] < Qmin) {
      Qf <- 0
    } else {
      Qf <- (Q[t-1] - Qmin)/(a + Q[t-1]- Qmin)
    }
    
    # Calculate K arrying capacity immediately following the disturbance
    K <- 100 + (400-100)*Qf
    
    # Function to calc. K as a function of time post-disturbance at a particular disturbance intensity
    K <- 100 + (400 - 100)*exp(-g*14)
    
    # Ricker model - pro = doesn't go negative
    #F_BAET <- F_BAET*exp(1.23*(1-(Total.N[t-1]/K)))
    
    #Ricker model from Mathmatica
    #F_BAET <- F_BAET*exp(-b * Total.N[t-1])
    
    # Beverton Holt from Mathmatica
    if (Total.N[t-1] < K){
      F_BAET <- F_BAET*((K - Total.N[t-1])/K)
    } else{
      F_BAET <- F_BAET*((K - (K-1))/K)
    }
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
    
 
    
#print(Total.N[t-1])
#print(F_BAET)
#print(Q)
#Baetidae

# if drought
if (Q[t-1] < 8000){
  BAET1 <- c(0, 0, F_BAET)
  BAET2 <- c(DG1, 0, 0)
  BAET3 <- c(DG2, 0, 0)
}

if (Q[t-1] >= 8000 & Q[t-1] < 10000){
  BAET1 <- c(0, 0, F_BAET)
  BAET2 <- c(G1, 0, 0)
  BAET3 <- c(0, G2, 0)
}

if (Q[t-1] >= 10000 & Q[t-1] < 15000){
  BAET1 <- c(0, 0, F_BAET)
  BAET2 <- c(SG1, 0, 0 )
  BAET3 <- c(0, SG2, 0)
}

if (Q[t-1] >= 15000){
  BAET1 <- c(0, 0, F_BAET)
  BAET2 <- c(LG1, 0, 0)
  BAET3 <- c(0, LG2, 0)
}

ABAET <- rbind( BAET1, BAET2, BAET3)

# if the water temp is warmer than the average water temp, then development is favored over growth
if (temps$Temperature[t-1] > mean_temp){
  ABAET[3,2] = ((temps$Temperature[t-1] - mean_temp)/mean_temp * ABAET[3,2]) + ABAET[3,2]
  ABAET[2,1] = ABAET[2,1] - ((temps$Temperature[t-1] - mean_temp)/mean_temp * ABAET[2,1]) 
}

# if the water temp is cooler than the average water temp, then growth is favored over development
if (temps$Temperature[t-1] < mean_temp){
  ABAET[2,1] = ((mean_temp - temps$Temperature[t-1])/mean_temp * ABAET[2,1]) + ABAET[2,1]
  ABAET[3,2] = ABAET[3,2] - ((temps$Temperature[t-1] - mean_temp)/mean_temp * ABAET[3,2])
  }

output.N.list[t, 1:3, 1] <- output.N.list[t-1, 1:3, 1] %*% ABAET
#replist[[1]][,,1] <- output.N.list[[1]]
Total.N[,iter] <- apply(output.N.list,1,sum)
}
}

# take a look at results
plot(timestep, Total.N[2:(length(timestep)+1)], type= "l")
plot(timestep, Q, type = "l")















repdf <- ldply(replist, function(x)
  adply(x, c(1,2,3))
})

names(repdf) <- c('Species', 'Timestep', 'Stage', 'Rep', 'Abund')
repdf$Timestep <- as.numeric(as.character(repdf$Timestep))

totn <- adply(Total.N, c(1,2))
names(totn) <- c('Timestep', 'Rep', 'Tot.abund')
totn$Timestep <- as.numeric(as.character(totn$Timestep))

## joining totn and repdf together
repdf <- left_join(totn, repdf)

## calculating relative abundance
repdf <- mutate(repdf, rel.abund = Abund/Tot.abund)

  }



F_BAET = rnorm(1, mean = 1104.5, sd = 42.75) * H_BAET #Baetidae egg minima and maxima from Degrange, 1960
if (Ntot < K){
  F_BAET = F_BAET*((K - Ntot)/K)
  } else {
    F_BAET = F_BAET*((K - 1)/K)
  }



output.N <- N %*% ABAET
Ntot = sum(output.N)

output.N.list[[sp]][t,1:3] <- output.N.list[[sp]][t-1, 1:3] %*% get(paste0('A', sp))
