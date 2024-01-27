##########################
# D sp model
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
# library(purrr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
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

source("1spFunctions.R")

#------------------------------------------------------------
# Set up location specific data
#-----------------------------------------------------------
#if looking at ColRiver temps read in temperature data from USGS gauge at Lees Ferry, AZ between _ to the end of last water year
# temp <- readNWISdv("09380000", "00010", "2007-10-01", "2021-09-30")
#-------------------------------------------------------------
# create a series of average 2 weekly data
#-------------------------------------------------------------
# calculate mean temperature data for each timestep (2 week interval)
# temps <- average.yearly.temp(temp, "X_00010_00003", "Date")

# Time <- c(1:1825)
# Date <- rep(c(1:365), times = 5)
# Day <- seq(as.Date("2022-01-01"), as.Date("2026-12-31"), by="days")
# Day <- Day[-which(Day == "2024-02-29")]
# 
# Temperature <-  -7.374528  * (cos(((2*pi)/365)*Date))  +  (-1.649263* sin(2*pi/(365)*Date)) + 10.956243
# 
# temp <- as.data.frame(cbind(Time, Day, Date, Temperature))
# peaklist <- 0 
# peakeach <- length(temp$Temperature)
# iteration <- 1
# baselineK <- 10000
# disturbanceK <- 40000

Dmodel <- function(flow.data, temp.data, baselineK, disturbanceK, Qmin, extinct, iteration, peaklist = NULL, peakeach = NULL, fecundity = 300, dds = 1500){
  
  # set up model
  source("NegExpSurv.R")

  Q <- as.numeric(flow.data)
  temps <- temp.data
  
  degreedays <- as.data.frame(cbind(temps$dts, temps$Temperature * 14))
  colnames(degreedays) <- c("dts", "DegreeDay")
  degreedays$dts <- as.Date(degreedays$dts, origin = "1970-01-01")
  
  # need to make ramped increasing hydropeaking index 
  hp <- c(rep(peaklist, each = peakeach))
  
  # specify iterations
  iterations <- iteration
  
  # baseline K in the absence of disturbance
  Kb <- as.numeric(baselineK)
  # max K after a big disturbance
  Kd <- as.numeric(disturbanceK)
  
  # specify baseline transition probabilities for each species at mean temps
  G1 =  0.109#move to Stage2 (subimago)
  G2 =  0.109 #move to Stage3 (adult)
  P1 =  0.81#stay in Stage1 (larvae)
  P2 =  0.81#stay in Stage2 (subimago)
  P3 = 0.6
  
  # want to run this for one year, in 14 day timesteps 
  timestep <- seq(2, (length(temps$Temperature) + 1), by = 1)
  
  # create an array to put our output into
  output.N.array <- array(0, dim = c(length(timestep) + 1))
  
  output.N.list <- list(output.N.array)
  
  # create array to put the total N of all species into
  Total.N <- array(0,
                   dim  <-c((length(timestep) +1 ), iterations),
                   dimnames <- list(1:(length(timestep) + 1), 1:iterations))
  
  # create list of arrays w/ abundance data for each spp
  reparray <- array(0,
                    
                    dim = c(length(timestep) + 1, 3, iterations),
                    dimnames = list(1:(length(timestep)+1), c("S1", "S2", "S3"), 1:iterations)
  )
  
  output.N.list <- reparray
  
  Qmin <- Qmin
  a <- 0.001
  g <- 1
  h <- low$m$getPars()[2]  
  k <- low$m$getPars()[1] 
  
  extinction <- extinct

  #-------------------------
  # Outer Loop of Iterations
  #--------------------------
  
  
  for (iter in c(1:iterations)) {
    K = Kb # need to reset K for each iteration
    
    # pull random values from a uniform distribution 
    output.N.list[1,1:3, iter]<- runif(3, min = 1, max = (0.3*K))
    
    # we often want to look at different parameter values after we run code, so we create some lists
    
    # list to input Ks
    Klist <- vector()
    Klist[1] <- 10000
    
    # list to imput flow morts
    flowmortlist <- vector()
    
    Flist <- vector()
    
    emergetime <- vector()
    
    sizelist <- vector()
    
    TempSurvival <- vector()
    for(c in temps$Temperature){
      
      b <- TempSurv(c)
      
      TempSurvival <- append(TempSurvival, b)
    }
    
    #-------------------------
    # Inner Loop of Timesteps
    #-------------------------
    for (t in timestep) {
      
      #----------------------------------------------------------
      # Calculate how many timesteps emerging adults have matured
      
      
      emergetime <- append(emergetime, back.count.degreedays(t, dds, degreedays)) # value from Sweeney et al 2017
      #---------------------------------------------------------
      # Calculate fecundity per adult
      
      # we start by pulling fecundities from normal distribution
      # assuming 50 50 sex ration, 0.22 of egg masses 'dissapearred', and 0.2 desiccation because of rock drying
      F3 = fecundity * hydropeaking.mortality(0.9, 1, h = hp[t-1]) 
      #F3 = rnorm(1, mean = 1104.5, sd = 42.75) * 0.5  #Baetidae egg minima and maxima from Degrange, 1960, assuming 1:1 sex ratio and 50% egg mortality
      
      # we can also relate fecundities to body mass.
      # Sweeney and Vannote 1980 have recorded dry body weight between 0.9 and 2.0 mg. 
      # That weight is related to fecundity Y = 614X - 300
      # we can "convert" emergetime to mg by multiplying to get dry weights between 0.9 - 2 mg, and then convert to fecunity
      # Issue: this data is for Ephemerella spp, not Baetidae spp
      # 
      
      x <- c(5,13)
      y <- c(fecundity*0.9, fecundity*1.1)
      mod <- lm(y~x)
      
      
      if (is.na(emergetime[t-1]) == F) { # will be erased in burn
        size <- emergetime[t-1]
        sizelist <- append(sizelist, size)
        F3 <- ((size*mod$coefficients[2])+mod$coefficients[1]) * hydropeaking.mortality(0.4, 0.6, h = hp[t-1])
        #F3 <- (57*size)+506 * 0.5 * hydropeaking.mortality(0.0, 0.2, h = hp[t-1]) * 0.78 * 0.65
      }
      #--------------------------------------------------
      # Calculate the disturbance magnitude-K relationship
      # Sets to 0 if below the Qmin
      Qf <- Qf.Function(Q[t-1], Qmin, a)
      
      #-------------------------------------------------------------------
      # Calculate K carrying capacity immediately following the disturbance
      K0 <- K + ((Kd-K)*Qf)
      
      # Calculate final K for timestep, including relationship between K and time since disturbance
      K <- post.dist.K(K0, Kb, g, t, Q, Qmin)
      
      Klist <- append(Klist, K)
      #---------------------------------------------
      # Calculate effect of density dependence on fecundity 
      
      # Logistic via Rogosch et al. Fish Model
      # no immediate egg mortality incorporated
      F3 <- Logistic.Dens.Dependence(F3, K, Total.N[t-1, iter])
      # 
      # add F_BAET to list
      Flist <- append(Flist, F3)
      #-----------------------------------------------
      # Calculate new transition probabilities based on temperature
      # This is the growth v development tradeoff
      
      # development measures
      # in this function, we assume that if below the min temp threshold (9) no maturation occurs (slow maturation, large growth)
      # if above the max temp threshold (15), no one remains more than 1 timestep in each stage (fast maturation, small growth)
      # if (5 > temps$Temperature[t-1]) {
      #   P1 <- (1-(1/8))*TempSurvival[t-1]
      #   P2 <- (1-(1/8))*TempSurvival[t-1]
      #   G1 <- (0.6/8)*TempSurvival[t-1]
      #   G2 <- (0.6/8)*TempSurvival[t-1]
      # }
      # 
      # if (temps$Temperature[t-1] > 20){
      #   P1 <- 1-(1/3)*TempSurvival[t-1]
      #   P2 <- 1-(1/3)*TempSurvival[t-1]
      #   G1 <- 0.6/3*TempSurvival[t-1]
      #   G2 <- 0.6/3*TempSurvival[t-1]
      # }
      # 
      
      if (is.na(emergetime[t-1]) == F){
        G1 <- (0.6/((emergetime[t-1]-1)/2))*TempSurvival[t-1]
        G2 <- (0.6/((emergetime[t-1]-1)/2))*TempSurvival[t-1]
        P1 <- (1-(1/((emergetime[t-1]-1)/2)))*TempSurvival[t-1]
        P2 <- (1-(1/((emergetime[t-1])/2)))*TempSurvival[t-1]
        P3 <- P3* TempSurvival[t-1]
      }
      
      if (is.na(emergetime[t-1]) == T) {
        G1 <- (0.6/((-0.588 * temps$Temperature[t-1]) + 18.764)) *TempSurvival[t-1]
        P1 <- (1-(1/((-0.588 * temps$Temperature[t-1]) + 18.764))) *TempSurvival[t-1]
        G2 <- (0.6/((-0.588 * temps$Temperature[t-1]) + 18.764))*TempSurvival[t-1]
        P2 <- (1-(1/((-0.588 * temps$Temperature[t-1]) + 18.764)))*TempSurvival[t-1]
        P3 <- P3*TempSurvival[t-1]
      }
      
      if (G1 > 1) G1 <- 1
      if (G1 < 0) G1 <- 0
      if (G2 > 1) G2 <- 1
      if (G2 < 0) G2 <- 0
      if (P1 > 1) P1 <- 1
      if (P2 < 0) P2 <- 0

      
    #-----------------------------------------------
      # Create Lefkovitch Matrix
      
      T1 <- c(P1, 0, F3)
      T2 <- c(G1, P2, 0)
      T3 <- c(0, G2, 0) 
      
      A <- rbind( T1, T2, T3)
      
      #--------------------------------------
      # Calculate abundances for each stage
      
      output.N.list[t, 1:3, iter] <- A %*% output.N.list[t-1, 1:3, iter] 
      
      #------------------------------------------
      #Calculate immediate mortality due to flows
      # mortality due to flooding follows N0 = Nz*e^-hQ
      
      #s1
      output.N.list[t, 1, iter] <- flood.mortality(output.N.list[t, 1, iter], k, h, Q[t-1], Qmin)
      #s2Q
      output.N.list[t,2,iter] <- flood.mortality(output.N.list[t,2,iter], k, h, Q[t-1], Qmin)
      
      output.N.list[t,3,iter] <- flood.mortality(output.N.list[t,3,iter], k, h, Q[t-1], Qmin)
      
      flowmortlist <- append(flowmortlist, flood.mortality(1, k, h, Q[t-1], Qmin))
      
      #------------------------------------------------------
      # check extinction threshold and if below set to 0
      Total.N[,iter] <- apply(output.N.list[,,iter],1,sum)
      if (Total.N[t, iter] < extinction){
        output.N.list[t,,iter] <- 0
        Total.N[t, iter] <- 0}
      
      
    } #-------------------------
    # End Inner Loop  
    #------------------------- 
  } #----------------------
  # End Outer Loop
  #----------------------
  return(output.N.list[ , 1:3, ])
}
