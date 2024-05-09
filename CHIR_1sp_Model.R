##########################
# CHIRis 1 sp model
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
# #read in flow data from USGS gauge at Lees Ferry, AZ between 1985 to the end of the last water year
discharge <- readNWISdv("09380000", "00060", "2007-10-01", "2023-05-01")
flow.magnitude <- TimestepDischarge(discharge, 85000)
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
temps <- TimestepTemperature(temp)

flow.data = flow.magnitude$Discharge
temp.data <- temps
baselineK <- 10000
disturbanceK <- 40000
Qmin <- 0.25
extinct <- 50
iteration <- 1
peaklist <- 0.13
peakeach<- length(temps$Temperature)

CHIRmodel <- function(flow.data, temp.data, baselineK, disturbanceK, Qmin, extinct, iteration, peaklist = NULL, peakeach = NULL){
  
  # set up model
  source("1spFunctions.R")
  source("CHIRSurvivorship.R")
  Q <- as.numeric(flow.data)
  temps <- temp.data
  
  degreedays <- as.data.frame(cbind(temps$dts, temps$Temperature * 14))
  colnames(degreedays) <- c("dts", "DegreeDay")
  degreedays$DegreeDay[degreedays$DegreeDay<0] <- 0
  degreedays$dts <- as.Date(degreedays$dts, origin = "1970-01-01")
  
  
  # need to make ramped increasing hydropeaking index 
  hp <- c(rep(peaklist, each = peakeach))
  
  # specify iterations
  iterations <- iteration
  
  # baseline K in the absence of disturbance
  Kb <- as.numeric(baselineK)
  # max K after a big disturbance
  Kd <- as.numeric(disturbanceK)
  
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
  h <- surv.fit.CHIR$m$getPars()[2]  
  k <- surv.fit.CHIR$m$getPars()[1] 
  
  extinction <- extinct
  
  
  #-------------------------
  # Outer Loop of Iterations
  #--------------------------
  
  
  for (iter in c(1:iterations)) {
    K = Kb # need to reset K for each iteration
    
    # pull random values from a uniform distribution 
    output.N.list[1,1:3, iter] <- runif(3, min = 1, max = (0.3*K))
    
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
      t <- t
      emergetime <- append(emergetime, back.count.degreedays(t, 600, degreedays)) # value from Ali et al 1985
      #   if (t < 7){
      #     development <- append(development, delta[t-1]*(MaturationRate(devtime(temps$Temperature[t-1])/14)))
      #   } else {
      #   development <- append(development, development[t-delta[t]]*((MaturationRate(devtime(temps$Temperature[t-1])/14))/(MaturationRate(devtime(temps$Temperature[t-delta[t]])/14))))
      #   }
      # }
      #---------------------------------------------------------
      # Calculate fecundity per adult
      
      # we start by pulling fecundities from normal distribution
      # assuming 50 50 sex ration
      F3 = 208 *0.5* hydropeaking.mortality(0.8, 1, h = hp[t-1]) * 0.698896
      #CHIR egg # and % mortality from Charles et al 2004
      # we can also relate fecundities to body size which is between 6 and 15 mm (also from Charles et al 2004)
      # we can "convert" emergetime to size by multiplying to get size between 6 and 15 mm and then convert to fecunity
      
      if (t > 19) {
        size <- 3*emergetime[t-1]-6
        sizelist <- append(sizelist, size)
        F3 <- (4.622*size)+159.468 *0.5* hydropeaking.mortality(0.8, 1, h = hp[t-1]) * 0.698896
        }
      # #--------------------------------------------------
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
      # add F_CHIR to list
      Flist <- append(Flist, F3)
      #-----------------------------------------------
      # Calculate new transition probabilities based on temperature
      # This is the growth v development tradeoff
      # using Birt et al 2009 calcs
      # using survivals from Charles et al. 2004
      
      # development measures
      # in this function, we assume that if below the min temp threshold (9) slow maturation
      # if above the max temp threshold (30), no one remains more than 1 timestep in each stage (fast maturation, small growth)

  if (11 > temps$Temperature[t-1]) {
  P1 <- (1-(1/3)) * TempSurvival[t-1]
  P2 <- P1 
  G1 <- 0 * TempSurvival[t-1]
  G2 <- 0 * TempSurvival[t-1]
  }
if (temps$Temperature[t-1] > 30){
  P1 <- 0
  P2 <- 0
  G1 <- 0.42 * TempSurvival[t-1]
  G2 <- G1
  }

      if (11 <= temps$Temperature[t-1] & temps$Temperature[t-1] <= 30 & is.na(emergetime[t-1] == F)){
        G1 <- (0.42/((emergetime[t-1])/2)) * TempSurvival[t-1]
        G2 <- G1
        P1 <- (1-(1/((emergetime[t-1])/2))) * TempSurvival[t-1]
        P2 <- P1
      }
      if (11 <= temps$Temperature[t-1] & temps$Temperature[t-1] <= 30 & is.na(emergetime[t-1] == T)) {
        G1 <- (0.42*((-0.136 * temps$Temperature[t-1]) + 5.088)) * TempSurvival[t-1]
        P1 <- (1-(1/((-0.136 * temps$Temperature[t-1]) + 5.088))) * TempSurvival[t-1]
        G2 <- G1
        P2 <- P1
      }
      
      
      # if (5 <= temps$Temperature[t-1] & temps$Temperature[t-1] <= 21){
      #   G1 <- 0.32/((delta[t-1]/2))
      #   G2 <- G1
      #   P1 <- 0.32*(1-1/(delta[t-1]/2))
      #   P2 <- P1
      # }
      
      
      # #if temp dependent surival curve is know, we can also use that, applying either a) constant mortality to all stages or b) giving specific weights of survival to different stages
      # if (7 <= temps$Temperature[t-1] & temps$Temperature[t-1] <= 25) G1 <- growth.development.tradeoff(temps$Temperature[t-1], 7, 25, 0.15, 0.25)
      # if (7 <= temps$Temperature[t-1] & temps$Temperature[t-1] <= 25) G2 <- growth.development.tradeoff(temps$Temperature[t-1], 7, 25, 0.15, 0.25)
      
      # growth (if below 10C, no growth can occur - everything basically freezes, if between 10 and 11, prob of remaining in same stage = 0.6395, if above 13, prob of transition to next stage is 0 )
      # P1 <- 0.55 - G1
      # P2 <- 0.55 - G2
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
      #s2Qt
      output.N.list[t,2,iter] <- flood.mortality(output.N.list[t,2,iter], k, h, Q[t-1], Qmin)
      
      output.N.list[t,3,iter] <- flood.mortality(output.N.list[t,3,iter], k, h, Q[t-1], Qmin)
      
      flowmortlist <- append(flowmortlist, flood.mortality(1, k, h, Q[t-1], Qmin))
      
      #------------------------------------------------------
      # check extinction threshold and if below set to 0
      Total.N[t,iter] <- sum(output.N.list[t,,iter])
      if (Total.N[t,iter] < extinction){
        output.N.list[t,,iter] <- 0
        Total.N[t, iter] <- 0}
      
      
    } #-------------------------
    # End Inner Loop  
    #------------------------- 
  } #----------------------
  # End Outer Loop
  #----------------------
  return(output.N.list[ , 1:2, ])
}
#------------------
# Analyzing Results
#-------------------
# summarizing iterations
# 
# out <- CHIRmodel(flow.data = discharge, temp.data = temp, disturbanceK = 40000, baselineK = 10000, Qmin = 0.25, extinct = 50, iteration = 3, peaklist = 0, peakeach = length(temp$Temperature))
# #
# adults <-as.data.frame(cbind((temp$dts), out[1:length(temp$dts),3,1]))
# colnames(adults) <- c("Time","Adult CHIRidae")
# adults$Time <- as.Date(as.POSIXct(adults$Time, origin = "1970-01-01"))
# #
# ## turning replist into a df
# means.list.CHIR <- mean.data.frame(out, burnin = 90, iteration = 3)
# means.list.CHIR <- cbind(means.list.CHIR[1:length(means.list.CHIR$mean.abund),], temp$dts[90:132])
# means.list.CHIR$`temp$dts` <- as.Date(means.list.CHIR$`temp$dts`)
# 
# # note how boom and bust this model is - K is set to be 10,000 not 100,000
# abund.trends.CHIR <- ggplot(data = adults[90:131,], aes(x = Time,
#                                        y = `Adult CHIRidae`/10000, group = 1)) +
#   geom_point()+
#   # geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
#   #                 ymax = mean.abund + 1.96 * se.abund),
#   #             colour = 'transparent',
#   #            alpha = .5,
#   #             show.legend = FALSE) +
#   geom_line(show.legend = FALSE) +
#   ylab('Adult CHIRis Abundance/ Baseline Reproductive Limit') +
#   xlab('Timestep')+
#   scale_x_date(date_labels="%B", date_breaks  ="4 months")+
# theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5),
#       axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# #
# plot(out[80:105, 2], type = "b", xlab = "Timesteps Jan 1 to Dec 31", ylab = "CHIR adult abundace (max on Aug 23 after prolonged warm period)" )
# 
# 
# plot(adults$`Adult CHIRidae`[90:130], adults$`Adult CHIRidae`[91:131], type = "b", xlab = "Adult CHIRids t", ylab = "Adult CHIRids t+1")
# #
# ggplot(data = means.list.CHIR, aes(x = `temp$dts`,
#                           y = mean.abund/10000, group = 1)) +
#   geom_point()+
#   # geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
#   #                 ymax = mean.abund + 1.96 * se.abund),
#   #             colour = 'transparent',
#   #            alpha = .5,
#   #             show.legend = FALSE) +
#   geom_line(show.legend = FALSE) +
#   geom_line(aes(y = discharge[89:131]*10), color = "blue")+
#   ylab('CHIRis spp. Abundance/Reproductive Limit') +
#   xlab('Timestep')+
#   scale_x_date(date_labels="%B", date_breaks  ="4 months")+
#   theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5),
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# 
# #
# # plot(means.list.CHIR$mean.abund[300:500], means.list.CHIR$mean.abund[301:501], type = "b", xlab = "Nt", ylab = "Nt+1")
# #
# # ggplot(data = NULL, mapping = aes(x = temps$dts, y = Total.N[2:2003]/10000))+
# #   geom_line(show.legend = FALSE) +
# #   ylab('CHIRis spp. Abundance/Reproductive Limit') +
# #   xlab(" ")
# plot(temp$dts[90:131], temp$Temperature[90:131], col = "red", type = "l", xlab = "Time", ylab = "Temperature C")
# 
