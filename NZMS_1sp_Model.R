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
source("1spFunctions.R")
#read in flow data from USGS gauge at Lees Ferry, AZ between 1985 to the end of the last water year
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2021-09-30")
temps <- average.yearly.temp(temp, "X_00010_00003", "Date")


n <- 77
# qr is the temp ramps I want to increase the average Lees Ferry temp by 
# how many years I want each temp ramp to last
qr <- 0
r <- 77

temps <- rep.avg.year(temps, n, change.in.temp = qr, years.at.temp = r)
discharge <- rep(0.1, times = length(temps$Temperature))
#temps$Temperature <- rep(12, times = length(temps$dts))

NZMSmodel <- function(flow.data, temp.data, baselineK, disturbanceK, Qmin, extinct, iteration, peaklist = NULL, peakeach = NULL){
  # source functions
  source("NZMSSurvivorship.R")
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
  
  # specify baseline transition probabilities
  # its speculated (Cross et al 2010) that survivorship is between 80 - 100% for NZMS in Grand Canyon - will say 90% survive, 10% baseline mortality
  # from timestep to timestep, we expect 90% to survive so for stage 1 (which lasts aproximately 14 timesteps, survival should be approx 09^14 = 0.2287679
  # from that, only 1/14th will transition out, 13/14 remain in stage (Birt et al 2009)
  
  # stage 1 G1 (prob transition to stage 2)
  # stage 1 P1 (prob remaining in stage 2) 
  # stage 2 P2 (prob remaining in stage 2)
  # stage 3 P3 (prob remaining in stage 3) 
  
  G1 = 0.9/14
  G2 = 0.9/7
  P1 = 6/7
  P2 = 13/14
  P3 = 6/7 
  
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
  a <- 0.1
  g <- 0.1
  h <- surv.fit.NZMS$m$getPars()[2]  
  k <- surv.fit.NZMS$m$getPars()[1] 
  
  extinction <- extinct
#-------------------------
# Outer Loop of Iterations
#--------------------------


for (iter in c(1:iterations)) {
  
  K = Kb # need to reset K for each iteration
  # we can pull random values from a uniform distribution 
  output.N.list[1,1:3, iter]<- runif(3, min = 1, max = (0.3*K))
  
  # we often want to look at different parameter values after we run code, so we create some lists
  # list to input Ks
  Klist <- vector()
  Klist[1] <- 10000
  Flist <- vector()
  # list to imput flow morts
  flowmortlist <- vector()
  Flist <- vector()
  emergetime <- vector()
  sizelist <- vector()
  
  #-------------------------
  # Inner Loop of Timesteps
  #-------------------------
  
  for (t in timestep) {
    
    #---------------------------------------------------------
    # Calculate starting fecundity per adult
    #temps
    # fecundities estimated from McKenzie et al. 2013 - reduced fecundity above 24 C and below 9 C. 
    # optimal temp between 16 and 19 C, but we don't really have parameterization for that
# 
      F2 <- 8.87473 * (-0.0001427 *  (temps$Temperature[t-1] - 17.5)^4 + 1)

      F3 <- 27.89665 * (-0.0001427 * (temps$Temperature[t-1] - 17.5)^4 + 1)
# #       
#       F2 <- 8.87473
#       F3 <- 27.89665
    
    
    #---------------------------------------------------
    # Calculate the disturbance magnitude-K relationship
    # Sets to 0 if below the Qmin
    Qf <- Qf.Function(Q[t-1], Qmin, a)
    
    #-------------------------------------------------------------------
    # Calculate K arrying capacity immediately following the disturbance
    K0 <- K + ((Kd-K)*Qf)
    
    # Calculate final K for timestep, including relationship between K and time since disturbance
    K <- post.dist.K(K0, Kb, g, t, Q, Qmin)
    Klist <- append(Klist, K)
    #---------------------------------------------
    # Calculate effect of density dependnce on fecundity 
    
    # Logistic Density Dependence on Fecundity via Rogosch et al. Fish Model
    # assume 97% Female poplation, and 30% instant mortality for eggs laid
    F2 <- Logistic.Dens.Dependence(F2, K, Total.N[t-1, iter]) * 0.97 * 0.7
    F3 <- Logistic.Dens.Dependence(F3, K, Total.N[t-1, iter]) * 0.97 * 0.7
    
    Flist <- append(Flist, (F2 + F3))
    #-----------------------------------------------
    # Create Lefkovitch Matrix
    
    S1 <- c(P1, F2, F3)
    S2 <- c(G1, P2, 0)
    S3 <- c(0, G2, P3) 
    
    A <- rbind( S1, S2, S3)
    
    #-----------------------------------------------
    # Calculate new transition probabilities based on temperature
    # This is the growth v development tradeoff
    #--------------------------------------
    # Calculate abundances for each stage
    
    output.N.list[t, 1:3, iter] <- A %*% output.N.list[t-1, 1:3, iter] 
    
    #------------------------------------------
    #Calculate immediate mortality due to flows

    #s1
    output.N.list[t, 1, iter] <- flood.mortality(output.N.list[t,1,iter], k, h, Q[t-1], Qmin)
    #s2
    output.N.list[t,2,iter] <- flood.mortality(output.N.list[t,2,iter], k, h, Q[t-1], Qmin)
    #3
    output.N.list[t,3,iter] <- flood.mortality(output.N.list[t,3,iter], k, h, Q[t-1], Qmin)
    
    # in case we want to look back on what the flow mortality rates were
    flowmortlist <- append(flowmortlist, flood.mortality(1, k, h, Q[t-1], Qmin))

    #-------------------------------------------------
    # Calculate sum of all stages (total population)
    Total.N[,iter] <- apply(output.N.list[,,iter],1,sum)
    # check extinction threshold
    if (Total.N[t, iter] < extinction){
      output.N.list[t,,iter] <- 0
      Total.N[t, iter] <- 0
    } #-------------------------
  # End Inner Loop  
  #------------------------- 
} #----------------------
# End Outer Loop
#----------------------
}
  return(output.N.list)
}

out <- NZMSmodel(flow.data = discharge, temp.data = temps, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 500, iteration = 1, peaklist = 0, peakeach = length(temps$Temperature))
#------------------
# Analyzing Results
#-------------------
# summarizing iterations

## turning replist into a df
means.list.NZMS <- mean.data.frame(out,burnin = 1, iteration= 1)
means.list.NZMS <- cbind(means.list.NZMS, temps$dts[1:length(means.list.NZMS$timesteps)])
means.list.NZMS$`temps$dts` <- as.Date(means.list.NZMS$`temps$dts`)
# plot abundance over time

arrows <- tibble(
  x1 = c("2005-01-07", "2007-01-07", "2009-01-07", "2011-01-07"),
  x2 = c("2005-01-07", "2007-01-07", "2009-01-07", "2011-01-07"),
  y1 = c(14500, 14500, 14500, 14500), 
  y2 = c(10000, 12500, 12500, 12500)
)

arrows$x1 <- as.Date(arrows$x1)
arrows$x2 <- as.Date(arrows$x2)

means.list.NZMS <- means.list.NZMS[1300:1501,]

abund.trends.NZMS <- ggplot(data = means.list.NZMS, aes(x = `temps$dts`,
                                              y = mean.abund/10000, group = 1)) +
  # geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
  #                 ymax = mean.abund + 1.96 * se.abund),
  #             colour = 'transparent',
  #             alpha = .5,
  #             show.legend = FALSE) +
  geom_line(show.legend = FALSE, linewidth = 0.7) +
  geom_point()+
  coord_cartesian(ylim = c(0,1.5)) +
  ylab('New Zealand Mudsnail Abundance/Recruitment Limit') +
  xlab(" ")+
  theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13))+
  scale_x_date(date_labels="%B", date_breaks  ="6 months")+
  annotate("segment", x = arrows$x1, y = arrows$y1, xend = arrows$x2, yend = arrows$y2,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")), color = "red")+
  annotate("text", x = arrows$x1[1], y = 15000, label = "+1°C", size = 5)+
  annotate("text", x = arrows$x1[2], y = 15000, label = "+2.5°C", size = 5)+
  annotate("text", x = arrows$x1[3], y = 15000, label = "+5°C", size = 5)+
  annotate("text", x = arrows$x1[4], y = 15000, label = "+7.5°C", size = 5 )




means.list.NZMS$`temps$dts` <- format(as.Date(means.list.NZMS$`temps$dts`), "%Y-%m")

# take a look at results
# 
# par(mfrow = c(1,1))
# plot(timestep[9:(length(timestep)+1)], output.N.list[9:(length(timestep)+1), 3, 1], type = "l", ylab = "Baetis spp. Adults", xlab = "Timestep (1 fortnight)")
plot(timestep[9:length(timestep)], Total.N[10:(length(timestep)+1)], type= "l", ylab = "New Zealand Mudsnails. Total N", xlab = "Timestep (1 fortnight)", ylim = c(0,18000))
abline(v = 16, col = "red")
abline(v = 30, col = "red")
abline(v = 48, col = "red")
abline(v = 274, col = "red")
abline(v = 298, col = "red")
abline(v = 316, col = "red")
abline(v = 334, col = "red")
abline(v = 381, col = "red")
abline(v = 586, col = "red")
abline(v = 670, col = "red")
abline(v = 683, col = "red")
lines(timestep[9:length(timestep)], Klist[10:(length(timestep)+1)], type = "l", col = "blue")
lines(timestep, Klist, type = "l", col = "blue")


ggplot(data = NULL, mapping = aes(x = temps$dts, y = Total.N[2:2003]/10000))+
  geom_line(show.legend = FALSE) +
  ylab('Hydrospyche spp. Abundance/Reproductive Limit') +
  xlab(" ")

os <- as.data.frame(cbind(means.list.NZMS$mean.abund[1:200], means.list.NZMS$mean.abund[2:201]))

plot(os$V1, os$V2, type = "b", xlab = "Nt", ylab = "Nt+1", main = "Snail Abundance")
