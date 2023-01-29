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
source("ColoradoRiverTempRamp.R")
# growth - according to _ growth is dependent on shell length, as is fecundity 
# we want 1 class of non-reproductive subadults (smaller than 3.2 mm) and two classes of reproductive adults (3.2 mm and larger)
# growth follows the equation growth per day = -0.006*length[t]+0.029 (Cross et al., 2010)
# this can be re-written as length[t+1] = 0.994*length[t]+0.029
shell.growth(0.994, 0.029, 0.5)
# on day 164, length ~ 3.2 (maturity) - so on week 24, maturity reached [stage 1 = 0.5mm - 3.2 mm] in about 24 weeks or 12 timesteps
# on day 266, length ~ 3.9538776 - [stage 2 = 3.2 - 3.9538776] in about 14 weeks or 7 timesteps
# that means [stage 3 = 3.9538776+] for about 14 weeks or 7 timesteps


#read in flow data from USGS gauge at Lees Ferry, AZ between 1985 to the end of the last water year
flow <- readNWISdv("09380000", "00060", "1985-10-01", "2021-09-30")
# read in temperature data from USGS gauge at Lees Ferry, AZ between _ to the end of last water year
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2021-09-30")

#flow.magnitude <- TimestepDischarge(flow, 85000) #discharge data from USGS is in cfs, bankfull dischage (pers comm TK) 85000 cfs
temps <- TimestepTemperature(temp, "Colorado River") # calculate mean temperature data for each timestep
degreedays <- TimestepDegreeDay(temp, "Colorado River")

## Uncomment if Using Colorado River Temp Ramp 
temps <- temp_seq
degreedays <- temps$Temperature * 14



# specify iterations
iterations <- 5

# baseline K in the absence of disturbance
Kb <- 10000
# max K after a big disturbance
Kd <- 40000

# specify baseline transition probabilities
# its speculated (Cross et al 2010) that survivorship is between 80 - 100% for NZMS in Grand Canyon - will say 90% survive, 10% baseline mortality
# from timestep to timestep, we expect 90% to survive so for stage 1 (which lasts aproximately 14 timesteps, survival should be approx 09^14 = 0.2287679
# from that, only 1/14th will transition out, 13/14 remain in stage

# stage 1 G1 (prob transition to stage 2) = 0.95^14/14
# stage 1 P1 (prob remaining in stage 2) = 0.95^14 -0.95^14/14
# stage 2 G2 (prob transition to stage 3) = 0.95^7/7
# stage 2 P2 (prob remaining in stage 2) = 0.95^7 - 0.95^7/7
# stage 3 P3 (prob remaining in stage 3) = 0.95 ^5 - 0.95^5/5

G1_NZMS = 0.1
G2_NZMS = 0.1
P1_NZMS = 0.4
P2_NZMS = 0.5
P3_NZMS = 0.8

# want to run this for one year, in 14 day timesteps 
#timestep <- seq(2, (length(flow.magnitude$Discharge) + 1), by = 1) # OR
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

# Q is equal to average discharge over 14 days
#Q <- flow.magnitude$Discharge

Q <- rep(0.1, length(temps$Temperature))

Qmin <- 0.25
a <- 0.1
g <- 0.1
h <- surv.fit.NZMS$m$getPars()[2]   
k <- surv.fit.NZMS$m$getPars()[1]
e = 2.71828
extinction <- 1

#-------------------------
# Outer Loop of Iterations
#--------------------------


for (iter in c(1:iterations)) {
  
  K = 10000 # need to reset K for each iteration
  
  # we can also create a random flow scenario by sampling flows
  #out_sample <- sample(out$Discharge,length(out$Discharge), replace=TRUE)
  #Q <- out_sample
  
  # another option is to keep flow the same each timestep (to test temp effects)
  #out_same <- rep(10000, length(out$Discharge))
  #Q <- out_same
  
  
  # need to assign starting value
  # we can start with 10 S1 individuals for each species
  #for (sp in species){
  #  output.N.list[[sp]][1,1] <- 10
  #}
  # or we can pull random values from a uniform distribution 
  output.N.list[1,1:3, iter]<- runif(3, min = 1, max = (0.5*K))
  
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
    
    # fecundities estimated from McKenzie et al. 2013; 
    if (temps$Temperature[t-1] <= 10) { 
      F2_NZMS <- 2
      F3_NZMS <- 2  
    } else {
      F2_NZMS <- 8.87473
      F3_NZMS <- 27.89665
    }
    
    #---------------------------------------------------
    # Calculate the disturbance magnitude-K relationship 

    # Calculate the disturbance magnitude-K relationship
    # Sets to 0 if below the Qmin
    Qf <- Qf.Function(Q[t-1], Qmin, a)
    
    #-------------------------------------------------------------------
    # Calculate K arrying capacity immediately following the disturbance
    K0 <- K + ((Kd-K)*Qf)
    
    # Calculate final K for timestep, including relationship between K and time since disturbance
    K <- post.dist.K(K0, Kb, g)
    Klist <- append(Klist, K)
    #---------------------------------------------
    # Calculate effect of density dependnce on fecundity 
    
    # Logistic Density Dependence on Fecundity via Rogosch et al. Fish Model \
    # assume 97% Female poplation, and 30% instant mortality for eggs laid
    F2_NZMS <- Logistic.Dens.Dependence(F2_NZMS, K, Total.N[t-1, iter]) * 0.97 * 0.7
    F3_NZMS <- Logistic.Dens.Dependence(F3_NZMS, K, Total.N[t-1, iter]) * 0.97 * 0.7
    
    Flist <- append(Flist, (F2_NZMS + F3_NZMS))
    #-----------------------------------------------
    # Create Lefkovitch Matrix
    
    NZMS1 <- c(P1_NZMS, F2_NZMS, F3_NZMS)
    NZMS2 <- c(G1_NZMS, P2_NZMS, 0)
    NZMS3 <- c(0, G2_NZMS, P3_NZMS) 
    
    ANZMS <- rbind( NZMS1, NZMS2, NZMS3)
    
    #-----------------------------------------------
    # Calculate new transition probabilities based on temperature
    # This is the growth v development tradeoff
    #--------------------------------------
    # Calculate abundances for each stage
    
    output.N.list[t, 1:3, iter] <- ANZMS %*% output.N.list[t-1, 1:3, iter] 
    
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
    extinction.threshold(extinction)
  
    } #-------------------------
  # End Inner Loop  
  #------------------------- 
} #----------------------
# End Outer Loop
#----------------------

#------------------
# Analyzing Results
#-------------------
# summarizing iterations

## turning replist into a df
means.list.NZMS <- mean.data.frame(output.N.list, stages = c(1,2,3), burnin = 25)
means.list.NZMS <- cbind(means.list.NZMS[2:339,], temps$dts)
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
abund.trends.NZMS <- ggplot(data = means.list.NZMS, aes(x = `temps$dts`,
                                              y = mean.abund, group = 1)) +
  # geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
  #                 ymax = mean.abund + 1.96 * se.abund),
  #             colour = 'transparent',
  #             alpha = .5,
  #             show.legend = FALSE) +
  geom_line(show.legend = FALSE, linewidth = 0.7) +
  coord_cartesian(ylim = c(0,15000)) +
  ylab('New Zealand Mudsnail Abundance') +
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



saveRDS(abund.trends, paste0('BAETplot', '.rds'))



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
