#####################################
# Ramped hydropeaking index
#####################################


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

source("HYOSSurvivorship.R")
source("1spFunctions.R")

#read in flow data from USGS gauge at Lees Ferry, AZ between 1985 to the end of the last water year
flow <- readNWISdv("09380000", "00060", "1985-10-01", "2021-09-30")
flow.magnitude <- TimestepDischarge(flow, 85000) #discharge data from USGS is in cfs, bankfull dischage (pers comm TK) 85000 cfs
# read in temp data
#temp <- temp <- readNWISdv("09380000", "00010", "2007-10-01", "2021-09-30")
temp <- read.delim("gcmrc20230123125915.tsv", header=T)

#---------------------------------------------------------------
# this chunk of code makes the repeating avg ColRiv temp series
colnames(temp) <- c("Date", "Temperature")
temps <- average.yearly.temp(temp, "Temperature", "Date")
n <- 13
# qr is the temp ramps I want to increase the average Lees Ferry temp by 
qr <- c(0, 0, 0, 0, 0)
# how many years I want each temp ramp to last
r <- c(5, 2, 2, 2, 2)
temps <- rep.avg.year(temps, n, change.in.temp = qr, years.at.temp = r)

degreedays <- as.data.frame(cbind(temps$dts, temps$Temperature * 14))
colnames(degreedays) <- c("dts", "DegreeDay")
degreedays$dts <- as.Date(degreedays$dts, origin = "1970-01-01")

#-------------------------------------------------------------
# need to make ramped increasing hydropeaking index 
hp <- c(rep(0, 130), rep(0.01, 52), rep(0.1, 52), rep(0.2, 52), rep(0.5, 52))
# specify iterations
iterations <- 50

# baseline K in the absence of disturbance
Kb <- 10000
# max K after a big disturbance
Kd <- 40000


# specify baseline transition probabilities for each species
# 3 stages - we have egg - larval instar V, pupae, and adult

G1_HYOS = 0.1 # move onto stage 2
G2_HYOS = 0.445 # move onto stage 3
P1_HYOS = 0.7 # remain in stage 1
P2_HYOS = 0.0 # remain in stage 2

# want to run this for one year, in 14 day timesteps 
#timestep <- seq(2, (length(flow.magnitude$Discharge) + 1), by = 1) # OR
timestep <- seq(2, (length(temps$Temperature) + 1), by = 1)
#timestep <- seq(2, (length(out_sample) + 1), by = 1)

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

output.N.list <- reparray
## Repeating the array 7 times 
#replist <- rep(list(reparray), 3)
#names(replist) <- species


# Q is equal to average discharge over 14 days
#Q <- flow.magnitude$Discharge
Q <- rep(0.1, length(temps$Temperature))
# Fall floods
Q[c(23, 49, 75, 101, 127, 153, 179)] <- 0.45
Q[c(189, 215, 241, 267, 293, 319)] <-  0.45
Qmin <- 0.25
a <- 0.1
g <- 0.1
h <- surv.fit.HYOS$m$getPars()[2]   
k <- surv.fit.HYOS$m$getPars()[1] 
extinction <- 500

#-------------------------
# Outer Loop of Iterations
#--------------------------


for (iter in c(1:iterations)) {
  K = 10000 # need to reset K for each iteration
  # we can also create a random flow scenario by sampleing flows
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
  
  # or we can pull randomw values from a uniform distribution 
  output.N.list[1,1:3, iter]<- runif(3, min = 1, max = (0.5*K))
  
  # we often want to look at different parameter values after we run code, so we create some lists
  
  # list to input Ks
  Klist <- vector()
  Klist[1] <- 10000
  
  # list to imput flow morts
  flowmortlist <- vector()
  
  Flist <- vector()
  
  emergetime <- vector()
  
  sizelist <- vector()
  #-------------------------
  # Inner Loop of Timesteps
  #-------------------------
  
  for (t in timestep) {
    
    #----------------------------------------------------------
    # Calculate how many timesteps emerging adults have matured
    # Calculate how many timesteps emerging adults have matured
    
    emergetime <- append(emergetime, back.count.degreedays(t, 1680))
    #---------------------------------------------------------
    # Calculate fecundity per adult
    
    F_HYOS = rnorm(1, mean = 235.6, sd = 11.05102 ) * 0.5 #* hydropeaking.mortality(lower = 0.1, upper = 0.4, h = hp[t-1])
    #from Willis Jr & Hendricks, sd calculated from 95% CI = 21.66 = 1.96*sd
    # * 0.5 assuming 50% female
    
    # we can scale fecundity based on the 95% CI of 21.66 (min = 213.94, max = 257.26) 
    if (t > 15) {
      size <- emergetime[t-1]
      sizelist <- append(sizelist, size)
      F_HYOS <- ((8.664 * size) + 127.3) * 0.5 #* hydropeaking.mortality(lower = 0.1, upper = 0.4, h = hp[t-1])
    }
    
    #---------------------------------------------------
    # Calculate the disturbance magnitude-K relationship 
    # Sets to 0 if below the Qmin
    # Calculate the disturbance magnitude-K relationship 
    # Sets to 0 if below the Qmin
    Qf <- Qf.Function(Q[t-1], Qmin, a)
    
    #-------------------------------------------------------------------
    # Calculate K arrying capacity immediately following the disturbance
    # Calculate K arrying capacity immediately following the disturbance
    K0 <- K + ((Kd-K)*Qf)
    
    # Calculate final K for timestep, including relationship between K and time since disturbance
    K <- post.dist.K(K0, Kb, g)
    Klist <- append(Klist, K)
    
    
    #---------------------------------------------
    # Calculate effect of density dependence on fecundity
    
    # Logistic via Rogosch et al. Fish Model
    F_HYOS <- Logistic.Dens.Dependence(F_HYOS, K, Total.N[t-1, iter]) * 0.5
    Flist <- append(Flist, F_HYOS)
    #-----------------------------------------------
    # Create Lefkovitch Matrix
    
    HYOS1 <- c(P1_HYOS, 0, F_HYOS)
    HYOS2 <- c(G1_HYOS, P2_HYOS, 0)
    HYOS3 <- c(0, G2_HYOS, 0) 
    
    AHYOS <- rbind( HYOS1, HYOS2, HYOS3)
    
    #-----------------------------------------------
    # Calculate new transition probabilities based on temperature
    # This is the growth v development tradeoff 
    # don't know if this exists for HYOS - they can emerge under a wide temp gradient (<5 - 25+ C) but relationship between growth and temp 
    
    # development measures (basically, if below 10C, no development, if between 10 and 12, follows a function, if above 12, prob of transition to next stage is 0.6395)
    # if (10 > temps$Temperature[t-1]) AHYOS[3,2] <- 0.001
    # if (temps$Temperature[t-1] > 13) AHYOS[3,2] <- 0.55
    # if (10 <= temps$Temperature[t-1] & temps$Temperature[t-1] <= 13) AHYOS[3,2] <- (0.183 * temps$Temperature[t-1]) -1.829 #(0.2 * temps$Temperature[t-1]) -2
    # AHYOS[2,1] <- AHYOS[3,2] 
    # # 
    # # # growth (if below 10C, no growth can occur - everything basically freezes, if between 10 and 11, prob of remaining in same stage = 0.6395, if above 13, prob of transition to next stage is 0 )
    # if (10 >= temps$Temperature[t-1]) AHYOS[2,2] <- 0.55
    # if (temps$Temperature[t-1] > 13) AHYOS[2,2] <- 0.001
    # if (10 < temps$Temperature[t-1] & temps$Temperature[t-1] <=  13) AHYOS[2,2] <- (-0.183 * temps$Temperature[t-1]) + 2.38 #(-0.1 temps$Temperature[t-1]) - 2.6
    # # 
    # AHYOS[1,1] <- AHYOS[2,2] 
    
    # 
    #Glist <-append(Glist, AHYOS[3,2])
    #Plist <- append(Plist, AHYOS[2,2])
    
    #--------------------------------------
    # Calculate abundances for each stage
    
    output.N.list[t, 1:3, iter] <- AHYOS %*% output.N.list[t-1, 1:3, iter] 
    
    #------------------------------------------
    #Calculate immediate mortality due to flows
    # mortality due to flooding follows N0 = Nz*e^-hQ
    # but what about using a sigmoidal logistic function so we can control the threshold point and rate of growth
    # following m = 1/1+e^-h*(x-xf)
    # where h is is shape value
    # x is Q, and xf is threshold point (100% of pop dies)
    #plot(x = Q, y = 1/(1+exp(-0.02*(Q-100000))))
    
    #s1
    output.N.list[t, 1, iter] <- flood.mortality(output.N.list[t, 1, iter], k, h, Q[t-1], Qmin)
    #s2
    output.N.list[t,2,iter] <- flood.mortality(output.N.list[t,2,iter], k, h, Q[t-1], Qmin)
    #3
    output.N.list[t,3,iter] <- flood.mortality(output.N.list[t,3,iter], k, h, Q[t-1], Qmin)
    
    flowmortlist <- append(flowmortlist, flood.mortality(1, k, h, Q[t-1], Qmin))
    #replist[[1]][,,1] <- output.N.list[[1]]
    Total.N[,iter] <- apply(output.N.list[,,iter],1,sum)
    #check extinction threshold
    if (Total.N[t, iter] < extinction){
      output.N.list[t,,iter] <- 0
      Total.N[t, iter] <- 0  } #-------------------------
  }  # End Inner Loop  
  #------------------------- 
} #----------------------

# End Outer Loop
#----------------------

#------------------
# Analyzing Results
#-------------------
# summarizing iterations
means.list.HYOS <- mean.data.frame(output.N.list, burnin = 27)
means.list.HYOS <- cbind(means.list.HYOS[27:339,], temps$dts[27:339])
means.list.HYOS$`temps$dts` <- as.Date(means.list.HYOS$`temps$dts`)
# plot abundance over time

falls <- tibble(
  x1 = c("2001-11-10", "2002-11-10", "2003-11-10", "2004-11-10", "2005-11-10", "2006-10-10"),
  x2 = c("2001-11-10", "2002-11-10", "2003-11-10", "2004-11-10", "2005-11-10", "2006-10-10"),
  y1 = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
  y2 = c(0.23, 0.23, 0.23, 0.23, 0.23, 0.23)
 )
springs <- tibble(
  x1 = c("2007-03-31", "2008-03-31", "2009-03-31", "2010-03-31", "2011-03-31", "2012-03-31"),
  x2 = c("2007-03-31", "2008-03-31", "2009-03-31", "2010-03-31", "2011-03-31", "2012-03-31"),
  y1 = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
  y2 = c(0.23, 0.23, 0.23, 0.23, 0.23, 0.23)
)


# arrows <- tibble (2013-03-31"
#   x1 = c("2003-01-07"),
#   x2 = c("2011-01-07"), 
#   y1 = c(1.0),
#   y2 = c(1.0)
# )
falls$x1 <- as.Date(falls$x1)
falls$x2 <- as.Date(falls$x2)
springs$x1 <- as.Date(springs$x1)
springs$x2 <- as.Date(springs$x2)

abund.trends.HYOS <- ggplot(data = means.list.HYOS, aes(x =  `temps$dts`,
                                                        y = mean.abund/10000, group = 1)) +
  #geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
  #                ymax = mean.abund + 1.96 * se.abund),
  #            colour = 'transparent',
  #            alpha = .5,
  #            show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  coord_cartesian(ylim = c(0,1.5)) +
  ylab('Hydrospyche spp. Abundance') +
  xlab(" ")+
  theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13))+
  scale_x_date(date_labels="%B", date_breaks  ="6 months", limits = as.Date(c("2001-01-27", "2012-12-04"
  )))+
  annotate("segment", x = falls$x1, y = falls$y1, xend = falls$x2, yend = falls$y2,
        arrow = arrow(type = "closed", length = unit(0.02, "npc")), color = "#a6611a")+
  annotate("text", x = falls$x1[1], y = 0, label = "Fall HFE (0.45 Bankflow)" ,hjust = 0, size = 5, color = "#a6611a")+
  annotate("segment", x = springs$x1, y = springs$y1, xend = springs$x2, yend = springs$y2,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")), color = "#018571")+
  annotate("text", x = springs$x1[1], y = 0, label = "Spring HFE (0.45 Bankflow)" ,hjust = 0, size = 5, color = "#018571")

