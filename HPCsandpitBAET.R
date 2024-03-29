##########################
# HPC Playground Baetis 1 sp model
###########################

# 
#library(purrr)
#library(tidyverse)
#library(lubridate)#
#library(plyr)
#library(dplyr)
#library(ggplot2)
# data retrieval tool from USGS
#library(dataRetrieval)

#Code for HPC - tidyverse has some issues on our HPC because one of the packages is deprecated
#We have to manually load all tidyverse packages
library(purrr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
library(tibble, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
library(tidyr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
library(readr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
library(stringr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
library(forcats, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
library(lubridate, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
library(plyr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
library(dplyr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
library(ggplot2, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
library(dataRetrieval, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")

source("/home/ib/kurthena/ColoradoRiverInverts/BAETSurvivorship.R")
source("/home/ib/kurthena/ColoradoRiverInverts/BAETSurvivorship.R/1spFunctions.R")

#read in flow data from USGS gauge at Lees Ferry, AZ between 1985 to the end of the last water year
#flow <- readNWISdv("09380000", "00060", "1985-10-01", "2021-09-30")
# read in temperature data from USGS gauge at Lees Ferry, AZ between _ to the end of last water year
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2021-09-30")

temps <- average.yearly.temp(temp, "X_00010_00003", "Date")

#flow.magnitude <- TimestepDischarge(flow, 85000) #discharge data from USGS is in cfs, bankfull dischage (pers comm TK) 85000 cfs
#temps <- TimestepTemperature(temp, "Colorado River") # calculate mean temperature data for each timestep
#degreedays <- TimestepDegreeDay(temp, "Colorado River")
n <- 13
# qr is the temp ramps I want to increase the average Lees Ferry temp by 
qr <- c(0, 0, 0, 0, 0)
# how many years I want each temp ramp to last
r <- c(5, 2, 2, 2, 2)
temps <- rep.avg.year(temps, n, change.in.temp = qr, years.at.temp = r)

# specify iterations
iterations <- 5
#species <- c("BAET", "SIMU", "CHIRO")

# baseline K in the absence of disturbance
Kb <- 10000
# max K after a big disturbance
Kd <- 40000


initial.dev.list <- seq(0.045, 0.6, length.out = 10)
initial.grow.list <- seq(0.025, 0.6, length.out = 10)
Fecundlist <- seq(0.089, 0.5, length.out = 10)
g1list <- seq(0.045,0.6, length.out = 10)
g2list <- seq(0.025, 0.6, length.out = 10)
p1list <- seq(0.045, 0.6, length.out = 10)
p2list <- seq(0.024, 0.6, length.out = 10)
for (initial.dev in initial.dev.list){
  for (initial.grow in initial.grow.list){
    for(Fecun in Fecundlist){
      for (g1 in g1list){
        for (g2 in g2list){
          for (p1 in p1list){
            for(p2 in p2list){
              


# specify baseline transition probabilities for each species
G1_BAET = initial.dev
G2_BAET = initial.dev
P1_BAET = initial.grow
P2_BAET = initial.grow


# want to run this for one year, in 14 day timesteps 
timestep <- seq(2, (length(temps$Temperature) + 1), by = 1) # OR
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
Q <- flow.magnitude$Discharge
Q <- rep(0.1, length(temps$Temperature))


Qmin <- 0.25
a <- 0.1
g <- 0.1
h <- surv.fit.BAET$m$getPars()[2]   
k <- surv.fit.BAET$m$getPars()[1]   
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
    
    emergetime <- append(emergetime, back.count.degreedays(t, 406)) # or 406? 266 degree days used from REF days above 10C. Add 266 + 10*14 
    #---------------------------------------------------------
    # Calculate fecundity per adult
    
    # we start by pulling fecundities from normal distribution
    F_BAET = rnorm(1, mean = 1104.5, sd = 42.75)*0.5*0.25 #* H_BAET #Baetidae egg minima and maxima from Degrange, 1960 *0.5 assuming 50% female and 0.8 because 80% survival of eggs according to McMullen 2019 supplemental materials
    
    # relate fecundities to temperature based on Sweeney et al., 2017  *0.5 assuming 50% female and * 0.5 assuming 60% mort.
    #F_BAET <- (-379.8021 * (temps$Temperature[t-1]) + 16.4664*(temps$Temperature[t-1]^2) - 0.2684* (temps$Temperature[t-1]^3) + 4196.8608) * 0.5 * 0.5
    
    # we can also relate fecundities to body mass. Sweeney and Vannote 1980 have recorded dry body weight between 0.9 and 2.0 mg. 
    # that weight is related to fecundity Y = 614X - 300
    # we can "convert" emergetime to mg by multiplying by 0.225 (to get dry weights between 0.9 - 2 mg)
    
    if (t > 5) {
      size <- (emergetime[t-1] * 0.55)-0.75
      sizelist <- append(sizelist, size)
      F_BAET <- ((614 * size) - 300)* 0.5*0.25 #* hydropeaking.mortality(lower = 0.0, upper = 0.2, h = hp[t-1])
    }
    #--------------------------------------------------
    # Calculate the disturbance magnitude-K relationship 
    #Klist[1] <- 10000
    
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
    
    # Logistic via Rogosch et al. Fish Model
    # no immediate egg mortality incorporated
    F_BAET <- Logistic.Dens.Dependence(F_BAET, K, Total.N[t-1, iter])
    
    Flist <- append(Flist, F_BAET)
    #-----------------------------------------------
    # Create Lefkovitch Matrix
    
    BAET1 <- c(P1_BAET, 0, F_BAET)
    BAET2 <- c(G1_BAET, P2_BAET, 0)
    BAET3 <- c(0, G2_BAET, 0) 
    
    ABAET <- rbind( BAET1, BAET2, BAET3)
    
    #-----------------------------------------------
    # Calculate new transition probabilities based on temperature
    # This is the growth v development tradeoff
    
    # development measures (basically, if below 10C, no development, if between 10 and 12, follows a function, if above 12, prob of transition to next stage is 0.6395
    ABAET[3, 2] <- growth.development.tradeoff(temps$Temperature[t-1],  10, 13, 0.001, g2)#G2
    ABAET[2,1] <- growth.development.tradeoff(temps$Temperature[t-1], 10, 13, 0.001, g1) #G1
    
    # growth (if below 10C, no growth can occur - everything basically freezes, if between 10 and 11, prob of remaining in same stage = 0.6395, if above 13, prob of transition to next stage is 0 )
    ABAET[2,2] <- growth.development.tradeoff(temps$Temperature[t-1], 10, 13, p2, 0.001)#P2
    ABAET[1,1] <- growth.development.tradeoff(temps$Temperature[t-1], 10, 13, p1, 0.001)#P1
    
    #--------------------------------------
    # Calculate abundances for each stage
    
    output.N.list[t, 1:3, iter] <- ABAET %*% output.N.list[t-1, 1:3, iter] 
    
    #------------------------------------------
    #Calculate immediate mortality due to flows
    # mortality due to flooding follows N0 = Nz*e^-hQ
    
    #s1
    output.N.list[t, 1, iter] <- flood.mortality(output.N.list[t, 1, iter], k, h, Q[t-1], Qmin)
    #s2Q
    output.N.list[t,2,iter] <- flood.mortality(output.N.list[t,2,iter], k, h, Q[t-1], Qmin)
    
    output.N.list[t,3,iter] <- flood.mortality(output.N.list[t,3,iter], k, h, Q[t-1], Qmin)
    
    flowmortlist <- append(flowmortlist, flood.mortality(1, k, h, Q[t-1], Qmin))
    #replist[[1]][,,1] <- output.N.list[[1]]
    # check extinction threshold
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

#------------------
# Analyzing Results
#-------------------
# summarizing iterations

## turning replist into a df
means.list.BAET <- mean.data.frame(output.N.list, burnin = 27)

abund.trends.BAET <- ggplot(data = means.list.BAET, aes(x = timesteps,
                                                        y = mean.abund, group = 1)) +
  geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
                  ymax = mean.abund + 1.96 * se.abund),
              colour = 'transparent',
              alpha = .5,
              show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  coord_cartesian(ylim = c(0,100000)) +
  ylab('Baetis Abundance') +
  xlab('Timestep')

if (is.na(means.list.BAET$mean.abund[length(means.list.BAET$mean.abund)]) == F & means.list.BAET$mean.abund[length(means.list.BAET$mean.abund)]> 0 ){
saveRDS(abund.trends, paste0('BAETplot', initial.dev, initial.grow, Fecun, g1, p1, g2, p2,'.rds'))}
            }
          }
        }
      }
    }
  }
}