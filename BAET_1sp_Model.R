##########################
# Baetis 1 sp model
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

source("BAETSurvivorship.R")
source("1spFunctions.R")

#------------------------------------------------------------
# Set up location specific data
#-----------------------------------------------------------
#if looking at ColRiver flows read in flow data from USGS gauge at Lees Ferry, AZ between 1985 to the end of the last water year
#flow <- readNWISdv("09380000", "00060", "1985-10-01", "2021-09-30")
#flow.magnitude <- TimestepDischarge(flow, 85000) #discharge data from USGS is in cfs, bankfull dischage (pers comm TK) 85000 cfs

#if looking at ColRiver temps read in temperature data from USGS gauge at Lees Ferry, AZ between _ to the end of last water year
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2021-09-30")
#temps <- TimestepTemperature(temp, "Colorado River") # can get time series data if desired
#degreedays <- TimestepDegreeDay(temp, "Colorado River")

#-------------------------------------------------------------
# create a series of average 2 weekly data
#-------------------------------------------------------------
# calculate mean temperature data for each timestep (2 week interval)
temps <- average.yearly.temp(temp, "X_00010_00003", "Date")

# if you want to augment ColRiver temps, set qr to different in C
qr <- 0
# if you want to augement how many years at temps, change r
r <- 13 # I want 13 years
temps <- rep.avg.year(temps, 13, change.in.temp = qr, years.at.temp = r)

# calculate accumulated degreedays, which is the days above critical threhold * temperature (degC)
# we assume accumulated degree days = 2-week average * 14
degreedays <- as.data.frame(cbind(temps$dts, temps$Temperature * 14))
colnames(degreedays) <- c("dts", "DegreeDay") # format
degreedays$dts <- as.Date(degreedays$DegreeDay, origin = "1970-01-01") # make date-time

# if desired, can pre-count how many degreedays each cohort should experience
#forward.count.degreedays(559)

#----------------------------------------------
# set model parameters
#----------------------------------------------
# specify iterations
iterations <- 50

# baseline K in the absence of disturbance
Kb <- 10000
# max K after a big disturbance
Kd <- 40000

# specify baseline transition probabilities for each species
G1_BAET = 0.25 #move to Stage2 (subimago)
G2_BAET = 0.25 #move to Stage3 (adult)
P1_BAET = 0.3 #stay in Stage1 (larvae)
P2_BAET = 0.3 #stay in Stage2 (subimago)

# want this to run as long as our temperature timeseries
timestep <- seq(2, (length(temps$Temperature) + 1), by = 1) 

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


# Q is equal to average discharge over 14 days/bankful discharge for the system
# in the Colorado River, bankful discharge = 85000 cfs (personal communication with Theodore Kennedy)
#Q <- flow.magnitude$Discharge 

# We want to create a model with no flow mortality, so set Q less than Qmin
Q <- rep(0.1, length(temps$Temperature))

Qmin <- 0.25 # Q min is the minimum flow required to impact mortality and carryin capactity (K)
a <- 0.1 #shape param for flow magnitude and K
g <- 0.1 #shape param for relationship between K and time since disturbance
h <- surv.fit.BAET$m$getPars()[2] # shape param for flood mortality
k <- surv.fit.BAET$m$getPars()[1] # shape param for flood mortality  
extinction <- 50 # extinction threshold - if Total abundance below this, population goes extinct

#-------------------------
# Outer Loop of Iterations
#--------------------------


for (iter in c(1:iterations)) {
  
  K = 10000 # need to reset K for each iteration
  
  # pull random values from a uniform distribution for starting pop
  output.N.list[1,1:3, iter]<- runif(3, min = 1, max = (0.3*K))
  
  # we often want to look at different parameter values after we run code, so we create some lists
  # list to input Ks
  Klist <- vector()
  Klist[1] <- 10000 #set first K
  
  # list to input flow mortality
  flowmortlist <- vector()
  
  # list to input fecundities
  Flist <- vector()
  
  # list to input back-looking emergence times
  emergetime <- vector()
  
  # list to input size
  sizelist <- vector()
  
  # list to input probs to remaining in same stage
  Glist <- vector()
  
  # list to input eigenvalue
  eigenlist <- vector()
  #-------------------------
  # Inner Loop of Timesteps
  #-------------------------
  
  for (t in timestep) {
    
    #----------------------------------------------------------
    # Calculate how many timesteps emerging adults have matured
    
    emergetime <- append(emergetime, back.count.degreedays(t, 559)) # value from Perry and Kennedy, 2016 
    #---------------------------------------------------------
    # Calculate fecundity per adult
    
    # we start by pulling fecundities from normal distribution
    F_BAET = rnorm(1, mean = 1104.5, sd = 42.75) * 0.5 * 0.5  #Baetidae egg minima and maxima from Degrange, 1960, assuming 1:1 sex ratio and 50% egg mortality
    
    # we can also relate fecundities to body mass.
    # Sweeney and Vannote 1980 have recorded dry body weight between 0.9 and 2.0 mg. 
    # That weight is related to fecundity Y = 614X - 300
    # we can "convert" emergetime to mg by multiplying to get dry weights between 0.9 - 2 mg, and then convert to fecunity
    # Issue: this data is for Ephemerella spp, not Baetidae spp
    
    if (t > 15) {
      size <- (emergetime[t-1] * 0.55)-0.75
      sizelist <- append(sizelist, size)
      F_BAET <- ((614 * size) - 300)* 0.5 * 0.5 
    }
    #--------------------------------------------------
    # Calculate the disturbance magnitude-K relationship
    # Sets to 0 if below the Qmin
    Qf <- Qf.Function(Q[t-1], Qmin, a)
    
    #-------------------------------------------------------------------
    # Calculate K carrying capacity immediately following the disturbance
    K0 <- K + ((Kd-K)*Qf)
    
    # Calculate final K for timestep, including relationship between K and time since disturbance
    K <- post.dist.K(K0, Kb, g)
    
    Klist <- append(Klist, K)
    #---------------------------------------------
    # Calculate effect of density dependence on fecundity 
    
    # Logistic via Rogosch et al. Fish Model
    # no immediate egg mortality incorporated
    F_BAET <- Logistic.Dens.Dependence(F_BAET, K, Total.N[t-1, iter])
    
    # add F_BAET to list
    Flist <- append(Flist, F_BAET)
    #-----------------------------------------------
    # Calculate new transition probabilities based on temperature
    # This is the growth v development tradeoff
    
    # development measures
    # in this function, we assume that if below the min temp threshold (9) no maturation occurs (slow maturation, large growth)
    # if above the max temp threshold (15), no one remains more than 1 timestep in each stage (fast maturation, small growth)
    
    # Probabilities of remaining in stages (when temps low, high prob of remaining)
    P1_BAET <- growth.development.tradeoff(temps$Temperature[t-1],  9, 13, 0.43, 0.0)
    P2_BAET <- growth.development.tradeoff(temps$Temperature[t-1], 9, 13, 0.43, 0)
    
    # growth (if below 10C, no growth can occur - everything basically freezes, if between 10 and 11, prob of remaining in same stage = 0.6395, if above 13, prob of transition to next stage is 0 )
    G1_BAET <- 0.43 - P1_BAET
    G2_BAET <- 0.43 - P2_BAET
    #-----------------------------------------------
    # Create Lefkovitch Matrix
    
    BAET1 <- c(P1_BAET, 0, F_BAET)
    BAET2 <- c(G1_BAET, P2_BAET, 0)
    BAET3 <- c(0, G2_BAET, 0) 
    
    ABAET <- rbind( BAET1, BAET2, BAET3)
    eigenlist <- append(eigenlist, eigen(ABAET)$values[1])

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
#------------------
# Analyzing Results
#-------------------
# summarizing iterations

## turning replist into a df
means.list.BAET <- mean.data.frame(output.N.list, burnin = 0)


# note how boom and bust this model is - K is set to be 10,000 not 100,000
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

