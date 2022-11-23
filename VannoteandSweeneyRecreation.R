##########################
# Baetis 1 sp model
###########################


library(purrr)
library(tidyverse)
library(lubridate)
library(plyr)
library(dplyr)
library(ggplot2)
library(scales)
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
# library(scales, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
# library(dataRetrieval, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")

## checkpos makes sure that the K-occupied term is positive, assigns 0 if not
checkpos <- function(x) {
  ifelse(x < 0, 0, x)
}

#read in flow data from USGS gauge at White Clay Creek tributary West Branch Brandywine Creek (39`57'42", 75`48'06") from Vannote and Sweeney, 1980
flow <- readNWISdv("01480617", "00060", "2018-10-01", "2022-09-05")

# so flow doesn't really matter because we are only looking at temp, so set all flows to 1000

out <- rep(1000, 27)
out <- rep(1000, 53)

# read in temperature data from USGS gauge at White Clay Creek at Neward, DE
temp <- readNWISdv("01480617", "00010", "2007-10-01", "2022-09-14")
temp$Date <- as_datetime(temp$Date)
temp$Date <- format(temp$Date, format = "%m-%d")


# option 1: for each date, calculate mean and standard error so we can create a new dataset
temp <- temp %>% group_by(Date) %>% dplyr::summarise(Temperature = mean(X_00010_00003), 
                                                     sd = sd(X_00010_00003))                                            
                                                    

# now starting at day 1 (jan 1) we can aggregate temps for every 2 weeks
# Make an index to be used for aggregating
ID <- as.numeric(as.factor(temp$Date))-1
# want it to be every 14 days, hence the 14
ID <- ID %/% 14
# aggregate over ID and TYPE for all numeric data.
outs <- aggregate(temp[sapply(temp,is.numeric)],
                  by=list(ID),
                  FUN=mean)


# format output
names(outs)[1:2] <-c("dts","Temperature")
# add the correct dates as the beginning of every period
outs$dts <-(temp$Date[(outs$dts*14)+1])

#warmWCC
outs <- rbind(outs[1:26,], outs)
#coolWCC
#outs$Temperature <- outs$Temperature - 4

degreeday <- aggregate(temp[sapply(temp,is.numeric)],
                       by=list(ID),
                       FUN=sum)
names(degreeday)[1:2] <-c("dts","DegreeDay")
# add the correct dates as the beginning of every period
# 
degreeday$dts <-(temp$Date[(degreeday$dts*14)+1])



degreeday$DegreeDay <- degreeday$DegreeDay 
# can't have negative numbers so turn those into 0s
degreeday$DegreeDay[which(degreeday$DegreeDay < 0)] <-  0

degreeday <- rbind(degreeday[1:26,], degreeday)
# specify iterations
iterations <- 5
#species <- c("BAET", "SIMU", "CHIRO")


# set carrying capacity
K = 10000
# baseline K in the absence of disturbance
Kb <- 10000
# max K after a big disturbance
Kd <- 40000

# specify baseline transition probabilities for each species
G1_EPSU = 0.6
G2_EPSU = 0.6
P1_EPSU = 0.2
P2_EPSU = 0.2

# want to run this for one year, in 14 day timesteps 
timestep <- seq(2, (length(outs$dts) + 1), by = 1) # OR
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


Qmin <- 40000
a <- 3000
g <- 0.01
h <- 0.000022 
b <- 1.221403 


#-------------------------
# Outer Loop of Iterations
#--------------------------


for (iter in c(1:iterations)) {
  
  Q <- out
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
    
    # for each timestep, we want to back calculate the number of degree days
    if (t == 1) {print("t = 1")}
    else {
      # create a sequence of time from last t to 1
      degseq <- seq(t-1, 1, by = -1)
      # create an empty vector to put the total number of degree days accumulated
      vec <- 0
      # for each value in that sequence, we will add the degree day values of 
      #the timestep prior and check if it adds up to our threshold to emergence
      for (s in degseq) {
        if(vec <= 525) { vec <- DDs$DegreeDay[s] + vec } # mean degree days from Higley 
        else {emerg <- t - s
        emergetime <- append(emergetime, emerg)
        break}
        # once we hit that threshold, we count the number of timesteps it took to reach that and add that to our emergetime vector  
      }
    }
    #---------------------------------------------------------
    # Calculate fecundity per adult
    
    # we start by pulling fecundities from normal distribution
    # relate fecundities to temperature based on Sweeney and Vannote 1981  *0.5 assuming 50% female and * 0.5 assuming 50% mort
    F_EPSU = 656.25 * 0.5 * 0.4 #mean female egg mass between 532.3 and 780.2

    # we can also relate fecundities to body mass. Sweeney and Vannote 1981 have recorded dry body weight between 6.6 and 10.8 mg. 
    # that weight is related to unit mass is related to emergetime via biomass = 0.467(emergetime) + 5.199
    # then we can convert that unit mass to fecundity using female fecunity = 143 + 59(biomass)
    if (t > 8) {
      size <- (emergetime[t-8]*0.467) + 5.199
      sizelist <- append(sizelist, size)
      F_EPSU <- ((59 * size) + 143) * 0.5 * 0.4
    }
    
    
    
    #---------------------------------------------------
    # Calculate the disturbance magnitude-K relationship 
    # Sets to 0 if below the Qmin
    Klist[1] <- 10000
    # Calculate the disturbance magnitude-K relationship 
    # Sets to 0 if below the Qmin
    if (Q[t-1] < Qmin) {
      Qf <- 0
    } else {
      Qf <- (Q[t-1] - Qmin)/(a + Q[t-1]- Qmin)
    }
    
    #-------------------------------------------------------------------
    # Calculate K arrying capacity immediately following the disturbance
    #-------------------------------------------------------------------
    # Calculate K arrying capacity immediately following the disturbance
    K <- K + ((Kd-K)*Qf)
    
    # Function to calc. K as a function of time post-disturbance at a particular disturbance intensity
    
    tau = (t-1) - (last(which(Q[1:t-1] > Qmin)))
    
    if (is.na(tau)==T) { tau <-  0}
    
    if (tau > 0) {
      
      K <- Kb + ((Klist[t-1] - Kb)*exp(-g*tau))}
    Klist <- append(Klist, K)
    
    #---------------------------------------------
    # Calculate effect of density dependnce on fecundity 
    
    # ricker model Nt+1= Nt* e^r(1-Nt/K)
    # overcompensatory, but could represent pops better? used in fisheries a lot
    # r is instrisice rate of increase when N is small 
    # r can be species specific, or the same for all species
    
    # Ricker model - pro = doesn't go negative
    #F_BAET <- F_BAET*exp(r*(1-(Total.N[t-1]/K)))
    
    #Ricker model from Mathmatica (after Recruitment = axe^-bx, see Bolker Ch 3 Deterministic Functions for
    #Ecological Modeling)
    
    #F_BAET <- F_BAET*exp(-b * Total.N[t-1, iter])
    
    # Ricker model reproductive output, from  (recruitment = axe^r-bx)
    #F_BAET <- F_BAET*exp(r-b * Total.N[t-1, iter])
    
    # beverton holt is Nt+1 = rNt/1-Nt(r-1)/K
    # it is supposed to be depensatory, so as t -> inf, Nt+1 -> K, BUT 
    # the discrete nature of this causes it overshoot by a lot, 
    # meaning it isn't any better or worse than traditional logistric growth models
    
    # Beverton Holt from Mathmatica - can't go negative
    #if (Total.N[t-1] < K){
    #  F_BAET <- F_BAET*((K - Total.N[t-1])/K)
    #} else{
    #  F_BAET <- F_BAET*((K - (K-1))/K)
    #}
    # Beverton Holt - issue, can go negative
    #if (Total.N[t-1] < K){
    #F_BAET <- F_BAET*((r*Total.N[t-1])/(1 - (Total.N[t-1]*(r-1)/K)))
    #} else{#
    #F_BAET <- F_BAET*((K - (K-1))/K)
    #}
    
    # Logistic via Rogosch et al. Fish Model
    F_EPSU <- F_EPSU * checkpos((K - Total.N[t-1, iter])/K)
    Flist <- append(Flist, F_EPSU)
    #-----------------------------------------------
    # Create Lefkovitch Matrix
    
    EPSU1 <- c(P1_EPSU, 0, F_EPSU)
    EPSU2 <- c(G1_EPSU, P2_EPSU, 0)
    EPSU3 <- c(0, G2_EPSU, 0) 
    
    AEPSU <- rbind( EPSU1, EPSU2, EPSU3)
    
    #-----------------------------------------------
    # Calculate new transition probabilities based on temperature
    # This is the growth v development tradeoff
    
    # development measures (basically, if below 10C, no development, if between 10 and 12, follows a function, if above 12, prob of transition to next stage is 0.6395)
    if (10 > outs$Temperature[t-1]) AEPSU[3,2] <- 0.001
    if (outs$Temperature[t-1] > 13) AEPSU[3,2] <- 0.55
    if (10 <= outs$Temperature[t-1] & outs$Temperature[t-1] <= 13) AEPSU[3,2] <- (0.183 * outs$Temperature[t-1]) -1.829 #(0.2 * temps$Temperature[t-1]) -2
    AEPSU[2,1] <- AEPSU[3,2] 
    
    # growth (if below 10C, no growth can occur - everything basically freezes, if between 10 and 11, prob of remaining in same stage = 0.6395, if above 13, prob of transition to next stage is 0 )
    if (10 > outs$Temperature[t-1]) AEPSU[2,2] <- 0.55
    if (outs$Temperature[t-1] > 13) AEPSU[2,2] <- 0.001
    if (10 <= outs$Temperature[t-1] & outs$Temperature[t-1] <=  13) AEPSU[2,2] <- (-0.183 * outs$Temperature[t-1]) + 2.38 #(-0.1 temps$Temperature[t-1]) - 2.6
    
    AEPSU[1,1] <- AEPSU[2,2] 
   
    #--------------------------------------
    # Calculate abundances for each stage
    
    output.N.list[t, 1:3, iter] <- AEPSU %*% output.N.list[t-1, 1:3, iter] 
    
    #------------------------------------------
    #Calculate immediate mortality due to flows
    # mortality due to flooding follows N0 = Nz*e^-hQ
    # but what about using a sigmoidal logistic function so we can control the threshold point and rate of growth
    # following m = 1/1+e^-h*(x-xf)
    # where h is is shape value
    # x is Q, and xf is threshold point (100% of pop dies)
    #plot(x = Q, y = 1/(1+exp(-0.02*(Q-100000))))
    
    #s1
    output.N.list[t, 1, iter] <- output.N.list[t, 1, iter]*b*exp(-h*Q[t-1])
    #s2
    output.N.list[t,2,iter] <- output.N.list[t,2,iter] *b*exp(-h*Q[t-1])
    #3
    output.N.list[t,3,iter] <- output.N.list[t,3,iter] *b*exp(-h*Q[t-1])
    
    flowmortlist <- append(flowmortlist, b*exp(-h*Q[t-1]))
    #replist[[1]][,,1] <- output.N.list[[1]]
    Total.N[,iter] <- apply(output.N.list[,,iter],1,sum)
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

df <- as.data.frame(cbind(sizelist, outs[8:53, 1], rep("Model", times = 46 ))) # max temp in summer is week of 7-29
colnames(df) <- c("Mean", "Date", "Group")
df$Mean <- as.numeric(df$Mean)
Vannote_Sweeney_Fig7 <- read.csv("~/ColoradoRiverInverts/Vannote_Sweeney_Fig7.csv")
Vannote_Sweeney_Fig7$Week <- as.Date(paste(2022, Vannote_Sweeney_Fig7$Week, 1, sep="-"), "%Y-%U-%u")
Vannote_Sweeney_Fig7$Week <- format(Vannote_Sweeney_Fig7$Week, format = "%m-%d")
# Vannote and Sweeney look at larvae and adults - we just look at adults
Vannote_Sweeney_adults <- Vannote_Sweeney_Fig7[which(Vannote_Sweeney_Fig7$X == "Adult"), 1:4]
Vannote_Sweeney_adults <- cbind(Vannote_Sweeney_adults, rep("Vannote & Sweeney, 1980", times = 4))
colnames(Vannote_Sweeney_adults) <- c("Date", "Max", "Mean", "Min", "Group")

Tempsize <- bind_rows(df[26:29, ], Vannote_Sweeney_adults)


abund.trends <- ggplot(data = Tempsize, aes(x = Date,
                                              y = Mean, group = Group, color = Group, fill = Group)) +
  geom_ribbon(data = subset(Tempsize, !is.na(Max)), aes(ymin =Min,
                  ymax = Max),
              colour = "NA",
              alpha = .25,
              show.legend = FALSE) +
  geom_line()+
  coord_cartesian(ylim = c(0,12)) +
  theme_bw()+
  ylab('E. subvaria size (mg)') +
  xlab('Date')

