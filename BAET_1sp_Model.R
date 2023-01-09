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

#read in flow data from USGS gauge at Lees Ferry, AZ between 1985 to the end of the last water year
flow <- readNWISdv("09380000", "00060", "1985-10-01", "2021-09-30")
# read in temperature data from USGS gauge at Lees Ferry, AZ between _ to the end of last water year
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2021-09-30")


flow.magnitude <- TimestepDischarge(flow, 85000) #discharge data from USGS is in cfs, bankfull dischage (pers comm TK) 85000 cfs
temps <- TimestepTemperature(temp, "Colorado River") # calculate mean temperature data for each timestep
degreedays <- TimestepDegreeDay(temp, "Colorado River")

# specify iterations
iterations <- 5
#species <- c("BAET", "SIMU", "CHIRO")

# baseline K in the absence of disturbance
Kb <- 10000
# max K after a big disturbance
Kd <- 40000

# specify baseline transition probabilities for each species
G1_BAET = 0.6
G2_BAET = 0.6
P1_BAET = 0.2
P2_BAET = 0.2

# want to run this for one year, in 14 day timesteps 
timestep <- seq(2, (length(out$Discharge) + 1), by = 1) # OR
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

Qmin <- 0.25
a <- 3000
g <- 0.1
h <- surv.fit.BAET$m$getPars()[2]   
k <- surv.fit.BAET$m$getPars()[1]   


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
    
    back.count.degreedays <- function(time, criticaldegreedays){
      # for each timestep, we want to back calculate the number of degree days
      if(time == 1) {print("t=1")
        emerg <- NA
      } else {
        # create a sequence of time from last t to 1
          degseq <- seq(time - 1, 1, by = -1)
          # create an empty vector to put the total number of degree days accumulated
          vec <- 0
          # for each value in that sequence, we will add the degree day values of 
          #the timestep prior and check if it adds up to our threshold to emergence
          for (s in degseq) {
            if(vec <= criticaldegreedays) { vec <- degreedays$DegreeDay[s] + vec}
            else {emerg <- time - s
            break}
            # once we hit that threshold, we count the number of timesteps it took to reach that and add that to our emergetime vector
          }
      }
      return(emerg)
    }
    
    emergetime <- append(emergetime, back.count.degreedays(t, 406)) # or 406? 266 degree days used from REF days above 10C. Add 266 + 10*14 
    #---------------------------------------------------------
    # Calculate fecundity per adult
    
    # we start by pulling fecundities from normal distribution
    F_BAET = rnorm(1, mean = 1104.5, sd = 42.75) * 0.5 *0.5  #* H_BAET #Baetidae egg minima and maxima from Degrange, 1960 *0.5 assuming 50% female and * 0.5 assuming 50% mort.
    
    # relate fecundities to temperature based on Sweeney et al., 2017  *0.5 assuming 50% female and * 0.5 assuming 60% mort.
    #F_BAET <- (-379.8021 * (temps$Temperature[t-1]) + 16.4664*(temps$Temperature[t-1]^2) - 0.2684* (temps$Temperature[t-1]^3) + 4196.8608) * 0.5 * 0.5
    
    # we can also relate fecundities to body mass. Sweeney and Vannote 1980 have recorded dry body weight between 0.9 and 2.0 mg. 
    # that weight is related to fecundity Y = 614X - 300
    # we can "convert" emergetime to mg by multiplying by 0.225 (to get dry weights between 0.9 - 2 mg)
    
    if (t > 4) {
      size <- emergetime[t-1] * 0.225
      sizelist <- append(sizelist, size)
      F_BAET <- ((614 * size) - 300) * 0.5 * 0.5
    }
    #--------------------------------------------------
    # Calculate the disturbance magnitude-K relationship 
    Klist[1] <- 10000
    
    # Calculate the disturbance magnitude-K relationship
    # Sets to 0 if below the Qmin
    Qf <- Qf.Function(Q[t-1], Qmin, a)
    
    #-------------------------------------------------------------------
    # Calculate K arrying capacity immediately following the disturbance
    K0 <- K + ((Kd-K)*Qf)
    
    # Calculate final K for timestep, including relationship between K and time since disturbance
    K <- post.dist.K(K0, Kb, g)
    
    
    #---------------------------------------------
    # Calculate effect of density dependnce on fecundity 
    
    # Logistic via Rogosch et al. Fish Model
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
    
    growth.development.tradeoff <- function(temp, thresholdtemp.min, thresholdtemp.max, min.rate, max.rate, m, b ){  # m and b from y = mx+b 
      if (thresholdtemp.min > temp) rate <- min.rate
      if (temp > thresholdtemp.max) rate <- max.rate
      if (thresholdtemp.min <= temp & temp <= thresholdtemp.max) rate <- (m * temp) + b
      return(rate)
    }
    # development measures (basically, if below 10C, no development, if between 10 and 12, follows a function, if above 12, prob of transition to next stage is 0.6395
    ABAET[3, 2] <- growth.development.tradeoff(temps$Temperature[t-1],  10, 13, 0.001, 0.55, 0.183, -1.829)
    ABAET[2,1] <- ABAET[3,2] 

    # growth (if below 10C, no growth can occur - everything basically freezes, if between 10 and 11, prob of remaining in same stage = 0.6395, if above 13, prob of transition to next stage is 0 )
    ABAET[2,2] <- growth.development.tradeoff(temps$Temperature[t-1], 10, 13, 0.55, 0.001, -0.183, 2.38)
    ABAET[1,1] <- ABAET[2,2]

    #--------------------------------------
    # Calculate abundances for each stage
    
    output.N.list[t, 1:3, iter] <- ABAET %*% output.N.list[t-1, 1:3, iter] 
    
    #------------------------------------------
    #Calculate immediate mortality due to flows
    # mortality due to flooding follows N0 = Nz*e^-hQ
    
    #s1
    output.N.list[t, 1, iter] <- flood.mortality(output.N.list[t, 1, iter], k, h, Q[t-1], Qmin)
    #s2
    output.N.list[t,2,iter] <- flood.mortality(output.N.list[t,2,iter], k, h, Q[t-1], Qmin)
    
    output.N.list[t,3,iter] <- flood.mortality(output.N.list[t,3,iter], k, h, Q[t-1], Qmin)
    
    flowmortlist <- append(flowmortlist, k* exp(-h*Q[t-1]))
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

## turning replist into a df
mean.data.frame(output.N.list, stages = c(1,2,3), burnin = 10)

abund.trends <- ggplot(data = burns.list, aes(x = timesteps,
                                       y = mean.abund, group = 1)) +
  geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
                  ymax = mean.abund + 1.96 * se.abund),
              colour = 'transparent',
              alpha = .5,
              show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  coord_cartesian(ylim = c(0,6000)) +
  ylab('Baetis Abundance') +
  xlab('Timestep')

saveRDS(abund.trends, paste0('BAETplot', '.rds'))



# take a look at results
# 
# par(mfrow = c(1,1))
plot(timestep[9:(length(timestep)+1)], output.N.list[9:(length(timestep)+1), 3, 1], type = "l", ylab = "Baetis spp. Adults", xlab = "Timestep (1 fortnight)")
plot(timestep[9:length(timestep)], Total.N[10:(length(timestep)+1)], type= "l", ylab = "Baetis spp. Total N", xlab = "Timestep (1 fortnight)")
# 
# 
# 
# #creating plots to analyze how temp relationship is working
# #create dataframe with timestemp, s1, s2, and s3 abundances, and tempterature
#  data <- as.data.frame(cbind(timestep, output.N.list[2:(length(timestep)+1) ,1, 1], output.N.list[2:(length(timestep)+1) ,2, 1], output.N.list[2:(length(timestep)+1) ,3, 1], temps$Temperature))
#  colnames(data) <- c("timestep", "Stage1", "Stage2", "Stage3", "Temperature")
# 
#  data <- data[10:60, ]
# ggplot(data = data, aes(x = timestep, y = Stage1, color = "Stage1"))+
#    geom_path()+
#    geom_path(aes(x = timestep, y = Stage2, color = "Stage2"))+
#    geom_path(aes(x = timestep, y = Stage3, color = "Stage3"))+
#    geom_path(aes(x = timestep, y = Temperature*200, color = "Temperature"))+
#    scale_y_continuous(
# 
#      # Features of the first axis
#      name = "Abundance",
# 
#      # Add a second axis and specify its features
#      sec.axis = sec_axis( ~.*0.005, name="Temperature C")
#    )
# 
#  # plot to show relationship between temp and fecundity
#  plot(temps$Temperature, Flist, ylab = "Fecundity per individual", xlab = "Temperature (C)", pch = 19)
# 
#  plot(timestep[10:60], output.N.list[10:60, 3,1] * 1, type = "l", ylim  = c(0, 500), ylab = " ")
#  lines(timestep[10:60], ((Flist[10:60]*output.N.list[10:60, 3,1])/10), col = "blue", ylim = c(0, 450))
#  lines(timestep[10:60], temps$Temperature[10:60]*15, col = "red")
# lines(timestep[10:60], emergetime[4:54]*18, col = "green")
#  lines(timestep[10:60], sizelist[4:54]*100, col = "magenta")
#  legend(25, 430, legend = c("Adult Abundance", "Realized Fecundity per female * 0.1", "Temperature (C) * 15", "Dry Weight (mg) * 100"),
#         col = c("black", "blue", "red", "magenta"), lty = 1, cex = 0.8)
# 
#  plot(timestep, temps$Temperature, type = "l", col = "blue")
# 
#  par(mfrow = c(1,2))
#  plot(timestep[10:37], Glist[10:37], type = "l", col = "blue")
#  lines(timestep[10:37], Plist[10:37], type = "l", col = "black")
#  legend(18, 0.1, legend=c("Growth (remain)", "Development (transition)"),
#         col=c("Black", "Blue"), lty=1, cex=0.8)
#  plot(timestep[10:37], temps$Temperature[10:37], type = "l", col = "red")
# 
#  plot(timestep[200:210], Total.N[201:211], type= "l", ylab = "Baetis spp. Total N", xlab = "Timestep (1 fortnight")
#  par(new=TRUE)
#  lines(timestep[200:210],temps$Temperature[201:211],col="green")
# 
#  Total.N
# 
#  r <-Total.N[2:(length(timestep)+1)]/Total.N[1:length(timestep)]
#  plot(timestep, r, type = "l")
# 
#  plot(Q, Klist[1:940])
# 
#  
#  
