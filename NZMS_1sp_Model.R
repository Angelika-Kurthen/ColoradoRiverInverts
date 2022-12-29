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


## checkpos makes sure that the K-occupied term is positive, assigns 0 if not
checkpos <- function(x) {
  ifelse(x < 0, 0, x)
}

## function to calculate daily shell size, to help determine stage transitions and fecundity at different stages
shell.growth <- function(m, b, start.size){
  l <- seq(1, 1000, by = 1)
  lengths <- vector(length = 1000)
  lengths[1] <- start.size
  for (i in l){
    lengths[i + 1] <- m*lengths[i]+b
  }
  return(lengths)
}

# function to index flow data, summarize as mean discharge per timestep, and relativize to flow magnitude (aka disturbance magnitude)
TimestepDischarge <- function(flow, bankfull_discharge){
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
  out$Discharge <- out$Discharge/bankfull_discharge # standardize to disturbance magnitude by taking discharge/bankfull_discharge
  return(out)
}

# function to index and summarize temperature data over timesteps length
TimestepTemperature <- function(temp, river){
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
  outs$dts <- as.POSIXct(temp$Date[(outs$dts*14)+1])
  # order by date in chronological order
  temps <- outs[order(outs$dts),]
  
  if (river == "Colorado River"){
    # there are less temperatures than discharge readings, so we will just repeat the same temp sequence for this exercise
    temps <- rbind(temps, temps, temps)
    temps <- temps[1:length(flow.magnitude$Discharge), ]
  }
  return(temps)
}

#function to calculate degree days accumulated every timestep
TimestepDegreeDay <- function(temp, river){
  # Make an index to be used for aggregating
  ID <- as.numeric(as.factor(temp$Date))-1
  # want it to be every 14 days, hence the 14
  ID <- ID %/% 14
  # aggregate over ID and T
  degreeday <- aggregate(temp[sapply(temp,is.numeric)],
                         by=list(ID),
                         FUN=sum)
  names(degreeday)[1:2] <-c("dts","DegreeDay")
  # add the correct dates as the beginning of every period
  degreeday$dts <- as.POSIXct(temp$Date[(degreeday$dts*14)+1])
  # order by date in chronological order
  degreeday <- degreeday[order(degreeday$dts),]
  # can't have negative numbers so turn those into 0s
  degreeday$DegreeDay[which(degreeday$DegreeDay < 0)] <-  0
  
  if (river == "Colorado River"){
    # there are less temperatures than discharge readings, so we will just repeat the same temp sequence for this exercise
    degreeday <- degreeday[1:363,] # last row isn't a full timestep
    DDs <- rbind(degreeday, degreeday, degreeday) 
    DDs <- DDs[1:length(flow.magnitude$Discharge), ]
  } 
  return(DDs)
}

# function to calculate Qf from McMullen et al 2017. Sets to 0 if below the Qmin
Qf.Function <- function(Q, Qmin, a){
  if (Q < Qmin) {
    Qf <- 0
  } else {
    Qf <- (Q - Qmin)/(a + Q- Qmin)
  }
  return(Qf)
}


# Function to calc. K as a function of time post-disturbance at a particular disturbance intensity
post.dist.K <- function(K0, Kb, g){
  #calculate tau (times since last distubance)
  tau = (t-1) - (last(which(Q[1:t-1] > Qmin)))
  if (is.na(tau)==T) { tau <-  0} 
  if (tau > 0) {
    K <- Kb + ((K0 - Kb)*exp(-g*tau))} # function from McMullen et al 2017, g is shape function
  Klist <- append(Klist, K)
  return(K)
}


# Function to calculate logistic density dependence on fecundity, after Rogosch et al 2019
Logistic.Dens.Dependence <- function(Fecundity, K, N){
  f.rate <- Fecundity * checkpos((K - N)/K) * 0.5
  return(f.rate)
}


#Ricker model (after Recruitment = axe^-bx, see Bolker Ch 3 Deterministic Functions for
#Ecological Modeling)
Ricker.Dens.Dependence <- function(b, N, fecundity){
  f.rate <- fecundity * exp(-b * N)
  return(f.rate)}
# b = 0.005
#F_NZMS <- Ricker.Dens.Dependence(b, Total.N[t-1, iter], F_NZMS) 

# beverton holt is Nt+1 = rNt/1-Nt(r-1)/K
# it is supposed to be depensatory, so as t -> inf, Nt+1 -> K, BUT 
# the discrete nature of this causes it overshoot by a lot, 
# meaning it isn't any better or worse than traditional logistric growth models
Bev.Holt.Dens.Dependence <- function(r, N, K, fecundity){
  if(N < K){
    f.rate <- fecundity * (K - N/K)
  } else {
    f.rate <- fecundity * (1/K)
  }
}

# F_NZMS <- Bev.Holt.Dens.Dependence(r, Total.N[t-1, iter], K, F_NZMS)




# mortality due to flooding follows N0 = Nz*e^-h
flood.mortality <- function(N, k, h, Q, Qmin){
  if (Q <= Qmin){
    newN <- N
  } else {
    newN <- N * k * exp(-h * Q)
  }
  return(newN)
}


# source functions
source("NZMSSurvivorship.R")


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

flow.magnitude <- TimestepDischarge(flow, 85000) #discharge data from USGS is in cfs, bankfull dischage (pers comm TK) 85000 cfs
temps <- TimestepTemperature(temp, "Colorado River") # calculate mean temperature data for each timestep
degreedays <- TimestepDegreeDay(temp, "Colorado River")

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
timestep <- seq(2, (length(flow.magnitude$Discharge) + 1), by = 1) # OR

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
Q <- flow.magnitude$Discharge

Qmin <- 0.25
a <- 
g <- 0.1
h <- surv.fit.NZMS$m$getPars()[2]   
k <- surv.fit.NZMS$m$getPars()[1]
e = 2.71828


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
    
    # # for each timestep, we want to back calculate the number of degree days
    # if (t == 1) {print("t = 1")}
    # else {
    #   # create a sequence of time from last t to 1
    #   degseq <- seq(t-1, 1, by = -1)
    #   # create an empty vector to put the total number of degree days accumulated
    #   vec <- 0
    #   # for each value in that sequence, we will add the degree day values of 
    #   #the timestep prior and check if it adds up to our threshold to emergence
    #   for (s in degseq) {
    #     if(vec <= 266) { vec <- DDs$DegreeDay[s] + vec }
    #     else {emerg <- t - s
    #     emergetime <- append(emergetime, emerg)
    #     break}
    #     # once we hit that threshold, we count the number of timesteps it took to reach that and add that to our emergetime vector  
    #   }
    # }
    #---------------------------------------------------------
    # Calculate starting fecundity per adult
    
    # fecundities estimated from McKenzie et al. 2013; 
    if (temps$Temperature[t-1] <= 10) { 
      F2_NZMS <- 0
      F3_NZMS <- 0  
    } else {
      F2_NZMS <- 8.87473
      F3_NZMS <- 27.89665
    }
    
    
    #---------------------------------------------------
    # Calculate the disturbance magnitude-K relationship 
    # Sets to 0 if below the Qmin
    Klist[1] <- 10000
    # Calculate the disturbance magnitude-K relationship 

    Qf <- Qf.Function(Q[t-1], Qmin, a)
    #-------------------------------------------------------------------
    # Calculate K arrying capacity immediately following the disturbance
    K0 <- K + ((Kd-K)*Qf)
    
    # Calculate final K for timestep, including relationship between K and time since disturbance
    K <- post.dist.K(K0, Kb, g)
    
    #---------------------------------------------
    # Calculate effect of density dependnce on fecundity 
    
    # Logistic Density Dependence on Fecundity via Rogosch et al. Fish Model
    F2_NZMS <- Logistic.Dens.Dependence(F2_NZMS, K, Total.N[t-1, iter])
    F3_NZMS <- Logistic.Dens.Dependence(F3_NZMS, K, Total.N[t-1, iter])
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
    flowmortlist <- append(flowmortlist, k* exp(-h*Q[t-1]))

    #-------------------------------------------------
    # calculate sum of all stages (total population)
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
repdf <- plyr::adply(output.N.list, c(1,2,3))

names(repdf) <- c('timesteps', 'stage', 'rep', 'abund')
repdf$timesteps <- as.numeric(as.character(repdf$timesteps))

totn <- adply(Total.N, c(1,2))
names(totn) <- c('timesteps', 'rep', 'tot.abund')
totn$timesteps <- as.numeric(as.character(totn$timesteps))

## joining totn and repdf together
repdf <- left_join(totn, repdf)

## calculating relative abundance
repdf <- mutate(repdf, rel.abund = abund/tot.abund)
repdf$timesteps <- as.factor(repdf$timesteps)
## Taking mean results to cf w/ observed data
means.list<- repdf %>%
  select(-tot.abund) %>%
  dplyr::group_by(timesteps, rep) %>% # combining stages
  dplyr::summarise(abund = sum(abund),
                   rel.abund = sum(rel.abund)) %>%
  ungroup() %>%
  dplyr::group_by(timesteps) %>%
  dplyr::summarise(mean.abund = mean(abund),
                   sd.abund = sd(abund),
                   se.abund = sd(abund)/sqrt(iterations),
                   mean.rel.abund = mean(rel.abund),
                   sd.rel.abund = sd(rel.abund),
                   se.rel.abund = sd(rel.abund)/sqrt(iterations)) %>%
  ungroup()
## Save the objects as .rds files - then use loadRDS in other file. 
saveRDS(means.list, paste0('modelresults', '.rds'))
#flowmeans.list[[flowmod]] <- ldply(flowlist, function(x) apply(x, 2, mean)) %>%
#  apply(2, mean)

# want to exclude first 10 timesteps for "burn in"
burns.list <- means.list[10:941, ]

abund.trends <- ggplot(data = burns.list, aes(x = timesteps,
                                              y = mean.abund, group = 1)) +
  geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
                  ymax = mean.abund + 1.96 * se.abund),
              colour = 'transparent',
              alpha = .5,
              show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  coord_cartesian(ylim = c(0,7000)) +
  ylab('NZMS Abundance') +
  xlab('Timestep')

saveRDS(abund.trends, paste0('BAETplot', '.rds'))



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


# 
lines(timestep, Klist, type = "l", col = "blue")


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
