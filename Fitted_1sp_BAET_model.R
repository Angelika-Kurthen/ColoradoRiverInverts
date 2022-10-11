
Vannote_Sweeney_temp <- read.csv("~/ColoradoRiverInverts/Vannote_Sweeney_temp.csv", header=FALSE)
Vannote_Sweeney_temp <- Vannote_Sweeney_temp[-1, ]
Vannote_Sweeney_temp <- as.data.frame(Vannote_Sweeney_temp)
Vannote_Sweeney_temp$V1 <- as.numeric(Vannote_Sweeney_temp$V1) * 7
Vannote_Sweeney_temp$V1 <- c(1:64)

fit <- lm(Vannote_Sweeney_temp$V4 ~ Vannote_Sweeney_temp$V1 + sin(2*pi/(365)*Vannote_Sweeney_temp$V1)+cos(2*pi/(365)*Vannote_Sweeney_temp$V1),data=Vannote_Sweeney_temp)
summary(fit)
Date <- c(1:365)

Temperature <-  -7.374528  * (cos(((2*pi)/365)*Date))  +  (-1.649263* sin(2*pi/(365)*Date)) + 10.956243

temp <- as.data.frame(cbind(Date, Temperature))

plot(Vannote_Sweeney_temp$V1, Vannote_Sweeney_temp$V4)
lines(Date, Temperature)

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


# so flow doesn't really matter because we are only looking at temp, so set all flows to 1000

out <- rep(1000, 26)

# #Make an index to be used for aggregating
# ID <- as.numeric(as.factor(flow$Date))-1
# # want it to be every 14 days, hence the 14
# ID <- ID %/% 14
# # aggregate over ID and TYPE for all numeric data.
# out <- aggregate(flow[sapply(flow,is.numeric)],
#                  by=list(ID,flow$X_00060_00003),
#                  FUN=mean)
# # format output
# names(out)[1:2] <-c("dts","Discharge")
# # add the correct dates as the beginning of every period
# out$dts <- as.POSIXct(flow$Date[(out$dts*14)+1])
# # order by date in chronological order
# out <- out[order(out$dts),]
# # get mean Discharge data for every 14 days
# out <- aggregate(out, by = list(out$dts), FUN = mean)
# 
# 
# temp$Date <- as_datetime(temp$Date)
# temp$Date <- yday(temp$Date)


# option 1: for each date, calculate mean and standard error so we can create a new dataset
# temp <- temp %>% group_by(Date) %>% dplyr::summarise(Temperature = mean(X_00010_00003), 
#                                                      sd = sd(X_00010_00003), 
#                                                      count = n(), 
#                                                      se = sd(X_00010_00003)/count)
#specify iterations
iterations <- 5000
# now we want to create some data based on this 
temparray <- array(0,
                   dim  <-c(length(temp$Date), iterations),
                   dimnames <- list(1:length(temp$Date), 1:iterations))
fortarray <- array(0, 
                   dim <- c(26, iterations), 
                   dimnames <- list(1:26, 1:iterations))
datearray <- fortarray

sizearray <- array(0, 
                   dim <- c(26, iterations), 
                   dimnames <- list(1:26, 1:iterations))

#species <- c("BAET", "SIMU", "CHIRO")


# set carrying capacity
K = 10000

# specify baseline transition probabilities for each species
G1_BAET = 0.8
G2_BAET = 0.8
P1_BAET = 0.8
P2_BAET = 0.8


# want to run this for one year, in 14 day timesteps 
timestep <- seq(2, (length(out) + 1), by = 1) # OR
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
Q <- out

Qmin <- 20000
a <- 100
g <- 0.1




# in this case, we will use the r that McMullen et al 2017 used for Beatis
#r = 1.23
e = 2.71828
b = 0.005

#-------------------------
# Outer Loop of Iterations
#--------------------------


for (iter in c(1:iterations)) {
  
  r = 1.23
  
  # # create array with 14 (or any multiple of 14) columns to pull from 
  # temps <- c()
  # for (day in 1:length(temp$Date)) {
  #   temps[day] <- rnorm(1, mean = temp$Temperature[day], sd = temp$se[day])
  # }
  # 
  # meantemp <- mean(temps)
  # # change amplitude
  # temps[which(temps > meantemp)] <- ((temps[which(temps > meantemp)] - meantemp)*2) + meantemp
  # temps[which(temps < meantemp)] <-  meantemp - ((meantemp - temps[which(temps < meantemp)])*2)
  # 
  temparray[,iter] <- temp$Temperature
  # temparray[,iter] <- temps
  
  
  
  # Make an index to be used for aggregating
  ID <- as.numeric(as.factor(temp$Date)) + iter 
  # want it to be every 14 days, hence the 14
  ID <- ID %/% 14
  
  ID[which(ID >=26)] <- ID[which(ID >= 26)] - 26
  # aggregate over ID and TYPEall numeric data.
  outs <- aggregate(temparray[,iter]
                    [sapply(temp,is.numeric)],
                    by=list(ID),
                    FUN=mean)
  
  
  
  # format output
  names(outs)[1:2] <-c("dts","Temperature")
  # add the correct dates as the beginning of every period
  outs$dts <- (outs$dts*14)+iter
  outs$dts[which(outs$dts >=365)] <- outs$dts[which(outs$dts >= 365)] - 364
  outs$dts <- strptime(outs$dts, "%j") ###Note need to subtract 365 if larger than 365
  
  # order by date in chronological order#
  #outs <- outs[order(outs$dts),]
  outs$dts <- as_date(outs$dts)
  
  fortarray[,iter] <- outs$Temperature
  datearray[,iter] <- outs$dts
  
  degreeday <- aggregate(temparray[,iter][sapply(temp,is.numeric)],
                         by=list(ID),
                         FUN=sum)
  names(degreeday)[1:2] <-c("dts","DegreeDay")
  # add the correct dates as the beginning of every period
  
  degreeday$dts <- (degreeday$dts*14)+iter
  degreeday$dts[which(degreeday$dts >=365)] <- degreeday$dts[which(degreeday$dts >= 365)] - 364
  degreeday$dts <- strptime(degreeday$dts, "%j") ###Note need to subtract 365 if larger than 365
  
  # order by date in chronological order
  #degreeday <- degreeday[order(degreeday$dts),]
  #degreeday <- degreeday[1:363,]
  
  degreeday$DegreeDay <- degreeday$DegreeDay - 100
  # can't have negative numbers so turn those into 0s
  degreeday$DegreeDay[which(degreeday$DegreeDay < 0)] <-  0
  
  
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
  
  # list of transitions to next change
  #Glist <- vector()
  
  # list of probability of remaining in stage
  #Plist <- vector()
  
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
        if(vec <= 266) { vec <- degreeday$DegreeDay[s] + vec
        }
        else {emerg <- t - s
        emergetime <- append(emergetime, emerg)
        break}
        # once we hit that threshold, we count the number of timesteps it took to reach that and add that to our emergetime vector  
      }
      if (vec  <= 266) {
        emerg = NA
        emergetime <- append(emergetime, emerg)}
    }
    #---------------------------------------------------------
    # Calculate fecundity per adult
    
    # we start by pulling fecundities from normal distribution
    F_BAET = rnorm(1, mean = 1104.5, sd = 42.75) * 0.5 *0.5  #* H_BAET #Baetidae egg minima and maxima from Degrange, 1960 *0.5 assuming 50% female and * 0.5 assuming 50% mort.
    
    # relate fecundities to temperature based on Sweeney et al., 2017  *0.5 assuming 50% female and * 0.5 assuming 50% mort.
    #F_BAET <- (-379.8021 * (temps$Temperature[t-1]) + 16.4664*(temps$Temperature[t-1]^2) - 0.2684* (temps$Temperature[t-1]^3) + 4196.8608) * 0.5 * 0.5
    
    # we can also relate fecundities to body mass. Sweeney and Vannote 1980 have recorded dry body weight between 0.9 and 2.0 mg. 
    # that weight is related to fecundity Y = 614X - 300
    # we can "convert" emergetime to mg by multiplying by 0.225 (to get dry weights between 0.9 - 2 mg)
    
    size <- emergetime[t-1] * 0.225
    sizelist <- append(sizelist, size)
    if (!is.na(size[t]) == T){
      F_BAET <- ((614 * size) - 300) * 0.5 * 0.5}
    
    
    Flist <- append(Flist, F_BAET)
    
    #---------------------------------------------------
    # Calculate the disturbance magnitude-K relationship 
    # Sets to 0 if below the Qmin
    if (Q[t-1] < Qmin) {
      Qf <- 0
    } else {
      Qf <- (Q[t-1] - Qmin)/(a + Q[t-1]- Qmin)
    }
    
    #-------------------------------------------------------------------
    # Calculate K arrying capacity immediately following the disturbance
    K <- 10000 + ((40000-10000)*Qf)
    
    # Function to calc. K as a function of time post-disturbance at a particular disturbance intensity
    K <- 10000 + ((K - 10000)*exp(-g*14))
    
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
    F_BAET <- F_BAET*exp(r-b * Total.N[t-1, iter])
    
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
    
    #-----------------------------------------------
    # Create Lefkovitch Matrix
    
    BAET1 <- c(P1_BAET, 0, F_BAET)
    BAET2 <- c(G1_BAET, P2_BAET, 0)
    BAET3 <- c(0, G2_BAET, 0) 
    
    ABAET <- rbind( BAET1, BAET2, BAET3)
    
    #-----------------------------------------------
    # Calculate new transition probabilities based on temperature
    # This is the growth v development tradeoff
    
    # development measures (basically, if below 10C, no development, if between 10 and 12, follows a function, if above 12, prob of transition to next stage is 0.6395)
    if (10 > outs$Temperature[t-1]) ABAET[3,2] <- 0.0
    if (10 > outs$Temperature[t-1]) ABAET[2,1] <- 0.0
    if (outs$Temperature[t-1] > 13) ABAET[3,2] <- 0.6  
    if (10 <= outs$Temperature[t-1] & outs$Temperature[t-1] <= 13) {ABAET[3,2] <- (0.2 * outs$Temperature[t-1]) -2
    ABAET[2,1] <- ABAET[3,2] }
    
    # growth (if below 10C, no growth can occur - everything basically freezes, if between 10 and 11, prob of remaining in same stage = 0.6395, if above 13, prob of transition to next stage is 0 )
    if (10 >= outs$Temperature[t-1]) ABAET[2,2] <- 0.6
    if (outs$Temperature[t-1] > 13) ABAET[1,1] <- 0.0
    if (outs$Temperature[t-1] > 13) ABAET[2,2] <- 0.0
    if (10 < outs$Temperature[t-1] & outs$Temperature[t-1] <=  13) {ABAET[2,2] <- (-0.2 * outs$Temperature[t-1]) + 2.6
    
    ABAET[1,1] <- ABAET[2,2] }
    
    # 
    #Glist <-append(Glist, ABAET[3,2])
    #Plist <- append(Plist, ABAET[2,2])
    
    #--------------------------------------
    # Calculate abundances for each stage
    
    output.N.list[t, 1:3, iter] <- ABAET %*% output.N.list[t-1, 1:3, iter] 
    
    #------------------------------------------
    #Calculate immediate mortality due to flows
    # mortality due to flooding follows N0 = Nz*e^-hQ
    # but what about using a sigmoidal logistic function so we can control the threshold point and rate of growth
    # following m = 1/1+e^-h*(x-xf)
    # where h is is shape value
    # x is Q, and xf is threshold point (100% of pop dies)
    #plot(x = Q, y = 1/(1+exp(-0.02*(Q-100000))))
    
    #s1
    output.N.list[t, 1, iter] <- output.N.list[t, 1, iter] - (Q[t-1] * 1/(1+exp(-0.02*(Q[t-1]-80000))))
    #s2
    output.N.list[t,2,iter] <- output.N.list[t,2,iter] - (Q[t-1] * 1/(1+exp(-0.02*(Q[t-1]-100000))))
    #3
    output.N.list[t,3,iter] <- output.N.list[t,3,iter] - (Q[t-1] * 1/(1+exp(-0.02*(Q[t-1]-100000))))
    
    flowmortlist <- append(flowmortlist, (Q[t-1] * 1/(1+exp(-0.02*(Q[t-1]-100000)))))
    #replist[[1]][,,1] <- output.N.list[[1]]
    Total.N[,iter] <- apply(output.N.list[,,iter],1,sum)
  } #-------------------------
  # End Inner Loop  
  #------------------------- 
  sizearray[, iter] <- sizelist
  
} #----------------------
out# End Outer Loop
#----------------------
#------------------
# Analyzing Results
#-------------------
# summarizing iterations

## turning replist into a df

repdf <- plyr::adply(output.N.list, c(1,2,3))
dates <- plyr::adply(datearray, c(1,2))
sizes <- plyr::adply(sizearray, c(1,2))



sizedf <- as.data.frame(cbind(dates$V1, sizes$V1))
names(sizedf) <- c('date', 'size')
sizedf$date <- as_date(sizedf$date, origin= "1970-01-01")
sizedf$date <- format(sizedf$date, format = "%m-%d")

sizedf <- sizedf[order(sizedf$date),]## Taking mean results to cf w/ observed data
sizedf <- na.omit(sizedf)
means.size<- sizedf %>%
  dplyr::group_by(date) %>% # 
  dplyr::summarise(mean.size = mean(size),
                   sd.size = sd(size),
                   count =  n(),
                   se.size = sd(size)/sqrt(count)) %>%
  ungroup()
means.size$date <- as.POSIXct(means.size$date, format = "%m-%d")

ggplot(data = means.size, aes(x = date,
                              y = mean.size, group = 1)) +
  geom_ribbon(aes(ymin = mean.size - 1.96 * se.size,
                  ymax = mean.size + 1.96 * se.size),
              colour = 'transparent',
              alpha = .5,
              show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  #coord_cartesian(xlim = ylim = c(0)) +
  ylab('Baetis Size (mg)') +
  xlab('Date') +
  scale_x_datetime(labels = date_format("%b"))

means_May_BV <- means.size[132:163, ]


ggplot(data = means_May_BV, aes(x = date,
                                y = mean.size, group = 1)) +
  geom_ribbon(aes(ymin = mean.size - 1.96 * se.size,
                  ymax = mean.size + 1.96 * se.size),
              colour = 'transparent',
              alpha = .5,
              show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  #coord_cartesian(xlim = ylim = c(0)) +
  ylab('Baetis Size (mg)') +
  xlab('Date') #+
#scale_x_datetime(labels = date_format("%b"))


# repdf$timesteps <- as.numeric(as.character(repdf$timesteps))

# totn <- adply(Total.N, c(1,2))
# names(totn) <- c('timesteps', 'rep', 'tot.abund')
# totn$timesteps <- as.numeric(as.character(totn$timesteps))

## joining totn and repdf together
#repdf <- left_join(totn, repdf)
repdf <- repdf[which(repdf$X2 == "S3"), ]
repdf <- repdf[which(repdf$X1 != 27), ]
repdf <- as.data.frame(cbind(repdf, dates$V1))

names(repdf) <- c('timestep', 'stage', 'rep', 'abund', "date")

repdf$date <- as_date(repdf$date, origin= "1970-01-01")

## Taking mean results to cf w/ observed data
S3.list<- repdf %>%
  dplyr::group_by(date) %>% # combining reps
  dplyr::summarise(abund = mean(abund))

plot(S3.list$date[-1], S3.list$abund[-1], type = "l")

# adult.lost <- repdf %>%
  # dplyr::group_by(timesteps, stage) %>%
  # dplyr::summarise(mean.adult = mean(abund)) %>%
  # dplyr::filter(stage == "S3")
