##########################
# HYOS 1 sp model
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

source("HYOSSurvivorship.R")
source("1spFunctions.R")
source("ColoradoRiverTempRamp.R")

#read in flow data from USGS gauge at Lees Ferry, AZ between 1985 to the end of the last water year
flow <- readNWISdv("09380000", "00060", "1985-10-01", "2021-09-30")
# read in temperature data from USGS gauge at Lees Ferry, AZ between _ to the end of last water year
temp <- read.delim("gcmrc20230123125915.tsv", header=T)

#---------------------------------------------------------------
# this chunk of code makes the repeating avg ColRiv temp series
colnames(temp) <- c("Date", "Temperature")
temp$Date <- as_datetime(temp$Date)
temp$Date <- yday(temp$Date)


# option 1: for each date, calculate mean and standard error so we can create a new dataset
temp <- temp %>% group_by(Date) %>% dplyr::summarise(Temperature = mean(Temperature))
temp

# Make an index to be used for aggregating
ID <- as.numeric(as.factor(temp$Date)) 
# want it to be every 14 days, hence the 14
ID <- ID %/% 14

ID[which(ID ==26)] <- 25
# aggregate over ID and TYPEall numeric data.
outs <- aggregate(temp[sapply(temp,is.numeric)],
                  by=list(ID),
                  FUN=mean)



# format output
names(outs)[2:3] <-c("dts","Temperature")
# add the correct dates as the beginning of every period
outs$dts <- strptime(round(outs$dts), "%j") ###Note need to subtract 365 if larger than 365

# order by date in chronological order#
#outs <- outs[order(outs$dts),]
outs$dts <- as_date(outs$dts)

n <- 13

# repeat this data frame for 13 years
temp_seq <- do.call("rbind", replicate(n, outs, simplify = FALSE))
temp_seq$dts <- as.Date(temp_seq$dts)
# now adjust the years so time can proceed in chronological order
year(temp_seq$dts[1:26]) <- 2000
year(temp_seq$dts[27:52]) <- 2001
year(temp_seq$dts[53:79]) <- 2002
year(temp_seq$dts[80:104]) <- 2003
year(temp_seq$dts[105:130]) <- 2004
year(temp_seq$dts[131:156]) <- 2005
year(temp_seq$dts[157:182]) <- 2006
year(temp_seq$dts[183:208]) <- 2007
year(temp_seq$dts[209:234]) <- 2008
year(temp_seq$dts[235:260]) <- 2009
year(temp_seq$dts[261:286]) <- 2010
year(temp_seq$dts[287:312]) <- 2011
year(temp_seq$dts[313:338]) <- 2012
#--------------------------------------
# years 6 and 7  represent a 1 C increase from baseline
temp_seq$Temperature[131:182] <- temp_seq$Temperature[131:182] + 1

# years 8 and 9 represent a 2.5 C increase from baseline
temp_seq$Temperature[183:234] <- temp_seq$Temperature[183:234] + 2.5

# years 10 and 11 represent at 5 C increase from baseline
temp_seq$Temperature[235:286] <- temp_seq$Temperature[235:286] + 5

# year 11 and 12 represent a 7.5 C increase from baseline
temp_seq$Temperature[287:338] <- temp_seq$Temperature[287:338] + 7.5

temps <- temp_seq
degreedays <- as.data.frame(cbind(temp_seq$dts, temps$Temperature * 14))
colnames(degreedays) <- c("dts", "DegreeDay")

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
Qmin <- 0.25
a <- 0.1
g <- 0.1
h <- surv.fit.HYOS$m$getPars()[2]   
k <- surv.fit.HYOS$m$getPars()[1] 

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
    
      F_HYOS = rnorm(1, mean = 235.6, sd = 11.05102 ) * 0.5  #from Willis Jr & Hendricks, sd calculated from 95% CI = 21.66 = 1.96*sd
    # * 0.5 assuming 50% female
      
      # we can scale fecundity based on the 95% CI of 21.66 (min = 213.94, max = 257.26) 
      if (t > 15) {
        size <- emergetime[t-1]
        sizelist <- append(sizelist, size)
        F_HYOS <- ((8.664 * size) + 127.3) * 0.5
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
    # 
    # # growth (if below 10C, no growth can occur - everything basically freezes, if between 10 and 11, prob of remaining in same stage = 0.6395, if above 13, prob of transition to next stage is 0 )
    # if (10 >= temps$Temperature[t-1]) AHYOS[2,2] <- 0.55
    # if (temps$Temperature[t-1] > 13) AHYOS[2,2] <- 0.001
    # if (10 < temps$Temperature[t-1] & temps$Temperature[t-1] <=  13) AHYOS[2,2] <- (-0.183 * temps$Temperature[t-1]) + 2.38 #(-0.1 temps$Temperature[t-1]) - 2.6
    # 
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
means.list.HYOS <- mean.data.frame(output.N.list, stages = c(1,2,3), burnin = 27)
means.list.HYOS <- cbind(means.list.HYOS[27:339,], temps$dts[27:339])
means.list.HYOS$`temps$dts` <- as.Date(means.list.HYOS$`temps$dts`)
# plot abundance over time

arrows <- tibble(
  x1 = c("2005-01-07", "2007-01-07", "2009-01-07", "2011-01-07"),
  x2 = c("2005-01-07", "2007-01-07", "2009-01-07", "2011-01-07"),
  y1 = c(1.4500, 1.4500, 1.4500, 1.4500), 
  y2 = c(0.9000, 0.9000, 0.9000, 0.9000)
)

arrows$x1 <- as.Date(arrows$x1)
arrows$x2 <- as.Date(arrows$x2)

abund.trends.HYOS <- ggplot(data = means.list.HYOS, aes(x =  `temps$dts`,
                                              y = mean.abund/10000, group = 1)) +
  geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
                  ymax = mean.abund + 1.96 * se.abund),
              colour = 'transparent',
              alpha = .5,
              show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  coord_cartesian(ylim = c(0.6,1.5000)) +
  ylab('Hydrospyche spp. Relative Abundance') +
  xlab(" ")+
  theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13))+
  scale_x_date(date_labels="%B", date_breaks  ="6 months")+
  annotate("segment", x = arrows$x1, y = arrows$y1, xend = arrows$x2, yend = arrows$y2,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")), color = "red")+
  annotate("text", x = arrows$x1[1], y = 1.5000, label = "+1°C", size = 5)+
  annotate("text", x = arrows$x1[2], y = 1.5000, label = "+2.5°C", size = 5)+
  annotate("text", x = arrows$x1[3], y = 1.5000, label = "+5°C", size = 5)+
  annotate("text", x = arrows$x1[4], y = 1.5000, label = "+7.5°C", size = 5 )




# take a look at results
# 
# par(mfrow = c(1,1))
plot(timestep[9:(length(timestep)+1)], output.N.list[9:(length(timestep)+1), 3, 1], type = "l", ylab = "Hydropsyche spp. Adults", xlab = "Timestep (1 fortnight)")
plot(timestep[9:length(timestep)], Total.N[10:(length(timestep)+1)], type= "l", ylab = "Hydropsyche spp. Total N", xlab = "Timestep (1 fortnight)")
# 
# 
# 
# #creating plots to analyze how temp relationship is working
#create dataframe with timestemp, s1, s2, and s3 abundances, and tempterature
 data <- as.data.frame(cbind(timestep, output.N.list[2:(length(timestep)+1) ,1, 1], output.N.list[2:(length(timestep)+1) ,2, 1], output.N.list[2:(length(timestep)+1) ,3, 1], temps$Temperature))
 colnames(data) <- c("timestep", "Stage1", "Stage2", "Stage3", "Temperature")

 data <- data[210:260, ]
ggplot(data = data, aes(x = timestep, y = Stage1, color = "Stage1"))+
   geom_path()+
   geom_path(aes(x = timestep, y = Stage2, color = "Stage2"))+
   geom_path(aes(x = timestep, y = Stage3, color = "Stage3"))+
   geom_path(aes(x = timestep, y = Temperature*200, color = "Temperature"))+
   scale_y_continuous(

     # Features of the first axis
     name = "Abundance",

     # Add a second axis and specify its features
     sec.axis = sec_axis( ~.*0.005, name="Temperature C")
   )
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