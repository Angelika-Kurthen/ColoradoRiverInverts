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
source("BAETSurvivorship.R")
source("HYOSSurvivorship.R")
source("1spFunctions.R")


#read in flow data from USGS gauge at Lees Ferry, AZ between 1985 to the end of the last water year
flow <- readNWISdv("09380000", "00060", "1985-10-01", "2021-09-30")
temp <- temp <- readNWISdv("09380000", "00010", "2007-10-01", "2021-09-30")
temps <- average.yearly.temp(temp, "X_00010_00003", "Date")
n <- 13
# qr is the temp ramps I want to increase the average Lees Ferry temp by 
qr <- c(templist[te])
# how many years I want each temp ramp to last
r <- c(13)
temps <- rep.avg.year(temps, n, change.in.temp = qr, years.at.temp = r)

# basically want to create a function where we can read in the spp, a temp list, a flow list, iterations, 
sppModel <- function(temp.data, flow.data, species, iterations, Kb = 10000, Kd = 40000, Qmin = 0.25, extinction, hydropeaking, h = NULL){
# set parameters for model
  Qmin <- Qmin
  extinction <- extinction
  Kb <- Kb
  Kd <- Kb
  temps <- temp.data
  Q <- flow.data
  a <- 0.1
  g <- 0.1
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
                    dimnames = list(1:(length(timestep)+1), c("S1", "S2", "S3"), 1:iterations))
  output.N.list <- reparray

  # set species specific parameters
  if (species == "BAET"){
    G1 = 0.25 #move to Stage2 (subimago)
    G2 = 0.25 #move to Stage3 (adult)
    P1 = 0.3 #stay in Stage1 (larvae)
    P2 = 0.3 #stay in Stage2 (subimago)
    F2 = 0 # no fecundity in Stage 2
    P3 = 0 # all adults die
    h <- surv.fit.BAET$m$getPars()[2]   
    k <- surv.fit.BAET$m$getPars()[1]
  }

  if (species == "NZMS"){
    G1 = 0.1
    G2 = 0.1
    P1 = 0.4
    P2 = 0.5
    P3 = 0.8
    h <- surv.fit.NZMS$m$getPars()[2]   
    k <- surv.fit.NZMS$m$getPars()[1]
  }
  
  if (species == "HYOS"){
    G1 = 0.1 # move onto stage 2
    G2 = 0.445 # move onto stage 3
    P1 = 0.7 # remain in stage 1
    P2 = 0.0 # remain in stage 2
    F2 = 0 # no fecundity in stage 2
    P3 = 0 # all adults die
    h <- surv.fit.HYOS$m$getPars()[2]   
    k <- surv.fit.HYOS$m$getPars()[1]
  }
  for (iter in c(1:iterations)) {
  
  K = Kb # need to reset K for each iteration
  # pull random starting values from a uniform distribution 
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
    #------------------------------------------------------------------
    # Calculate fecundity per adult and 
    # Calculate effect of density dependnce on fecundity 
    
    # Logistic Density Dependence on Fecundity via Rogosch et al. Fish Model \
    # assume 97% Female poplation, and 30% instant mortality for eggs laid
    
    if (species == "BAET"){
      emergetime <- append(emergetime, back.count.degreedays(t, 559)) # value from Perry and Kennedy, 2016 
      # we start by pulling fecundities from normal distribution
    F3 = rnorm(1, mean = 1104.5, sd = 42.75) * 0.5 *0.089 * 0.15  #Baetidae egg minima and maxima from Degrange, 1960, assuming 1:1 sex ratio and 50% egg mortality
      
    # we can also relate fecundities to body mass.
    # Sweeney and Vannote 1980 have recorded dry body weight between 0.9 and 2.0 mg. 
    # That weight is related to fecundity Y = 614X - 300
    # we can "convert" emergetime to mg by multiplying to get dry weights between 0.9 - 2 mg, and then convert to fecunity
    # Issue: this data is for Ephemerella spp, not Baetidae spp
      
    if (t > 15) {
      size <- (emergetime[t-1] * 0.55)-0.75
      sizelist <- append(sizelist, size)
      F3 <- ((614 * size) - 300)* 0.5 * 0.089 * 0.15
    }
    if (hydropeaking == T){
    F3 <- F3 * hydropeaking.mortality(lower = 0.0, upper = 0.2, h = hp[t-1])
    }
    F3 <- Logistic.Dens.Dependence(F3, K, Total.N[t-1, iter])
    
    # add F_BAET to list
    Flist <- append(Flist, F3)}
    
    if(species == "NZMS"){
    # fecundities estimated from McKenzie et al. 2013; 
    if (temps$Temperature[t-1] <= 10) { 
      F2 <- 2
      F3 <- 2  
    } else {
      F2 <- 8.87473
      F3 <- 27.89665
    }    

    
    F2 <- Logistic.Dens.Dependence(F2, K, Total.N[t-1, iter]) * 0.97 * 0.7
    F3 <- Logistic.Dens.Dependence(F3, K, Total.N[t-1, iter]) * 0.97 * 0.7
    
    Flist <- append(Flist, (F2_NZMS + F3_NZMS))
    }
    
    if(species == "HYOS"){
      # Calculate how many timesteps emerging adults have matured
      emergetime <- append(emergetime, back.count.degreedays(t, 1680))
      # Calculate fecundity per adult
      
      F3 = rnorm(1, mean = 235.6, sd = 11.05102 )
      #from Willis Jr & Hendricks, sd calculated from 95% CI = 21.66 = 1.96*sd
      # * 0.5 assuming 50% female
      
      # we can scale fecundity based on the 95% CI of 21.66 (min = 213.94, max = 257.26) 
      if (t > 15) {
        size <- emergetime[t-1]
        sizelist <- append(sizelist, size)
        F3 <- ((8.664 * size) + 127.3) * 0.5
      }
    # add hydropeaking mortality
    if (hydropeaking == T){
      F3 <- F3 *hydropeaking.mortality(lower = 0.2, upper = 0.4, h = hp[t-1])
    }
      # denisty dependence
      F3 <- Logistic.Dens.Dependence(F3, K, Total.N[t-1, iter])
      }
    
    
    #-----------------------------------------------
    # Create Lefkovitch Matrix
    
    R1 <- c(P1, F2, F3)
    R2 <- c(G1, P2, 0)
    R3 <- c(0, G2, P3) 
    
    A <- rbind(R1, R2, R3)
    
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
    # Check extinction threshold 
    if (Total.N[t, iter] < extinction){
      output.N.list[t,,iter] <- 0
      Total.N[t, iter] <- 0    
  }} #-------------------------
  # End Inner Loop  
  #------------------------- 
} #----------------------
# End Outer Loop
#----------------------
  return(output.N.list)
}

sppModel(temps, flow.data = rep(0.1, times = length(temps$Temperature)), species = "HYOS", iterations = 50, extinction = 500, hydropeaking = F)
#------------------
# Analyzing Results
#-------------------
# summarizing iterations

## turning replist into a df
means.list.NZMS <- mean.data.frame(output.N.list[,2:3,], burnin = 25)
means.list.NZMS <- cbind(means.list.NZMS[27:339,], temps$dts[27:339])
means.list.NZMS$`temps$dts` <- as.Date(means.list.NZMS$`temps$dts`)
# plot abundance over time

falls <- tibble(
  x1 = c("2001-11-10", "2002-11-10", "2003-11-10", "2004-11-10", "2005-11-10", "2006-10-10"),
  x2 = c("2001-11-10", "2002-11-10", "2003-11-10", "2004-11-10", "2005-11-10", "2006-10-10"),
  y1 = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05),
  y2 = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2)
)
springs <- tibble(
  x1 = c("2007-03-31", "2008-03-31", "2009-03-31", "2010-03-31", "2011-03-31", "2012-03-31"),
  x2 = c("2007-03-31", "2008-03-31", "2009-03-31", "2010-03-31", "2011-03-31", "2012-03-31"),
  y1 = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05),
  y2 = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2)
)
falls$x1 <- as.Date(falls$x1)
falls$x2 <- as.Date(falls$x2)
springs$x1 <- as.Date(springs$x1)
springs$x2 <- as.Date(springs$x2)

abund.trends.NZMS <- ggplot(data = means.list.NZMS, aes(x = `temps$dts`,
                                                        y = mean.abund/10000, group = 1)) +
  # geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
  #                 ymax = mean.abund + 1.96 * se.abund),
  #             colour = 'transparent',
  #             alpha = .5,
  #             show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  coord_cartesian(ylim = c(0,2)) +
  ylab(paste0('New Zealand Mudsnail Abundance at ', templist[te],"°C")) +
  xlab(" ")+
  theme(text = element_text(size = 13), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 12.5))+
  scale_x_date(date_labels="%B", date_breaks  ="6 months")+
  annotate("segment", x = falls$x1, y = falls$y1, xend = falls$x2, yend = falls$y2,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")), color = "#a6611a")+
  annotate("text", x = falls$x1[1], y = 0, label = "Fall HFE (0.45 Bankflow)" ,hjust = 0, size = 5, color = "#a6611a")+
  annotate("segment", x = springs$x1, y = springs$y1, xend = springs$x2, yend = springs$y2,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")), color = "#018571")+
  annotate("text", x = springs$x1[1], y = 0, label = "Spring HFE (0.45 Bankflow)" ,hjust = 0, size = 5, color = "#018571")
ggsave(abund.trends.NZMS, filename = paste0("NZMSTempFlow_", templist[te],".png"))
plotlist <- append(plotlist, paste0("NZMSTempFlow_", templist[te],".png"))


}
library(png)
library(grid)
library(gridExtra)
for (te in 1:length(templist)){
  assign(paste0("p",te), readPNG(plotlist[te]))
}
grid.arrange(rasterGrob(p1), rasterGrob(p2), rasterGrob(p3), rasterGrob(p4), rasterGrob(p5), ncol = 3 )
