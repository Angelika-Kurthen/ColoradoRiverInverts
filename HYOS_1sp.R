#####################################
# Ramped hydropeaking index, temperature, and flows
#####################################


library(purrr)
library(tidyverse)
library(lubridate)
library(plyr)
library(dplyr)
library(ggplot2)
# data retrieval tool from USGS
library(dataRetrieval)
source("1spFunctions.R")

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



#read in flow data from USGS gauge at Lees Ferry, AZ between 1985 to the end of the last water year
#flow <- readNWISdv("09380000", "00060", "1985-10-01", "2021-09-30")
#flow.magnitude <- TimestepDischarge(flow, 85000) #discharge data from USGS is in cfs, bankfull dischage (pers comm TK) 85000 cfs
# read in temp data
#temp <- temp <- readNWISdv("09380000", "00010", "2007-10-01", "2021-09-30")
temp <- read.delim("gcmrc20230123125915.tsv", header=T)
colnames(temp) <- c("Date", "Temperature")
# peaklist <- c(0.01, 0.1, 0.2, 0.5)
# tempslist <- c(0, 0.5, 1, 2.5, 5, 7.5)

n <- 13

# qr is the temp ramps I want to increase the average temp by 
qr <- 0
# how many years I want each temp ramp to last
r <- 13

temps <- average.yearly.temp(temp, "Temperature", "Date")

temps <- rep.avg.year(temps, n, change.in.temp = qr, years.at.temp = r)
discharge <- rep(0.1, times = length(temps$Temperature))
HYOSmodel <- function(flow.data, temp.data, baselineK, disturbanceK, Qmin, extinct, iteration, peaklist = NULL, peakeach = NULL){
#---------------------------------------------------------------
# set up model
source("HYOSSurvivorship.R")

Q <- as.numeric(flow.data)
temps <- temp.data

degreedays <- as.data.frame(cbind(temps$dts, temps$Temperature * 14))
colnames(degreedays) <- c("dts", "DegreeDay")
degreedays$dts <- as.Date(degreedays$dts, origin = "1970-01-01")

# need to make ramped increasing hydropeaking index 
hp <- c(rep(peaklist, each = peakeach))

# specify iterations
iterations <- iteration

# baseline K in the absence of disturbance
Kb <- as.numeric(baselineK)
# max K after a big disturbance
Kd <- as.numeric(disturbanceK)

# specify baseline transition probabilities for each species
# 3 stages - we have egg - larval instar V, pupae, and adult

G1 = 0.1 # move onto stage 2
G2 = 0.445 # move onto stage 3
P1 = 0.7 # remain in stage 1
P2 = 0.0 # remain in stage 2

# want to run this for one year, in 14 day timesteps 
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

Qmin <- Qmin
a <- 0.1
g <- 0.1
h <- surv.fit.HYOS$m$getPars()[2]  
k <- surv.fit.HYOS$m$getPars()[1] 

extinction <- extinct
#-------------------------
# Outer Loop of Iterations
#--------------------------


for (iter in c(1:iterations)) {
  K = Kb # need to reset K for each iteration
  
  # pull random values from a uniform distribution 
  output.N.list[1,1:3, iter]<- runif(3, min = 1, max = (0.3*K))
  
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
    
    emergetime <- append(emergetime, back.count.degreedays(t, 1680, degreedays))
    #---------------------------------------------------------
    # Calculate fecundity per adult
    
    F3 = rnorm(1, mean = 235.6, sd = 11.05102 ) * 0.5 * hydropeaking.mortality(lower = 0.4, upper = 0.6, h = hp[t-1])
    #from Willis Jr & Hendricks, sd calculated from 95% CI = 21.66 = 1.96*sd
    # * 0.5 assuming 50% female
    
    # we can scale fecundity based on the 95% CI of 21.66 (min = 213.94, max = 257.26) 
    if (t > 15) {
      size <- emergetime[t-1]
      sizelist <- append(sizelist, size)
      F3 <- ((8.664 * size) + 127.3) * 0.5 * hydropeaking.mortality(lower = 0.4, upper = 0.6, h = hp[t-1])
    }
    
    #---------------------------------------------------
    # Calculate the disturbance magnitude-K relationship 
    # Sets to 0 if below the Qmin
    Qf <- as.numeric(Qf.Function(Q[t-1], Qmin, a))
    
    #-------------------------------------------------------------------
    # Calculate K arrying capacity immediately following the disturbance
    # Calculate K arrying capacity immediately following the disturbance
    K0 <- as.numeric(K + ((Kd-K)*Qf))
    
    # Calculate final K for timestep, including relationship between K and time since disturbance
    K <- post.dist.K(K0, Kb, g, t, Q, Qmin)
    Klist <- append(Klist, K)
    
    
    #---------------------------------------------
    # Calculate effect of density dependence on fecundity
    
    # Logistic via Rogosch et al. Fish Model
    F3 <- Logistic.Dens.Dependence(F3, K, Total.N[t-1, iter]) * 0.5
    Flist <- append(Flist, F3)
    
    # Calculate new transition probabilities based on temperature
    # This is the growth v development tradeoff 
    # don't know if this exists for HYOS - they can emerge under a wide temp gradient (<5 - 25+ C) but relationship between growth and temp 
    
    #development measures (basically, if below 10C, no development, if between 10 and 12, follows a function, if above 12, prob of transition to next stage is 0.6395)
    if (5 > temps$Temperature[t-1]) {
      G1 <- 0.001
      G2 <- 0.001}
    
    if (temps$Temperature[t-1] > 25){
      G1 <-0.799
      G2 <-0.444}
    if (5 <= temps$Temperature[t-1] & temps$Temperature[t-1] <= 25) G1 <- growth.development.tradeoff(temps$Temperature[t-1], 5, 25, 0.001, 0.799)
    if (5 <= temps$Temperature[t-1] & temps$Temperature[t-1] <= 25) G2 <- growth.development.tradeoff(temps$Temperature[t-1], 5, 25, 0.001, 0.444)
    P1 <- 0.8 - G1
    P2 <- 0.445 - G2
    
   
    #-----------------------------------------------
    # Create Lefkovitch Matrix
    
    H1 <- c(P1, 0, F3)
    H2 <- c(G1, P2, 0)
    H3 <- c(0, G2, 0) 
    
    A <- rbind( H1, H2, H3)
    
    #-----------------------------------------------
 
    # # growth (if below 10C, no growth can occur - everything basically freezes, if between 10 and 11, prob of remaining in same stage = 0.6395, if above 13, prob of transition to next stage is 0 
    
    #Glist <-append(Glist, AHYOS[3,2])
    #Plist <- append(Plist, AHYOS[2,2])
    
    #--------------------------------------
    # Calculate abundances for each stage
    
    output.N.list[t, 1:3, iter] <- A %*% output.N.list[t-1, 1:3, iter] 
    
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
return(output.N.list)
}


out <- HYOSmodel(flow.data = discharge, temp.data = temps, disturbanceK = 40000, baselineK = 10000, Qmin = 0.24, extinct = 500, iteration = 5, peaklist = 0, peakeach = length(temps$Temperature))
#------------------
# Analyzing Results
#-------------------
# summarizing iterations
means.list.HYOS <- mean.data.frame(out, burnin = 27, iteration = 5)
means.list.HYOS <- cbind(means.list.HYOS[27:length(means.list.HYOS$mean.abund),], temps$dts[27:length(means.list.HYOS$mean.abund)])
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
# arrows <- tibble(
#   x1 = c("2005-01-07", "2007-01-07", "2009-01-07", "2011-01-07"),
#   x2 = c("2005-01-07", "2007-01-07", "2009-01-07", "2011-01-07"),
#   y1 = c(1.45, 1.45, 1.45, 1.45), 
#   y2 = c(1, 1, 1, 1)
# )
# 
# arrows$x1 <- as.Date(arrows$x1)
# arrows$x2 <- as.Date(arrows$x2)

falls$x1 <- as.Date(falls$x1)
falls$x2 <- as.Date(falls$x2)
springs$x1 <- as.Date(springs$x1)
springs$x2 <- as.Date(springs$x2)

abund.trends.HYOS <- ggplot(data = means.list.HYOS, aes(x =  `temps$dts`,
                                                        y = mean.abund/10000, group = 1)) +
  geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
                 ymax = mean.abund + 1.96 * se.abund),
             colour = 'transparent',
             alpha = .5,
             show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  coord_cartesian(ylim = c(0,2)) +
  ylab(paste0('Hydrospyche spp. Abundance adding ', tempslist[te], "°C")) +
  xlab(" ")+
  theme(text = element_text(size = 13), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13))+
  scale_x_date(date_labels="%B", date_breaks  ="6 months", limits = as.Date(c("2001-01-27", "2012-12-04"
  )))+
  
  # annotate("segment", x = arrows$x1, y = arrows$y1, xend = arrows$x2, yend = arrows$y2,
  #          arrow = arrow(type = "closed", length = unit(0.02, "npc")), color = "red")+
  # annotate("text", x = arrows$x1[1], y = 1.5, label = "HI = 0.01", size = 4)+
  # annotate("text", x = arrows$x1[2], y = 1.5, label = "HI = 0.1", size = 4)+
  # annotate("text", x = arrows$x1[3], y = 1.5, label = "HI = 0.2", size = 4)+
  # annotate("text", x = arrows$x1[4], y = 1.5, label = "HI = 0.5", size = 4 )

  annotate("segment", x = falls$x1, y = falls$y1, xend = falls$x2, yend = falls$y2,
        arrow = arrow(type = "closed", length = unit(0.02, "npc")), color = "#a6611a")+
  annotate("text", x = falls$x1[1], y = 0, label = "Fall HFE (0.45 Bankflow)" ,hjust = 0, size = 5, color = "#a6611a")+
  annotate("segment", x = springs$x1, y = springs$y1, xend = springs$x2, yend = springs$y2,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")), color = "#018571")+
  annotate("text", x = springs$x1[1], y = 0, label = "Spring HFE (0.45 Bankflow)" ,hjust = 0, size = 5, color = "#018571")
  # 
  ggsave(abund.trends.HYOS, filename = paste0("HYOSTempFlowHI", tempslist[te], peaklist[pe],".png"))
  #ggsave(abund.trends.HYOS, filename = paste0("HYOSTempFlowHI", tempslist[te],".png"))
  plotlist <- append(plotlist, paste0("HYOSTempFlowHI", tempslist[te], peaklist[pe],".png"))

  
}
plots <- lapply(ll <- plotlist ,function(x){
  img <- as.raster(readPNG(x))
  rasterGrob(img, interpolate = FALSE)
})
  ggsave(filename = paste0("HYOSTempFlowHI",peaklist[pe],".pdf"),width=8.5, height=11, 
       marrangeGrob(grobs = plots, nrow = 3, ncol=2))
  plotlist <- NULL
}
#   for (te in 1:length(tempslist)){
#   assign(paste0("p",te), readPNG(plotlist[te]))}
#   grid.arrange(rasterGrob(p1), rasterGrob(p2), rasterGrob(p3), rasterGrob(p4), rasterGrob(p5), rasterGrob(p6), ncol = 3 )
# 
# # 
# # 
# pdf(paste0("HYOSTemp_", tempslist[te], "_Flood_HI_", peaklist[pe]), width = 8.27, height = 11.69)
# gr <- grid.arrange(rasterGrob(p1), rasterGrob(p2), rasterGrob(p3), rasterGrob(p4), rasterGrob(p5), rasterGrob(p6), ncol = 3 )
# 
# library(patchwork) 
# 
# 
# plotlist <- append(plotlist, paste0("HYOSTempFlowHI", tempslist[te], peaklist[pe],".png"))
# 
# 
# ggsave(abund.trends.HYOS, filename = paste0("HYOSTempFlowHI", tempslist[te], peaklist[pe],".png"))
# plotlist <- append(plotlist, paste0("HYOSTempFlowHI", tempslist[te], peaklist[pe],".png"))
