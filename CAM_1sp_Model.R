source("BAETSurvivorship.R")
source("1spFunctions.R")

#read in flow data from USGS gauge at Lees Ferry, AZ between 1985 to the end of the last water year
#flow <- readNWISdv("09380000", "00060", "1985-10-01", "2021-09-30")
# read in temperature data from USGS gauge at Lees Ferry, AZ between _ to the end of last water year
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2021-09-30")
temps <- average.yearly.temp(temp, "X_00010_00003", "Date")
n <- 13
# qr is the temp ramps I want to increase the average Lees Ferry temp by 
qr <- c(0,0, 0, 0, 0)
# how many years I want each temp ramp to last
r <- c(5, 2, 2, 2, 2)
temps <- rep.avg.year(temps, n, change.in.temp = qr, years.at.temp = r)

degreedays <- as.data.frame(cbind(temps$dts, temps$Temperature * 14))
colnames(degreedays) <- c("dts", "DegreeDay")
degreedays$dts <- as.Date(degreedays$DegreeDay, origin = "1970-01-01")

#emergetime <- back.count.degreedays(525)



#flow.magnitude <- TimestepDischarge(flow, 85000) #discharge data from USGS is in cfs, bankfull dischage (pers comm TK) 85000 cfs
#temps <- TimestepTemperature(temp, "Colorado River") # calculate mean temperature data for each timestep
#degreedays <- TimestepDegreeDay(temp, "Colorado River")

# specify iterations
iterations <- 50
#species <- c("EPSU", "SIMU", "CHIRO")

# baseline K in the absence of disturbance
Kb <- 10000
# max K after a big disturbance
Kd <- 40000

# specify baseline transition probabilities for each species
G1_EPSU = 0.3
G2_EPSU = 0.3
P1_EPSU = 0.0
P2_EPSU = 0.0

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
#Q <- flow.magnitude$Discharge
Q <- rep(0.1, length(temps$Temperature))


Qmin <- 0.25
a <- 0.1
g <- 0.1
h <- surv.fit.BAET$m$getPars()[2]   
k <- surv.fit.BAET$m$getPars()[1]   
extinction <- 50


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
  
  Glist <- vector()
  Plist <- vector()
  #-------------------------
  # Inner Loop of Timesteps
  #-------------------------
  
  for (t in timestep) {
    
    #----------------------------------------------------------
    # Calculate how many timesteps emerging adults have matured
    
    emergetime <- append(emergetime, back.count.degreedays(t, 559)) 
    #---------------------------------------------------------
    # Calculate fecundity per adult
    # Calculate fecundity per adult
    
    # we start by pulling fecundities from normal distribution
    # relate fecundities to temperature based on Sweeney and Vannote 1981  *0.5 assuming 50% female and * 0.5 assuming 50% mort
    F_EPSU = 281.05  #mean female egg mass between 532.3 and 780.2
    
    # we can also relate time in larval stage to fecundity. Wagner (1995) estimated 225 eggs per female and Perry and Kennedy () estimated 337.7. We can scale time spent (emergetime) to fit those values
    if (t > 10) {
      F_EPSU <- ((emergetime[t-1]*56.35) -0.4) }
    #----------------------------------------------
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
    F_EPSU <- Logistic.Dens.Dependence(F_EPSU, K, Total.N[t-1, iter])
    
    Flist <- append(Flist, F_EPSU)
    
    
    # Probabilities of remaining in stages (when temps low, high prob of remaining)
    P1_BAET <- growth.development.tradeoff(temps$Temperature[t-1],  9, 13, 0.43, 0.0)
    P2_BAET <- growth.development.tradeoff(temps$Temperature[t-1], 9, 13, 0.43, 0)
    
    # growth (if below 10C, no growth can occur - everything basically freezes, if between 10 and 11, prob of remaining in same stage = 0.6395, if above 13, prob of transition to next stage is 0 )
    G1_BAET <- 0.43 - P1_BAET
    G2_BAET <- 0.43 - P2_BAET
    #-----------------------------------------------
    # Create Lefkovitch Matrix
    
    EPSU1 <- c(P1_EPSU, 0, F_EPSU)
    EPSU2 <- c(G1_EPSU, P2_EPSU, 0)
    EPSU3 <- c(0, G2_EPSU, 0) 
    
    AEPSU <- rbind( EPSU1, EPSU2, EPSU3)
    
    #-----------------------------------------------
    # Calculate new transition probabilities based on temperature
    # This is the growth v development tradeoff

    # # # development measures (basically, if below 10C, no development, if between 10 and 12, follows a function, if above 12, prob of transition to next stage is 0.6395
    # AEPSU[3, 2] <- growth.development.tradeoff(temps$Temperature[t-1],  10, 13, 0.001, 0.6)
    # AEPSU[2,1] <- growth.development.tradeoff(temps$Temperature[t-1], 10, 13, 0.001, 0.6)
    # # #
    # # # # growth (if below 10C, no growth can occur - everything basically freezes, if between 10 and 11, prob of remaining in same stage = 0.6395, if above 13, prob of transition to next stage is 0 )
    # AEPSU[2,2] <- growth.development.tradeoff(temps$Temperature[t-1], 10, 13, 0.6, 0.001)
    # AEPSU[1,1] <- growth.development.tradeoff(temps$Temperature[t-1], 10, 13, 0.6, 0.001)
    # Glist <- append(Glist, AEPSU[2,1])
    # # Plist <- append(Plist, AEPSU[1,1])
    # # development measures (basically, if below 10C, no development, if between 10 and 12, follows a function, if above 12, prob of transition to next stage is 0.6395
    # ABAET[3, 2] <- growth.development.tradeoff(temps$Temperature[t-1],  min(temps$Temperature), max(temps$Temperature), 0.001, 0.5)
    # ABAET[2,1] <- growth.development.tradeoff(temps$Temperature[t-1], min(temps$Temperature), max(temps$Temperature), 0.001, 0.44)
    # 
    # # growth (if below 10C, no growth can occur - everything basically freezes, if between 10 and 11, prob of remaining in same stage = 0.6395, if above 13, prob of transition to next stage is 0 )
    # ABAET[2,2] <- 0.44 - ABAET[3,2]
    # ABAET[1,1] <-0.44 - ABAET[2,1]
    #--------------------------------------
    # Calculate abundances for each stage
    
    output.N.list[t, 1:3, iter] <- AEPSU %*% output.N.list[t-1, 1:3, iter] 
    
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
means.list.EPSU <- mean.data.frame(output.N.list, burnin = 27)

abund.trends.EPSU <- ggplot(data = means.list.EPSU, aes(x = timesteps,y = mean.abund, group = 1)) +
      geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
                      ymax = mean.abund + 1.96 * se.abund),
                  colour = 'transparent',
                  alpha = .5,
                  show.legend = FALSE) +
      geom_line(show.legend = FALSE) +
      coord_cartesian(ylim = c(0,100000)) +
      ylab('EPSU Abundance') +
      xlab('Timestep')


