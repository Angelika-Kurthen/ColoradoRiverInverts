#############################
## 1 species Model Functions
#############################

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
      if(vec <= criticaldegreedays) {vec <- degreedays$DegreeDay[s] + vec
      emerg <- NA}
      else {emerg <- time - s
      break}
      # once we hit that threshold, we count the number of timesteps it took to reach that and add that to our emergetime vector
    }
  }
  return(emerg)
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
  if (is.na(tau)==T | tau == 0) { tau <-  0
  K <- K0} 
  if (tau > 0) {
    K <- Kb + ((K0 - Kb)*exp(-g*tau))} # function from McMullen et al 2017, g is shape function
  return(K)
}


# Function to calculate logistic density dependence on fecundity, after Rogosch et al 2019
Logistic.Dens.Dependence <- function(Fecundity, K, N){
  f.rate <- Fecundity * checkpos((K - N)/K) 
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

#function to summarize code into mean population abundance over iterations
mean.data.frame <- function(data, stages, burnin){
  repdf <- plyr::adply(data, stages)
  names(repdf) <- c('timesteps', 'stage', 'rep', 'abund')
  repdf$timesteps <- as.numeric(as.character(repdf$timesteps))
  
  totn <- plyr::adply(Total.N, c(1,2))
  names(totn) <- c('timesteps', 'rep', 'tot.abund')
  totn$timesteps <- as.numeric(as.character(totn$timesteps))
  
  # joining totn and repdf together
  repdf <- dplyr::left_join(totn, repdf)
  
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
  
  if (burnin == T){
    means.list <- means.list[ , burnin:means.list$timesteps]
  }
  return(means.list)
}

