####################
## Multispecies? 
####################
library(dataRetrieval)

source("1spFunctions.R")
source("HYOSSurvivorship.R")
source("BAETSurvivorship.R")
source("NZMS_shell_length_fecundity.R")
source("NZMSSurvivorship.R")
source("CHIRSurvivorship.R")
source("GAMMSurvivorship.R")
# temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
# temps <- TimestepTemperature(temp)
# discharge <- readNWISdv("09380000", "00060", "2007-10-01", "2023-05-01")
# flow.magnitude <- TimestepDischarge(discharge, 85000)
# flow.data <- flow.magnitude$Discharge
# temp.data <- temps
# baselineK <- 10000
# disturbanceK <- 40000
# Qmin <- 0.25
# extinct <- 50
# iteration <- 1
# peakeach <- length(temps$Temperature)
# peaklist <- 0
# stage_output <- "size"
# discharge <temps# discharge <- rep(0.1, times = length(temps$Temperature))
#modify_paramter <- NULL
Multispp <- function(flow.data, temp.data, baselineK, disturbanceK, Qmin,
                     extinct, iteration, peaklist = NULL, peakeach = NULL, 
                     stage_output = "all", modify_parameter = NULL, increment = NULL ){

  Q <- as.numeric(flow.data)
  temps <- temp.data
  
  degreedays <- as.data.frame(cbind(temps$dts, temps$Temperature * 14))
  colnames(degreedays) <- c("dts", "DegreeDay")
  degreedays$DegreeDay[degreedays$DegreeDay<0] <- 0
  degreedays$dts <- as.POSIXct(degreedays$dts, origin = "1970-01-01")
  
  # could do ramped hydropeaking, we do not - all
  hp <- c(rep(peaklist, each = peakeach))
  
  # Compute hydropeaking mortality for different taxa
  hydro.mort_HYOS <- hydropeaking.mortality(lower = 0.4, upper = 0.6, h = unique(hp))
  hydro.mort_BAET <- hydropeaking.mortality(lower = 0.0, upper = 0.15, h = unique(hp))
  hydro.mort_GAMM <- hydropeaking.mortality(lower = 0, upper = 1, h = unique(hp))
  hydro.mort_NZMS <- hydropeaking.mortality(lower = 0, upper = 1, h = unique(hp))
  hydro.mort_CHIR <- hydropeaking.mortality(lower = 0, upper = 0.6, h = unique(hp))  # Fixed h = unique(hp)
  
  
  if (!is.null(modify_parameter) && nzchar(modify_parameter)) {
    # List of valid parameter names
    valid_parameters <- c("hydro.mort_HYOS", "hydro.mort_BAET", "hydro.mort_GAMM", "hydro.mort_NZMS", "hydro.mort_CHIR")
    
    # Check if modify_parameter is a valid parameter name
    if (modify_parameter %in% valid_parameters) {
      current_value <- get(modify_parameter)
      sens <- checkpos(current_value + increment)
      sens <- ifelse(sens > 1, 1, sens)
      assign(modify_parameter, sens, envir = .GlobalEnv)  
      
    }}
  # mortality_lookup <- sapply(unique_h, function(hp, lower, upper) {
  #   hydropeaking.mortality(lower = X, upper = Y, h = h)
  # })
  # names(mortality_lookup) <- unique_h  # Name the vector for fast lookup
  # specify iterations
  iterations <- iteration
  
  # baseline K in the absence of disturbance
  Kb <- as.numeric(baselineK)
  # max K after a big disturbance
  Kd <- as.numeric(disturbanceK)
  
  # want to run this for one year, in 14 day timesteps 
  timestep <- seq(2, (length(temps$Temperature) + 1), by = 1)
  
  # create array to put the total N of all species into
  Total.Biomass <- array(0,
                   dim  <-c((length(timestep) +1 ), iterations),
                   dimnames <- list(1:(length(timestep) + 1), 1:iterations))
  
  # create list of arrays w/ abundance data for each spp
  reparray <- array(0,
                    dim = c(length(timestep) + 1, 3, iterations, 5),
                    dimnames = list(1:(length(timestep)+1), c("S1", "S2", "S3"), 1:iterations, c("HYOS", "BAET", "NZMS", "CHIR", "GAMM")))
  
  # one for Ns 
  output.N.list <- reparray
  # one for abundances
  output.Biomass.list <- reparray
  
  Qmin <- Qmin
  
  a <- 0.001 # the same for all
  g <- 1 # the same for all 
  
  # load h and k estimates for Negative Exp relationship between 
  HYOS_h <- unname(surv.fit.HYOS$m$getPars()[2])
  HYOS_k <- unname(surv.fit.HYOS$m$getPars()[1])
  
  BAET_h <- unname(surv.fit.BAET$m$getPars()[2])
  BAET_k <- unname(surv.fit.BAET$m$getPars()[2])
  
  NZMS_h <- unname(surv.fit.NZMS$m$getPars()[2])
  NZMS_k <- unname(surv.fit.NZMS$m$getPars()[1])
  
  CHIR_h <- unname(surv.fit.CHIR$m$getPars()[2])
  CHIR_k <- unname(surv.fit.CHIR$m$getPars()[1])
  
  GAMM_h <- unname(surv.fit.GAMM$m$getPars()[2])
  GAMM_k <- unname(surv.fit.GAMM$m$getPars()[1])
  
  extinction <- extinct # the same for all
  
  flood.mort_HYOS <- sapply(Q, flood.mortality, N = 1, k = HYOS_k, h = HYOS_h, Qmin = Qmin)
  flood.mort_BAET <- sapply(Q, flood.mortality, N = 1, k = BAET_k, h = BAET_h, Qmin = Qmin)
  flood.mort_CHIR <- sapply(Q, flood.mortality, N = 1, k = CHIR_k, h = CHIR_h, Qmin = Qmin)
  flood.mort_GAMM <- sapply(Q, flood.mortality, N = 1, k = GAMM_k, h = GAMM_h, Qmin = Qmin)
  flood.mort_NZMS <- sapply(Q, flood.mortality, N = 1, k = NZMS_k, h = NZMS_h, Qmin = Qmin)
  
  if (!is.null(modify_parameter) && nzchar(modify_parameter)) {
    # List of valid parameter names
  valid_parameters <- c("flood.mort_HYOS", "flood.mort_BAET", "flood.mort_GAMM", "flood.mort_NZMS", "flood.mort_CHIR")
  
  # Check if modify_parameter is a valid parameter name
  if (modify_parameter %in% valid_parameters) {
        current_value <- get(modify_parameter)
        sens <- checkpos(current_value + increment)
       sens <- ifelse(sens > 1, 1, sens)
        assign(modify_parameter, sens, envir = .GlobalEnv)  
    
  }}
  
  # temperature dependent survival
  TempSurvival_HYOS <- sapply(temps$Temperature, TempSurv_HYOS)
  TempSurvival_BAET <- sapply(temps$Temperature, TempSurv_BAET)
  TempSurvival_NZMS <- sapply(temps$Temperature, TempSurv_NZMS)
  TempSurvival_CHIR <- sapply(temps$Temperature, TempSurv_CHIR)
  TempSurvival_GAMM <- sapply(temps$Temperature, TempSurv_GAMM)
  
  if (!is.null(modify_parameter) && nzchar(modify_parameter)) {
    # List of valid parameter names
    valid_parameters <- c("TempSurvival_HYOS", "TempSurvival_BAET", "TempSurvival_NZMS", "TempSurvival_CHIR", "TempSurvival_GAMM")
    
    # Check if modify_parameter is a valid parameter name
    if (modify_parameter %in% valid_parameters) {
      current_value <- get(modify_parameter)
      sens <- checkpos(current_value + increment)
      sens <- ifelse(sens > 1, 1, sens)
      assign(modify_parameter, sens, envir = .GlobalEnv)  
    }}
  # Calculate how many timesteps emerging adults have matured
  emergetime_HYOS <- sapply(timestep, back.count.degreedays, criticaldegreedays = 1680, degreedays)
  emergetime_BAET <- sapply(timestep, back.count.degreedays, criticaldegreedays = 250, degreedays)
  emergetime_CHIR <- sapply(timestep, back.count.degreedays, criticaldegreedays = 600, degreedays) # value from Ali et al 1985
  emergetime_GAMM <- sapply(timestep, back.count.degreedays, criticaldegreedays = 1000, degreedays) # value from Sweeney et al 2017
  
  # how big are those adults?
  sizelist_HYOS <- emergetime_HYOS 
  sizelist_BAET <- emergetime_BAET
  sizelist_CHIR <- 3*emergetime_CHIR-6
  sizelist_GAMM <- emergetime_GAMM
  
  #----------------------------------------------------------
  # Outer Loop of Iterations
  #--------------------------
  
  for (iter in c(1:iterations)) {
    K = Kb # need to reset K for each iteration
    
    # pull random values from a uniform distribution 
    output.N.list[1,1:3, iter, ]<- runif(15, min = 1, max = (0.3*K))
    # ok can also start everyone the same 
    # output.N.list[1,1:3, iter, ] <- c(1000, 100,10)

    # calculate biomass at t = 1, we don't know size, so we assume mean
    output.Biomass.list[1,1:3, iter, "HYOS"] <- output.N.list[1,1:3, iter, "HYOS"] * (0.0046 * 12.73^2.926) #mean size is around 12.73846
    output.Biomass.list[1,1:3, iter, "BAET"] <- output.N.list[1,1:3, iter, "BAET"] * (0.0053 * 3.183^2.875) # mean size is around 3.186528
    output.Biomass.list[1,1:3, iter, "CHIR"] <- output.N.list[1,1:3, iter, "CHIR"] * (0.0018 * 10.795^2.617)
    
    # size strucured so we know specific values
    output.Biomass.list[1,1, iter, "NZMS"] <- output.N.list[1,1, iter, "NZMS"]*(0.02 * mean(c(0.5, 3.2))^2.4315)
    output.Biomass.list[1,2, iter, "NZMS"] <- output.N.list[1,2, iter, "NZMS"]*(0.02 * mean(c(3.2, 4))^2.4315)
    output.Biomass.list[1,3, iter, "NZMS"] <- output.N.list[1,3, iter, "NZMS"]*(0.02 * mean(c(4, 5.5))^2.4315)
    
    output.Biomass.list[1,1, iter, "GAMM"] <- output.N.list[1,1, iter, "GAMM"] *(0.063 * mean(c(2.5, 7))^2.46)
    output.Biomass.list[1,2, iter, "GAMM"] <- output.N.list[1,1, iter, "GAMM"] *(0.063 * mean(c(7, 9))^2.46)
    output.Biomass.list[1,3, iter, "GAMM"] <- output.N.list[1,1, iter, "GAMM"] *(0.063 * mean(c(9, 12))^2.46)
    
    # fill in first total.biomass
    Total.Biomass[1, iter] <- sum(output.Biomass.list[1,1:3, iter, ])
    # list to input Ks
    Klist <- vector()
    Klist[1] <- K
    
    # fecundity
    Flist_HYOS <- numeric(length(timestep) +1)
    Flist_BAET <- numeric(length(timestep) +1)
    Flist_NZMS <- numeric(length(timestep) +1)
    Flist_CHIR <- numeric(length(timestep) +1)
    Flist_GAMM <- numeric(length(timestep) +1)

    #-------------------------
    # Inner Loop of Timesteps
    #-------------------------
    
    for (t in timestep) {
      t <- t
      #---------------------------------------------------------
      # Calculate fecundity per adult
      # Hydropsyche spp.
      F3_HYOS = 235.6 * hydro.mort_HYOS
      #from Willis Jr & Hendricks* 0.5 assuming 50% female.
      # we can scale fecundity to size (based on emergetime) to fecundities
      if (t > 15) {
        F3_HYOS <- ((7.219 * sizelist_HYOS[t-1]) + 180.4) * hydro.mort_HYOS
      }
      
      #Baetidae spp.
      F3_BAET = 1104.4 * hydro.mort_BAET
      #Baetidae egg maxima * 0.5 from Degrange, 1960, assuming 1:1 sex ratio and 50% egg mortality
      # fecundities estimated from McKenzie et al. 2013 - reduced fecundity above 24 C and below 9 C. 
      # optimal temp between 16 and 19 C, but we don't really have parameterization for that
      
      if (t > 19) {
        F3_BAET <- ((200*sizelist_BAET[t-1])+200) * hydro.mort_BAET
      }
      #NZMS
      F2_NZMS <- 8.87473 * (-0.0001427 *(temps$Temperature[t-1] - 17.5)^4 + 1) *  hydro.mort_NZMS

      F3_NZMS <-  27.89665 *(-0.0001427 * (temps$Temperature[t-1] - 17.5)^4 + 1) * hydro.mort_NZMS
      # equation sometimes goes below 0, so make sure it doesn't
      if (F2_NZMS < 0){
        F2_NZMS <- 0
      }
      if (F3_NZMS < 0){
        F3_NZMS <- 0
      }
  
      # CHIR
      F3_CHIR = 300 * 0.5 * hydro.mort_CHIR
      #CHIR egg # and % mortality from Charles et al 
      # we can also relate fecundities to body size which is between 6 and 15 mm (also from Charles et al 2004)
      # we can "convert" emergetime to size by multiplying to get size between 6 and 15 mm and then convert to fecunity
      
      if (t > 19) {
        F3_CHIR <- (13.33*sizelist_CHIR[t-1]+180) * 0.5* hydro.mort_CHIR
      }
      
      #GAMM
      F2_GAMM <- 10.3 * 0.5 * hydro.mort_GAMM
      
      F3_GAMM <-  23.3 * 0.5 * hydro.mort_GAMM
      # # #   
      if (F2_GAMM < 0){
        F2_GAMM <- 0
      }
      if (F3_GAMM < 0){
        F3_GAMM <- 0
      }
      
      if (t > 19) {
        F2 <-  (1.934*(sizelist_GAMM[t-1]) - 4.67) *0.5* hydro.mort_GAMM
        F3 <-  ( 2.812*(sizelist_GAMM[t-1]) + 2.94)*0.5* hydro.mort_GAMM #* 0.78 * 0.65
      }    
      #---------------------------------------------------
      # Calculate the disturbance magnitude-K relationship- this will be the same
      
      # Sets to 0 if below the Qmin
      Qf <- as.numeric(Qf.Function(Q[t-1], Qmin, a))
      
      #-------------------------------------------------------------------
      # Calculate K arrying capacity immediately following the disturbance
      K0 <- as.numeric(K + ((Kd-K)*Qf))
      
      # Calculate final K for timestep, including relationship between K and time since disturbance
      K <- post.dist.K(K0, Kb, g, t, Q, Qmin)
      Klist <- append(Klist, K)
      
      #---------------------------------------------
      # Calculate effect of community density dependence on fecundity
      
      # Logistic via Rogosch et al. Fish Model
      F3_HYOS <- Logistic.Dens.Dependence(F3_HYOS, K, Total.Biomass[t-1, iter])
      Flist_HYOS[t] <- F3_HYOS
      
      F3_BAET <- Logistic.Dens.Dependence(F3_BAET, K, Total.Biomass[t-1, iter]) 
      Flist_BAET[t] <-  F3_BAET
      
      
      F2_NZMS <- Logistic.Dens.Dependence(F2_NZMS, K, Total.Biomass[t-1, iter])  
      F3_NZMS <- Logistic.Dens.Dependence(F3_NZMS, K, Total.Biomass[t-1, iter]) 
      #once again, sometimes they goes below 0, so make sure not negative
      if (F2_NZMS < 0){
        F2_NZMS <- 0
      }
      if (F3_NZMS < 0){
        F3_NZMS <- 0
      }
      
      Flist_NZMS[t] <- F2_NZMS + F3_NZMS
      
      # CHIR
      F3_CHIR <- Logistic.Dens.Dependence(F3_CHIR, K, Total.Biomass[t-1, iter])
      # 
      # add F_CHIR to list
      Flist_CHIR[t] <- F3_CHIR
      
      #GAMM
      F2_GAMM <- Logistic.Dens.Dependence(F2_GAMM, K, Total.Biomass[t-1, iter])
      F3_GAMM <- Logistic.Dens.Dependence(F3_GAMM, K, Total.Biomass[t-1, iter]) 
      #
      if (F2_GAMM < 0){
        F2_GAMM <- 0
      }
      if (F3_GAMM < 0){
        F3_GAMM <- 0
      }
      # add F_BAET to list
      Flist_GAMM[t] <- F3_GAMM+F2_GAMM
      #------------------------------------------------
       
      # Calculate new transition probabilities based on temperature
      # This is the growth v development tradeoff 
      # at cold temps, takes 20 timesteps to complete Stage 1
      # at warm temps, take 1 timestep to complete Stage 1
      
      # ------------------------------------------------
      # for Hydropsyche spp. 
      #include temperature dependent mortality
      # section 1 is for extreme temperatures 
      if (5 > temps$Temperature[t-1])  {
        G1_HYOS <- 0.36/20 *TempSurvival_HYOS[t-1]
        P1_HYOS <- 1-(1/20)* TempSurvival_HYOS[t-1]
      }
      if (30 > temps$Temperature[t-1]){ # Spend about 20 timesteps in Stage 2 if less than 30 C
        G2_HYOS <- 0.04/20 * TempSurvival_HYOS[t-1] 
        P2_HYOS <- 1-(1/20) *TempSurvival_HYOS[t-1]
      } else {
        G2_HYOS = 0.04/3 *TempSurvival_HYOS[t-1]  # make them grow fast when super warm
        P2_HYOS = 1 - (1/3) *TempSurvival_HYOS[t-1]
      }
      
      if (temps$Temperature[t-1] > 30){ # immediately move onto next stage when warm
        G1_HYOS <- 0.36 *TempSurvival_HYOS[t-1]
        P1_HYOS <- 0
      }
      
      # if within standard water temps, use stage duration
      # Stage 1 can vary from a few timesteps to many (saw Instar 2 last from 5 to 250 days)
      if (5 <= temps$Temperature[t-1] & temps$Temperature[t-1] <=30 & is.na(emergetime_HYOS[t-1] == F)){
        G1_HYOS <- 0.36/(emergetime_HYOS[t-1]) *TempSurvival_HYOS[t-1]
        P1_HYOS <- 1-(1/(emergetime_HYOS[t-1])) *TempSurvival_HYOS[t-1]
      }
      
      # if stage duration not available, approximate linearly 
      if (5 <= temps$Temperature[t-1] & temps$Temperature[t-1] <= 30 & is.na(emergetime_HYOS[t-1] == T)) {
        G1_HYOS <- 0.36/((-0.95 * temps$Temperature[t-1]) + 24.75) *TempSurvival_HYOS[t-1]
        P1_HYOS <- 1-(1/((-0.95 * temps$Temperature[t-1]) + 24.75))*TempSurvival_HYOS[t-1]
      }
      
      #-----------------------------------------------
      # for Baetidae spp.
      
      
      # development measures
      # in this function, we assume that if below the min temp threshold (9) no maturation occurs (slow maturation, large growth)
      # if above the max temp threshold (30), no one remains more than 1 timestep in each stage (fast maturation, small growth)
      
      if (9 > temps$Temperature[t-1]) {
        P1_BAET <- (1-(1/9)) *TempSurvival_BAET[t-1]
        P2_BAET <- P1_BAET
        G1_BAET <- 0.3/9 *TempSurvival_BAET[t-1]
        G2_BAET <- G1_BAET
      }
      if (temps$Temperature[t-1] > 30){
        P1_BAET <- (1-(1/1.5)) *TempSurvival_BAET[t-1]
        P2_BAET <- P1_BAET
        G1_BAET <- 0.3/1.5 *TempSurvival_BAET[t-1]
        G2_BAET <- G1_BAET
      }
      
      if (9 <= temps$Temperature[t-1] & temps$Temperature[t-1] <= 30 & is.na(emergetime_BAET[t-1] == F)){
        G1_BAET <- 0.3/((emergetime_BAET[t-1])/2) *TempSurvival_BAET[t-1]
        G2_BAET <- G1_BAET
        P1_BAET <- (1-(1/((emergetime_BAET[t-1])/2))) *TempSurvival_BAET[t-1]
        P2_BAET <- P1_BAET
      }
      if (9 <= temps$Temperature[t-1] & temps$Temperature[t-1] <= 30 & is.na(emergetime_BAET[t-1] == T)) {
        G1_BAET <- 0.3/((-0.786 * temps$Temperature[t-1]) + 18) *TempSurvival_BAET[t-1]
        P1_BAET <- (1-(1/((-0.786 * temps$Temperature[t-1]) + 18))) *TempSurvival_BAET[t-1]
        G2_BAET <- G1_BAET
        P2_BAET <- P1_BAET
      }
      #-----------------------------------------------
      # NZMS
      # its speculated (Cross et al 2010) that survivorship is between 80 - 100% for NZMS in Grand Canyon - will say 80% survival to be on the conservative side
      # from timestep to timestep, we expect 90% to survive so for stage 1 (which lasts aproximately 14 timesteps, survival should be approx 09^14 = 0.2287679
      # from that, only 1/14th will transition out, 13/14 remain in stage (Birt et al 2009)
      
      stageduration1_NZMS <- timestep_to_mat(temps$Temperature[t-1])[[1]]
      stageduration2_NZMS <- timestep_to_mat(temps$Temperature[t-1])[[2]]
      
      G1_NZMS = (0.91/stageduration1_NZMS) *TempSurvival_NZMS[t-1]
      G2_NZMS = (0.91/stageduration2_NZMS)*TempSurvival_NZMS[t-1]
      P1_NZMS = (1-(1/stageduration1_NZMS)) *TempSurvival_NZMS[t-1]
      P2_NZMS = (1-(1/stageduration2_NZMS)) *TempSurvival_NZMS[t-1]
      P3_NZMS = (1-(1/7))*TempSurvival_NZMS[t-1]
      
      #-----------------------------------------------
      # Chironomidae spp. 
      if (5 > temps$Temperature[t-1]) {
        P1_CHIR <- (1-(1/20)) * TempSurvival_CHIR[t-1]
        P2_CHIR <- P1_CHIR
        G1_CHIR <- (0.74/20) * TempSurvival_CHIR[t-1]
        G2_CHIR <- (0.66/20) * TempSurvival_CHIR[t-1]
      }
      if (temps$Temperature[t-1] > 30){
        P1_CHIR <- 0
        P2_CHIR <- 0
        G1_CHIR <- 0.74 * TempSurvival_CHIR[t-1]
        G2_CHIR <- 0.66 * TempSurvival_CHIR[t-1]
      }
      
      if (5 <= temps$Temperature[t-1] & temps$Temperature[t-1] <= 30 & (is.na(emergetime_CHIR[t-1]) == F)){
        G1_CHIR <- (0.74/((emergetime_CHIR[t-1])/2)) * TempSurvival_CHIR[t-1]
        G2_CHIR <- (0.66/((emergetime_CHIR[t-1])/2)) * TempSurvival_CHIR[t-1]
        P1_CHIR <- (1-(1/((emergetime_CHIR[t-1])/2))) * TempSurvival_CHIR[t-1]
        P2_CHIR <- P1_CHIR
      }
      if (5 <= temps$Temperature[t-1] & temps$Temperature[t-1] <= 30 & (is.na(emergetime_CHIR[t-1] == T))) {
        G1_CHIR <- (0.74*((-0.136 * temps$Temperature[t-1]) + 5.088)) * TempSurvival_CHIR[t-1]
        P1_CHIR <- (1-(1/((-0.136 * temps$Temperature[t-1]) + 5.088))) * TempSurvival_CHIR[t-1]
        G2_CHIR <- (0.66*((-0.136 * temps$Temperature[t-1]) + 5.088)) * TempSurvival_CHIR[t-1]
        P2_CHIR <- P1_CHIR
      }
      
      # For G. lacustris
      if (temps$Temperature[t-1] < 10){
        G1_GAMM = 0
        G2_GAMM = 0
      }
      
      if (is.na(emergetime_GAMM[t-1])== F){
        G1_GAMM = (0.95/(emergetime_GAMM[t-1]/2)) *TempSurvival_GAMM[t-1]
        G2_GAMM = (0.99/(emergetime_GAMM[t-1]/2))*TempSurvival_GAMM[t-1]
        P1_GAMM = (1-(1/(emergetime_GAMM[t-1]/2))) *TempSurvival_GAMM[t-1]
        P2_GAMM = (1-(1/(emergetime_GAMM[t-1]/2))) *TempSurvival_GAMM[t-1]
        P3_GAMM = (1-(0.64))*TempSurvival_GAMM[t-1] }#survival for last stage 
      
      if(is.na(emergetime_GAMM[t-1])==T){
        G1_GAMM = (0.95/( ((-0.233 * temps$Temperature[t-1]) + 11)/2)) *TempSurvival_GAMM[t-1]
        G2_GAMM = (0.99/( ((-0.233 * temps$Temperature[t-1]) + 11)/2))*TempSurvival_GAMM[t-1]
        P1_GAMM = (1-(1/( ((-0.233 * temps$Temperature[t-1]) + 11)/2))) *TempSurvival_GAMM[t-1]
        P2_GAMM = (1-(1/( ((-0.233 * temps$Temperature[t-1]) + 11)/2))) *TempSurvival_GAMM[t-1]
        P3_GAMM = (1-(0.64))*TempSurvival_GAMM[t-1] }
      
      #-----------------------------------------------
      # Create Lefkovitch Matrix
      if (!is.null(modify_parameter) && nzchar(modify_parameter)) {
      
      # List of valid parameter names
      valid_parameters <- c("G1_HYOS", "G2_HYOS", "G1_BAET", "G2_BAET", 
                            "G1_CHIR", "G2_CHIR", "G1_GAMM", "G2_GAMM", 
                            "G1_NZMS", "G2_NZMS")
      
      # Check if modify_parameter is a valid parameter name
      if (modify_parameter %in% valid_parameters) {
        current_value <- get(modify_parameter)
        sens <- checkpos(current_value + increment)
        sens <- ifelse(sens > 1, 1, sens)
        assign(modify_parameter, sens)  
        
      }}
      

      # Compile matrices
      HYOS1 <- c(P1_HYOS, 0, F3_HYOS)
      HYOS2 <- c(G1_HYOS, P2_HYOS, 0)
      HYOS3 <- c(0, G2_HYOS, 0) 
      
      A_HYOS <- rbind( HYOS1, HYOS2, HYOS3)
      
      BAET1 <- c(P1_BAET, 0, F3_BAET)
      BAET2 <- c(G1_BAET, P2_BAET, 0)
      BAET3 <- c(0, G2_BAET, 0) 
      
      A_BAET <- rbind( BAET1, BAET2, BAET3)
      
      NZMS1 <- c(P1_NZMS, F2_NZMS, F3_NZMS)
      NZMS2 <- c(G1_NZMS, P2_NZMS, 0)
      NZMS3 <- c(0, G2_NZMS, P3_NZMS) 
      
      A_NZMS <- rbind(NZMS1, NZMS2, NZMS3)
      
      
      CHIR1 <- c(P1_CHIR, 0, F3_CHIR)
      CHIR2 <- c(G1_CHIR, P2_CHIR, 0)
      CHIR3 <- c(0, G2_CHIR, 0)
      
      A_CHIR <- rbind(CHIR1, CHIR2, CHIR3)
      
      GAMM1 <- c(P1_GAMM, F2_GAMM, F3_GAMM)
      GAMM2 <- c(G1_GAMM, P2_GAMM, 0)
      GAMM3 <- c(0, G2_GAMM, P3_GAMM)
      
      A_GAMM <- rbind(GAMM1, GAMM2, GAMM3)
      
      
      
      #--------------------------------------
      # Calculate abundances for each stage
      
      output.N.list[t, 1:3, iter, "HYOS"] <- A_HYOS %*% output.N.list[t-1, 1:3, iter, "HYOS"] 
      
      # Calculate abundances for each stage
      
      output.N.list[t, 1:3, iter, "BAET"] <- A_BAET %*% output.N.list[t-1, 1:3, iter, "BAET"] 
      
      # Calculate abudances for each stage
      output.N.list[t, 1:3, iter, "NZMS"] <- A_NZMS %*% output.N.list[t-1, 1:3, iter, "NZMS"] 
      
      # Calculate abundances for each stage
      output.N.list[t, 1:3, iter, "CHIR"] <- A_CHIR %*% output.N.list[t-1, 1:3, iter, "CHIR"] 
      
      # Calculate abundances for each stage
      output.N.list[t, 1:3, iter, "GAMM"] <- A_GAMM %*% output.N.list[t-1, 1:3, iter, "GAMM"] 
      
      #------------------------------------------
      #Calculate immediate mortality due to flows
      # mortality due to flooding follows N0 = Nz*e^-hQ
      # but what about using a sigmoidal logistic function so we can control the threshold point and rate of growth
      # following m = 1/1+e^-h*(x-xf)
      # where h is is shape value
      # x is Q, and xf is threshold point (100% of pop dies)
      #plot(x = Q, y = 1/(1+exp(-0.02*(Q-100000))))
      #------------------------------------------------------------
      # Hydropsyche spp.
      #s1
      output.N.list[t, 1, iter, "HYOS"] <- output.N.list[t, 1, iter, "HYOS"] * flood.mort_HYOS[t-1]
      #s2
      output.N.list[t,2,iter, "HYOS"] <- output.N.list[t, 2, iter, "HYOS"] * flood.mort_HYOS[t-1]
      #3
      #output.N.list[t,3,iter, "HYOS"] <- flood.mortality(output.N.list[t,3,iter, "HYOS"], HYOS_k, HYOS_h, Q[t-1], Qmin)
      
      #Baetidae spp.
      #s1
      output.N.list[t, 1, iter, "BAET"] <-  output.N.list[t, 1, iter, "BAET"] * flood.mort_BAET[t-1]
      #s2
      output.N.list[t,2,iter, "BAET"] <- output.N.list[t,2,iter, "BAET"] * flood.mort_BAET[t-1]

      #NZMS
      #s1
      output.N.list[t, 1, iter, "NZMS"] <- output.N.list[t, 1, iter, "NZMS"] * flood.mort_NZMS[t-1]
      #s2
      output.N.list[t,2,iter, "NZMS"] <- output.N.list[t,2,iter, "NZMS"] * flood.mort_NZMS[t-1]
      #3
      output.N.list[t,3,iter, "NZMS"] <- output.N.list[t,3,iter, "NZMS"] * flood.mort_NZMS[t-1]
      
      # CHIR
      #s1
      output.N.list[t, 1, iter, "CHIR"] <- output.N.list[t, 1, iter, "CHIR"] * flood.mort_CHIR[t-1]
      #s2
      output.N.list[t,2,iter, "CHIR"] <- output.N.list[t,2,iter, "CHIR"] * flood.mort_CHIR[t-1]

      #GAMM
      #s1
      output.N.list[t, 1, iter, "GAMM"] <- output.N.list[t, 1, iter, "GAMM"] * flood.mort_GAMM[t-1]
      #s2Qt
      output.N.list[t,2,iter, "GAMM"] <- output.N.list[t,2,iter, "GAMM"] * flood.mort_GAMM[t-1]
      
      output.N.list[t,3,iter, "GAMM"] <- output.N.list[t,3,iter, "GAMM"] * flood.mort_GAMM[t-1]
      
      #check extinction threshold for each spp + add rescue effect
      # Hydrospyche spp. 
      if (sum(output.N.list[t, ,iter, "HYOS"]) < extinction){
        output.N.list[t,1:3,iter, "HYOS"] <- c(50,10,5) #rescue effect
      } 
      
      #Baetidae spp
      if (sum(output.N.list[t, ,iter, "BAET"]) < extinction){
        output.N.list[t,1:3,iter, "BAET"] <- c(50,10,5)
      } 
      #NZMS
      if (sum(output.N.list[t, ,iter, "NZMS"]) < extinction){
        output.N.list[t,1:3,iter, "NZMS"] <- c(50,10,5)
      } 
      #CHIR
      if (sum(output.N.list[t, ,iter, "CHIR"]) < extinction){
        output.N.list[t,1:3,iter, "CHIR"] <- c(50,10,5)
      } 
      
      if (sum(output.N.list[t, ,iter, "GAMM"]) < extinction){
        output.N.list[t,1:3,iter, "GAMM"] <- c(50,10,5)
      } 
      # convert to biomass
      # for timesteps before we have backcounted stage duration  we use mean
      # these will be removed anyway due to burnin
      
      # for hydropsyche spp.
      output.Biomass.list[t,1:3, iter, "HYOS"] <- output.N.list[t,1:3, iter, "HYOS"] * (0.0046 * 12.73^2.926) #mean size is around 12.73846
      if (t > 19) {
        #s1
        output.Biomass.list[t, 1, iter, "HYOS"] <- output.N.list[t,1, iter, "HYOS"] * (0.0046 * emergetime_HYOS[t-1]^2.926)
        #s2
        output.Biomass.list[t, 2, iter, "HYOS"] <- output.N.list[t,2, iter, "HYOS"] * (0.0046 * (emergetime_HYOS[t-1]+3)^2.926)
        #s3
        output.Biomass.list[t, 3, iter, "HYOS"] <- output.N.list[t,3, iter, "HYOS"] * (0.0046 * (emergetime_HYOS[t-1]+3)^2.926)
      }
      # for Baetidae spp.
      output.Biomass.list[t,1:3, iter, "BAET"] <- output.N.list[t,1:3, iter, "BAET"] * (0.0053 * 3.183^2.875) # mean size is around 3.186528
      
      if (t > 19 ){
        #S1
        output.Biomass.list[t,1, iter, "BAET"] <- output.N.list[t,1, iter, "BAET"] * (0.0053 * (emergetime_BAET[t-1]/2)^2.875) 
        #S2
        output.Biomass.list[t,2, iter, "BAET"] <- output.N.list[t,2, iter, "BAET"] * (0.0053 * emergetime_BAET[t-1]^2.875) 
        #S3
        output.Biomass.list[t,3, iter, "BAET"] <- output.N.list[t,3, iter, "BAET"] * (0.0053 * emergetime_BAET[t-1]^2.875) 
      }
      # for NZMS
        #S1
        output.Biomass.list[t,1, iter, "NZMS"] <- output.N.list[t,1, iter, "NZMS"] * (0.02 * mean(c(0.5, 3.2))^2.4315)
        #S2
        output.Biomass.list[t,2, iter, "NZMS"] <- output.N.list[t,2, iter, "NZMS"] * (0.02 * mean(c(3.2, 4))^2.4315)
        #S3
        output.Biomass.list[t,3, iter, "NZMS"] <- output.N.list[t,3, iter, "NZMS"] * (0.02 * mean(c(4, 5.5))^2.4315)
       
        # for CHIR
        output.Biomass.list[t,1:3, iter, "CHIR"] <- output.N.list[1,1:3, iter, "CHIR"] * (0.0018 * 10.795^2.617)
      
        if (t > 19){
        output.Biomass.list[t,1, iter, "CHIR"] <- output.N.list[t,1, iter, "CHIR"] * (0.0018 * (sizelist_CHIR[t-1]/2)^2.617)
        #S2
        output.Biomass.list[t,2, iter, "CHIR"] <- output.N.list[t,2, iter, "CHIR"] * (0.0018 * (sizelist_CHIR[t-1])^2.617)
        #S3
        output.Biomass.list[t,3, iter, "CHIR"] <- output.N.list[t,3, iter, "CHIR"] * (0.0018 * (sizelist_CHIR[t-1])^2.617)
      }
        
        output.Biomass.list[t,1, iter, "GAMM"] <- output.N.list[t,1, iter, "GAMM"]* (0.063 * mean(c(2.5, 7))^2.46)
        output.Biomass.list[t,2, iter, "GAMM"] <- output.N.list[t,2, iter, "GAMM"]* (0.063 * mean(c(7, 9))^2.46)
        output.Biomass.list[t,3, iter, "GAMM"] <- output.N.list[t,3, iter, "GAMM"]* (0.063 * mean(c(9, 12))^2.46)
        
      # calculate total biomass
      Total.Biomass[t, iter] <- sum(output.Biomass.list[t, ,iter,])
      }
      #-------------------------
    }  # End Inner Loop  
    #------------------------- 
  #----------------------
  # End Outer Loop
  #----------------------
  results <- list()
  
  if ("larvae" %in% stage_output) {
    results$N <- output.N.list[ ,1:2, ,]
  }  
  if ("all" %in% stage_output) {
    results$N <- output.N.list[ , 1:3, ,]
  }  
  if ("3" %in% stage_output) {
    results$N <- output.N.list[ , 3, ,]
  }  
  if ("size" %in% stage_output) {
    results$size <- cbind(sizelist_HYOS, sizelist_BAET, sizelist_CHIR)
  }  
  if ("biomass" %in% stage_output) {
    results$biomass <- output.Biomass.list[ ,1:3, ,]
  }  
  
  # Return the list
  return(results)
}


# out <- Multispp(flow.data = flow.data, temp.data = temp.data,
#                            baselineK =10000 , disturbanceK = 100000, Qmin = 0.25,
#                            extinct = 50, iteration = 2, peaklist = 0.17,
#                            peakeach = length(temps$Temperature), stage_output = c("biomass", "size"), modify_parameter = "TempSurvival_HYOS", increment = 0.001)



                     
