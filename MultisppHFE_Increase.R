##########
# Multispecies Hydropeaking Simulation
#############

# Load required scripts with custom functions and models
source("Multispp.R")
source("MultisppFunctions.R")

# Read in Lees Ferry temperature and discharge data from 2007 to 2023
#temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")  # Water temperature data
#discharge <- readNWISdv("09380000", "00060", "2007-10-01", "2023-05-01")  # River discharge data

# read in from csv to save time
temp <- read.csv2("LFtemp2007to2023.csv", header= T)
discharge <- read.csv2("LFdischarge2007to2023.csv", header = T)

# Calculate average yearly flows from discharge data
flow <- average.yearly.flows(flowdata = discharge, flow.column_name = "X_00060_00003", date.column_name = "Date")

# Calculate average yearly temperatures from temperature data
temps <- average.yearly.temp(tempdata = temp, temp.column_name = "X_00010_00003", date.column_name = "Date")

# Create a time series of average temperatures extending for 100 years with no temperature change
temps <- rep.avg.year(temps, n = 100, change.in.temp = 0, years.at.temp = 0)

# Create a time series of average flows for 100 years by repeating the original yearly data
flows <- do.call("rbind", replicate(100, flow, simplify = FALSE))

# Match the dates between temperature and flow datasets
flows$dts <- as.Date(temps$dts)

# Normalize discharge values by dividing by the bankfull discharge (85,000 cfs)
flows$Discharge <- flows$Discharge / 85000

# Create a sequence of hydropeaking intensity levels from 0.00 to 0.70 in increments of 0.05
HFEs <- seq(0, 2, by = 0.25)
HFEs <- 1
# Get row indices for high flow events
HFE_rows <- which(month(flows$dts) == 3 & day(flows$dts) == 31)
HFE_rows <- HFE_rows[seq(1, length(HFE_rows), by = 3)]  # Keep every 3rd occurrence

# Initialize lists to store results for abundance, biomass, and annual biomass calculations
Multispp_temp_abund <- data.frame(timesteps= numeric(), taxa = factor(), abundance=numeric(), temperature=factor())
Multispp_temp_biomass <- data.frame(timesteps= numeric(), taxa = factor(), biomass=numeric(), temperature=factor())

# Loop through each hfe intensity level
results <- mclapply(HFEs, function(HFE) {
  set.seed(123)
  flows$Discharge[HFE_rows] <- HFE
  
  out <- Multispp(flow.data = flows$Discharge, temp.data = temps, baselineK = 10000, 
                  disturbanceK = 100000, Qmin = 0.25, extinct = 50, iteration = 1000, 
                  peaklist = 0, peakeach = length(temps$Temperature), stage_output = c("all", "biomass"))
  
  means.abund <- multispp.data.frame(out$N, burnin = 260, iteration = 1000, value = "abund")
  means.biomass <- multispp.data.frame(out$biomass, burnin = 260, iteration = 1000, value = "biomass")

  # Convert to data frame instead of cbind matrix
  average_means <- data.frame(means.abund, temperature = rep(HFE, nrow(means.abund)))
  colnames(average_means) <- colnames(Multispp_temp_abund)
  
  average_biomass <- data.frame(means.biomass, temperature = rep(HFE, nrow(means.biomass)))
  colnames(average_biomass) <- colnames(Multispp_temp_biomass)
  
  # Ensure factor consistency
  average_means$taxa <- as.factor(average_means$taxa)
  average_means$temperature <- as.factor(average_means$temperature)
  
  average_biomass$taxa <- as.factor(average_biomass$taxa)
  average_biomass$temperature <- as.factor(average_biomass$temperature)
  
  return(list(Multispp_temp_abund = average_means, 
              Multispp_temp_biomass = average_biomass))
 # }, mc.cores = 1)
}, mc.cores = detectCores() - 1)

Multispp_temp_abund <- do.call(rbind, lapply(results, `[[`, "Multispp_temp_abund"))
Multispp_temp_biomass <- do.call(rbind, lapply(results, `[[`, "Multispp_temp_biomass"))

# Write results to CSV files
write.csv(Multispp_temp_abund, "Multispp_temp_HFE_abund.csv", row.names = FALSE)
write.csv(Multispp_temp_biomass, "Multispp_temp_HFE_biomass.csv", row.names = FALSE)


