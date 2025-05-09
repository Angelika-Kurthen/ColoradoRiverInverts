############################
# Make partial eta files for multispp temp + hfe --> run on HPC
################################

source("Multispp.R")
source("MultisppFunctions.R")


# Load necessary libraries
library(doParallel)
library(foreach)
library(dataRetrieval)
# Read in Lees Ferry temperature and discharge data from 2007 to 2023
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")  # Water temperature data
discharge <- readNWISdv("09380000", "00060", "2007-10-01", "2023-05-01")  # River discharge data

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

# Get row indices for high flow events
HFE_rows <- which(month(flows$dts) == 3 & day(flows$dts) == 31)
HFE_rows <- HFE_rows[seq(1, length(HFE_rows), by = 3)]  # Keep every 3rd occurrence


# Update Discharge to 0.45 for those dates
flows$Discharge[HFE_rows] <- 0.45

# Create a sequences
temp_seq <- c(1,1.1,1.2,1.5)
z_seq <- c(1, 2, 4, 8)

# Initialize final storage data frames

mclapply(z_seq, function(z_val) {
  set.seed(123)
  zs <- list(HYOS = z_val, BAET = z_val, NZMS = z_val, CHIR = z_val, GAMM = z_val)
  
  results_for_z <- lapply(temp_seq, function(te) {
    temps$Temperature <- temps$Temperature * te
    
    out <- Multispp(flow.data = flows$Discharge, temp.data = temps,
                    baselineK = 10000, disturbanceK = 100000,
                    Qmin = 0.25, extinct = 50, iteration = 1000,
                    peaklist = 0.17, peakeach = length(temps$Temperature),
                    stage_output = c("all", "biomass"), z = zs)
    
    means.biomass <- multispp.data.frame(out$biomass, burnin = 260, iteration = 1000, value = "biomass")
    means.biomass$temperature <- as.factor(te)
    means.biomass$q <- as.factor(z_val)
    
    temps$Temperature <- temps$Temperature / te
    return(means.biomass)
  })
  
  combined_df <- do.call(rbind, results_for_z)
  
  fname <- file.path(paste0("Multispp_Temp_Hyd_HFE", "_Z", z_val, ".rds"))
  dir.create(dirname(fname), recursive = TRUE, showWarnings = FALSE)
  saveRDS(combined_df, file = fname)
  
  return(NULL)
}, mc.cores = detectCores() - 1, mc.preschedule = FALSE)

## Now we add the summer spike to temperatures
# read in LF temp and discharge data from 2007 to 2023
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
# calculate average yearly temperatures
temps <- average.yearly.temp(tempdata = temp, temp.column_name = "X_00010_00003", date.column_name = "Date")
# create summertime spike (up to 21 C, then scale from there)
temps$Temperature[16:23] <- c(14, 16, 18, 21, 21, 18, 16, 14)

# create a timeseries of average temperatures 100 years long
temps <- rep.avg.year(temps, n = 100, change.in.temp = 0, years.at.temp = 0)
mclapply(z_seq, function(z_val) {
  set.seed(123)
  zs <- list(HYOS = z_val, BAET = z_val, NZMS = z_val, CHIR = z_val, GAMM = z_val)
  
  results_for_z <- lapply(temp_seq, function(te) {
    temps$Temperature <- temps$Temperature * te
    
    out <- Multispp(flow.data = flows$Discharge, temp.data = temps,
                    baselineK = 10000, disturbanceK = 100000,
                    Qmin = 0.25, extinct = 50, iteration = 1000,
                    peaklist = 0.17, peakeach = length(temps$Temperature),
                    stage_output = c("all", "biomass"), z = zs)
    
    means.biomass <- multispp.data.frame(out$biomass, burnin = 260, iteration = 1000, value = "biomass")
    means.biomass$temperature <- as.factor(te)
    means.biomass$q <- as.factor(z_val)
    
    temps$Temperature <- temps$Temperature / te
    return(means.biomass)
  })
  
  combined_df <- do.call(rbind, results_for_z)

  fname <- file.path(paste0("Multispp_Temp_Hyd_HFE_Spike", "_Z", z_val, ".rds"))
  dir.create(dirname(fname), recursive = TRUE, showWarnings = FALSE)
  saveRDS(combined_df, file = fname)
  
  return(NULL)
}, mc.cores = detectCores() - 1, mc.preschedule = FALSE)