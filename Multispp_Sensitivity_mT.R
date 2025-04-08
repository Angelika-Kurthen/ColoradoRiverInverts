###################
## Multispp Sensitivity Temperature Mortality
##################
library(ggplot2)
library(dplyr)
library(dataRetrieval)
library(parallel)

source("1spFunctions.R")
source("HYOSSurvivorship.R")
source("BAETSurvivorship.R")
source("NZMS_shell_length_fecundity.R")
source("NZMSSurvivorship.R")
source("CHIRSurvivorship.R")
source("GAMMSurvivorship.R")
source("MultisppFunctions.R")
source("Multispp.R")

# Read in Lees Ferry temperature and discharge data from 2007 to 2023
# temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")  # Water temperature data
# discharge <- readNWISdv("09380000", "00060", "2007-10-01", "2023-05-01")  # River discharge data
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

# Define your species-stage parameters
parameters <- c("TempSurvival_HYOS", "TempSurvival_BAET", "TempSurvival_CHIR", "TempSurvival_NZMS", "TempSurvival_GAMM")
# Sensitivity increments from -0.01 to +0.01 in 0.001 steps
increments <- seq(-0.001, 0.001, by = 0.0001)

temp_seq <- c(1, 1.1, 1.2, 1.5)
# # Create a grid of all parameter/increment combinations
param_grid <- expand.grid(Parameter = parameters, Increment = increments, 
                          stringsAsFactors = FALSE)
# 
# # # Function to modify the parameter and run the model
run_model_with_modification <- function(param, inc) {
  # Modify the model based on param and increment (Replace this with your actual model function)
  modified_model <- Multispp(flow.data = flows$Discharge, temp.data = temps,
                            baselineK =10000 , disturbanceK = 100000, Qmin = 0.25,
                            extinct = 50, iteration = 1000, peaklist = 0,
                            peakeach = length(temps$Temperature), stage_output = c("all", "biomass"),
                            modify_parameter = param, increment =inc)

    average_abundance <- apply(modified_model$N[-c(1:259),,,], c(2, 4), function(x) mean(x, na.rm = TRUE))
    average_biomass <- apply(modified_model$N[-c(1:259),,,], c(2, 4), function(x) mean(x, na.rm = TRUE))
    # Create a data frame from the average abundances + biomass
    average_abundance_df <- as.data.frame(as.table(average_abundance))
    average_biomass_df <- as.data.frame(as.table(average_biomass))
    output <- cbind(average_abundance_df, average_biomass_df$Freq)
    colnames(output) <- c("Var1", "Var2", "Abundance", "Biomass")
    # Combine stage and taxon into a single column
    output$StageGroup <- paste(output$Var1,"_",output$Var2, sep = "")
    # Rename columns for clarity
    output <- as.data.frame(cbind(output[ , c("StageGroup", "Abundance", "Biomass")],
                                                rep(param, times = length(output$Var1)),
                                                rep(inc,times = length(output$Var1))))
    colnames(output) <- c("StageGroup", "Abundance", "Biomass", "Parameter", "SensitivityIncrement")
    # Extract relevant values (ensure your model function returns these)
    data.frame(
      output
  )
}
# 
# # Determine the number of cores to use (leave one core free)
# numCores <- detectCores() - 1
# 
# # Use mclapply to loop over each row in the grid in parallel.
# # Note: mclapply returns a list that we then combine into a single data frame.
# sensitivity_results_list <- mclapply(1:nrow(param_grid), function(i) {
#   param <- as.character(param_grid$Parameter[i])
#   inc <- param_grid$Increment[i]
#   run_model_with_modification(param, inc)
# }, mc.cores = numCores)
# 
# numCores
# # Combine the list into a data frame
# sensitivity_results <- do.call(rbind, sensitivity_results_list)
# 
# # Save the data frame as a CSV file
# write.csv(sensitivity_results, "Multispp_sensitivity_vitalrates.csv", row.names = FALSE)
# 
# 

# Create a grid of all temperature factors and parameter modifications
scenario_grid <- expand.grid(temp_factor = temp_seq, param_idx = 1:nrow(param_grid))

run_scenario <- function(idx) {
  temp_factor <- scenario_grid$temp_factor[idx]
  param_idx <- scenario_grid$param_idx[idx]  # Ensure single index extraction

  # Modify temperature dataset
  temps_modified <- temps
  temps_modified$Temperature <- temps_modified$Temperature * temp_factor

  # Extract parameter modification details (Ensure single values!)
  param <- as.character(param_grid$Parameter[[param_idx]])
  inc <- param_grid$Increment[[param_idx]]

  # Run the model
  result <- run_model_with_modification(param, inc)

  # Add metadata to results
  result$TemperatureFactor <- temp_factor
  result$Parameter <- param
  result$Increment <- inc
  return(result)
}

# Run all scenarios in parallel
library(parallel)
numCores <- detectCores() - 1  # Use available cores minus 1 for safety
all_results <- mclapply(1:nrow(scenario_grid), run_scenario, mc.cores = numCores)
#all_results <- mclapply(1:3, run_scenario, mc.cores = 1)

# Combine results into a single dataframe
final_results <- do.call(rbind, all_results)

# Save to CSV
write.csv(final_results, "Multispp_mT_sens_temps.csv", row.names = FALSE)

