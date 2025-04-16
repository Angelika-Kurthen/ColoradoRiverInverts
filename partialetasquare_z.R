# ## Multispecies Partial Eta Square
library(data.table)
library(purrr)
library(effectsize)
library(dplyr)
library(parallel)

# List of input files and their metadata
files <- list(
  "Multispp_temp_biomass_z_hyos_full.csv" = list(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  "Multispp_temp_biomass_spike_z_hyos_full.csv" = list(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  "Multispp_temp_biomass_hyd_z_hyos_full.csv" = list(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  "Multispp_temp_hyd_biomass_spike_z_hyos_full.csv" = list(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  "Multispp_temp_biomass_hfe_z_hyos_full.csv" = list(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
  "Multispp_temp_HFE_biomass_spike_z_hyos_full.csv" = list(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
  "Multispp_temp_hyd_HFE_biomass_z_hyos_full.csv" = list(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
  "Multispp_temp_hyd_HFE_biomass_spike_z_hyos_full.csv" = list(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
)

# List of taxa
taxa_list <- c("BAET", "HYOS", "CHIR", "GAMM", "NZMS")

# Loop through each taxon separately
walk(taxa_list, function(taxon) {
  message("Processing taxon: ", taxon)
  
  # Combine results across files for this taxon
  eta_df <- map_df(names(files), function(fname) {
    meta <- files[[fname]]
    
    # Read in just this taxon's data, avoid memory mapping
    dt <- tryCatch({
      fread(fname, select = c("taxa", "biomass", "temperature", "q"))[taxa == taxon]
    }, error = function(e) {
      message("Error reading ", fname, ": ", conditionMessage(e))
      return(NULL)
    })
    
    if (is.null(dt) || nrow(dt) == 0) return(NULL)
    
    dt[, `:=`(
      biomass = as.numeric(biomass),
      temperature = as.factor(temperature),
      q = as.factor(q),
      source = as.factor(meta$source),
      temp_spike = as.factor(meta$temp_spike),
      hydropeaking = as.factor(meta$hydropeaking),
      HFE = as.factor(meta$HFE)
    )]
    
    # Fit simplified model
    model <-  aov(biomass ~ (temperature + temp_spike + hydropeaking + HFE + q)^2, data = dt)

    
    # Extract partial eta squared
    eta <- eta_squared(model, partial = TRUE) %>%
      as.data.frame() %>%
      arrange(desc(Eta2_partial))
    
    eta$scenario <- meta$source
    eta$temp_spike <- meta$temp_spike
    eta$hydropeaking <- meta$hydropeaking
    eta$HFE <- meta$HFE
    
    return(eta)
  })
  
  # Save result per taxon
  if (nrow(eta_df) > 0) {
    eta_df$taxon <- taxon
    outname <- paste0("partialetasq_biomass_hyos_", tolower(taxon), ".csv")
    fwrite(eta_df, outname)
    message("Saved: ", outname)
  } else {
    message("No data processed for taxon: ", taxon)
  }
})


#baetidae

# # List of input files and their metadata
# files <- list(
#   "Multispp_temp_biomass_z_baet_full.csv" = list(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
#   "Multispp_temp_biomass_spike_z_baet_full.csv" = list(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
#   "Multispp_temp_biomass_hyd_z_hyos_full.csv" = list(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
#   "Multispp_temp_hyd_biomass_spike_z_baet_full.csv" = list(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
#   "Multispp_temp_biomass_hfe_z_baet_full.csv" = list(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
#   "Multispp_temp_HFE_biomass_spike_z_baet_full.csv" = list(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
#   "Multispp_temp_hyd_HFE_biomass_z_baet_full.csv" = list(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
#   "Multispp_temp_hyd_HFE_biomass_spike_z_baet_full.csv" = list(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
# )
# 
# # List of taxa
# taxa_list <- c("BAET", "HYOS", "CHIR", "GAMM", "NZMS")
# 
# # Run per-taxon
# all_eta <- map_df(taxa_list, function(taxon) {
#   # Read and filter just rows for this taxon
#   data_list <- map(names(files), function(fname) {
#     meta <- files[[fname]]
#     dt <- fread(fname)[taxa == taxon]
#     dt[, `:=`(
#       biomass = biomass,  # Rename for clarity
#       temperature = as.factor(temperature),
#       q = as.factor(q),
#       source = meta$source,
#       temp_spike = meta$temp_spike,
#       hydropeaking = meta$hydropeaking,
#       HFE = meta$HFE
#     )]
#     return(dt)
#   })
#   
#   # Combine all rows for this taxon
#   taxon_df <- rbindlist(data_list)
#   
#   # Fit model and compute partial eta squared
#   model <- aov(biomass ~ (temperature + temp_spike + hydropeaking + HFE + q)^2, data = taxon_df)
#   eta_sq_df <- eta_squared(model, partial = TRUE) %>%
#     as.data.frame() %>%
#     arrange(desc(Eta2_partial))
#   
#   eta_sq_df$taxon <- taxon
#   return(eta_sq_df)
# })
# 
# # Save result
# fwrite(all_eta, "Multispp_partialetasquare_biomass_z_baet.csv")
# rm(data_list, model, eta_sq_df)
# 
# # chironomidae
# 
# # List of input files and their metadata
# files <- list(
#   "Multispp_temp_biomass_z_chir_full.csv" = list(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
#   "Multispp_temp_biomass_spike_z_chir_full.csv" = list(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
#   "Multispp_temp_biomass_hyd_z_chir_full.csv" = list(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
#   "Multispp_temp_hyd_biomass_spike_z_chir_full.csv" = list(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
#   "Multispp_temp_biomass_hfe_z_chir_full.csv" = list(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
#   "Multispp_temp_HFE_biomass_spike_z_chir_full.csv" = list(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
#   "Multispp_temp_hyd_HFE_biomass_z_chir_full.csv" = list(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
#   "Multispp_temp_hyd_HFE_biomass_spike_z_chir_full.csv" = list(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
# )
# 
# # List of taxa
# taxa_list <- c("BAET", "HYOS", "CHIR", "GAMM", "NZMS")
# 
# # Run per-taxon
# all_eta <- map_df(taxa_list, function(taxon) {
#   # Read and filter just rows for this taxon
#   data_list <- map(names(files), function(fname) {
#     meta <- files[[fname]]
#     dt <- fread(fname)[taxa == taxon]
#     dt[, `:=`(
#       biomass = biomass,  # Rename for clarity
#       temperature = as.factor(temperature),
#       q = as.factor(q),
#       source = meta$source,
#       temp_spike = meta$temp_spike,
#       hydropeaking = meta$hydropeaking,
#       HFE = meta$HFE
#     )]
#     return(dt)
#   })
#   
#   # Combine all rows for this taxon
#   taxon_df <- rbindlist(data_list)
#   
#   # Fit model and compute partial eta squared
#   model <- aov(biomass ~ (temperature + temp_spike + hydropeaking + HFE + q)^2, data = taxon_df)
#   eta_sq_df <- eta_squared(model, partial = TRUE) %>%
#     as.data.frame() %>%
#     arrange(desc(Eta2_partial))
#   
#   eta_sq_df$taxon <- taxon
#   return(eta_sq_df)
# })
# 
# # Save result
# fwrite(all_eta, "Multispp_partialetasquare_biomass_z_chir.csv")
# rm(data_list, model, eta_sq_df)
# 
# #gammarus
# 
# # List of input files and their metadata
# files <- list(
#   "Multispp_temp_biomass_z_gamm_full.csv" = list(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
#   "Multispp_temp_biomass_spike_z_gamm_full.csv" = list(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
#   "Multispp_temp_biomass_hyd_z_gamm_full.csv" = list(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
#   "Multispp_temp_hyd_biomass_spike_z_gamm_full.csv" = list(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
#   "Multispp_temp_biomass_hfe_z_gamm_full.csv" = list(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
#   "Multispp_temp_HFE_biomass_spike_z_gamm_full.csv" = list(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
#   "Multispp_temp_hyd_HFE_biomass_z_gamm_full.csv" = list(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
#   "Multispp_temp_hyd_HFE_biomass_spike_z_gamm_full.csv" = list(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
# )
# 
# # List of taxa
# taxa_list <- c("BAET", "HYOS", "CHIR", "GAMM", "NZMS")
# 
# # Run per-taxon
# all_eta <- map_df(taxa_list, function(taxon) {
#   # Read and filter just rows for this taxon
#   data_list <- map(names(files), function(fname) {
#     meta <- files[[fname]]
#     dt <- fread(fname)[taxa == taxon]
#     dt[, `:=`(
#       biomass = biomass,  # Rename for clarity
#       temperature = as.factor(temperature),
#       q = as.factor(q),
#       source = meta$source,
#       temp_spike = meta$temp_spike,
#       hydropeaking = meta$hydropeaking,
#       HFE = meta$HFE
#     )]
#     return(dt)
#   })
#   
#   # Combine all rows for this taxon
#   taxon_df <- rbindlist(data_list)
#   
#   # Fit model and compute partial eta squared
#   model <- aov(biomass ~ (temperature + temp_spike + hydropeaking + HFE + q)^2, data = taxon_df)
#   eta_sq_df <- eta_squared(model, partial = TRUE) %>%
#     as.data.frame() %>%
#     arrange(desc(Eta2_partial))
#   
#   eta_sq_df$taxon <- taxon
#   return(eta_sq_df)
# })
# 
# # Save result
# fwrite(all_eta, "Multispp_partialetasquare_biomass_z_gamm.csv")
# rm(data_list, model, eta_sq_df)
# 
# 
# library(data.table)
# library(purrr)
# library(effectsize)
# library(dplyr)
# library(parallel)
# 
# # List of input files and their metadata
# files <- list(
#   "Multispp_temp_biomass_z_nzms_full.csv" = list(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
#   "Multispp_temp_biomass_spike_z_nzms_full.csv" = list(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
#   "Multispp_temp_biomass_hyd_z_nzms_full.csv" = list(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
#   "Multispp_temp_hyd_biomass_spike_z_nzms_full.csv" = list(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
#   "Multispp_temp_biomass_hfe_z_nzms_full.csv" = list(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
#   "Multispp_temp_HFE_biomass_spike_z_nzms_full.csv" = list(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
#   "Multispp_temp_hyd_HFE_biomass_z_nzms_full.csv" = list(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
#   "Multispp_temp_hyd_HFE_biomass_spike_z_nzms_full.csv" = list(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
# )
# 
# # List of taxa
# taxa_list <- c("BAET", "HYOS", "CHIR", "GAMM", "NZMS")
# 
# # Run per-taxon
# all_eta <- map_df(taxa_list, function(taxon) {
#   # Read and filter just rows for this taxon
#   data_list <- map(names(files), function(fname) {
#     meta <- files[[fname]]
#     dt <- fread(fname)[taxa == taxon]
#     dt[, `:=`(
#       biomass = biomass,  # Rename for clarity
#       temperature = as.factor(temperature),
#       q = as.factor(q),
#       source = meta$source,
#       temp_spike = meta$temp_spike,
#       hydropeaking = meta$hydropeaking,
#       HFE = meta$HFE
#     )]
#     return(dt)
#   })
#   
#   # Combine all rows for this taxon
#   taxon_df <- rbindlist(data_list)
#   
#   # Fit model and compute partial eta squared
#   model <- aov(biomass ~ (temperature + temp_spike + hydropeaking + HFE + q)^2, data = taxon_df)
#   eta_sq_df <- eta_squared(model, partial = TRUE) %>%
#     as.data.frame() %>%
#     arrange(desc(Eta2_partial))
#   
#   eta_sq_df$taxon <- taxon
#   return(eta_sq_df)
# })
# 
# # Save result
# fwrite(all_eta, "Multispp_partialetasquare_biomass_z_nzms.csv")
# rm(data_list, model, eta_sq_df)
# 
# #all 
# 
# # List of input files and their metadata
# files <- list(
#   "Multispp_temp_biomass_z_full.csv" = list(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
#    "Multispp_temp_biomass_spike_z_full.csv" = list(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
#   "Multispp_temp_biomass_hyd_z_full.csv" = list(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
#   "Multispp_temp_hyd_biomass_spike_z_full.csv" = list(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
#   "Multispp_temp_HFE_biomass_spike_z_full.csv" = list(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
#   "Multispp_temp_hyd_HFE_biomass_z_full.csv" = list(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
#   "Multispp_temp_hyd_HFE_biomass_spike_z_full.csv" = list(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
# )
# 
# # List of taxa
# taxa_list <- c("BAET", "HYOS", "CHIR", "GAMM", "NZMS")
# 
# # Run per-taxon
# all_eta <- map_df(taxa_list, function(taxon) {
#   # Read and filter just rows for this taxon
#   data_list <- map(names(files), function(fname) {
#     meta <- files[[fname]]
#     dt <- fread(fname)[taxa == taxon]
#     dt[, `:=`(
#       biomass = biomass,  # Rename for clarity
#       temperature = as.factor(temperature),
#       q = as.factor(q),
#       source = meta$source,
#       temp_spike = meta$temp_spike,
#       hydropeaking = meta$hydropeaking,
#       HFE = meta$HFE
#     )]
#     return(dt)
#   })
#   
#   # Combine all rows for this taxon
#   taxon_df <- rbindlist(data_list)
#   
#   # Fit model and compute partial eta squared
#   model <- aov(biomass ~ (temperature + temp_spike + hydropeaking + HFE + q)^2, data = taxon_df)
#   eta_sq_df <- eta_squared(model, partial = TRUE) %>%
#     as.data.frame() %>%
#     arrange(desc(Eta2_partial))
#   
#   eta_sq_df$taxon <- taxon
#   return(eta_sq_df)
# })
# 
# # Save result
# fwrite(all_eta, "Multispp_partialetasquare_biomass_z_hyos.csv")
# rm(data_list, model, eta_sq_df)
# 


# 
# library(data.table)
# library(purrr)
# library(effectsize)
# library(dplyr)
# library(parallel)
# 
# # List of input files and their metadata
# files <- list(
#   "Multispp_temp_biomass_z_hyos_full.csv" = list(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
#   "Multispp_temp_biomass_spike_z_hyos_full.csv" = list(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
#   "Multispp_temp_biomass_hyd_z_hyos_full.csv" = list(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
#   "Multispp_temp_biomass_spike_hyd_z_hyos_full.csv" = list(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
#   "Multispp_temp_biomass_hfe_z_hyos_full.csv" = list(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
#   "Multispp_temp_HFE_biomass_spike_z_hyos_full.csv" = list(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
#   "Multispp_temp_hyd_HFE_biomass_z_hyos_full.csv" = list(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
#   "Multispp_temp_hyd_HFE_biomass_spike_z_hyos_full.csv" = list(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
# )
# 
# # List of taxa
# taxa_list <- c("BAET", "HYOS", "CHIR", "GAMM", "NZMS")
# 
# # Run per-taxon
# all_eta <- map_df(taxa_list, function(taxon) {
#   # Read and filter just rows for this taxon
#   data_list <- map(names(files), function(fname) {
#     meta <- files[[fname]]
#     dt <- fread(fname)[taxa == taxon]
#     dt[, `:=`(
#       source = meta$source,
#       temp_spike = meta$temp_spike,
#       hydropeaking = meta$hydropeaking,
#       HFE = meta$HFE,
#       q = as.factor(meta$q)
#     )]
#     return(dt)
#   })
# 
#   # Combine all rows for this taxon
#   taxon_df <- rbindlist(data_list)
# 
#   # Fit model and compute partial eta squared
#   model <- aov(biomass ~ (temperature + temp_spike + hydropeaking + HFE + q)^5, data = taxon_df)
#   eta_sq_df <- eta_squared(model, partial = TRUE) %>%
#     as.data.frame() %>%
#     arrange(desc(Eta2_partial))
# 
#   eta_sq_df$taxon <- taxon
#   return(eta_sq_df)
# })
# 
# # Save result
# fwrite(all_eta, "Multispp_partialetasquare_biomass_z_hyos.csv")
# rm(data_list, model, eta_sq_df)
# 
# 
# 
# # List of input files and their metadata
# files <- list(
#   "Multispp_temp_biomass_z_baet_full.csv" = list(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
#   "Multispp_temp_biomass_spike_z_baet_full.csv" = list(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
#   "Multispp_temp_biomass_hyd_z_baet_full.csv" = list(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
#   "Multispp_temp_hyd_biomass_spike_z_nzms_full.csv" = list(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
#   "Multispp_temp_biomass_hfe_z_baet_full.csv" = list(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
#   "Multispp_temp_HFE_biomass_spike_z_baet_full.csv" = list(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
#   "Multispp_temp_hyd_HFE_biomass_z_baet_full.csv" = list(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
#   "Multispp_temp_hyd_HFE_biomass_spike_z_baet_full.csv" = list(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
# )
# 
# # List of taxa
# taxa_list <- c("BAET", "HYOS", "CHIR", "GAMM", "NZMS")
# 
# # Run per-taxon
# all_eta <- map_df(taxa_list, function(taxon) {
#   # Read and filter just rows for this taxon
#   data_list <- map(names(files), function(fname) {
#     meta <- files[[fname]]
#     dt <- fread(fname)[taxa == taxon]
#     dt[, `:=`(
#       source = meta$source,
#       temp_spike = meta$temp_spike,
#       hydropeaking = meta$hydropeaking,
#       HFE = meta$HFE,
#       q = as.factor(meta$q)
#     )]
#     return(dt)
#   })
#   
#   # Combine all rows for this taxon
#   taxon_df <- rbindlist(data_list)
#   
#   # Fit model and compute partial eta squared
#   model <- aov(biomass ~ (temperature + temp_spike + hydropeaking + HFE + q)^5, data = taxon_df)
#   eta_sq_df <- eta_squared(model, partial = TRUE) %>%
#     as.data.frame() %>%
#     arrange(desc(Eta2_partial))
#   
#   eta_sq_df$taxon <- taxon
#   return(eta_sq_df)
# })
# 
# # Save result
# fwrite(all_eta, "Multispp_partialetasquare_biomass_z_baet.csv")
# rm(data_list, model, eta_sq_df)
# 
# 
# 
# # List of input files and their metadata
# files <- list(
#   "Multispp_temp_biomass_z_chir_full.csv" = list(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
#   "Multispp_temp_biomass_spike_z_chir_full.csv" = list(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
#   "Multispp_temp_biomass_hyd_z_chir_full.csv" = list(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
#   "Multispp_temp_hyd_biomass_spike_z_nzms_full.csv" = list(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
#   "Multispp_temp_biomass_hfe_z_chir_full.csv" = list(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
#   "Multispp_temp_HFE_biomass_spike_z_chir_full.csv" = list(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
#   "Multispp_temp_hyd_HFE_biomass_z_chir_full.csv" = list(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
#   "Multispp_temp_hyd_HFE_biomass_spike_z_chir_full.csv" = list(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
# )
# 
# # List of taxa
# taxa_list <- c("BAET", "HYOS", "CHIR", "GAMM", "NZMS")
# 
# # Run per-taxon
# all_eta <- map_df(taxa_list, function(taxon) {
#   # Read and filter just rows for this taxon
#   data_list <- map(names(files), function(fname) {
#     meta <- files[[fname]]
#     dt <- fread(fname)[taxa == taxon]
#     dt[, `:=`(
#       source = meta$source,
#       temp_spike = meta$temp_spike,
#       hydropeaking = meta$hydropeaking,
#       HFE = meta$HFE,
#       q = as.factor(meta$q)
#     )]
#     return(dt)
#   })
#   
#   # Combine all rows for this taxon
#   taxon_df <- rbindlist(data_list)
#   
#   # Fit model and compute partial eta squared
#   model <- aov(biomass ~ (temperature + temp_spike + hydropeaking + HFE + q)^5, data = taxon_df)
#   eta_sq_df <- eta_squared(model, partial = TRUE) %>%
#     as.data.frame() %>%
#     arrange(desc(Eta2_partial))
#   
#   eta_sq_df$taxon <- taxon
#   return(eta_sq_df)
# })
# 
# # Save result
# fwrite(all_eta, "Multispp_partialetasquare_biomass_z_chir.csv")
# rm(data_list, model, eta_sq_df)
# gc()
# 
# # List of input files and their metadata
# files <- list(
#   "Multispp_temp_biomass_z_gamm_full.csv" = list(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
#   "Multispp_temp_biomass_spike_z_gamm_full.csv" = list(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
#   "Multispp_temp_biomass_hyd_z_gamm_full.csv" = list(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
#   "Multispp_temp_hyd_biomass_spike_z_nzms_full.csv" = list(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
#   "Multispp_temp_biomass_hfe_z_gamm_full.csv" = list(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
#   "Multispp_temp_HFE_biomass_spike_z_gamm_full.csv" = list(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
#   "Multispp_temp_hyd_HFE_biomass_z_gamm_full.csv" = list(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
#   "Multispp_temp_hyd_HFE_biomass_spike_z_gamm_full.csv" = list(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
# )
# 
# # List of taxa
# taxa_list <- c("BAET", "HYOS", "CHIR", "GAMM", "NZMS")
# 
# # Run per-taxon
# all_eta <- map_df(taxa_list, function(taxon) {
#   # Read and filter just rows for this taxon
#   data_list <- map(names(files), function(fname) {
#     meta <- files[[fname]]
#     dt <- fread(fname)[taxa == taxon]
#     dt[, `:=`(
#       source = meta$source,
#       temp_spike = meta$temp_spike,
#       hydropeaking = meta$hydropeaking,
#       HFE = meta$HFE,
#       q = as.factor(meta$q)
#     )]
#     return(dt)
#   })
#   
#   # Combine all rows for this taxon
#   taxon_df <- rbindlist(data_list)
#   
#   # Fit model and compute partial eta squared
#   model <- aov(biomass ~ (temperature + temp_spike + hydropeaking + HFE + q)^5, data = taxon_df)
#   eta_sq_df <- eta_squared(model, partial = TRUE) %>%
#     as.data.frame() %>%
#     arrange(desc(Eta2_partial))
#   
#   eta_sq_df$taxon <- taxon
#   return(eta_sq_df)
# })
# 
# # Save result
# fwrite(all_eta, "Multispp_partialetasquare_biomass_z_gamm.csv")
# rm(data_list, model, eta_sq_df)
# gc()
# 
# # List of input files and their metadata
# files <- list(
#   "Multispp_temp_biomass_z_nzms_full.csv" = list(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
#   "Multispp_temp_biomass_spike_z_nzms_full.csv" = list(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
#   "Multispp_temp_biomass_hyd_z_nzms_full.csv" = list(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
#   "Multispp_temp_hyd_biomass_spike_z_nzms_full.csv" = list(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
#   "Multispp_temp_biomass_hfe_z_nzms_full.csv" = list(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
#   "Multispp_temp_HFE_biomass_spike_z_nzms_full.csv" = list(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
#   "Multispp_temp_hyd_HFE_biomass_z_nzms_full.csv" = list(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
#   "Multispp_temp_hyd_HFE_biomass_spike_z_nzms_full.csv" = list(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
# )
# 
# # List of taxa
# taxa_list <- c("BAET", "HYOS", "CHIR", "GAMM", "NZMS")
# 
# # Run per-taxon
# all_eta <- map_df(taxa_list, function(taxon) {
#   # Read and filter just rows for this taxon
#   data_list <- map(names(files), function(fname) {
#     meta <- files[[fname]]
#     dt <- fread(fname)[taxa == taxon]
#     dt[, `:=`(
#       source = meta$source,
#       temp_spike = meta$temp_spike,
#       hydropeaking = meta$hydropeaking,
#       HFE = meta$HFE,
#       q = as.factor(meta$q)
#     )]
#     return(dt)
#   })
#   
#   # Combine all rows for this taxon
#   taxon_df <- rbindlist(data_list)
#   
#   # Fit model and compute partial eta squared
#   model <- aov(biomass ~ (temperature + temp_spike + hydropeaking + HFE + q)^5, data = taxon_df)
#   eta_sq_df <- eta_squared(model, partial = TRUE) %>%
#     as.data.frame() %>%
#     arrange(desc(Eta2_partial))
#   
#   eta_sq_df$taxon <- taxon
#   return(eta_sq_df)
# })
# 
# # Save result
# fwrite(all_eta, "Multispp_partialetasquare_biomass_z_nzms.csv")
# rm(data_list, model, eta_sq_df)
# gc()
# 
# 
