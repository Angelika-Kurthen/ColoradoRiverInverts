# ## Multispecies Partial Eta Square
library(data.table)
library(purrr)
library(effectsize)
library(dplyr)
library(parallel)

# List of input files and their metadata
hyos <- bind_rows(
  "Multispp_temp_biomass_z_hyos_full_HYOS.csv",
  "Multispp_temp_biomass_spike_z_hyos_full_HYOS.csv", 
  "Multispp_temp_biomass_hyd_z_hyos_full_HYOS.csv",
  "Multispp_temp_hyd_biomass_spike_z_hyos_full_HYOS.csv",
  "Multispp_temp_biomass_hfe_z_hyos_full_HYOS.csv",
  "Multispp_temp_HFE_biomass_spike_z_hyos_full_HYOS.csv", 
  "Multispp_temp_hyd_HFE_biomass_spike_z_hyos_full_HYOS.csv", 
)

label.names <- c(
  "Temperature", 
  "Temperature +\nSpike", 
  "Temperature +\nHydropeaking",
  "Temperature +\nSpike + Hydropeaking",
  "Temperature +\nHFE",
  "Temperature +\nSpike + HFE",
  "Temperature +\nHydropeaking + HFE",
  "Temperature +\nSpike + Hydropeaking + HFE"
)

hyos$source <- ordered(hyos$source, levels = c(
  "Temperature", 
  "Temperature & spike",
  "Temperature & hydropeaking",
  "Temperature & spike & hydropeaking",
  "Temperature & HFE",
  "Temperature & spike & HFE",
  "Temperature & hydropeaking & HFE",
  "Temperature & spike & hydropeaking & HFE"
))  

# Function to fit ANOVA and extract partial eta squared
compute_eta_sq <- function(df, response_var) {
  model <- aov(as.formula(paste(response_var, "~ (temperature + temp_spike + hydropeaking + HFE + q)^2")), data = df)
  eta_squared(model, partial = TRUE) %>%
    as.data.frame()
}

# Run analysis per taxon for abundance
eta_sq_taxon_abundance <- hyos %>%
  map_df(~ compute_eta_sq(.x, "biomass"))
write.csv(eta_sq_df.multi.b, "Multispp_partialetasquare_biomass_hyos_hyos.csv", row.names = FALSE)
rm(hyos)


# List of input files and their metadata
baet <- bind_rows(
  "Multispp_temp_biomass_z_hyos_full_BAET.csv",
  "Multispp_temp_biomass_spike_z_hyos_full_BAET.csv", 
  "Multispp_temp_biomass_hyd_z_hyos_full_BAET.csv",
  "Multispp_temp_hyd_biomass_spike_z_hyos_full_BAET.csv",
  "Multispp_temp_biomass_hfe_z_hyos_full_BAET.csv",
  "Multispp_temp_HFE_biomass_spike_z_hyos_full_BAET.csv", 
  "Multispp_temp_hyd_HFE_biomass_spike_z_hyos_full_BAET.csv", 
)

label.names <- c(
  "Temperature", 
  "Temperature +\nSpike", 
  "Temperature +\nHydropeaking",
  "Temperature +\nSpike + Hydropeaking",
  "Temperature +\nHFE",
  "Temperature +\nSpike + HFE",
  "Temperature +\nHydropeaking + HFE",
  "Temperature +\nSpike + Hydropeaking + HFE"
)

baet$source <- ordered(baet$source, levels = c(
  "Temperature", 
  "Temperature & spike",
  "Temperature & hydropeaking",
  "Temperature & spike & hydropeaking",
  "Temperature & HFE",
  "Temperature & spike & HFE",
  "Temperature & hydropeaking & HFE",
  "Temperature & spike & hydropeaking & HFE"
))  

# Function to fit ANOVA and extract partial eta squared
compute_eta_sq <- function(df, response_var) {
  model <- aov(as.formula(paste(response_var, "~ (temperature + temp_spike + hydropeaking + HFE + q)^2")), data = df)
  eta_squared(model, partial = TRUE) %>%
    as.data.frame()
}

# Run analysis per taxon for abundance
eta_sq_taxon_abundance <- baet %>%
  map_df(~ compute_eta_sq(.x, "biomass"))
write.csv(eta_sq_df.multi.b, "Multispp_partialetasquare_biomass_hyos_baet.csv", row.names = FALSE)
rm(baet)


# List of input files and their metadata
chir <- bind_rows(
  "Multispp_temp_biomass_z_hyos_full_CHIR.csv",
  "Multispp_temp_biomass_spike_z_hyos_full_CHIR.csv", 
  "Multispp_temp_biomass_hyd_z_hyos_full_CHIR.csv",
  "Multispp_temp_hyd_biomass_spike_z_hyos_full_CHIR.csv",
  "Multispp_temp_biomass_hfe_z_hyos_full_CHIR.csv",
  "Multispp_temp_HFE_biomass_spike_z_hyos_full_CHIR.csv", 
  "Multispp_temp_hyd_HFE_biomass_spike_z_hyos_full_CHIR.csv", 
)

label.names <- c(
  "Temperature", 
  "Temperature +\nSpike", 
  "Temperature +\nHydropeaking",
  "Temperature +\nSpike + Hydropeaking",
  "Temperature +\nHFE",
  "Temperature +\nSpike + HFE",
  "Temperature +\nHydropeaking + HFE",
  "Temperature +\nSpike + Hydropeaking + HFE"
)

chir$source <- ordered(chir$source, levels = c(
  "Temperature", 
  "Temperature & spike",
  "Temperature & hydropeaking",
  "Temperature & spike & hydropeaking",
  "Temperature & HFE",
  "Temperature & spike & HFE",
  "Temperature & hydropeaking & HFE",
  "Temperature & spike & hydropeaking & HFE"
))  

# Function to fit ANOVA and extract partial eta squared
compute_eta_sq <- function(df, response_var) {
  model <- aov(as.formula(paste(response_var, "~ (temperature + temp_spike + hydropeaking + HFE + q)^2")), data = df)
  eta_squared(model, partial = TRUE) %>%
    as.data.frame()
}

# Run analysis per taxon for abundance
eta_sq_taxon_abundance <- chir %>%
  map_df(~ compute_eta_sq(.x, "biomass"))
write.csv(eta_sq_df.multi.b, "Multispp_partialetasquare_biomass_hyos_chir.csv", row.names = FALSE)
rm(chir)


# List of input files and their metadata
gamm <- bind_rows(
  "Multispp_temp_biomass_z_hyos_full_GAMM.csv",
  "Multispp_temp_biomass_spike_z_hyos_full_GAMM.csv", 
  "Multispp_temp_biomass_hyd_z_hyos_full_GAMM.csv",
  "Multispp_temp_hyd_biomass_spike_z_hyos_full_GAMM.csv",
  "Multispp_temp_biomass_hfe_z_hyos_full_GAMM.csv",
  "Multispp_temp_HFE_biomass_spike_z_hyos_full_GAMM.csv", 
  "Multispp_temp_hyd_HFE_biomass_spike_z_hyos_full_GAMM.csv", 
)

label.names <- c(
  "Temperature", 
  "Temperature +\nSpike", 
  "Temperature +\nHydropeaking",
  "Temperature +\nSpike + Hydropeaking",
  "Temperature +\nHFE",
  "Temperature +\nSpike + HFE",
  "Temperature +\nHydropeaking + HFE",
  "Temperature +\nSpike + Hydropeaking + HFE"
)

gamm$source <- ordered(gamm$source, levels = c(
  "Temperature", 
  "Temperature & spike",
  "Temperature & hydropeaking",
  "Temperature & spike & hydropeaking",
  "Temperature & HFE",
  "Temperature & spike & HFE",
  "Temperature & hydropeaking & HFE",
  "Temperature & spike & hydropeaking & HFE"
))  

# Function to fit ANOVA and extract partial eta squared
compute_eta_sq <- function(df, response_var) {
  model <- aov(as.formula(paste(response_var, "~ (temperature + temp_spike + hydropeaking + HFE + q)^2")), data = df)
  eta_squared(model, partial = TRUE) %>%
    as.data.frame()
}

# Run analysis per taxon for abundance
eta_sq_taxon_abundance <- gamm %>%
  map_df(~ compute_eta_sq(.x, "biomass"))
write.csv(eta_sq_df.multi.b, "Multispp_partialetasquare_biomass_hyos_gamm.csv", row.names = FALSE)
rm(gamm)


# List of input files and their metadata
nzms <- bind_rows(
  "Multispp_temp_biomass_z_hyos_full_NZMS.csv",
  "Multispp_temp_biomass_spike_z_hyos_full_NZMS.csv", 
  "Multispp_temp_biomass_hyd_z_hyos_full_NZMS.csv",
  "Multispp_temp_hyd_biomass_spike_z_hyos_full_NZMS.csv",
  "Multispp_temp_biomass_hfe_z_hyos_full_NZMS.csv",
  "Multispp_temp_HFE_biomass_spike_z_hyos_full_NZMS.csv", 
  "Multispp_temp_hyd_HFE_biomass_spike_z_hyos_full_NZMS.csv", 
)

label.names <- c(
  "Temperature", 
  "Temperature +\nSpike", 
  "Temperature +\nHydropeaking",
  "Temperature +\nSpike + Hydropeaking",
  "Temperature +\nHFE",
  "Temperature +\nSpike + HFE",
  "Temperature +\nHydropeaking + HFE",
  "Temperature +\nSpike + Hydropeaking + HFE"
)

nzms$source <- ordered(nzms$source, levels = c(
  "Temperature", 
  "Temperature & spike",
  "Temperature & hydropeaking",
  "Temperature & spike & hydropeaking",
  "Temperature & HFE",
  "Temperature & spike & HFE",
  "Temperature & hydropeaking & HFE",
  "Temperature & spike & hydropeaking & HFE"
))  

# Function to fit ANOVA and extract partial eta squared
compute_eta_sq <- function(df, response_var) {
  model <- aov(as.formula(paste(response_var, "~ (temperature + temp_spike + hydropeaking + HFE + q)^2")), data = df)
  eta_squared(model, partial = TRUE) %>%
    as.data.frame()
}

# Run analysis per taxon for abundance
eta_sq_taxon_abundance <- nzms %>%
  map_df(~ compute_eta_sq(.x, "biomass"))
write.csv(eta_sq_df.multi.b, "Multispp_partialetasquare_biomass_hyos_nzms.csv", row.names = FALSE)
rm(nzms)