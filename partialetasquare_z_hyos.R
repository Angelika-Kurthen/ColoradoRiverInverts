# ## Multispecies Partial Eta Square
library(data.table)
library(purrr)
library(effectsize)
library(dplyr)
library(parallel)

# List of input Hydropsyche files when zHydropsyche is changed
biomass <- read.csv2("Multispp_temp_biomass_z_hyos_full_HYOS.csv")
biomass_spike <- read.csv2("Multispp_temp_biomass_spike_z_hyos_full_HYOS.csv") 
hyd <- read.csv2("Multispp_temp_biomass_hyd_z_hyos_full_HYOS.csv")
hyd_spike <- read.csv2("Multispp_temp_hyd_biomass_spike_z_hyos_full_HYOS.csv")
hfe <- read.csv2("Multispp_temp_biomass_hfe_z_hyos_full_HYOS.csv")
hfe_spike <- read.csv2("Multispp_temp_HFE_biomass_spike_z_hyos_full_HYOS.csv")
hyd_hfe <- read.csv2("Multispp_temp_hyd_HFE_biomass_z_hyos_full_HYOS.csv")
hyd_hfe_spike <- read.csv2("Multispp_temp_hyd_HFE_biomass_spike_z_hyos_full_HYOS.csv")


hyos <- bind_rows(
  biomass, 
  biomass_spike, 
  hyd, 
  hyd_spike, 
  hfe, 
  hfe_spike, 
  hyd_hfe, 
  hyd_hfe_spike
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


# List of input Baetidae files when zHydropsyche is changed
biomass <- read.csv2("Multispp_temp_biomass_z_hyos_full_BAET.csv")
biomass_spike <- read.csv2("Multispp_temp_biomass_spike_z_hyos_full_BAET.csv") 
hyd <- read.csv2("Multispp_temp_biomass_hyd_z_hyos_full_BAET.csv")
hyd_spike <- read.csv2("Multispp_temp_hyd_biomass_spike_z_hyos_full_BAET.csv")
hfe <- read.csv2("Multispp_temp_biomass_hfe_z_hyos_full_BAET.csv")
hfe_spike <- read.csv2("Multispp_temp_HFE_biomass_spike_z_hyos_full_BAET.csv")
hyd_hfe <- read.csv2("Multispp_temp_hyd_HFE_biomass_z_hyos_full_BAET.csv")
hyd_hfe_spike <- read.csv2("Multispp_temp_hyd_HFE_biomass_spike_z_hyos_full_BAET.csv")


baet <- bind_rows(
  biomass, 
  biomass_spike, 
  hyd, 
  hyd_spike, 
  hfe, 
  hfe_spike, 
  hyd_hfe, 
  hyd_hfe_spike
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


# List of input Chironomidae files when zHydropsyche is changed
biomass <- read.csv2("Multispp_temp_biomass_z_hyos_full_CHIR.csv")
biomass_spike <- read.csv2("Multispp_temp_biomass_spike_z_hyos_full_CHIR.csv") 
hyd <- read.csv2("Multispp_temp_biomass_hyd_z_hyos_full_CHIR.csv")
hyd_spike <- read.csv2("Multispp_temp_hyd_biomass_spike_z_hyos_full_CHIR.csv")
hfe <- read.csv2("Multispp_temp_biomass_hfe_z_hyos_full_CHIR.csv")
hfe_spike <- read.csv2("Multispp_temp_HFE_biomass_spike_z_hyos_full_CHIR.csv")
hyd_hfe <- read.csv2("Multispp_temp_hyd_HFE_biomass_z_hyos_full_CHIR.csv")
hyd_hfe_spike <- read.csv2("Multispp_temp_hyd_HFE_biomass_spike_z_hyos_full_HYOS.csv")


chir <- bind_rows(
  biomass, 
  biomass_spike, 
  hyd, 
  hyd_spike, 
  hfe, 
  hfe_spike, 
  hyd_hfe, 
  hyd_hfe_spike
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


# List of input Gammarus files when zHydropsyche is changed
# List of input Chironomidae files when zHydropsyche is changed
biomass <- read.csv2("Multispp_temp_biomass_z_hyos_full_GAMM.csv")
biomass_spike <- read.csv2("Multispp_temp_biomass_spike_z_hyos_full_GAMM.csv") 
hyd <- read.csv2("Multispp_temp_biomass_hyd_z_hyos_full_GAMM.csv")
hyd_spike <- read.csv2("Multispp_temp_hyd_biomass_spike_z_hyos_full_GAMM.csv")
hfe <- read.csv2("Multispp_temp_biomass_hfe_z_hyos_full_GAMM.csv")
hfe_spike <- read.csv2("Multispp_temp_HFE_biomass_spike_z_hyos_full_GAMM.csv")
hyd_hfe <- read.csv2("Multispp_temp_hyd_HFE_biomass_z_hyos_full_GAMM.csv")
hyd_hfe_spike <- read.csv2("Multispp_temp_hyd_HFE_biomass_spike_z_hyos_full_GAMM.csv")


gamm <- bind_rows(
  biomass, 
  biomass_spike, 
  hyd, 
  hyd_spike, 
  hfe, 
  hfe_spike, 
  hyd_hfe, 
  hyd_hfe_spike
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


# List of input NZMS files when zHydropsyche is changed
biomass <- read.csv2("Multispp_temp_biomass_z_hyos_full_NZMS.csv")
biomass_spike <- read.csv2("Multispp_temp_biomass_spike_z_hyos_full_NZMS.csv") 
hyd <- read.csv2("Multispp_temp_biomass_hyd_z_hyos_full_NZMS.csv")
hyd_spike <- read.csv2("Multispp_temp_hyd_biomass_spike_z_hyos_full_NZMS.csv")
hfe <- read.csv2("Multispp_temp_biomass_hfe_z_hyos_full_NZMS.csv")
hfe_spike <- read.csv2("Multispp_temp_HFE_biomass_spike_z_hyos_full_NZMS.csv")
hyd_hfe <- read.csv2("Multispp_temp_hyd_HFE_biomass_z_hyos_full_NZMS.csv")
hyd_hfe_spike <- read.csv2("Multispp_temp_hyd_HFE_biomass_spike_z_hyos_full_NZMS.csv")


nzms <- bind_rows(
  biomass, 
  biomass_spike, 
  hyd, 
  hyd_spike, 
  hfe, 
  hfe_spike, 
  hyd_hfe, 
  hyd_hfe_spike
)
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