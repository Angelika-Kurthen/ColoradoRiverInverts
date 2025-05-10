# ## Multispecies Partial Eta Square
library(data.table)
library(purrr)
library(effectsize)
library(dplyr)
library(parallel)

# # List of input Hydropsyche files when zi is changed
# biomass <- read.csv2("Multispp_temp_biomass_z_full_HYOS.csv")
# biomass_spike <- read.csv2("Multispp_temp_biomass_spike_z_full_HYOS.csv") 
# hyd <- read.csv2("Multispp_temp_biomass_hyd_z_full_HYOS.csv")
# hyd_spike <- read.csv2("Multispp_temp_hyd_biomass_spike_z_full_HYOS.csv")
# hfe <- read.csv2("Multispp_temp_biomass_hfe_z_full_HYOS.csv")
# hfe_spike <- read.csv2( "Multispp_temp_HFE_biomass_spike_z_full_HYOS.csv")
# hyd_hfe <- read.csv2("Multispp_temp_hyd_HFE_biomass_z_full_HYOS.csv")
# hyd_hfe_spike <- read.csv2("Multispp_temp_hyd_HFE_biomass_spike_z_full_HYOS.csv")
# 
# 
# hyos <- bind_rows(
#   biomass, 
#   biomass_spike, 
#   hyd, 
#   hyd_spike, 
#   hfe, 
#   hfe_spike, 
#   hyd_hfe, 
#   hyd_hfe_spike
# )
# 
# label.names <- c(
#   "Temperature", 
#   "Temperature +\nSpike", 
#   "Temperature +\nHydropeaking",
#   "Temperature +\nSpike + Hydropeaking",
#   "Temperature +\nHFE",
#   "Temperature +\nSpike + HFE",
#   "Temperature +\nHydropeaking + HFE",
#   "Temperature +\nSpike + Hydropeaking + HFE"
# )
# 
# hyos$source <- ordered(hyos$source, levels = c(
#   "Temperature", 
#   "Temperature & spike",
#   "Temperature & hydropeaking",
#   "Temperature & spike & hydropeaking",
#   "Temperature & HFE",
#   "Temperature & spike & HFE",
#   "Temperature & hydropeaking & HFE",
#   "Temperature & spike & hydropeaking & HFE"
# ))  
# 
# # Function to fit ANOVA and extract partial eta squared
# compute_eta_sq <- function(df, response_var) {
#   model <- aov(as.formula(paste(response_var, "~ (temperature + temp_spike + hydropeaking + HFE + q)^2")), data = df)
#   eta_squared(model, partial = TRUE) %>%
#     as.data.frame()
# }
# 
# # Run analysis per taxon for abundance
# eta_sq_taxon_abundance <- hyos %>%
#   map_df(~ compute_eta_sq(.x, "biomass"))
# write.csv(eta_sq_df.multi.b, "Multispp_partialetasquare_biomass_hyos_hyos.csv", row.names = FALSE)
# rm(hyos)
# 
# 
# # List of input Baetidae files when zi is changed
# biomass <- read.csv2("Multispp_temp_biomass_z_full_BAET.csv")
# biomass_spike <- read.csv2("Multispp_temp_biomass_spike_z_full_BAET.csv") 
# hyd <- read.csv2("Multispp_temp_biomass_hyd_z_full_BAET.csv")
# hyd_spike <- read.csv2("Multispp_temp_hyd_biomass_spike_z_full_BAET.csv")
# hfe <- read.csv2("Multispp_temp_biomass_hfe_z_full_BAET.csv")
# hfe_spike <- read.csv2( "Multispp_temp_HFE_biomass_spike_z_full_BAET.csv")
# hyd_hfe <- read.csv2("Multispp_temp_hyd_HFE_biomass_z_full_BAET.csv")
# hyd_hfe_spike <- read.csv2("Multispp_temp_hyd_HFE_biomass_spike_z_full_BAET.csv")
# 
# 
# 
# baet <- bind_rows(
#   biomass, 
#   biomass_spike, 
#   hyd, 
#   hyd_spike, 
#   hfe, 
#   hfe_spike, 
#   hyd_hfe, 
#   hyd_hfe_spike
# )
# 
# label.names <- c(
#   "Temperature", 
#   "Temperature +\nSpike", 
#   "Temperature +\nHydropeaking",
#   "Temperature +\nSpike + Hydropeaking",
#   "Temperature +\nHFE",
#   "Temperature +\nSpike + HFE",
#   "Temperature +\nHydropeaking + HFE",
#   "Temperature +\nSpike + Hydropeaking + HFE"
# )
# 
# baet$source <- ordered(baet$source, levels = c(
#   "Temperature", 
#   "Temperature & spike",
#   "Temperature & hydropeaking",
#   "Temperature & spike & hydropeaking",
#   "Temperature & HFE",
#   "Temperature & spike & HFE",
#   "Temperature & hydropeaking & HFE",
#   "Temperature & spike & hydropeaking & HFE"
# ))  
# 
# # Function to fit ANOVA and extract partial eta squared
# compute_eta_sq <- function(df, response_var) {
#   model <- aov(as.formula(paste(response_var, "~ (temperature + temp_spike + hydropeaking + HFE + q)^2")), data = df)
#   eta_squared(model, partial = TRUE) %>%
#     as.data.frame()
# }
# 
# # Run analysis per taxon for abundance
# eta_sq_taxon_abundance <- baet %>%
#   map_df(~ compute_eta_sq(.x, "biomass"))
# write.csv(eta_sq_df.multi.b, "Multispp_partialetasquare_biomass_hyos_baet.csv", row.names = FALSE)
# rm(baet)
# 
# 
# # List of input Chironomidae files when zi is changed
# biomass <- read.csv2("Multispp_temp_biomass_z_full_CHIR.csv")
# biomass_spike <- read.csv2("Multispp_temp_biomass_spike_z_full_CHIR.csv") 
# hyd <- read.csv2("Multispp_temp_biomass_hyd_z_full_CHIR.csv")
# hyd_spike <- read.csv2("Multispp_temp_hyd_biomass_spike_z_full_CHIR.csv")
# hfe <- read.csv2("Multispp_temp_biomass_hfe_z_full_CHIR.csv")
# hfe_spike <- read.csv2( "Multispp_temp_HFE_biomass_spike_z_full_CHIR.csv")
# hyd_hfe <- read.csv2("Multispp_temp_hyd_HFE_biomass_z_full_CHIR.csv")
# hyd_hfe_spike <- read.csv2("Multispp_temp_hyd_HFE_biomass_spike_z_full_CHIR.csv")
# 
# 
# 
# chir <- bind_rows(
#   biomass, 
#   biomass_spike, 
#   hyd, 
#   hyd_spike, 
#   hfe, 
#   hfe_spike, 
#   hyd_hfe, 
#   hyd_hfe_spike
# )
# 
# label.names <- c(
#   "Temperature", 
#   "Temperature +\nSpike", 
#   "Temperature +\nHydropeaking",
#   "Temperature +\nSpike + Hydropeaking",
#   "Temperature +\nHFE",
#   "Temperature +\nSpike + HFE",
#   "Temperature +\nHydropeaking + HFE",
#   "Temperature +\nSpike + Hydropeaking + HFE"
# )
# 
# chir$source <- ordered(chir$source, levels = c(
#   "Temperature", 
#   "Temperature & spike",
#   "Temperature & hydropeaking",
#   "Temperature & spike & hydropeaking",
#   "Temperature & HFE",
#   "Temperature & spike & HFE",
#   "Temperature & hydropeaking & HFE",
#   "Temperature & spike & hydropeaking & HFE"
# ))  
# 
# # Function to fit ANOVA and extract partial eta squared
# compute_eta_sq <- function(df, response_var) {
#   model <- aov(as.formula(paste(response_var, "~ (temperature + temp_spike + hydropeaking + HFE + q)^2")), data = df)
#   eta_squared(model, partial = TRUE) %>%
#     as.data.frame()
# }
# 
# # Run analysis per taxon for abundance
# eta_sq_taxon_abundance <- chir %>%
#   map_df(~ compute_eta_sq(.x, "biomass"))
# write.csv(eta_sq_df.multi.b, "Multispp_partialetasquare_biomass_hyos_chir.csv", row.names = FALSE)
# rm(chir)
# 
# 
# # List of input Gammarus files when zi is changed
# biomass <- read.csv2("Multispp_temp_biomass_z_full_GAMM.csv")
# biomass_spike <- read.csv2("Multispp_temp_biomass_spike_z_full_GAMM.csv") 
# hyd <- read.csv2("Multispp_temp_biomass_hyd_z_full_GAMM.csv")
# hyd_spike <- read.csv2("Multispp_temp_hyd_biomass_spike_z_full_GAMM.csv")
# hfe <- read.csv2("Multispp_temp_biomass_hfe_z_full_GAMM.csv")
# hfe_spike <- read.csv2( "Multispp_temp_HFE_biomass_spike_z_full_GAMM.csv")
# hyd_hfe <- read.csv2("Multispp_temp_hyd_HFE_biomass_z_full_GAMM.csv")
# hyd_hfe_spike <- read.csv2("Multispp_temp_hyd_HFE_biomass_spike_z_full_GAMM.csv")
# 
# 
# 
# gamm <- bind_rows(
#   biomass, 
#   biomass_spike, 
#   hyd, 
#   hyd_spike, 
#   hfe, 
#   hfe_spike, 
#   hyd_hfe, 
#   hyd_hfe_spike
# )
# 
# label.names <- c(
#   "Temperature", 
#   "Temperature +\nSpike", 
#   "Temperature +\nHydropeaking",
#   "Temperature +\nSpike + Hydropeaking",
#   "Temperature +\nHFE",
#   "Temperature +\nSpike + HFE",
#   "Temperature +\nHydropeaking + HFE",
#   "Temperature +\nSpike + Hydropeaking + HFE"
# )
# 
# gamm$source <- ordered(gamm$source, levels = c(
#   "Temperature", 
#   "Temperature & spike",
#   "Temperature & hydropeaking",
#   "Temperature & spike & hydropeaking",
#   "Temperature & HFE",
#   "Temperature & spike & HFE",
#   "Temperature & hydropeaking & HFE",
#   "Temperature & spike & hydropeaking & HFE"
# ))  
# 
# # Function to fit ANOVA and extract partial eta squared
# compute_eta_sq <- function(df, response_var) {
#   model <- aov(as.formula(paste(response_var, "~ (temperature + temp_spike + hydropeaking + HFE + q)^2")), data = df)
#   eta_squared(model, partial = TRUE) %>%
#     as.data.frame()
# }
# 
# # Run analysis per taxon for abundance
# eta_sq_taxon_abundance <- gamm %>%
#   map_df(~ compute_eta_sq(.x, "biomass"))
# write.csv(eta_sq_df.multi.b, "Multispp_partialetasquare_biomass_hyos_gamm.csv", row.names = FALSE)
# rm(gamm)
# 
# 
# # List of input NZMS files when zi is changed
# biomass <- read.csv2("Multispp_temp_biomass_z_full_NZMS.csv")
# biomass_spike <- read.csv2("Multispp_temp_biomass_spike_z_full_NZMS.csv") 
# hyd <- read.csv2("Multispp_temp_biomass_hyd_z_full_NZMS.csv")
# hyd_spike <- read.csv2("Multispp_temp_hyd_biomass_spike_z_full_NZMS.csv")
# hfe <- read.csv2("Multispp_temp_biomass_hfe_z_full_NZMS.csv")
# hfe_spike <- read.csv2( "Multispp_temp_HFE_biomass_spike_z_full_NZMS.csv")
# hyd_hfe <- read.csv2("Multispp_temp_hyd_HFE_biomass_z_full_NZMS.csv")
# hyd_hfe_spike <- read.csv2("Multispp_temp_hyd_HFE_biomass_spike_z_full_NZMS.csv")
# 
# 
# nzms <- bind_rows(
#   biomass, 
#   biomass_spike, 
#   hyd, 
#   hyd_spike, 
#   hfe, 
#   hfe_spike, 
#   hyd_hfe, 
#   hyd_hfe_spike
# )
# nzms <- bind_rows(
#   "Multispp_temp_biomass_z_hyos_full_NZMS.csv",
#   "Multispp_temp_biomass_spike_z_hyos_full_NZMS.csv", 
#   "Multispp_temp_biomass_hyd_z_hyos_full_NZMS.csv",
#   "Multispp_temp_hyd_biomass_spike_z_hyos_full_NZMS.csv",
#   "Multispp_temp_biomass_hfe_z_hyos_full_NZMS.csv",
#   "Multispp_temp_HFE_biomass_spike_z_hyos_full_NZMS.csv", 
#   "Multispp_temp_hyd_HFE_biomass_spike_z_hyos_full_NZMS.csv", 
# )
# 
# label.names <- c(
#   "Temperature", 
#   "Temperature +\nSpike", 
#   "Temperature +\nHydropeaking",
#   "Temperature +\nSpike + Hydropeaking",
#   "Temperature +\nHFE",
#   "Temperature +\nSpike + HFE",
#   "Temperature +\nHydropeaking + HFE",
#   "Temperature +\nSpike + Hydropeaking + HFE"
# )
# 
# nzms$source <- ordered(nzms$source, levels = c(
#   "Temperature", 
#   "Temperature & spike",
#   "Temperature & hydropeaking",
#   "Temperature & spike & hydropeaking",
#   "Temperature & HFE",
#   "Temperature & spike & HFE",
#   "Temperature & hydropeaking & HFE",
#   "Temperature & spike & hydropeaking & HFE"
# ))  
# 
# # Function to fit ANOVA and extract partial eta squared
# compute_eta_sq <- function(df, response_var) {
#   model <- aov(as.formula(paste(response_var, "~ (temperature + temp_spike + hydropeaking + HFE + q)^2")), data = df)
#   eta_squared(model, partial = TRUE) %>%
#     as.data.frame()
# }
# 
# # Run analysis per taxon for abundance
# eta_sq_taxon_abundance <- nzms %>%
#   map_df(~ compute_eta_sq(.x, "biomass"))
# write.csv(eta_sq_df.multi.b, "Multispp_partialetasquare_biomass_hyos_nzms.csv", row.names = FALSE)
# rm(nzms)


# List of input NZMS files when zi is changed

# zi = 1 
temp1 <- readRDS("/home/kurthena/ColoradoRiverInverts/Multispp_Temp_Z1.rds")
temp_spike1 <- readRDS("/home/kurthena/ColoradoRiverInverts/Multispp_Temp_Spike_Z1.rds")

temp1 <- readRDS("C:/Users/jelly/Documents/ColoradoRiverInverts/Multispp_Temp_Z1.rds")
temp_spike1 <- readRDS("C:/Users/jelly/Documents/ColoradoRiverInverts/Multispp_Temp_Spike_Z1.rds")




hyd1 <- readRDS("/home/kurthena/ColoradoRiverInverts/Multispp_Temp_Hyd_Z1.rds")
hyd_spike1 <- readRDS("/home/kurthena/ColoradoRiverInverts/Multispp_Temp_Hyd_Spike_Z1.rds")
hfe1 <- readRDS("/home/kurthena/ColoradoRiverInverts/Multispp_Temp_HFE_Z1.rds")
hfe_spike1 <- readRDS("/home/kurthena/ColoradoRiverInverts/Multispp_Temp_HFE_Spike_Z1.rds")
hyd_hfe1 <- readRDS("/home/kurthena/ColoradoRiverInverts/Multispp_Temp_Hyd_HFE_Z1.rds")
hyd_hfe_spike1 <-readRDS("/home/kurthena/ColoradoRiverInverts/Multispp_Temp_Z1.rds")


partialeta <- bind_rows(
    temp1  %>%  mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
    temp_spike1 %>% mutate(source = "Temperature", temp_spike = 1, hydropeaking = 0, HFE = 0),
    hyd1 %>%  mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
    hyd_spike1 %>% mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
    hfe1 %>% mutate(source = "Temperature & HFE",  temp_spike = 0, hydropeaking = 0, HFE = 1),
    hfe_spike1 %>% mutate(source = "Temperature & spike & HFE",  temp_spike = 1, hydropeaking = 0, HFE = 1),
    hyd_hfe1 %>% mutate(source = "Temperature & hydropeaking & HFE",  temp_spike = 0, hydropeaking = 1, HFE = 1),
    hyd_hfe_spike1 %>% mutate(source = "Temperature & spike & hydropeaking & HFE",  temp_spike = 1, hydropeaking = 1, HFE = 1)
)
  
rm(temp1)
rm(temp_spike1)
rm(hyd1)
rm(hyd_spike1)
rm(hfe1)
rm(hfe_spike1)
rm(hyd_hfe1)
rm(hyd_hfe_spike1)

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

partialeta$source <- ordered(partialeta$source, levels = c(
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
  model <- aov(as.formula(paste(response_var, "~ (temperature + temp_spike + hydropeaking + HFE)^4")), data = df)
  eta_squared(model, partial = TRUE) %>%
    as.data.frame()
}

# Run analysis per taxon for biomass
eta_sq_taxon_biomass <- partialeta %>%
  group_by(taxa) %>%
  group_split() %>%
  map_df(~ compute_eta_sq(.x, "rel_biomass"), .id = "taxon")

write.csv(eta_sq_taxon_biomass, "Multispp_partialetasquare_biomass_z1.csv")
