## Multispecies Partial Eta Square
library(tidyverse)
library(effectsize)

# read in the csv's
Multispp_temp_biomass <- read_csv("ColoradoRiverInverts/Multispp_temp_biomass_z_full.csv")
Multispp_temp_biomass_spike <- read_csv("ColoradoRiverInverts/Multispp_temp_biomass_spike_z_full.csv")
Multispp_temp_hyd_biomass_spike <- read_csv("ColoradoRiverInverts/Multispp_temp_hyd_biomass_spike_z_full.csv")

Multispp_temp_HFE_biomass <- read_csv("ColoradoRiverInverts/Multispp_temp_biomass_hfe_z_full.csv")
Multispp_temp_HFE_biomass_spike <- read_csv("ColoradoRiverInverts/Multispp_temp_biomass_hfe_z_spike_full.csv")
Multispp_temp_hyd_HFE_biomass <- read_csv("ColoradoRiverInverts/Multispp_temp_hyd_HFE_biomass_z_full.csv")
Multispp_temp_hyd_HFE_biomas_spike <- read_csv("ColoradoRiverInverts/Multispp_temp_hyd_HFE_biomass_spike_z_full.csv")

multi_biomass_combo <- bind_rows(
  Multispp_temp_biomass %>%  mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  Multispp_temp_biomass_spike %>% mutate(source = "Temperature", temp_spike = 1, hydropeaking = 0, HFE = 0),
  Multispp_temp_hyd_biomass %>%  mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  Multispp_temp_hyd_biomass_spike %>% mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  Multispp_temp_HFE_biomass %>% mutate(source = "Temperature & HFE",  temp_spike = 0, hydropeaking = 0, HFE = 1),
  Multispp_temp_HFE_biomass_spike %>% mutate(source = "Temperature & spike & HFE",  temp_spike = 1, hydropeaking = 0, HFE = 1),
  Multispp_temp_hyd_HFE_biomass %>% mutate(source = "Temperature & hydropeaking & HFE",  temp_spike = 0, hydropeaking = 1, HFE = 1),
  Multispp_temp_hyd_HFE_biomass_spike %>% mutate(source = "Temperature & spike & hydropeaking & HFE",  temp_spike = 1, hydropeaking = 1, HFE = 1)
)

multi_biomass_combo <- bind_rows(
  Multispp_temp_biomass %>%  mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  Multispp_temp_biomass_spike %>% mutate(source = "Temperature", temp_spike = 1, hydropeaking = 0, HFE = 0),
  Multispp_temp_hyd_biomass %>%  mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  Multispp_temp_hyd_biomass_spike %>% mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  Multispp_temp_HFE_biomass %>% mutate(source = "Temperature & HFE",  temp_spike = 0, hydropeaking = 0, HFE = 1),
  Multispp_temp_HFE_biomass_spike %>% mutate(source = "Temperature & spike & HFE",  temp_spike = 1, hydropeaking = 0, HFE = 1),
  Multispp_temp_hyd_HFE_biomass %>% mutate(source = "Temperature & hydropeaking & HFE",  temp_spike = 0, hydropeaking = 1, HFE = 1),
  Multispp_temp_hyd_HFE_biomass_spike %>% mutate(source = "Temperature & spike & hydropeaking & HFE",  temp_spike = 1, hydropeaking = 1, HFE = 1)
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

multi_biomass_combo$source <- ordered(multi_biomass_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking",
  "Temperature & spike & hydropeaking",
  "Temperature & HFE",
  "Temperature & spike & HFE",
  "Temperature & hydropeaking & HFE",
  "Temperature & spike & hydropeaking & HFE"
))  
multi_biomass_combo$taxa <- ordered(multi_biomass_combo$taxa, levels = c("BAET", "HYOS", "CHIR", "GAMM", "NZMS"))
# 
# # run this on HPC since laptop cannot allocate this much memory
# # because if categorical, fit manova
# multi_biomass <- aov(biomass ~ (temperature + temp_spike + hydropeaking + HFE)^4, data = multi_biomass_combo)
# multi_abund <- aov(abundance ~ (temperature + temp_spike + hydropeaking + HFE)^4, data = multi_abund_combo)
# 
# # Extract partial eta squared values
# eta_sq_df.multi.b <- eta_squared(multi_biomass, partial = TRUE) %>%
#   as.data.frame() %>%
#   arrange(desc(Eta2_partial))  # Order from largest to smallest effect
# 
# eta_sq_df.multi.a <- eta_squared(multi_abund, partial = TRUE) %>%
#   as.data.frame() %>%
#   arrange(desc(Eta2_partial))  # Order from largest to smallest effect
# 


# Function to fit ANOVA and extract partial eta squared
compute_eta_sq <- function(df, response_var) {
  model <- aov(as.formula(paste(response_var, "~ (temperature + temp_spike + hydropeaking + HFE + q)^5")), data = df)
  eta_squared(model, partial = TRUE) %>%
    as.data.frame() %>%
    arrange(desc(Eta2_partial))
}

# Run analysis per taxon for biomass
eta_sq_taxon_biomass <- multi_biomass_combo %>%
  group_by(taxa) %>%
  group_split() %>%
  map_df(~ compute_eta_sq(.x, "biomass"), .id = "taxon")


write.csv(eta_sq_taxon_biomass, "Multispp_partialetasquare_biomass_z.csv", row.names = FALSE)
