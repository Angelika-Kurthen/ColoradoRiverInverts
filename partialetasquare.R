## Multispecies Partial Eta Square
library(tidyverse)
library(effectsize)

# read in the csv's
Multispp_temp_biomass <- read_csv("Multispp_temp_biomass.csv")
Multispp_temp_abund <- read_csv("Multispp_temp_abund.csv")
Multispp_temp_abund_spike <- read_csv("Multispp_temp_abund_spike.csv")
Multispp_temp_biomass_spike <- read_csv("Multispp_temp_biomass_spike.csv")
Multispp_temp_hyd_biomass <- read_csv("Multispp_temp_biomass_hyd.csv")
Multispp_temp_hyd_biomass_spike <- read_csv("Multispp_temp_hyd_biomass_spike.csv")
Multispp_temp_hyd_HFE_biomass_spike <- read_csv("Multispp_temp_hyd_HFE_biomass_spike.csv")
Multispp_temp_hyd_HFE_biomass <- read_csv("Multispp_temp_hyd_HFE_biomass.csv")
Multispp_temp_hyd_HFE_abund_spike <- read_csv("Multispp_temp_hyd_HFE_abund_spike.csv")
Multispp_temp_HFE_abund_spike <- read_csv("Multispp_temp_HFE_abund_spike.csv")
Multispp_temp_HFE_abund <- read_csv("Multispp_temp_abund_HFE.csv")
Multispp_temp_HFE_biomass <- read_csv("Multispp_temp_biomass_HFE.csv")
Multispp_temp_hyd_abund <- read_csv("Multispp_temp_abund_hyd.csv")
Multispp_temp_HFE_biomass_spike <- read_csv("Multispp_temp_HFE_biomass_spike.csv")
Multispp_temp_hyd_abund_spike <- read_csv("Multispp_temp_hyd_abund_spike.csv")
Multispp_temp_hyd_HFE_abund <- read_csv("Multispp_temp_hyd_HFE_abund.csv")
Multispp_temp_S3 <- read_csv("MultisppS3_temp_biomass.csv")
Multispp_temp_S3_spike <- read_csv("MultisppS3_temp_biomass_spike.csv")

multi_abund_combo <- bind_rows(
  Multispp_temp_abund %>%  mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  Multispp_temp_abund_spike %>% mutate(source = "Temperature", temp_spike = 1, hydropeaking = 0, HFE = 0),
  Multispp_temp_hyd_abund %>%  mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  Multispp_temp_hyd_abund_spike %>% mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  Multispp_temp_HFE_abund %>% mutate(source = "Temperature & HFE",  temp_spike = 0, hydropeaking = 0, HFE = 1),
  Multispp_temp_HFE_abund_spike %>% mutate(source = "Temperature & spike & HFE",  temp_spike = 1, hydropeaking = 0, HFE = 1),
  Multispp_temp_hyd_HFE_abund %>% mutate(source = "Temperature & hydropeaking & HFE",  temp_spike = 0, hydropeaking = 1, HFE = 1),
  Multispp_temp_hyd_HFE_abund_spike %>% mutate(source = "Temperature & spike & hydropeaking & HFE",  temp_spike = 1, hydropeaking = 1, HFE = 1)
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

multi_abund_combo$source <- ordered(multi_abund_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike",
  "Temperature & hydropeaking",
  "Temperature & spike & hydropeaking",
  "Temperature & HFE",
  "Temperature & spike & HFE",
  "Temperature & hydropeaking & HFE",
  "Temperature & spike & hydropeaking & HFE"
  ))  
multi_abund_combo$taxa <- ordered(multi_abund_combo$taxa, levels = c("BAET", "HYOS", "CHIR", "GAMM", "NZMS"))

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
  model <- aov(as.formula(paste(response_var, "~ (temperature + temp_spike + hydropeaking + HFE)^4")), data = df)
  eta_squared(model, partial = TRUE) %>%
    as.data.frame() %>%
    arrange(desc(Eta2_partial))
}

# Run analysis per taxon for biomass
eta_sq_taxon_biomass <- multi_biomass_combo %>%
  group_by(taxa) %>%
  group_split() %>%
  map_df(~ compute_eta_sq(.x, "biomass"), .id = "taxon")

# Run analysis per taxon for abundance
eta_sq_taxon_abundance <- multi_abund_combo %>%
  group_by(taxa) %>%
  group_split() %>%
  map_df(~ compute_eta_sq(.x, "abundance"), .id = "taxon")
write.csv(eta_sq_df.multi.b, "Multispp_partialetasquare_biomass.csv", row.names = FALSE)
write.csv(eta_sq_df.multi.a, "Multispp_partialetasquare_abund.csv", row.names = FALSE)
