#####################
## Code to make sensitivity plots
#######################
library(tidyverse)
Multispp_mF_sens_temps_hfe <- read_csv("Multispp_mF_sens_temps_hfe.csv")

# we did this for all temperautres, but mainly want to look at baseline right now
sens_hfe <- Multispp_mF_sens_temps_hfe %>% 
  filter(TemperatureFactor == 1)

mod  <- lm(sens_hfe$Abundance ~ sens_hfe$SensitivityIncrement)
mod_sum <- summary(mod)

#coef(mod_sum)[8] = p value
# mod_sum[9] = adj r squared
# plot grid 

p1 <- ggplot(sens_hfe, aes(x = SensitivityIncrement,
               y = Abundance,
               color = StageGroup)) +
  geom_point(alpha = 0.5) +
  facet_grid(rows = vars(StageGroup), cols = vars(Parameter), scales = "free_x")+
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    x = "Sensitivity Increment to Flood Mortality",
    y = "Species Stage Abundance"
  )

x11()
p1


sens_hfe_slopes <- sens_hfe %>%
  group_by(StageGroup, Parameter) %>%
  summarise(
    scale_abund_slope = coef(lm(scale(Abundance) ~ scale(SensitivityIncrement)))[2],
    R_abund = coef(summary(lm((Abundance) ~ (SensitivityIncrement))))[8],
    p_abund = summary(lm((Abundance) ~ (SensitivityIncrement)))[9],
    scale_biomass_slope = coef(lm(scale(Biomass) ~ scale(SensitivityIncrement)))[2],
    R_biomass = coef(summary(lm((Biomass) ~ (SensitivityIncrement))))[8],
    p_biomass = summary(lm((Biomass) ~ (SensitivityIncrement)))[9],
    .groups = "drop"
  )
sens_hfe_slopes <- sens_hfe_slopes %>%
  mutate(Parameter = forcats::fct_reorder(
    Parameter, abs(scale_abund_slope), .fun = max
  ))

sens_hfe_slopes <- sens_hfe_slopes %>%
  mutate(
    abund_effect = case_when(
      abs(scale_abund_slope) < 0.3 ~ "small",
      abs(scale_abund_slope) < 0.5 ~ "moderate",
      TRUE ~ "large"
    ),
    biomass_effect = case_when(
      abs(scale_biomass_slope) < 0.3 ~ "small",
      abs(scale_biomass_slope) < 0.5 ~ "moderate",
      TRUE ~ "large"
    )
  )

Multispp_mh_sens_hyd <- read_csv("Multispp_mh_sens_hyd.csv")

ggplot(sens_hfe_slopes, aes(x = Parameter,
                            y = StageGroup,
                            fill = scale_abund_slope)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",      # negative slopes
    mid = "white",     # zero
    high = "red",      # positive slopes
    midpoint = 0,      # center at 0
    limits = c(-1, 1)  # fix min and max
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(sens_hfe_slopes, aes(x = Parameter,
                            y = StageGroup,
                            fill = scale_biomass_slope)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",      # negative slopes
    mid = "white",     # zero
    high = "red",      # positive slopes
    midpoint = 0,      # center at 0
    limits = c(-1, 1)  # fix min and max
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 
Multispp_mh_sens_temps <- read_csv("Multispp_mh_sens_temps.csv")
# we did this for all temperautres, but mainly want to look at baseline right now
sens_hyd <- Multispp_mh_sens_temps %>% 
  filter(TemperatureFactor == 1)

p2 <- ggplot(sens_hyd, aes(x = SensitivityIncrement,
                           y = Abundance,
                           color = StageGroup)) +
  geom_point(alpha = 0.5) +
  facet_grid(rows = vars(StageGroup), cols = vars(Parameter), scales = "free_x")+
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    x = "Sensitivity Increment to Hydropeaking Mortality",
    y = "Species Stage Abundance"
  )

sens_hyd_scaled <- sens_hyd %>%
  mutate(
    SensitivityIncrement_s = scale(SensitivityIncrement),
    Abundance_s = scale(Abundance),
    Biomass_s = scale(Biomass)
  )

sens_hyd_slopes <- sens_hyd_scaled %>%
  group_by(StageGroup, Parameter) %>%
  summarise(
    scale_abund_slope =
      coef(lm(Abundance_s ~ SensitivityIncrement_s))[2],
    scale_biomass_slope =
      coef(lm(Biomass_s ~ SensitivityIncrement_s))[2],
    .groups = "drop"
  )


sens_hyd_slopes <- sens_hyd_slopes %>%
  mutate(
    abund_effect = case_when(
      abs(scale_abund_slope) < 0.3 ~ "small",
      abs(scale_abund_slope) < 0.5 ~ "moderate",
      TRUE ~ "large"
    ),
    biomass_effect = case_when(
      abs(scale_biomass_slope) < 0.3 ~ "small",
      abs(scale_biomass_slope) < 0.5 ~ "moderate",
      TRUE ~ "large"
    )
  )


ggplot(sens_hyd_slopes, aes(x = Parameter,
                            y = StageGroup,
                            fill = scale_abund_slope)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",      # negative slopes
    mid = "white",     # zero
    high = "red",      # positive slopes
    midpoint = 0,      # center at 0
    limits = c(-1, 1)  # fix min and max
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(sens_hyd_slopes, aes(x = Parameter,
                            y = StageGroup,
                            fill = scale_biomass_slope)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",      # negative slopes
    mid = "white",     # zero
    high = "red",      # positive slopes
    midpoint = 0,      # center at 0
    limits = c(-1, 1)  # fix min and max
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

Multispp_mT_sens_temps <- read_csv("Multispp_mT_sens_temps.csv")
sens_temp <- Multispp_mT_sens_temps %>% 
  filter(TemperatureFactor == 1)

p2 <- ggplot(sens_temp, aes(x = SensitivityIncrement,
                           y = Abundance,
                           color = StageGroup)) +
  geom_point(alpha = 0.5) +
  facet_grid(rows = vars(StageGroup), cols = vars(Parameter), scales = "free_x")+
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    x = "Sensitivity Increment to Temperature Mortality",
    y = "Species Stage Abundance"
  )

sens_temp_scaled <- sens_temp %>%
  mutate(
    SensitivityIncrement_s = scale(SensitivityIncrement),
    Abundance_s = scale(Abundance),
    Biomass_s = scale(Biomass)
  )

sens_temp_slopes <- sens_temp_scaled %>%
  group_by(StageGroup, Parameter) %>%
  summarise(
    scale_abund_slope =
      coef(lm(Abundance_s ~ SensitivityIncrement_s))[2],
    scale_biomass_slope =
      coef(lm(Biomass_s ~ SensitivityIncrement_s))[2],
    .groups = "drop"
  )


sens_temp_slopes <- sens_temp_slopes %>%
  mutate(
    abund_effect = case_when(
      abs(scale_abund_slope) < 0.3 ~ "small",
      abs(scale_abund_slope) < 0.5 ~ "moderate",
      TRUE ~ "large"
    ),
    biomass_effect = case_when(
      abs(scale_biomass_slope) < 0.3 ~ "small",
      abs(scale_biomass_slope) < 0.5 ~ "moderate",
      TRUE ~ "large"
    )
  )


ggplot(sens_temp_slopes, aes(x = Parameter,
                            y = StageGroup,
                            fill = scale_abund_slope)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",      # negative slopes
    mid = "white",     # zero
    high = "red",      # positive slopes
    midpoint = 0,      # center at 0
    limits = c(-1, 1)  # fix min and max
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(sens_temp_slopes, aes(x = Parameter,
                            y = StageGroup,
                            fill = scale_biomass_slope)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",      # negative slopes
    mid = "white",     # zero
    high = "red",      # positive slopes
    midpoint = 0,      # center at 0
    limits = c(-1, 1)  # fix min and max
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

