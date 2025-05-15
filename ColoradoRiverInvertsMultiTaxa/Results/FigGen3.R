###########
## Figure generation for 3rd manuscript
##########
library(ggpubr)
library(patchwork)
library(tidyverse)
library(svglite)
library(ggraph)
library(dplyr)
library(broom) 
library(readr)
library(lm.beta)
library(effectsize)

# increasing HFEs: 
taxa.colors =c("#66CCEE", "#AA3377", "#228833", "#CCBB44", "#4477AA")

# 8 different scenarios zi 
# Load and combine all datasets, tagging them hyos
multi_biomass_combo_hyos <- bind_rows(
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_biomass_z_hyos.csv") %>%
    mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_biomass_spike_z_hyos.csv") %>%
    mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_biomass_hyd_z_hyos.csv") %>%
    mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_hyd_biomass_spike_z_hyos.csv") %>%
    mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_biomass_hfe_z_hyos.csv") %>%
    mutate(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_HFE_biomass_spike_z_hyos.csv") %>%
    mutate(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_hyd_HFE_biomass_z_hyos.csv") %>%
    mutate(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_hyd_HFE_biomass_spike_z_hyos.csv") %>%
    mutate(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
)


taxa.colors =c("#66CCEE", "#AA3377", "#228833", "#CCBB44", "#4477AA")

multi_biomass_combo_hyos$q <- factor(as.character(multi_biomass_combo_hyos$q), levels = c("1", "2", "4", "8"))

multi_biomass_combo_hyos <- filter(multi_biomass_combo_hyos, q != "1")
# Define the labels for q as character strings (without expression)
q_labels <- as_labeller(c(
  "2" = "z[Hydropsyche]==2",
  "4" = "z[Hydropsyche]==4",
  "8" = "z[Hydropsyche]==8"), default = label_parsed
)
label.names <- c(
  "Temperature", 
  "Temperature +\nSpike" , 
  "Temperature +\nHydropeaking", 
  "Temperature +\nSpike + Hydropeaking",
  "Temperature +\nHFE", 
  "Temperature +\nSpike + HFE",
  "Temperature +\nHydropeaking + HFE", 
  "Temperature +\nSpike + Hydropeaking + HFE"
)
multi_biomass_combo_hyos$source <- ordered(multi_biomass_combo_hyos$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE",
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"
))  
multi_biomass_combo_hyos$taxa <- ordered(multi_biomass_combo_hyos$taxa, levels = c("BAET", "HYOS", "CHIR", "GAMM", "NZMS"))
# multi_biomass_combo <- multi_biomass_combo %>% 
#   filter(q == 1 | q == 2 | q == 4)

multi_biomass_hyos <- ggplot(data = multi_biomass_combo_hyos, 
                      aes(x = as.factor(temperature), y = (mean_biomass), color = taxa)) +
  geom_errorbar(aes(ymin = (mean_biomass - sd_biomass), ymax = (mean_biomass + sd_biomass)), 
                linewidth = 0.75, width = 0.2, position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%")) +
  facet_wrap(. ~ source + q, nrow = 8,
             labeller = labeller(q = q_labels))+
  theme_bw() +
  scale_color_manual(name = " ", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Hydropsyche"), " spp.")),expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")),  expression(italic("P. antipodarum"))), values=taxa.colors) +
  xlab("Temperature") +
  labs(y= "Relative Biomass") +
  theme(text = element_text(size = 13), 
        axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), 
        legend.key = element_rect(fill = "transparent"))

ggsave("multi_z_hyos.png", multi_biomass_hyos, device = "png", width = 14.5, height = 12, units = "in", dpi = "retina")


# Load and combine all datasets, tagging them baet
multi_biomass_combo_baet <- bind_rows(
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_biomass_z_baet.csv") %>%
    mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_biomass_spike_z_baet.csv") %>%
    mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_biomass_hyd_z_baet.csv") %>%
    mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_hyd_biomass_spike_z_baet.csv") %>%
    mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_biomass_hfe_z_baet.csv") %>%
    mutate(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_HFE_biomass_spike_z_baet.csv") %>%
    mutate(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_hyd_HFE_biomass_z_baet.csv") %>%
    mutate(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_hyd_HFE_biomass_spike_z_baet.csv") %>%
    mutate(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
)


taxa.colors =c("#66CCEE", "#AA3377", "#228833", "#CCBB44", "#4477AA")

q_labels <- as_labeller(c(
  "2" = "z[Baetidae]==2",
  "4" = "z[Baetidae]==4",
  "8" = "z[Baetidae]==8"), default = label_parsed
)

label.names <- c(
  "Temperature", 
  "Temperature +\nSpike" , 
  "Temperature +\nHydropeaking", 
  "Temperature +\nSpike + Hydropeaking",
  "Temperature +\nHFE", 
  "Temperature +\nSpike + HFE",
  "Temperature +\nHydropeaking + HFE", 
  "Temperature +\nSpike + Hydropeaking + HFE"
)

multi_biomass_combo_baet$source <- ordered(multi_biomass_combo_baet$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE",
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"
))  
multi_biomass_combo_baet$taxa <- ordered(multi_biomass_combo_baet$taxa, levels = c("BAET", "HYOS", "CHIR", "GAMM", "NZMS"))
multi_biomass_combo_baet$q <- as.factor(multi_biomass_combo_baet$q)

multi_biomass_combo_baet <- filter(multi_biomass_combo_baet, q != "1")
multi_biomass_baet <- ggplot(data = multi_biomass_combo_baet, 
                        aes(x = as.factor(temperature), y = (mean_biomass), color = taxa)) +
  geom_errorbar(aes(ymin = (mean_biomass - sd_biomass), ymax = (mean_biomass + sd_biomass)), 
                linewidth = 0.75, width = 0.2, position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%")) +
  facet_wrap(. ~ source + q, nrow = 8,
             labeller = labeller(q = q_labels))+
  theme_bw() +
  scale_color_manual(name = " ", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Hydropsyche"), " spp.")),expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")),  expression(italic("P. antipodarum"))), values=taxa.colors) +
  xlab("Temperature") +
  labs(y= "Relative Biomass") +
  theme(text = element_text(size = 13), 
        axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), 
        legend.key = element_rect(fill = "transparent"))
ggsave("multi_z_baet.png", multi_biomass_baet, device = "png", width = 14.5, height = 12, units = "in", dpi = "retina")

# Load and combine all datasets, tagging them chir
multi_biomass_combo_chir <- bind_rows(
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_biomass_z_chir.csv") %>%
    mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_biomass_spike_z_chir.csv") %>%
    mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_biomass_hyd_z_chir.csv") %>%
    mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_hyd_biomass_spike_z_chir.csv") %>%
    mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_biomass_hfe_z_chir.csv") %>%
    mutate(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_HFE_biomass_spike_z_chir.csv") %>%
    mutate(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_hyd_HFE_biomass_z_chir.csv") %>%
    mutate(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_hyd_HFE_biomass_spike_z_chir.csv") %>%
    mutate(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
)


taxa.colors =c("#66CCEE", "#AA3377", "#228833", "#CCBB44", "#4477AA")

q_labels <- as_labeller(c(
  "2" = "z[Chironomidae]==2",
  "4" = "z[Chironomidae]==4",
  "8" = "z[Chironomidae]==8"), default = label_parsed
)

label.names <- c(
  "Temperature", 
  "Temperature +\nSpike" , 
  "Temperature +\nHydropeaking", 
  "Temperature +\nSpike + Hydropeaking",
  "Temperature +\nHFE", 
  "Temperature +\nSpike + HFE",
  "Temperature +\nHydropeaking + HFE", 
  "Temperature +\nSpike + Hydropeaking + HFE"
)

multi_biomass_combo_chir$source <- ordered(multi_biomass_combo_chir$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE",
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"
))  
multi_biomass_combo_chir$taxa <- ordered(multi_biomass_combo_chir$taxa, levels = c("BAET", "HYOS", "CHIR", "GAMM", "NZMS"))
multi_biomass_combo_chir$q <- as.factor(multi_biomass_combo_chir$q)

multi_biomass_combo_chir <- filter(multi_biomass_combo_chir, q != "1")

multi_biomass_chir <- ggplot(data = multi_biomass_combo_chir, 
                        aes(x = as.factor(temperature), y = (mean_biomass), color = taxa)) +
  geom_errorbar(aes(ymin = (mean_biomass - sd_biomass), ymax = (mean_biomass + sd_biomass)), 
                linewidth = 0.75, width = 0.2, position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%")) +
  facet_wrap(. ~ source + q, nrow = 8,
             labeller = labeller(q = q_labels))+
  theme_bw() +
  scale_color_manual(name = " ", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Hydropsyche"), " spp.")),expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")),  expression(italic("P. antipodarum"))), values=taxa.colors) +
  xlab("Temperature") +
  labs(y= "Relative Biomass") +
  theme(text = element_text(size = 13), 
        axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), 
        legend.key = element_rect(fill = "transparent"))

ggsave("multi_z_chir.png", multi_biomass_chir, device = "png", width = 14.5, height = 12, units = "in", dpi = "retina")

# Load and combine all datasets, tagging them gamm
multi_biomass_combo_gamm <- bind_rows(
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_biomass_z_gamm.csv") %>%
    mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_biomass_spike_z_gamm.csv") %>%
    mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_biomass_hyd_z_gamm.csv") %>%
    mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_hyd_biomass_spike_z_gamm.csv") %>%
    mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_biomass_hfe_z_gamm.csv") %>%
    mutate(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_HFE_biomass_spike_z_gamm.csv") %>%
    mutate(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_hyd_HFE_biomass_z_gamm.csv") %>%
    mutate(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_hyd_HFE_biomass_spike_z_gamm.csv") %>%
    mutate(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
)

q_labels <- as_labeller(c(
  "2" = "z[G.lacustris]==2",
  "4" = "z[G.lacustris]==4",
  "8" = "z[G.lacustris]==8"), default = label_parsed
)

taxa.colors =c("#66CCEE", "#AA3377", "#228833", "#CCBB44", "#4477AA")

label.names <- c(
  "Temperature", 
  "Temperature +\nSpike" , 
  "Temperature +\nHydropeaking", 
  "Temperature +\nSpike + Hydropeaking",
  "Temperature +\nHFE", 
  "Temperature +\nSpike + HFE",
  "Temperature +\nHydropeaking + HFE", 
  "Temperature +\nSpike + Hydropeaking + HFE"
)

multi_biomass_combo_gamm$source <- ordered(multi_biomass_combo_gamm$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE",
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"
))  
multi_biomass_combo_gamm$taxa <- ordered(multi_biomass_combo_gamm$taxa, levels = c("BAET", "HYOS", "CHIR", "GAMM", "NZMS"))
multi_biomass_combo_gamm$q <- as.factor(multi_biomass_combo_gamm$q)

multi_biomass_combo_gamm <- filter(multi_biomass_combo_gamm, q != "1")

multi_biomass_gamm <- ggplot(data = multi_biomass_combo_gamm, 
                        aes(x = as.factor(temperature), y = (mean_biomass), color = taxa)) +
  geom_errorbar(aes(ymin = (mean_biomass - sd_biomass), ymax = (mean_biomass + sd_biomass)), 
                linewidth = 0.75, width = 0.2, position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%")) +
  facet_wrap(. ~ source + q, nrow = 8,
             labeller = labeller(q = q_labels))+
  theme_bw() +
  scale_color_manual(name = " ", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Hydropsyche"), " spp.")),expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")),  expression(italic("P. antipodarum"))), values=taxa.colors) +
  xlab("Temperature") +
  labs(y= "Relative Biomass") +
  theme(text = element_text(size = 13), 
        axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), 
        legend.key = element_rect(fill = "transparent"))
ggsave("multi_z_gamm.png", multi_biomass_gamm, device = "png", width = 14.5, height = 12, units = "in", dpi = "retina")

# Load and combine all datasets, tagging them nzms
multi_biomass_combo_nzms <- bind_rows(
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_biomass_z_nzms.csv") %>%
    mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_biomass_spike_z_nzms.csv") %>%
    mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_biomass_hyd_z_nzms.csv") %>%
    mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_hyd_biomass_spike_z_nzms.csv") %>%
    mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_biomass_hfe_z_nzms.csv") %>%
    mutate(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_HFE_biomass_spike_z_nzms.csv") %>%
    mutate(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_hyd_HFE_biomass_z_nzms.csv") %>%
    mutate(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_hyd_HFE_biomass_spike_z_nzms.csv") %>%
    mutate(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
)


taxa.colors =c("#66CCEE", "#AA3377", "#228833", "#CCBB44", "#4477AA")

q_labels <- as_labeller(c(
  "1" = "z[P.antipodarum]==1",
  "2" = "z[P.antipodarum]==2",
  "4" = "z[P.antipodarum]==4",
  "8" = "z[P.antipodarum]==8"), default = label_parsed
)

label.names <- c(
  "Temperature", 
  "Temperature +\nSpike" , 
  "Temperature +\nHydropeaking", 
  "Temperature +\nSpike + Hydropeaking",
  "Temperature +\nHFE", 
  "Temperature +\nSpike + HFE",
  "Temperature +\nHydropeaking + HFE", 
  "Temperature +\nSpike + Hydropeaking + HFE"
)

multi_biomass_combo_nzms$source <- ordered(multi_biomass_combo_nzms$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE",
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"
))  
multi_biomass_combo_nzms$taxa <- ordered(multi_biomass_combo_nzms$taxa, levels = c("BAET", "HYOS", "CHIR", "GAMM", "NZMS"))
z_nzms <- multi_biomass_combo_nzms %>% 
  filter(source == "Temperature" | source == "Temperature & spike & hydropeaking & HFE") %>%
  filter(temperature == 1 | temperature == 1.5)
multi_z_nzms <- ggplot(data = z_nzms, 
                       aes(x = q, y = (mean_biomass), color = taxa)) +
  geom_errorbar(aes(ymin = (mean_biomass - sd_biomass), ymax = (mean_biomass + sd_biomass)), 
                linewidth = 0.75, width = 0.2, position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_line(linewidth = 0.8, position = position_dodge(0.5))+
  facet_wrap(. ~ source + temperature, nrow = 2, labeller = labeller(temperature = temp_labels))+
  theme_bw() +
  scale_color_manual(name = " ", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Hydropsyche"), " spp.")),expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")),  expression(italic("P. antipodarum"))), values=taxa.colors) +
  scale_x_continuous(breaks = c(1, 2, 4, 8))+
  xlab(expression(z[P.antipodarum]))+
  labs(y= "Relative Biomass") +
  theme(text = element_text(size = 13), 
        axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), 
        legend.key = element_rect(fill = "transparent"))

ggsave("z_plot_nzms.png", multi_z_nzms, device = "png", width = 9, height = 7, units = c("in"), dpi = "retina")

multi_biomass_combo_nzms$q <- as.factor(multi_biomass_combo_nzms$q)

#multi_biomass_combo_nzms <- filter(multi_biomass_combo_nzms, q != "1")

multi_biomass_nzms <- ggplot(data = multi_biomass_combo_nzms, 
                        aes(x = as.factor(temperature), y = (mean_biomass), color = taxa)) +
  geom_errorbar(aes(ymin = (mean_biomass - sd_biomass), ymax = (mean_biomass + sd_biomass)), 
                linewidth = 0.75, width = 0.2, position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%")) +
  facet_wrap(. ~ source + q, nrow = 8,
             labeller = labeller(q = q_labels))+
  theme_bw() +
  scale_color_manual(name = " ", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Hydropsyche"), " spp.")),expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")),  expression(italic("P. antipodarum"))), values=taxa.colors) +
  xlab("Temperature") +
  labs(y= "Relative Biomass") +
  theme(text = element_text(size = 13), 
        axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), 
        legend.key = element_rect(fill = "transparent"))

ggsave("multi_z_nzms.png", multi_biomass_nzms, device = "png", width = 14.5, height = 12, units = "in", dpi = "retina")


# Load and combine all datasets, tagging them
# now when we alter all zis 
multi_biomass_combo <- bind_rows(
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_biomass_z.csv") %>%
    mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_biomass_spike_z.csv") %>%
    mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_biomass_hyd_z.csv") %>%
    mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_hyd_biomass_spike_z.csv") %>%
    mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_biomass_hfe_z.csv") %>%
    mutate(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_HFE_biomass_spike_z.csv") %>%
    mutate(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_hyd_HFE_biomass_z.csv") %>%
    mutate(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
  read_csv("ColoradoRiverInvertsMultiTaxa/Data/Multispp_temp_hyd_HFE_biomass_spike_z.csv") %>%
    mutate(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
)


taxa.colors =c("#66CCEE", "#AA3377", "#228833", "#CCBB44", "#4477AA")

temp_labels <- as_labeller(c(
  "1" = "Baseline",
  "1.5" = "+5~degree*C"
), default = label_parsed)


q_labels <- as_labeller(c(
  "1" = "z[i]==1",
  "2" = "z[i]==2",
  "4" = "z[i]==4",
  "8" = "z[i]==8"), default = label_parsed
)
label.names <- c(
  "Temperature", 
  "Temperature +\nSpike" , 
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

z_df <- multi_biomass_combo %>% 
  filter(source == "Temperature" | source == "Temperature & spike & hydropeaking & HFE") %>%
  filter(temperature == 1 | temperature == 1.5)

multi_z <- ggplot(data = z_df, 
                        aes(x = q, y = (mean_biomass), color = taxa)) +
  geom_errorbar(aes(ymin = (mean_biomass - sd_biomass), ymax = (mean_biomass + sd_biomass)), 
                linewidth = 0.75, width = 0.2, position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_line(linewidth = 0.8, position = position_dodge(0.5))+
  facet_wrap(. ~ source + temperature, nrow = 2, labeller = labeller(temperature = temp_labels))+
  theme_bw() +
  scale_color_manual(name = " ", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Hydropsyche"), " spp.")),expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")),  expression(italic("P. antipodarum"))), values=taxa.colors) +
  scale_x_continuous(breaks = c(1, 2, 4, 8))+
  xlab(expression(z[i]))+
  labs(y= "Relative Biomass") +
  theme(text = element_text(size = 13), 
        axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), 
        legend.key = element_rect(fill = "transparent"))

ggsave("z_plot.png", multi_z, device = "png", width = 9, height = 7, units = c("in"), dpi = "retina")

multi_biomass_combo$q <- as.factor(multi_biomass_combo$q)
multi_biomass_1s <- filter(multi_biomass_combo, q == "1")


multi_biomass_combo <- subset(multi_biomass_combo, q == "1" | q == "2")


multi_biomass <- ggplot(data = multi_biomass_combo, 
                        aes(x = as.factor(temperature), y = (mean_biomass), color = taxa)) +
  geom_errorbar(aes(ymin = (mean_biomass - sd_biomass), ymax = (mean_biomass + sd_biomass)), 
                linewidth = 0.75, width = 0.2, position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%")) +
  facet_wrap(. ~ source * q, ncol = 4,
             labeller = labeller(q = q_labels))+
  #facet_grid(q ~ source, labeller = labeller(q = q_labels))+
  theme_bw() +
  scale_color_manual(name = " ", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Hydropsyche"), " spp.")),expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")),  expression(italic("P. antipodarum"))), values=taxa.colors) +
  xlab("Temperature") +
  labs(y= "Relative Biomass") +
  theme(text = element_text(size = 13), 
        axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), 
        legend.key = element_rect(fill = "transparent"))

multi_biomass1 <- ggplot(data = multi_biomass_1s, 
                        aes(x = as.factor(temperature), y = (mean_biomass), color = taxa)) +
  geom_errorbar(aes(ymin = (mean_biomass - sd_biomass), ymax = (mean_biomass + sd_biomass)), 
                linewidth = 0.75, width = 0.2, position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%")) +
  facet_wrap(. ~ source + q, nrow = 4,
             labeller = labeller(q = q_labels))+
  theme_bw() +
  scale_color_manual(name = " ", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Hydropsyche"), " spp.")),expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")),  expression(italic("P. antipodarum"))), values=taxa.colors) +
  xlab("Temperature") +
  labs(y= "Relative Biomass") +
  theme(text = element_text(size = 13), 
        axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), 
        legend.key = element_rect(fill = "transparent"))

ggsave("multi_z.png", multi_biomass, device = "png", width = 15, height = 12, units = "in", dpi = "retina")

ggsave("multi_z1.png", multi_biomass1, device = "png", width = 10, height = 10, units = "in", dpi = "retina")

