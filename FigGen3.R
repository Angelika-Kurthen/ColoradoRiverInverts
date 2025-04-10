###########
## Figure generation for 3rd manuscript
##########
library(ggpubr)
library(patchwork)
library(tidyverse)
#install.packages("svglite")
library(svglite)
library(igraph)
library(ggraph)
library(dplyr)
library(broom) 
library(readr)
library(lm.beta)
library(effectsize)

# 
# #source("Multispp.R")
# source("MultisppHydropeak.R")
# 
# # make sure in the correct format
# hydropeak.abunds$`rep(hydropeak[hyd], times = 5)` <- as.numeric(hydropeak.abunds$`rep(hydropeak[hyd], times = 5)`)
# hydropeak.abunds$avg.abund<- as.numeric(hydropeak.abunds$avg.abund)
# hydropeak.abunds$taxa <- as.factor(hydropeak.abunds$taxa)
# 
# #plot
# hyd_abund_multi <- ggplot(data = hydropeak.abunds, aes(`rep(hydropeak[hyd], times = 5)`, avg.abund, group = taxa, color = taxa))+ 
#   geom_line(linewidth = 1, alpha = 0.8)+
#   # geom_ribbon(data = Hydro_Abund, aes(ymin = scale(means - sd),
#   #                 ymax = scale(means + sd), fill = V4),
#   #                        alpha = .15,
#   #                        show.legend = T) +
#   scale_color_manual(name = "Taxon", labels=c("Hydropsyche spp.", "Baetidae spp.", "P. anitpodarum", "Chironomidae spp.", "G. lacustris"), values=c("#AA3377","#66CCEE","#4477AA","#228833", "#CCBB44"))+
#   geom_vline(aes(xintercept = 0.01), linetype = "dotted", 
#              linewidth=1)+
#   geom_vline(aes(xintercept = 0.17 ), linetype="dashed", 
#              linewidth=1)+
#   geom_vline(aes(xintercept = 0.55 ), linetype = "dotdash", 
#              linewidth=1)+
#   scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
#   theme_bw()+
#   xlab("Hydropeaking Intensity")+
#   ylab("Relative Abundance")+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# 
# # make sure in correct format
# hydropeak.biomass$`rep(hydropeak[hyd], times = 5)` <- as.numeric(hydropeak.biomass$`rep(hydropeak[hyd]`)
# hydropeak.biomass$avg.biomass <- as.numeric(hydropeak.biomass$avg.biomass)
# hydropeak.biomass$taxa <- as.factor(hydropeak.biomass$taxa) 
# #plot
# hyd_ts_multi <- ggplot(data = hydropeak.biomass, aes(`rep(hydropeak[hyd], times = 5)`,avg.biomass , color = taxa)) + 
#   # geom_ribbon(aes(ymin = sizemeans - sizesd,
#   #                   ymax = sizemeans + sizesd),
#   #               colour = V3,
#   #               alpha = .15,
#   #               show.legend = T) +
#   geom_line(linewidth = 1, alpha = 0.8)+
#   scale_color_manual(name = "Taxon", labels=c("Hydropsyche spp.", "Baetidae spp.", "P. anitpodarum", "Chironomidae spp.", "G. lacustris"), values=c("#AA3377","#66CCEE","#4477AA","#228833", "#CCBB44"))+
#   geom_vline(aes(xintercept = 0.01), linewidth = 1, linetype = "dotted")+
#   geom_vline(aes(xintercept = 0.17 ), linewidth = 1, linetype="dashed")+ 
#   geom_vline(aes(xintercept = 0.55 ), linewidth = 1, linetype = "dotdash")+
#   theme_bw()+
#   xlab("Hydropeaking Intensity")+
#   ylab("Relative Biomass (mg)")+
#   scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# 
# #make sure in correct format
# # can also make this emergent by only looking at HYOS, CHIR, and BAET
# hydropeak.avg.annual$`rep(hydropeak[hyd], times = 5)` <- as.numeric(hydropeak.avg.annual$`rep(hydropeak[hyd], times = 5)`)
# hydropeak.avg.annual$taxa <- as.factor(hydropeak.avg.annual$taxa)
# hydropeak.avg.annual$`means.s3.biomass$mean.S3.biomass` <- as.numeric(hydropeak.avg.annual$`means.s3.biomass$mean.S3.biomass`)
# hyd_yr_multi <- ggplot(data = hydropeak.avg.annual, aes(`rep(hydropeak[hyd], times = 5)`, log(`means.s3.biomass$mean.S3.biomass`), group = taxa, color = taxa)) + 
#   # geom_ribbon(aes(ymin = sizemeans - sizesd,
#   #                   ymax = sizemeans + sizesd),
#   #               colour = V3,
#   #               alpha = .15,
#   #               show.legend = T) +
#   geom_line(linewidth = 1, alpha = 0.8)+
#   scale_color_manual(name = "Taxon", labels=c("Hydropsyche spp.", "Baetidae spp.", "P. anitpodarum", "Chironomidae spp.", "G. lacustris"), values=c("#AA3377","#66CCEE","#4477AA","#228833", "#CCBB44"))+
#   geom_vline(aes(xintercept = 0.01), linewidth = 1, linetype = "dotted")+
#   geom_vline(aes(xintercept = 0.17 ), linewidth = 1, linetype="dashed")+ 
#   geom_vline(aes(xintercept = 0.55 ), linewidth = 1, linetype = "dotdash")+
#   theme_bw()+
#   scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
#   xlab("Hydropeaking Intensity")+
#   ylab("Log Annual S3 Biomass (mg)")+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# 
# ggarrange(hyd_abund_multi, hyd_ts_multi , hyd_yr_multi,
#           labels = c("a", "b", "c"),
#           ncol = 3, nrow = 1, common.legend = T)
# 

# increasing HFEs: 
taxa.colors =c("#66CCEE", "#AA3377", "#228833", "#CCBB44", "#4477AA")
Multispp_temp_HFE_abund <- read_csv("Multispp_temp_HFE_abund.csv")
# note that I mislabeled HFE magnitude as temperature
Multispp_temp_HFE_abund$taxa <- as.factor(Multispp_temp_HFE_abund$taxa )
Multispp_temp_HFE_abund$taxa <- ordered(Multispp_temp_HFE_abund$taxa, levels = c("BAET", "HYOS", "CHIR", "GAMM", "NZMS"))
Multispp_temp_HFE_abund$temperature <- as.factor(Multispp_temp_HFE_abund$temperature)
Multispp_HFE_summary <- Multispp_temp_HFE_abund %>%
    group_by(temperature, taxa) %>%
    summarise(mean_abund = mean(abundance, na.rm = TRUE),
              sd_abund = sd(abundance, na.rm = TRUE))


hfe_abund <- ggplot(data = Multispp_HFE_summary, 
                      aes(x = as.factor(temperature), y = mean_abund, color = taxa, group = taxa)) +
  geom_errorbar(aes(ymin = mean_abund - sd_abund, ymax = mean_abund + sd_abund), 
                linewidth = 0.75, width = 0.2, position = position_dodge(0.5)) +
  geom_path(linewidth = 0.75, alpha = 0.7)+
  geom_point(size = 3, position = position_dodge(0.5)) +
  theme_bw() +
  ylim(c(0, 1))+
  scale_color_manual(name = " ", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Hydropsyche"), " spp.")),expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")),  expression(italic("P. antipodarum"))), values=taxa.colors) +
  xlab("Spring HFE Magnitude") +
  labs(y= "Relative Abundance") +
  theme(text = element_text(size = 13), 
        axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), 
        legend.key = element_rect(fill = "transparent"))

abund.hfe <- Multispp_temp_HFE_abund %>%
  group_by(taxa) %>%
  group_modify(~ {
    model <- eta_squared(aov(abundance~ temperature, data = .x))
  })

Multispp_temp_HFE_biomass <- read_csv("Multispp_temp_HFE_biomass.csv")
# note that I mislabeled HFE magnitude as temperature
Multispp_temp_HFE_biomass$taxa <- as.factor(Multispp_temp_HFE_biomass$taxa)
Multispp_temp_HFE_biomass$taxa <- ordered(Multispp_temp_HFE_biomass$taxa, levels = c("BAET", "HYOS", "CHIR", "GAMM", "NZMS"))
Multispp_temp_HFE_biomass$temperature <- as.factor(Multispp_temp_HFE_biomass$temperature)


Multispp_HFE_bio_summary <- Multispp_temp_HFE_biomass %>%
  group_by(temperature, taxa) %>%
  summarise(mean_biomass = mean(biomass, na.rm = TRUE),
            sd_biomass = sd(biomass, na.rm = TRUE), .groups = "drop")

hfe_biomass <- ggplot(data = Multispp_HFE_bio_summary, 
                    aes(x = as.factor(temperature), y = mean_biomass, color = taxa, group = taxa)) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), 
                linewidth = 0.75, width = 0.2, position = position_dodge(0.5)) +
  geom_path(linewidth = 0.75, alpha = 0.7)+
  geom_point(size = 3, position = position_dodge(0.5)) +
  theme_bw() +
  ylim(c(0,1))+
  scale_color_manual(name = " ", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Hydropsyche"), " spp.")),expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")),  expression(italic("P. antipodarum"))), values=taxa.colors) +
  xlab("Spring HFE Magnitude") +
  labs(y= "Relative Biomass") +
  theme(text = element_text(size = 13), 
        axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), 
        legend.key = element_rect(fill = "transparent"))

biomass.hfe <- Multispp_temp_HFE_biomass %>%
  group_by(taxa) %>%
  group_modify(~ {
    model <- eta_squared(aov(biomass ~ temperature, data = .x))
  })

multispp_hfe <- ggarrange(hfe_abund, hfe_biomass, 
          labels = c("a", "b"), 
        common.legend = T)
ggsave("multispp_hfe.png", plot = multispp_hfe, device = "png", width = 8, height = 6,  dpi = "retina")


# hyd
Multispp_temp_hyd_abund <- read_csv("Multispp_temp_hyd_abund.csv")
# note that I mislabeled HFE magnitude as temperature
Multispp_temp_hyd_abund$taxa <- as.factor(Multispp_temp_hyd_abund$taxa )
Multispp_temp_hyd_abund$taxa <- ordered(Multispp_temp_hyd_abund$taxa, levels = c("BAET", "HYOS", "CHIR", "GAMM", "NZMS"))
Multispp_temp_hyd_abund$hi <- as.factor(Multispp_temp_hyd_abund$hi)
Multispp_hyd_summary <- Multispp_temp_hyd_abund %>%
  group_by(hi, taxa) %>%
  summarise(mean_abund = mean(abundance, na.rm = TRUE),
            sd_abund = sd(abundance, na.rm = TRUE))


hyd_abund <- ggplot(data = Multispp_hyd_summary, 
                    aes(x = as.factor(hi), y = mean_abund, color = taxa, group = taxa)) +
  geom_errorbar(aes(ymin = mean_abund - sd_abund, ymax = mean_abund + sd_abund), 
                linewidth = 0.75, width = 0.2, position = position_dodge(0.5)) +
  geom_path(linewidth = 0.75, alpha = 0.7)+
  geom_point(size = 3, position = position_dodge(0.5)) +
  theme_bw() +
  ylim(c(0, 1))+
  scale_color_manual(name = " ", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Hydropsyche"), " spp.")),expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")),  expression(italic("P. antipodarum"))), values=taxa.colors) +
  xlab("Hydropeaking Index") +
  labs(y= "Relative Abundance") +
  theme(text = element_text(size = 13), 
        axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), 
        legend.key = element_rect(fill = "transparent"))

abund.hyd <- Multispp_temp_hyd_abund %>%
  group_by(taxa) %>%
  group_modify(~ {
    model <- eta_squared(aov(abundance~ hi, data = .x))
  })

Multispp_temp_hyd_abundno1 <- filter(Multispp_temp_hyd_abund, hi != 1)
abund.hydno1 <- Multispp_temp_hyd_abundno1 %>%
  group_by(taxa) %>%
  group_modify(~ {
    model <- eta_squared(aov(abundance~ hi, data = .x))
  })

Multispp_temp_hyd_biomass <- read_csv("Multispp_temp_hyd_biomass.csv")
# note that I mislabeled hyd magnitude as temperature
Multispp_temp_hyd_biomass$taxa <- as.factor(Multispp_temp_hyd_biomass$taxa)
Multispp_temp_hyd_biomass$taxa <- ordered(Multispp_temp_hyd_biomass$taxa, levels = c("BAET", "HYOS", "CHIR", "GAMM", "NZMS"))
Multispp_temp_hyd_biomass$hi <- as.factor(Multispp_temp_hyd_biomass$hi)


Multispp_hyd_bio_summary <- Multispp_temp_hyd_biomass %>%
  group_by(hi, taxa) %>%
  summarise(mean_biomass = mean(biomass, na.rm = TRUE),
            sd_biomass = sd(biomass, na.rm = TRUE), .groups = "drop")

hyd_biomass <- ggplot(data = Multispp_hyd_bio_summary, 
                      aes(x = as.factor(hi), y = mean_biomass, color = taxa, group = taxa)) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), 
                linewidth = 0.75, width = 0.2, position = position_dodge(0.5)) +
  geom_path(linewidth = 0.75, alpha = 0.7)+
  geom_point(size = 3, position = position_dodge(0.5)) +
  theme_bw() +
  ylim(c(0,1))+
  scale_color_manual(name = " ", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Hydropsyche"), " spp.")),expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")),  expression(italic("P. antipodarum"))), values=taxa.colors) +
  xlab("Hydropeaking Index") +
  labs(y= "Relative Biomass") +
  theme(text = element_text(size = 13), 
        axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), 
        legend.key = element_rect(fill = "transparent"))

biomass.hyd <- Multispp_temp_hyd_biomass %>%
  group_by(taxa) %>%
  group_modify(~ {
    model <- eta_squared(aov(biomass ~ hi, data = .x))
  })

Multispp_temp_hyd_biomassno1 <- filter(Multispp_temp_hyd_biomass, hi != 1)
biomass.hydno1 <- Multispp_temp_hyd_biomassno1 %>%
  group_by(taxa) %>%
  group_modify(~ {
    model <- eta_squared(aov(biomass ~ hi, data = .x))
  })

multispp_hyd <- ggarrange(hyd_abund, hyd_biomass, 
                          labels = c("a", "b"), 
                          common.legend = T)
ggsave("multispp_hyd.png", plot = multispp_hyd, device = "png", width = 9, height = 6,  dpi = "retina")


# 8 different scenarios (proportions)

# read in the csv's
Multispp_temp_biomass <- read_csv("Multispp_temp_biomass1.2.csv")
Multispp_temp_abund <- read_csv("Multispp_temp_abund1.2.csv")
Multispp_temp_abund_spike <- read_csv("Multispp_temp_abund_spike1.2.csv")
Multispp_temp_biomass_spike <- read_csv("Multispp_temp_biomass_spike1.2.csv")
Multispp_temp_hyd_biomass <- read_csv("Multispp_temp_biomass_hyd.csv")
Multispp_temp_hyd_biomass_spike <- read_csv("Multispp_temp_hyd_biomass_spike.csv")
Multispp_temp_hyd_HFE_biomass_spike <- read_csv("Multispp_temp_hyd_HFE_biomass_spike.csv")
Multispp_temp_hyd_HFE_biomass <- read_csv("Multispp_temp_hyd_HFE_biomass.csv")
Multispp_temp_hyd_HFE_abund_spike <- read_csv("Multispp_temp_hyd_HFE_abund_spike.csv")
Multispp_temp_HFE_abund_spike <- read_csv("Multispp_temp_HFE_abund_spike1.2.csv")
Multispp_temp_HFE_abund <- read_csv("Multispp_temp_HFE_abund1.2.csv")
Multispp_temp_HFE_biomass <- read_csv("Multispp_temp_HFE_biomass1.2.csv")
Multispp_temp_hyd_abund <- read_csv("Multispp_temp_abund_hyd.csv")
Multispp_temp_HFE_biomass_spike <- read_csv("Multispp_temp_HFE_biomass_spike1.2.csv")
Multispp_temp_hyd_abund_spike <- read_csv("Multispp_temp_hyd_abund_spike.csv")
Multispp_temp_hyd_HFE_abund <- read_csv("Multispp_temp_hyd_HFE_abund.csv")
Multispp_temp_S3 <- read_csv("MultisppS3_temp_biomass.csv")
Multispp_temp_S3_spike <- read_csv("MultisppS3_temp_biomass_spike.csv")

multi_abund_combo <- bind_rows(
  Multispp_temp_abund %>%  mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  Multispp_temp_abund_spike %>% mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  Multispp_temp_hyd_abund %>%  mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  Multispp_temp_hyd_abund_spike %>% mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  Multispp_temp_HFE_abund %>% mutate(source = "Temperature & HFE",  temp_spike = 0, hydropeaking = 0, HFE = 1),
  Multispp_temp_HFE_abund_spike %>% mutate(source = "Temperature & spike & HFE",  temp_spike = 1, hydropeaking = 0, HFE = 1),
  Multispp_temp_hyd_HFE_abund %>% mutate(source = "Temperature & hydropeaking & HFE",  temp_spike = 0, hydropeaking = 1, HFE = 1),
  Multispp_temp_hyd_HFE_abund_spike %>% mutate(source = "Temperature & spike & hydropeaking & HFE",  temp_spike = 1, hydropeaking = 1, HFE = 1)
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


multi_abund_summary <- multi_abund_combo %>%
  group_by(temperature, taxa, source) %>%
  summarise(mean_abund = mean(abundance, na.rm = TRUE),
            sd_abund = sd(abundance, na.rm = TRUE), .groups = "drop")

head(multi_abund_summary)

multi_abund <- ggplot(data = multi_abund_summary, 
                     aes(x = as.factor(temperature), y = log(mean_abund), color = taxa)) +
  geom_errorbar(aes(ymin = log(mean_abund - sd_abund), ymax = log(mean_abund + sd_abund)), 
                linewidth = 0.75, width = 0.2, position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%")) +
  facet_wrap(.~source, nrow = 4) +
  theme_bw() +
  scale_color_manual(name = " ", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Hydropsyche"), " spp.")),expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")),  expression(italic("P. antipodarum"))), values=taxa.colors) +
  xlab("Temperature") +
  labs(y= "Log Relative Abundance") +
  theme(text = element_text(size = 13), 
        axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), 
        legend.key = element_rect(fill = "transparent"))

multi_biomass_combo <- bind_rows(
  Multispp_temp_biomass %>%  mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  Multispp_temp_biomass_spike %>% mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  #Multispp_temp_hyd_biomass %>%  mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  #Multispp_temp_hyd_biomass_spike %>% mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  Multispp_temp_HFE_biomass %>% mutate(source = "Temperature & HFE",  temp_spike = 0, hydropeaking = 0, HFE = 1),
  Multispp_temp_HFE_biomass_spike %>% mutate(source = "Temperature & spike & HFE",  temp_spike = 1, hydropeaking = 0, HFE = 1)
  #Multispp_temp_hyd_HFE_biomass %>% mutate(source = "Temperature & hydropeaking & HFE",  temp_spike = 0, hydropeaking = 1, HFE = 1),
  #Multispp_temp_hyd_HFE_biomass_spike %>% mutate(source = "Temperature & spike & hydropeaking & HFE",  temp_spike = 1, hydropeaking = 1, HFE = 1)
 )


label.names <- c(
  "Temperature", 
  "Temperature +\nSpike", 
  #"Temperature +\nHydropeaking", 
  #"Temperature +\nSpike + Hydropeaking",
  "Temperature +\nHFE", 
  "Temperature +\nSpike + HFE"
  #"Temperature +\nHydropeaking + HFE", 
  #"Temperature +\nSpike + Hydropeaking + HFE"
)

multi_biomass_combo$source <- ordered(multi_biomass_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  #"Temperature & hydropeaking", 
  #"Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE" 
  #"Temperature & hydropeaking & HFE", 
  #"Temperature & spike & hydropeaking & HFE"
  ))  
multi_biomass_combo$taxa <- ordered(multi_biomass_combo$taxa, levels = c("BAET", "HYOS", "CHIR", "GAMM", "NZMS"))

write.csv2(multi_abund_combo, file = "multi_abund_combo.csv")
write.csv2(multi_biomass_combo, file = "multi_biomass_combo.csv")


multi_biomass_summary <- multi_biomass_combo %>%
    group_by(temperature, taxa, source) %>%
    summarise(mean_biomass = mean(biomass, na.rm = TRUE),
              sd_biomass = sd(biomass, na.rm = TRUE), .groups = "drop")
  
head(multi_biomass_summary)
multi_biomass <- ggplot(data = multi_biomass_summary, 
                      aes(x = as.factor(temperature), y = log(mean_biomass), color = taxa)) +
  geom_errorbar(aes(ymin = log(mean_biomass - sd_biomass), ymax = log(mean_biomass + sd_biomass)), 
                linewidth = 0.75, width = 0.2, position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%")) +
  facet_wrap(.~source, nrow= 4) +
  theme_bw() +
  scale_color_manual(name = " ", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Hydropsyche"), " spp.")),expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")),  expression(italic("P. antipodarum"))), values=taxa.colors) +
  xlab("Temperature") +
  labs(y= "Relative Biomass") +
  theme(text = element_text(size = 13), 
        axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), 
        legend.key = element_rect(fill = "transparent"))



###
# 8 different scenarios zi 

Multispp_temp_biomass <- read_csv("Multispp_temp_biomass_z.csv")
Multispp_temp_abund <- read_csv("Multispp_temp_abund_z.csv")
Multispp_temp_hyd_HFE_abund <- read_csv("Multispp_temp_hyd_HFE_abund_z.csv")
Multispp_temp_hyd_HFE_biomass <- read_csv("Multispp_temp_hyd_HFE_biomass_z.csv")
Multispp_temp_abund_hfe <- read_csv("Multispp_temp_abund_hfe_z.csv")
Multispp_temp_biomass_hfe <- read_csv("Multispp_temp_biomass_hfe_z.csv")

multi_biomass_combo <- bind_rows(
  Multispp_temp_biomass %>%  mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  #Multispp_temp_biomass_spike %>% mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  #Multispp_temp_hyd_biomass %>%  mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  #Multispp_temp_hyd_biomass_spike %>% mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  Multispp_temp_biomass_hfe %>% mutate(source = "Temperature & HFE",  temp_spike = 0, hydropeaking = 0, HFE = 1),
 # Multispp_temp_HFE_biomass_spike %>% mutate(source = "Temperature & spike & HFE",  temp_spike = 1, hydropeaking = 0, HFE = 1),
 Multispp_temp_hyd_HFE_biomass %>% mutate(source = "Temperature & hydropeaking & HFE",  temp_spike = 0, hydropeaking = 1, HFE = 1),
 # Multispp_temp_hyd_HFE_biomass_spike %>% mutate(source = "Temperature & spike & hydropeaking & HFE",  temp_spike = 1, hydropeaking = 1, HFE = 1)
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
multi_biomass_combo$q <- as.factor(multi_biomass_combo$q)
multi_biomass_combo <- multi_biomass_combo %>% 
  filter(q == 1 | q == 2 | q == 4)

multi_biomass <- ggplot(data = multi_biomass_combo, 
                      aes(x = as.factor(temperature), y = (mean_biomass), color = taxa)) +
  geom_errorbar(aes(ymin = (mean_biomass - sd_biomass), ymax = (mean_biomass + sd_biomass)), 
                linewidth = 0.75, width = 0.2, position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%")) +
  facet_wrap(.~source + q, nrow = 4) +
  theme_bw() +
  scale_color_manual(name = " ", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Hydropsyche"), " spp.")),expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")),  expression(italic("P. antipodarum"))), values=taxa.colors) +
  xlab("Temperature") +
  labs(y= "Relative Biomass") +
  theme(text = element_text(size = 13), 
        axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), 
        legend.key = element_rect(fill = "transparent"))

# #----------------------------------------------------------------
# Multispp_mh_sens_temps <- read_csv("Multispp_mh_sens_temps.csv")
# Multispp_mh_sens_temps$TemperatureFactor <- as.factor(Multispp_mh_sens_temps$TemperatureFactor)
# sensgrid4 <- ggplot(Multispp_mh_sens_temps, aes(x = SensitivityIncrement, y = Biomass, color = StageGroup)) +
#   geom_point(size = 1) +  # Alternatively, you can add geom_line() for smoother trends.
#   facet_grid(StageGroup ~ Parameter + TemperatureFactor, scales = "free_y") +  # Facet by both StageGroup & Parameter
#   theme_bw() +
#   labs(
#     x = "Sensitivity Increment",
#     y = "Biomass",
#     title = "Sensitivity Analysis by StageGroup and Parameter"
#   ) +
#   theme(
#     strip.text = element_text(size = 8, face = "bold"),  # Adjust facet label size
#     axis.text.x = element_text(angle = 90, hjust = 1),
#     legend.position = "none"  
#   )
# 
# ggsave("sensgrid4.png", sensgrid4, device = "png", width = 8, height = 11.5, units= c("in"), dpi = "retina")
# # Fit a linear model for each StageGroup × Parameter combination
# interaction_results_mh <- Multispp_mh_sens_temps %>%
#   group_by(StageGroup, Parameter, TemperatureFactor) %>%
#   do({
#     model <- lm(Biomass ~ SensitivityIncrement, data = .)  # Fit linear model
#     model_summary <- summary(model)
#     
#     # Extract key values
#     tibble(
#       Slope = coef(model)[2],  # Extract slope
#       P_value = coef(summary(model))["SensitivityIncrement", "Pr(>|t|)"],  # Get p-value for slope
#       R2 = model_summary$r.squared  # Get R
#     )
#   }) %>%
#   ungroup() %>%
#   # Filter results based on criteria
#   filter(R2 >= 0.13, P_value < 0.05) #moderate plus?
# 
# # Rename parameters: Change "G1" to "S1", "G2" to "S2"
# interaction_results_mh <- interaction_results_mh %>%
#   mutate(Parameter = case_when(
#     grepl("hydro.mort_NZMS", Parameter) ~ "NZMS",
#     grepl("hydro.mort_GAMM", Parameter) ~ "GAMM",
#     grepl("hydro.mort_CHIR", Parameter) ~ "CHIR",
#     grepl("hydro.mort_HYOS", Parameter) ~ "HYOS",
#     grepl("hydro.mort_BAET", Parameter) ~ "BAET",
#     TRUE ~ Parameter
#   ))
# 
# interaction_results_mh <- interaction_results_mh %>%
#   mutate(StageGroup = case_when(
#     grepl("S3_NZMS", StageGroup) ~ "NZMS",
#     grepl("S3_GAMM", StageGroup) ~ "GAMM",
#     grepl("S3_CHIR", StageGroup) ~ "CHIR",
#     grepl("S3_HYOS", StageGroup) ~ "HYOS",
#     grepl("S3_BAET", StageGroup) ~ "BAET",
#     grepl("S2_NZMS", StageGroup) ~ "NZMS",
#     grepl("S2_GAMM", StageGroup) ~ "GAMM",
#     grepl("S2_CHIR", StageGroup) ~ "CHIR",
#     grepl("S2_HYOS", StageGroup) ~ "HYOS",
#     grepl("S2_BAET", StageGroup) ~ "BAET",
#     grepl("S1_NZMS", StageGroup) ~ "NZMS",
#     grepl("S1_GAMM", StageGroup) ~ "GAMM",
#     grepl("S1_CHIR", StageGroup) ~ "CHIR",
#     grepl("S1_HYOS", StageGroup) ~ "HYOS",
#     grepl("S1_BAET", StageGroup) ~ "BAET",
#     TRUE ~ StageGroup
#   ))
# 
# 
# # Extract the taxon type from parameter names
# interaction_results_mh <- interaction_results_mh %>%
#   mutate(Taxon = case_when(
#     grepl("HYOS", Parameter) ~ "HYOS",
#     grepl("GAMM", Parameter) ~ "GAMM",
#     grepl("NZMS", Parameter) ~ "NZMS",
#     grepl("CHIR", Parameter) ~ "CHIR",
#     grepl("BAET", Parameter) ~ "BAET",
#     TRUE ~ "Other"
#   ))
# 
# # Create an edge list
# edges_mh <- interaction_results_mh%>%
#   select(Parameter, StageGroup, Slope, TemperatureFactor, Taxon)
# 
# colnames(edges_mh) <- c("from", "to", "weight", "TemperatureFactor", "taxon")  # Rename for igraph format
# 
# # Create nodes from unique 'from' and 'to' values
# # Create nodes for both networks based on the union of nodes from both graphs
# all_nodes <- unique(c(edges_mh$from, edges_mh$to, edges_mh$from, edges_mh$to))
# 
# # Create nodes with Taxon and type for both graphs
# nodes_all <- tibble(name = all_nodes) %>%
#   mutate(Taxon = str_extract(name, "(HYOS|GAMM|NZMS|CHIR|BAET)"),
#          type = ifelse(name %in% edges_mh$from | name %in% edges_mh$from, "Source", "Target"))
# 
# 
# # Create igraph object
# network_graph_mh <- graph_from_data_frame(d = edges_mh, vertices = nodes_all, directed = TRUE)
# 
# 
# # Assign a layout manually for better visualization
# network_layout_mh <- create_layout(network_graph_mh, layout = "circle")
# 
# # Define taxon colors
# taxon_colors <- c(
#   "HYOS" = "#AA3377",
#   "GAMM" = "#CCBB44",
#   "NZMS" = "#4477AA",
#   "CHIR" = "#228833",
#   "BAET" = "#66CCEE"
# )
# 
# # Plot the network
# hyd_mort_net <- ggraph(network_layout_mh, layout = "circle") +
#   geom_edge_arc(aes(edge_width = abs(weight), color = taxon),
#                 arrow = arrow(length = unit(4, "mm")),
#                 end_cap = circle(4, "mm"),
#                 start_cap = circle(4, "mm"),
#                 alpha = 0.8,  
#                 curvature = 0.3, 
#                 show.legend = FALSE) +
#   geom_node_point(size = 8, aes(color = Taxon)) +
#   geom_node_text(aes(label = " "), repel = F, size = 4, fontface = "bold") +
#   scale_edge_width(range = c(0.5, 3)) +
#   scale_edge_color_manual(values = taxon_colors) +
#   scale_color_manual(values = taxon_colors, 
#                      labels=c(expression(paste(italic("Baetidae")," spp.")), 
#                               expression(paste(italic("Chironomidae"), " spp.")), 
#                               expression(italic("G. lacustris")), 
#                               expression(paste(italic("Hydropsyche"), " spp.")), 
#                               expression(italic("P. antipodarum")))) +
#   theme(
#     panel.background = element_rect(fill = "white", color = "black"), 
#     panel.grid = element_blank(),  
#     strip.background = element_rect(fill = "white", color = "black"),  
#     strip.text = element_text(color = "black"),  
#     legend.key = element_rect(fill = "white", color = NA),
#     plot.title = element_text(hjust = 0.5)
#   ) +
#   labs(title = "Hydropeaking Mortality Sensitivity Analysis", color = "Taxon", edge_width = "Slope")+
#   facet_wrap(~TemperatureFactor, labeller = labeller(TemperatureFactor = c(
#     "1" = "Baseline", 
#     "1.1" = "+1°C", 
#     "1.2" = "+2°C", 
#     "1.5" = "+5°C"
#   )))
# #-------------------------------------------------------------------------------
# # sensitivity network analysis for flood mortality
# 
# nodes_all <- tibble(name = all_nodes) %>%
#   mutate(Taxon = str_extract(name, "(HYOS|GAMM|NZMS|CHIR|BAET)"),
#          type = ifelse(name %in% edges_mh$from | name %in% edges_mh$from, "Source", "Target"))
# 
# 
# Multispp_mF_temps <- read_csv("Multispp_mF_sens_temps.csv")
# 
# sensgrid2 <- ggplot(Multispp_mF_temps, aes(x = SensitivityIncrement, y = Biomass, color = StageGroup)) +
#   geom_point(size = 1) +  # Alternatively, you can add geom_line() for smoother trends.
#   facet_grid(StageGroup ~ Parameter + TemperatureFactor, scales = "free_y") +  # Facet by both StageGroup & Parameter
#   theme_bw() +
#   labs(
#     x = "Sensitivity Increment",
#     y = "Biomass",
#     title = "Sensitivity Analysis by StageGroup and Parameter"
#   ) +
#   theme(
#     strip.text = element_text(size = 8, face = "bold"),  # Adjust facet label size
#     axis.text.x = element_text(angle = 90, hjust = 1),
#     legend.position = "none"  # Hide legend (since color is redundant with facet_grid)
#   )
# 
# ggsave("sensgrid2.png", sensgrid2, device = "png", width = 8, height = 11.5, units= c("in"), dpi = "retina")
# # Fit a linear model for each StageGroup × Parameter combination
# interaction_results_mf <- Multispp_mF_temps %>%
#   group_by(StageGroup, Parameter, TemperatureFactor) %>%
#   do({
#     model <- lm(Biomass ~ SensitivityIncrement, data = .)  # Fit linear model
#     model_summary <- summary(model)
#     
#     # Extract key values
#     tibble(
#       Slope = coef(model)[2],  # Extract slope
#       P_value = coef(summary(model))["SensitivityIncrement", "Pr(>|t|)"],  # Get p-value for slope
#       R2 = model_summary$r.squared  # Get R
#     )
#   }) %>%
#   ungroup() %>%
#   # Filter results based on criteria
#   filter(R2 >= 0.13, P_value < 0.05) #moderate plus?
# 
# # Rename parameters: Change "G1" to "S1", "G2" to "S2"
# interaction_results_mf <- interaction_results_mf %>%
#   mutate(Parameter = case_when(
#     grepl("flood.mort_NZMS", Parameter) ~ "NZMS",
#     grepl("flood.mort_GAMM", Parameter) ~ "GAMM",
#     grepl("flood.mort_CHIR", Parameter) ~ "CHIR",
#     grepl("flood.mort_HYOS", Parameter) ~ "HYOS",
#     grepl("flood.mort_BAET", Parameter) ~ "BAET",
#     TRUE ~ Parameter
#   ))
# 
# interaction_results_mf <- interaction_results_mf %>%
#   mutate(StageGroup = case_when(
#     grepl("S3_NZMS", StageGroup) ~ "NZMS",
#     grepl("S3_GAMM", StageGroup) ~ "GAMM",
#     grepl("S3_CHIR", StageGroup) ~ "CHIR",
#     grepl("S3_HYOS", StageGroup) ~ "HYOS",
#     grepl("S3_BAET", StageGroup) ~ "BAET",
#     grepl("S2_NZMS", StageGroup) ~ "NZMS",
#     grepl("S2_GAMM", StageGroup) ~ "GAMM",
#     grepl("S2_CHIR", StageGroup) ~ "CHIR",
#     grepl("S2_HYOS", StageGroup) ~ "HYOS",
#     grepl("S2_BAET", StageGroup) ~ "BAET",
#     grepl("S1_NZMS", StageGroup) ~ "NZMS",
#     grepl("S1_GAMM", StageGroup) ~ "GAMM",
#     grepl("S1_CHIR", StageGroup) ~ "CHIR",
#     grepl("S1_HYOS", StageGroup) ~ "HYOS",
#     grepl("S1_BAET", StageGroup) ~ "BAET",
#     TRUE ~ StageGroup
#   ))
# 
# 
# # Extract the taxon type from parameter names
# interaction_results_mf <- interaction_results_mf %>%
#   mutate(Taxon = case_when(
#     grepl("HYOS", Parameter) ~ "HYOS",
#     grepl("GAMM", Parameter) ~ "GAMM",
#     grepl("NZMS", Parameter) ~ "NZMS",
#     grepl("CHIR", Parameter) ~ "CHIR",
#     grepl("BAET", Parameter) ~ "BAET",
#     TRUE ~ "Other"
#   ))
# 
# # Create an edge list
# edges_mf <- interaction_results_mf %>%
#   select(Parameter, StageGroup, Slope, TemperatureFactor, Taxon)
# 
# colnames(edges_mf) <- c("from", "to", "weight", "TemperatureFactor", "taxon")  # Rename for igraph format
# 
# 
# # Create igraph object
# network_graph_mf <- graph_from_data_frame(d = edges_mf, vertices = nodes_all, directed = TRUE)
# network_graph_mf <- add_edges(network_graph_mf, c("BAET", "BAET"))
# E(network_graph_mf)$TemperatureFactor[18] <- c(1.2)
# 
# 
# # Assign a layout manually for better visualization
# network_layout_mf <- create_layout(network_graph_mf, layout = "circle")
# 
# # Define taxon colors
# taxon_colors <- c(
#   "HYOS" = "#AA3377",
#   "GAMM" = "#CCBB44",
#   "NZMS" = "#4477AA",
#   "CHIR" = "#228833",
#   "BAET" = "#66CCEE"
# )
# 
# # Plot the network
# flood_mort_net <- ggraph(network_layout_mf, layout = "circle") +
#   geom_edge_arc(aes(edge_width = abs(weight), color = taxon),
#                 arrow = arrow(length = unit(4, "mm")),
#                 end_cap = circle(4, "mm"),
#                 start_cap = circle(4, "mm"),
#                 alpha = 0.8,  
#                 curvature = 0.3, 
#                 show.legend = FALSE) +
#   geom_node_point(size = 8, aes(color = Taxon)) +
#   geom_node_text(aes(label = " "), repel = F, size = 4, fontface = "bold") +
#   scale_edge_width(range = c(0.5, 3)) +
#   scale_edge_color_manual(values = taxon_colors) +
#   scale_color_manual(values = taxon_colors, 
#                      labels=c(expression(paste(italic("Baetidae")," spp.")), 
#                               expression(paste(italic("Chironomidae"), " spp.")), 
#                               expression(italic("G. lacustris")), 
#                               expression(paste(italic("Hydropsyche"), " spp.")), 
#                               expression(italic("P. antipodarum")))) +
#   theme(
#     panel.background = element_rect(fill = "white", color = "black"),
#     panel.grid = element_blank(),
#     strip.background = element_rect(fill = "white", color = "black"),
#     strip.text = element_text(color = "black"),
#     legend.key = element_rect(fill = "white", color = NA),
#     plot.title = element_text(hjust = 0.5)
#   ) +
#   labs(title = "Flood Mortality Sensitivity Analysis", color = "Taxon", edge_width = "Slope")+
#   facet_wrap(~TemperatureFactor, labeller = labeller(TemperatureFactor = c(
#     "1" = "Baseline", 
#     "1.1" = "+1°C", 
#     "1.2" = "+2°C", 
#     "1.5" = "+5°C"
#   )))
# 
# #----------------------------------------------------------
# # sensitivity network analysis for flood mortality with HFE
# Multispp_mF_sens_temps_hfe <- read_csv("Multispp_mF_sens_temps_hfe.csv")
# Multispp_mF_sens_temps_hfe$TemperatureFactor <- as.factor(Multispp_mF_sens_temps_hfe$TemperatureFactor)
# sensgrid3 <- ggplot(Multispp_mF_sens_temps_hfe, aes(x = SensitivityIncrement, y = Biomass, color = StageGroup)) +
#   geom_point(size = 1) +  # Alternatively, you can add geom_line() for smoother trends.
#   facet_grid(StageGroup ~ Parameter + TemperatureFactor, scales = "free_y") +  # Facet by both StageGroup & Parameter
#   theme_bw() +
#   labs(
#     x = "Sensitivity Increment",
#     y = "Biomass",
#     title = "Sensitivity Analysis by StageGroup and Parameter"
#   ) +
#   theme(
#     strip.text = element_text(size = 8, face = "bold"),  # Adjust facet label size
#     axis.text.x = element_text(angle = 90, hjust = 1),
#     legend.position = "none"  
#   )
# 
# ggsave("sensgrid3.png", sensgrid3, device = "png", width = 8, height = 11.5, units= c("in"), dpi = "retina")
# # Fit a linear model for each StageGroup × Parameter combination
# interaction_results_mf_hfe <- Multispp_mF_sens_temps_hfe %>%
#   group_by(StageGroup, Parameter, TemperatureFactor) %>%
#   do({
#     model <- lm(Biomass ~ SensitivityIncrement, data = .)  # Fit linear model
#     model_summary <- summary(model)
#     
#     # Extract key values
#     tibble(
#       Slope = coef(model)[2],  # Extract slope
#       P_value = coef(summary(model))["SensitivityIncrement", "Pr(>|t|)"],  # Get p-value for slope
#       R2 = model_summary$r.squared  # Get R
#     )
#   }) %>%
#   ungroup() %>%
#   # Filter results based on criteria
#   filter(R2 >= 0.13, P_value < 0.05) #moderate plus?
# 
# # Rename parameters: Change "G1" to "S1", "G2" to "S2"
# interaction_results_mf_hfe <- interaction_results_mf_hfe %>%
#   mutate(Parameter = case_when(
#     grepl("flood.mort_NZMS", Parameter) ~ "NZMS",
#     grepl("flood.mort_GAMM", Parameter) ~ "GAMM",
#     grepl("flood.mort_CHIR", Parameter) ~ "CHIR",
#     grepl("flood.mort_HYOS", Parameter) ~ "HYOS",
#     grepl("flood.mort_BAET", Parameter) ~ "BAET",
#     TRUE ~ Parameter
#   ))
# 
# interaction_results_mf_hfe <- interaction_results_mf_hfe %>%
#   mutate(StageGroup = case_when(
#     grepl("S3_NZMS", StageGroup) ~ "NZMS",
#     grepl("S3_GAMM", StageGroup) ~ "GAMM",
#     grepl("S3_CHIR", StageGroup) ~ "CHIR",
#     grepl("S3_HYOS", StageGroup) ~ "HYOS",
#     grepl("S3_BAET", StageGroup) ~ "BAET",
#     grepl("S2_NZMS", StageGroup) ~ "NZMS",
#     grepl("S2_GAMM", StageGroup) ~ "GAMM",
#     grepl("S2_CHIR", StageGroup) ~ "CHIR",
#     grepl("S2_HYOS", StageGroup) ~ "HYOS",
#     grepl("S2_BAET", StageGroup) ~ "BAET",
#     grepl("S1_NZMS", StageGroup) ~ "NZMS",
#     grepl("S1_GAMM", StageGroup) ~ "GAMM",
#     grepl("S1_CHIR", StageGroup) ~ "CHIR",
#     grepl("S1_HYOS", StageGroup) ~ "HYOS",
#     grepl("S1_BAET", StageGroup) ~ "BAET",
#     TRUE ~ StageGroup
#   ))
# 
# 
# # Extract the taxon type from parameter names
# interaction_results_mf_hfe <- interaction_results_mf_hfe %>%
#   mutate(Taxon = case_when(
#     grepl("HYOS", Parameter) ~ "HYOS",
#     grepl("GAMM", Parameter) ~ "GAMM",
#     grepl("NZMS", Parameter) ~ "NZMS",
#     grepl("CHIR", Parameter) ~ "CHIR",
#     grepl("BAET", Parameter) ~ "BAET",
#     TRUE ~ "Other"
#   ))
# 
# # Create an edge list
# edges_mf_hfe <- interaction_results_mf_hfe %>%
#   select(Parameter, StageGroup, Slope, TemperatureFactor, Taxon)
# 
# colnames(edges_mf_hfe) <- c("from", "to", "weight", "TemperatureFactor", "taxon")  # Rename for igraph format
# 
# # Create igraph object
# network_graph_mf_hfe <- graph_from_data_frame(d = edges_mf_hfe, vertices = nodes_all, directed = TRUE)
# network_graph_mf_hfe <- add_edges(network_graph_mf_hfe, c("BAET", "CHIR", "BAET", "CHIR"))
# E(network_graph_mf_hfe)$TemperatureFactor[10:11] <- c(1.2, 1.5)
# # Assign a layout manually for better visualization
# network_layout_mf_hfe <- create_layout(network_graph_mf_hfe, layout = "circle")
# 
# # Define taxon colors
# taxon_colors <- c(
#   "HYOS" = "#AA3377",
#   "GAMM" = "#CCBB44",
#   "NZMS" = "#4477AA",
#   "CHIR" = "#228833",
#   "BAET" = "#66CCEE"
# )
# 
# # Plot the network
# flood_mort_net_hfe <- ggraph(network_layout_mf_hfe, layout = "circle") +
#   geom_edge_arc(aes(edge_width = abs(weight), color = taxon),
#                 arrow = arrow(length = unit(4, "mm")),
#                 end_cap = circle(4, "mm"),
#                 start_cap = circle(4, "mm"),
#                 alpha = 0.8,  
#                 curvature = 0.3, 
#                 show.legend = FALSE) +
#   geom_node_point(size = 8, aes(color = Taxon)) +
#   geom_node_text(aes(label = " "), repel = F, size = 4, fontface = "bold") +
#   scale_edge_width(range = c(0.5, 3)) +
#   scale_edge_color_manual(values = taxon_colors) +
#   scale_color_manual(values = taxon_colors, 
#                      labels=c(expression(paste(italic("Baetidae")," spp.")), 
#                               expression(paste(italic("Chironomidae"), " spp.")), 
#                               expression(italic("G. lacustris")), 
#                               expression(paste(italic("Hydropsyche"), " spp.")), 
#                               expression(italic("P. antipodarum")))) +
#   theme(
#     panel.background = element_rect(fill = "white", color = "black"),  
#     panel.grid = element_blank(), 
#     strip.background = element_rect(fill = "white", color = "black"), 
#     strip.text = element_text(color = "black"),  
#     legend.key = element_rect(fill = "white", color = NA),
#     plot.title = element_text(hjust = 0.5)
#   ) +
#   labs(title = "Flood Mortality Sensitivity Analysis with HFE", color = "Taxon", edge_width = "Slope")+
#   facet_wrap(~TemperatureFactor, labeller = labeller(TemperatureFactor = c(
#     "1" = "Baseline", 
#     "1.1" = "+1°C", 
#     "1.2" = "+2°C", 
#     "1.5" = "+5°C"
#   )))
# 
# 
# #-------------------------------------------------------------------------------
# Multispp_mh_sens_hyd <- read_csv("Multispp_mh_sens_hyd.csv")
# Multispp_mh_sens_hyd$TemperatureFactor <- as.factor(Multispp_mh_sens_hyd$TemperatureFactor)
# sensgrid4 <- ggplot(Multispp_mh_sens_hyd, aes(x = SensitivityIncrement, y = Biomass, color = StageGroup)) +
#   geom_point(size = 1) +  # Alternatively, you can add geom_line() for smoother trends.
#   facet_grid(StageGroup ~ Parameter + TemperatureFactor, scales = "free_y") +  # Facet by both StageGroup & Parameter
#   theme_bw() +
#   labs(
#     x = "Sensitivity Increment",
#     y = "Biomass",
#     title = "Sensitivity Analysis by StageGroup and Parameter"
#   ) +
#   theme(
#     strip.text = element_text(size = 8, face = "bold"),  # Adjust facet label size
#     axis.text.x = element_text(angle = 90, hjust = 1),
#     legend.position = "none"  
#   )
# 
# ggsave("sensgrid4.png", sensgrid4, device = "png", width = 8, height = 11.5, units= c("in"), dpi = "retina")
# # Fit a linear model for each StageGroup × Parameter combination
# interaction_results_mh_hyd <- Multispp_mh_sens_hyd %>%
#   group_by(StageGroup, Parameter, TemperatureFactor) %>%
#   do({
#     model <- lm(Biomass ~ SensitivityIncrement, data = .)  # Fit linear model
#     model_summary <- summary(model)
#     
#     # Extract key values
#     tibble(
#       Slope = coef(model)[2],  # Extract slope
#       P_value = coef(summary(model))["SensitivityIncrement", "Pr(>|t|)"],  # Get p-value for slope
#       R2 = model_summary$r.squared  # Get R
#     )
#   }) %>%
#   ungroup() %>%
#   # Filter results based on criteria
#   filter(R2 >= 0.13, P_value < 0.05) #moderate plus?
# 
# # Rename parameters: Change "G1" to "S1", "G2" to "S2"
# interaction_results_mh_hyd <- interaction_results_mh_hyd %>%
#   mutate(Parameter = case_when(
#     grepl("hydro.mort_NZMS", Parameter) ~ "NZMS",
#     grepl("hydro.mort_GAMM", Parameter) ~ "GAMM",
#     grepl("hydro.mort_CHIR", Parameter) ~ "CHIR",
#     grepl("hydro.mort_HYOS", Parameter) ~ "HYOS",
#     grepl("hydro.mort_BAET", Parameter) ~ "BAET",
#     TRUE ~ Parameter
#   ))
# 
# interaction_results_mh_hyd <- interaction_results_mh_hyd %>%
#   mutate(StageGroup = case_when(
#     grepl("S3_NZMS", StageGroup) ~ "NZMS",
#     grepl("S3_GAMM", StageGroup) ~ "GAMM",
#     grepl("S3_CHIR", StageGroup) ~ "CHIR",
#     grepl("S3_HYOS", StageGroup) ~ "HYOS",
#     grepl("S3_BAET", StageGroup) ~ "BAET",
#     grepl("S2_NZMS", StageGroup) ~ "NZMS",
#     grepl("S2_GAMM", StageGroup) ~ "GAMM",
#     grepl("S2_CHIR", StageGroup) ~ "CHIR",
#     grepl("S2_HYOS", StageGroup) ~ "HYOS",
#     grepl("S2_BAET", StageGroup) ~ "BAET",
#     grepl("S1_NZMS", StageGroup) ~ "NZMS",
#     grepl("S1_GAMM", StageGroup) ~ "GAMM",
#     grepl("S1_CHIR", StageGroup) ~ "CHIR",
#     grepl("S1_HYOS", StageGroup) ~ "HYOS",
#     grepl("S1_BAET", StageGroup) ~ "BAET",
#     TRUE ~ StageGroup
#   ))
# 
# 
# # Extract the taxon type from parameter names
# interaction_results_mh_hyd <- interaction_results_mh_hyd %>%
#   mutate(Taxon = case_when(
#     grepl("HYOS", Parameter) ~ "HYOS",
#     grepl("GAMM", Parameter) ~ "GAMM",
#     grepl("NZMS", Parameter) ~ "NZMS",
#     grepl("CHIR", Parameter) ~ "CHIR",
#     grepl("BAET", Parameter) ~ "BAET",
#     TRUE ~ "Other"
#   ))
# 
# # Create an edge list
# edges_mh_hyd <- interaction_results_mh_hyd %>%
#   select(Parameter, StageGroup, Slope, TemperatureFactor, Taxon)
# 
# colnames(edges_mh_hyd) <- c("from", "to", "weight", "TemperatureFactor", "taxon")  # Rename for igraph format
# 
# # Create igraph object
# network_graph_mh_hyd <- graph_from_data_frame(d = edges_mh_hyd, vertices = nodes_all, directed = TRUE)
# 
# # Create new edges with TemperatureFactor values (assuming BAET->CHIR is an example edge)
# network_graph_mh_hyd <- add_edges(network_graph_mh_hyd, c("BAET", "CHIR", "BAET", "CHIR", "BAET", "CHIR"))
# E(network_graph_mh_hyd)$TemperatureFactor[6:8] <- c(1.1, 1.2, 1.5)
# 
# # Assign a layout manually for better visualization
# network_layout_mh_hyd <- create_layout(network_graph_mh_hyd, layout = "circle")
# 
# # Define taxon colors
# taxon_colors <- c(
#   "HYOS" = "#AA3377",
#   "GAMM" = "#CCBB44",
#   "NZMS" = "#4477AA",
#   "CHIR" = "#228833",
#   "BAET" = "#66CCEE"
# )
# 
# # Plot the network
# hyd_mort_net_hyd <- ggraph(network_layout_mh_hyd, layout = "circle") +
#   geom_edge_arc(aes(edge_width = abs(weight), color = taxon),
#                 arrow = arrow(length = unit(4, "mm")),
#                 end_cap = circle(4, "mm"),
#                 start_cap = circle(4, "mm"),
#                 alpha = 0.8,  
#                 curvature = 0.3, 
#                 show.legend = FALSE) +
#   geom_node_point(size = 8, aes(color = Taxon)) +
#   geom_node_text(aes(label = " "), repel = F, size = 4, fontface = "bold") +
#   scale_edge_width(range = c(0.5, 3)) +
#   scale_edge_color_manual(values = taxon_colors) +
#   scale_color_manual(values = taxon_colors, 
#                      labels=c(expression(paste(italic("Baetidae")," spp.")), 
#                               expression(paste(italic("Chironomidae"), " spp.")), 
#                               expression(italic("G. lacustris")), 
#                               expression(paste(italic("Hydropsyche"), " spp.")), 
#                               expression(italic("P. antipodarum")))) +
#   theme(
#     panel.background = element_rect(fill = "white", color = "black"),
#     panel.grid = element_blank(), 
#     strip.background = element_rect(fill = "white", color = "black"), 
#     strip.text = element_text(color = "black"), 
#     legend.key = element_rect(fill = "white", color = NA), 
#     plot.title = element_text(hjust = 0.5)
#   ) +
#   labs(title = "Hydropeaking Mortality Sensitivity Analysis with Hydropeaking", color = "Taxon", edge_width = "Slope")+
#   facet_wrap(~TemperatureFactor, labeller = labeller(TemperatureFactor = c(
#     "1" = "Baseline", 
#     "1.1" = "+1°C", 
#     "1.2" = "+2°C", 
#     "1.5" = "+5°C"
#   )))
# 
# #--------------------------------------------------------------------------------
# # temperature mortality sensitivity network
# Multispp_mt_temps <- read_csv("Multispp_mT_sens_temps.csv")
# 
# 
# sensgrid1 <- ggplot(Multispp_mt_temps, aes(x = SensitivityIncrement, y = Biomass, color = StageGroup)) +
#   geom_point(size = 1) +  # Alternatively, you can add geom_line() for smoother trends.
#   facet_grid(StageGroup ~ TemperatureFactor * Parameter, scales = "free_y") +  # Facet by both StageGroup & Parameter
#   theme_bw() +
#   labs(
#     x = "Sensitivity Increment",
#     y = "Biomass",
#     title = "Sensitivity Analysis by StageGroup and Parameter"
#   ) +
#   theme(
#     strip.text = element_text(size = 8, face = "bold"),  # Adjust facet label size
#     axis.text.x = element_text(angle = 90, hjust = 1),
#     legend.position = "none"  # Hide legend (since color is redundant with facet_grid)
#   )
# 
# ggsave("sensgrid1.png", sensgrid1, device = "png", width = 8, height = 11.5, units= c("in"), dpi = "retina")
# # Fit a linear model for each StageGroup × Parameter combination
# interaction_results_mt <- Multispp_mt_temps %>%
#   group_by(StageGroup, Parameter, TemperatureFactor) %>%
#   do({
#     model <- lm(Abundance ~ SensitivityIncrement, data = .)  # Fit linear model
#     model_summary <- summary(model)
#     
#     # Extract key values
#     tibble(
#       Slope = coef(model)[2],  # Extract slope
#       P_value = coef(summary(model))["SensitivityIncrement", "Pr(>|t|)"],  # Get p-value for slope
#       R2 = model_summary$r.squared  # Get R²
#     )
#   }) %>%
#   ungroup() %>%
#   # Filter results based on criteria
#   filter(R2 >= 0.13) #Cohen J. (1988). Statistical Power Analysis for the Behavioral Sciences, 2nd Ed. Hillsdale, NJ: Laurence Erlbaum Associates p 413-414 
# 
# # Rename parameters: Change "G1" to "S1", "G2" to "S2"
# interaction_results_mt <- interaction_results_mt %>%
#   mutate(Parameter = case_when(
#     grepl("TempSurvival_NZMS", Parameter) ~ "NZMS",
#     grepl("TempSurvival_GAMM", Parameter) ~ "GAMM",
#     grepl("TempSurvival_CHIR", Parameter) ~ "CHIR",
#     grepl("TempSurvival_HYOS", Parameter) ~ "HYOS",
#     grepl("TempSurvival_BAET", Parameter) ~ "BAET",
#     TRUE ~ Parameter
#   ))
# 
# interaction_results_mt <- interaction_results_mt %>%
#   mutate(StageGroup = case_when(
#     grepl("S3_NZMS", StageGroup) ~ "NZMS",
#     grepl("S3_GAMM", StageGroup) ~ "GAMM",
#     grepl("S3_CHIR", StageGroup) ~ "CHIR",
#     grepl("S3_HYOS", StageGroup) ~ "HYOS",
#     grepl("S3_BAET", StageGroup) ~ "BAET",
#     grepl("S2_GAMM", StageGroup) ~ "GAMM",
#     grepl("S2_HYOS", StageGroup) ~ "HYOS",
#     TRUE ~ StageGroup
#   ))
# 
# 
# # Extract the taxon type from parameter names
# interaction_results_mt <- interaction_results_mt %>%
#   mutate(Taxon = case_when(
#     grepl("HYOS", Parameter) ~ "HYOS",
#     grepl("GAMM", Parameter) ~ "GAMM",
#     grepl("NZMS", Parameter) ~ "NZMS",
#     grepl("CHIR", Parameter) ~ "CHIR",
#     grepl("BAET", Parameter) ~ "BAET",
#     TRUE ~ "Other"
#   ))
# 
# # Create an edge list
# edges_mt <- interaction_results_mt %>%
#   select(Parameter, StageGroup, Slope, TemperatureFactor, Taxon)
# 
# colnames(edges_mt) <- c("from", "to", "weight", "TemperatureFactor", "taxon")  # Rename for igraph format
# 
# # Create igraph object
# network_graph_mt <- graph_from_data_frame(d = edges_mt, vertices = nodes_all, directed = TRUE)
# 
# # Assign a layout manually for better visualization
# network_layout_mt <- create_layout(network_graph_mt, layout = "circle")
# 
# # Define taxon colors
# taxon_colors <- c(
#   "HYOS" = "#AA3377",
#   "GAMM" = "#CCBB44",
#   "NZMS" = "#4477AA",
#   "CHIR" = "#228833",
#   "BAET" = "#66CCEE"
# )
# 
# # Plot the network
# temp_mort_net <- ggraph(network_layout_mt, layout = "circle") +
#   geom_edge_arc(aes(edge_width = abs(weight), color = taxon),
#                 arrow = arrow(length = unit(4, "mm")),
#                 end_cap = circle(4, "mm"),
#                 start_cap = circle(4, "mm"),
#                 alpha = 0.8,  
#                 curvature = 0.3, 
#                 show.legend = FALSE) +
#   geom_node_point(size = 8, aes(color = Taxon)) +
#   geom_node_text(aes(label = " "), repel = F, size = 4, fontface = "bold") +
#   scale_edge_width(range = c(0.5, 3)) +
#   scale_edge_color_manual(values = taxon_colors) +
#   scale_color_manual(values = taxon_colors, 
#                      labels=c(expression(paste(italic("Baetidae")," spp.")), 
#                               expression(paste(italic("Chironomidae"), " spp.")), 
#                               expression(italic("G. lacustris")), 
#                               expression(paste(italic("Hydropsyche"), " spp.")), 
#                               expression(italic("P. antipodarum")))) +
#   theme(
#     panel.background = element_rect(fill = "white", color = "black"), 
#     panel.grid = element_blank(),  
#     strip.background = element_rect(fill = "white", color = "black"), 
#     strip.text = element_text(color = "black"), 
#     legend.key = element_rect(fill = "white", color = NA),
#     plot.title = element_text(hjust = 0.5)
#       ) +
#   labs(title = "Temperature Mortality Sensitivity Analysis", color = "Taxon", edge_width = "Slope")+
#   facet_wrap(~TemperatureFactor, labeller = labeller(TemperatureFactor = c(
#     "1" = "Baseline", 
#     "1.1" = "+1°C", 
#     "1.2" = "+2°C", 
#     "1.5" = "+5°C"
#   )))
# # no sig relationships for abundance or biomass
# #---------------------------------------------------------------------------
# Multispp_mt_sens_temps_sp <- read_csv("Multispp_mT_sens_temps_spike.csv")
# Multispp_mt_sens_temps_sp$TemperatureFactor <- as.factor(Multispp_mt_sens_temps_sp$TemperatureFactor)
# sensgrid4 <- ggplot(Multispp_mt_sens_temps_sp, aes(x = SensitivityIncrement, y = Biomass, color = StageGroup)) +
#   geom_point(size = 1) +  # Alternatively, you can add geom_line() for smoother trends.
#   facet_grid(StageGroup ~ Parameter + TemperatureFactor, scales = "free_y") +  # Facet by both StageGroup & Parameter
#   theme_bw() +
#   labs(
#     x = "Sensitivity Increment",
#     y = "Biomass",
#     title = "Sensitivity Analysis by StageGroup and Parameter"
#   ) +
#   theme(
#     strip.text = element_text(size = 8, face = "bold"),  # Adjust facet label size
#     axis.text.x = element_text(angle = 90, hjust = 1),
#     legend.position = "none"  
#   )
# 
# ggsave("sensgrid4.png", sensgrid4, device = "png", width = 8, height = 11.5, units= c("in"), dpi = "retina")
# # Fit a linear model for each StageGroup × Parameter combination
# interaction_results_mt_sp <- Multispp_mt_sens_temps_sp %>%
#   group_by(StageGroup, Parameter, TemperatureFactor) %>%
#   do({
#     model <- lm(Biomass ~ SensitivityIncrement, data = .)  # Fit linear model
#     model_summary <- summary(model)
#     
#     # Extract key values
#     tibble(
#       Slope = coef(model)[2],  # Extract slope
#       P_value = coef(summary(model))["SensitivityIncrement", "Pr(>|t|)"],  # Get p-value for slope
#       R2 = model_summary$r.squared  # Get R
#     )
#   }) %>%
#   ungroup() %>%
#   # Filter results based on criteria
#   filter(R2 >= 0.13) #moderate plus?
# 
# # Rename parameters: Change "G1" to "S1", "G2" to "S2"
# interaction_results_mt_sp <- interaction_results_mt_sp %>%
#   mutate(Parameter = case_when(
#     grepl("TempSurvival_NZMS", Parameter) ~ "NZMS",
#     grepl("TempSurvival_GAMM", Parameter) ~ "GAMM",
#     grepl("TempSurvival_CHIR", Parameter) ~ "CHIR",
#     grepl("TempSurvival_HYOS", Parameter) ~ "HYOS",
#     grepl("TempSurvival_BAET", Parameter) ~ "BAET",
#     TRUE ~ Parameter
#   ))
# 
# interaction_results_mt_sp <- interaction_results_mt_sp %>%
#   mutate(StageGroup = case_when(
#     grepl("S3_NZMS", StageGroup) ~ "NZMS",
#     grepl("S3_GAMM", StageGroup) ~ "GAMM",
#     grepl("S3_CHIR", StageGroup) ~ "CHIR",
#     grepl("S3_HYOS", StageGroup) ~ "HYOS",
#     grepl("S3_BAET", StageGroup) ~ "BAET",
#     grepl("S2_NZMS", StageGroup) ~ "NZMS",
#     grepl("S2_GAMM", StageGroup) ~ "GAMM",
#     grepl("S2_CHIR", StageGroup) ~ "CHIR",
#     grepl("S2_HYOS", StageGroup) ~ "HYOS",
#     grepl("S2_BAET", StageGroup) ~ "BAET",
#     grepl("S1_NZMS", StageGroup) ~ "NZMS",
#     grepl("S1_GAMM", StageGroup) ~ "GAMM",
#     grepl("S1_CHIR", StageGroup) ~ "CHIR",
#     grepl("S1_HYOS", StageGroup) ~ "HYOS",
#     grepl("S1_BAET", StageGroup) ~ "BAET",
#     TRUE ~ StageGroup
#   ))
# 
# 
# # Extract the taxon type from parameter names
# interaction_results_mt_sp <- interaction_results_mt_sp %>%
#   mutate(Taxon = case_when(
#     grepl("HYOS", Parameter) ~ "HYOS",
#     grepl("GAMM", Parameter) ~ "GAMM",
#     grepl("NZMS", Parameter) ~ "NZMS",
#     grepl("CHIR", Parameter) ~ "CHIR",
#     grepl("BAET", Parameter) ~ "BAET",
#     TRUE ~ "Other"
#   ))
# 
# # Create an edge list
# edges_mt_sp <- interaction_results_mt_sp%>%
#   select(Parameter, StageGroup, Slope, TemperatureFactor, Taxon)
# 
# colnames(edges_mt_sp) <- c("from", "to", "weight", "TemperatureFactor", "taxon")  # Rename for igraph format
# 
# 
# # Create igraph object
# network_graph_mt_sp <- graph_from_data_frame(d = edges_mt_sp, vertices = nodes_all, directed = TRUE)
# 
# 
# # Assign a layout manually for better visualization
# network_layout_mt_sp <- create_layout(network_graph_mt_sp, layout = "circle")
# 
# # Define taxon colors
# taxon_colors <- c(
#   "HYOS" = "#AA3377",
#   "GAMM" = "#CCBB44",
#   "NZMS" = "#4477AA",
#   "CHIR" = "#228833",
#   "BAET" = "#66CCEE"
# )
# 
# # Plot the network
# temp_mort_net_spike <- ggraph(network_layout_mt_sp, layout = "circle") +
#   geom_edge_arc(aes(edge_width = abs(weight), color = taxon),
#                 arrow = arrow(length = unit(4, "mm")),
#                 end_cap = circle(4, "mm"),
#                 start_cap = circle(4, "mm"),
#                 alpha = 0.8,  
#                 curvature = 0.3, 
#                 show.legend = FALSE) +
#   geom_node_point(size = 8, aes(color = Taxon)) +
#   geom_node_text(aes(label = " "), repel = F, size = 4, fontface = "bold") +
#   scale_edge_width(range = c(0.5, 3)) +
#   scale_edge_color_manual(values = taxon_colors) +
#   scale_color_manual(values = taxon_colors, 
#                      labels=c(expression(paste(italic("Baetidae")," spp.")), 
#                               expression(paste(italic("Chironomidae"), " spp.")), 
#                               expression(italic("G. lacustris")), 
#                               expression(paste(italic("Hydropsyche"), " spp.")), 
#                               expression(italic("P. antipodarum")))) +
#   theme(
#     panel.background = element_rect(fill = "white", color = "black"),  
#     panel.grid = element_blank(),  
#     strip.background = element_rect(fill = "white", color = "black"),  
#     strip.text = element_text(color = "black"), 
#     legend.key = element_rect(fill = "white", color = NA), 
#     plot.title = element_text(hjust = 0.5)
#   ) +
#   labs(title = "Temperature Mortality Sensitivity Analysis with Temperature Spike", color = "Taxon", edge_width = "Slope")+
#   facet_wrap(~TemperatureFactor, labeller = labeller(TemperatureFactor = c(
#     "1" = "Baseline", 
#     "1.1" = "+1°C", 
#     "1.2" = "+2°C", 
#     "1.5" = "+5°C"
#   )))
# #-----------------------------------------------------------------------------
# 
# mort_network <- ggarrange(flood_mort_net_hfe, hyd_mort_net_hyd,
#           labels = c("a", "b"),
#           ncol= 2,
#           common.legend = T)
# ggsave("mort_network.png", plot = flood_mort_network, device = "png", dpi = "retina", height = 8, width  = 14, units = "in")
# hyd_mort_network <- ggarrange(hyd_mort_net, hyd_mort_net_hyd, 
#           labels = c("a", "b"),
#           ncol= 2,
#           common.legend = T)
# ggsave("hyd_mort_network.png", plot = hyd_mort_network, device = "png", dpi = "retina", height = 8, width  = 14, units = "in")
# 
#----------------------------------------------------------
# sensitivity network analysis for flood mortality with HFE
all_nodes <- unique(c("S1_NZMS", "S2_NZMS", "S3_NZMS","S1_GAMM","S2_GAMM","S3_GAMM",
                    "S1_CHIR","S2_CHIR", "S3_CHIR","S1_HYOS", "S2_HYOS","S3_HYOS",
                    "S1_BAET","S2_BAET","S3_BAET"))
nodes_all <- tibble(name = all_nodes) %>%
  mutate(Taxon = case_when(
    str_detect(name, "S[123]_NZMS") ~ "NZMS",
    str_detect(name, "S[123]_GAMM") ~ "GAMM",
    str_detect(name, "S[123]_CHIR") ~ "CHIR",
    str_detect(name, "S[123]_HYOS") ~ "HYOS",
    str_detect(name, "S[123]_BAET") ~ "BAET",
    TRUE ~ NA_character_
  ),
  type = ifelse(name %in% edges_mh$from | name %in% edges_mh$to, "Source", "Target"))



Multispp_mF_sens_temps_hfe <- read_csv("Multispp_mF_sens_temps_hfe.csv")
Multispp_mF_sens_temps_hfe$TemperatureFactor <- as.factor(Multispp_mF_sens_temps_hfe$TemperatureFactor)
sensgrid3 <- ggplot(Multispp_mF_sens_temps_hfe, aes(x = SensitivityIncrement, y = Biomass, color = StageGroup)) +
  geom_point(size = 1) +  # Alternatively, you can add geom_line() for smoother trends.
  facet_grid(StageGroup ~ Parameter + TemperatureFactor, scales = "free_y") +  # Facet by both StageGroup & Parameter
  theme_bw() +
  labs(
    x = "Sensitivity Increment",
    y = "Biomass",
    title = "Sensitivity Analysis by StageGroup and Parameter"
  ) +
  theme(
    strip.text = element_text(size = 8, face = "bold"),  # Adjust facet label size
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none"  
  )

ggsave("sensgrid3.png", sensgrid3, device = "png", width = 8, height = 11.5, units= c("in"), dpi = "retina")
# Fit a linear model for each StageGroup × Parameter combination
interaction_results_mf_hfe <- Multispp_mF_sens_temps_hfe %>%
  group_by(StageGroup, Parameter, TemperatureFactor) %>%
  do({
    model <- lm(Biomass ~ SensitivityIncrement, data = .)  # Fit linear model
    model_summary <- summary(model)
    
    # Extract key values
    tibble(
      Slope = coef(model)[2],  # Extract slope
      P_value = coef(summary(model))["SensitivityIncrement", "Pr(>|t|)"],  # Get p-value for slope
      R2 = model_summary$r.squared  # Get R
    )
  }) %>%
  ungroup() %>%
  # Filter results based on criteria
  filter(R2 >= 0.13, P_value < 0.05) #moderate plus?

meanBiomass <- Multispp_mF_sens_temps_hfe %>%
  group_by(StageGroup, TemperatureFactor, Parameter) %>%
  dplyr::summarize(meanBio = mean(Biomass, na.rm = T))

interaction_results_mf_hfe = merge(interaction_results_mf_hfe, meanBiomass, by=c("StageGroup", "Parameter", "TemperatureFactor"))

# Rename parameters: Change "G1" to "S1", "G2" to "S2"
interaction_results_mf_hfe <- interaction_results_mf_hfe %>%
  mutate(Parameter = case_when(
    grepl("flood.mort_NZMS", Parameter) ~ "S2_NZMS",
    grepl("flood.mort_GAMM", Parameter) ~ "S2_GAMM",
    grepl("flood.mort_CHIR", Parameter) ~ "S2_CHIR",
    grepl("flood.mort_HYOS", Parameter) ~ "S2_HYOS",
    grepl("flood.mort_BAET", Parameter) ~ "S2_BAET",
    TRUE ~ Parameter
  ))

# interaction_results_mf_hfe <- interaction_results_mf_hfe %>%
#   mutate(StageGroup = case_when(
#     grepl("S3_NZMS", StageGroup) ~ "NZMS",
#     grepl("S3_GAMM", StageGroup) ~ "GAMM",
#     grepl("S3_CHIR", StageGroup) ~ "CHIR",
#     grepl("S3_HYOS", StageGroup) ~ "HYOS",
#     grepl("S3_BAET", StageGroup) ~ "BAET",
#     grepl("S2_NZMS", StageGroup) ~ "NZMS",
#     grepl("S2_GAMM", StageGroup) ~ "GAMM",
#     grepl("S2_CHIR", StageGroup) ~ "CHIR",
#     grepl("S2_HYOS", StageGroup) ~ "HYOS",
#     grepl("S2_BAET", StageGroup) ~ "BAET",
#     grepl("S1_NZMS", StageGroup) ~ "NZMS",
#     grepl("S1_GAMM", StageGroup) ~ "GAMM",
#     grepl("S1_CHIR", StageGroup) ~ "CHIR",
#     grepl("S1_HYOS", StageGroup) ~ "HYOS",
#     grepl("S1_BAET", StageGroup) ~ "BAET",
#     TRUE ~ StageGroup
#   ))


# Extract the taxon type from parameter names
interaction_results_mf_hfe <- interaction_results_mf_hfe %>%
  mutate(Taxon = case_when(
    grepl("HYOS", Parameter) ~ "HYOS",
    grepl("GAMM", Parameter) ~ "GAMM",
    grepl("NZMS", Parameter) ~ "NZMS",
    grepl("CHIR", Parameter) ~ "CHIR",
    grepl("BAET", Parameter) ~ "BAET",
    TRUE ~ "Other"
  ))


# Create an edge list
edges_mf_hfe <- interaction_results_mf_hfe %>%
  select(Parameter, StageGroup, Slope, meanBio, TemperatureFactor, Taxon)

colnames(edges_mf_hfe) <- c("from", "to", "weight", "size", "TemperatureFactor", "taxon")  # Rename for igraph format

# Create igraph object
network_graph_mf_hfe <- graph_from_data_frame(d = edges_mf_hfe, vertices = nodes_all, directed = TRUE)
network_graph_mf_hfe <- add_edges(network_graph_mf_hfe, c("S2_BAET", "S2_CHIR", "S2_BAET", "S2_CHIR"))
E(network_graph_mf_hfe)$TemperatureFactor[10:11] <- c(1.2, 1.5)
E(network_graph_mf_hfe)$size[which(is.na(E(network_graph_mf_hfe)$size)==T)] <- 0
# vertex size is scaled by the mean abundance of the vertex spp found in 'mean_ab_l[[i]]' 
# to do this i'm using the match function
# Assign vertex size based on meanBio from the edge list
V(network_graph_mf_hfe)$size <- edges_mf_hfe$size[match(V(network_graph_mf_hfe)$name, edges_mf_hfe$to)] * 1.5
V(network_graph_mf_hfe)$size[which(is.na(V(network_graph_mf_hfe)$size)==T)] <- min(V(network_graph_mf_hfe)$size, na.rm = T)/2
# Assign a layout manually for better visualization
network_layout_mf_hfe <- create_layout(network_graph_mf_hfe, layout = "circle")

# Define taxon colors
taxon_colors <- c(
  "HYOS" = "#AA3377",
  "GAMM" = "#CCBB44",
  "NZMS" = "#4477AA",
  "CHIR" = "#228833",
  "BAET" = "#66CCEE"
)

# Create a vector of labels for the nodes
node_labels <- rep(c("S1", "S2", "S3"), length.out = 60)

# Assign these labels to the vertices of the graph
V(network_graph_mf_hfe)$label <- node_labels * 3

# Plot the network
flood_mort_net_hfe <- ggraph(network_layout_mf_hfe, layout = "circle") +
  geom_edge_arc(aes(edge_width = abs(weight), color = taxon),
                arrow = arrow(length = unit(4, "mm")),
                end_cap = circle(4, "mm"),
                start_cap = circle(4, "mm"),
                alpha = 0.8,  
                curvature = 0.3, 
                show.legend = FALSE) +
  geom_node_point(aes(color = Taxon, size = size)) +
  geom_node_text(aes(label = node_labels, x = 1.2 * x,
                 y = 1.2 * y)) +  
  scale_edge_width(range = c(0.5, 3)) +
  scale_edge_color_manual(values = taxon_colors) +
  scale_color_manual(values = taxon_colors, 
                     labels=c(expression(paste(italic("Baetidae")," spp.")), 
                              expression(paste(italic("Chironomidae"), " spp.")), 
                              expression(italic("G. lacustris")), 
                              expression(paste(italic("Hydropsyche"), " spp.")), 
                              expression(italic("P. antipodarum")))) +
  scale_size_continuous(range = c(2, 10)) +  # Scale node size
  theme(
    panel.background = element_rect(fill = "white", color = "black"),  
    panel.grid = element_blank(), 
    strip.background = element_rect(fill = "white", color = "black"), 
    strip.text = element_text(color = "black"),  
    legend.key = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "Flood Mortality Sensitivity Analysis with HFE", 
       color = "Taxon", 
       edge_width = "Slope") +
  guides(size = "none") +  # Remove size legend
  facet_wrap(~TemperatureFactor, labeller = labeller(TemperatureFactor = c(
    "1" = "Baseline", 
    "1.1" = "+1°C", 
    "1.2" = "+2°C", 
    "1.5" = "+5°C"
  )))


# Print the plot
print(flood_mort_net_hfe)
ggsave(filename = "flow_net.png", flood_mort_net_hfe, device = "png", 
       width = 6.25 , height = 5.5, units = c("in"), dpi= "retina")

#-------------------------------------------------------------------------------
Multispp_mh_sens_hyd <- read_csv("Multispp_mh_sens_hyd.csv")
Multispp_mh_sens_hyd$TemperatureFactor <- as.factor(Multispp_mh_sens_hyd$TemperatureFactor)
sensgrid4 <- ggplot(Multispp_mh_sens_hyd, aes(x = SensitivityIncrement, y = Biomass, color = StageGroup)) +
  geom_point(size = 1) +  # Alternatively, you can add geom_line() for smoother trends.
  facet_grid(StageGroup ~ Parameter + TemperatureFactor, scales = "free_y") +  # Facet by both StageGroup & Parameter
  theme_bw() +
  labs(
    x = "Sensitivity Increment",
    y = "Biomass",
    title = "Sensitivity Analysis by StageGroup and Parameter"
  ) +
  theme(
    strip.text = element_text(size = 8, face = "bold"),  # Adjust facet label size
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none"  
  )

ggsave("sensgrid4.png", sensgrid4, device = "png", width = 8, height = 11.5, units= c("in"), dpi = "retina")
# Fit a linear model for each StageGroup × Parameter combination
interaction_results_mh_hyd <- Multispp_mh_sens_hyd %>%
  group_by(StageGroup, Parameter, TemperatureFactor) %>%
  do({
    model <- lm(Biomass ~ SensitivityIncrement, data = .)  # Fit linear model
    model_summary <- summary(model)
    
    # Extract key values
    tibble(
      Slope = coef(model)[2],  # Extract slope
      P_value = coef(summary(model))["SensitivityIncrement", "Pr(>|t|)"],  # Get p-value for slope
      R2 = model_summary$r.squared  # Get R
    )
  }) %>%
  ungroup() %>%
  # Filter results based on criteria
  filter(R2 >= 0.13, P_value < 0.05) #moderate plus?

meanBiomass <- Multispp_mh_sens_hyd %>%
  group_by(StageGroup, TemperatureFactor, Parameter) %>%
  dplyr::summarize(meanBio = mean(Biomass, na.rm = T))

interaction_results_mh_hyd = merge(interaction_results_mh_hyd, meanBiomass, by=c("StageGroup", "Parameter", "TemperatureFactor"))

# Rename parameters: Change "G1" to "S1", "G2" to "S2"
interaction_results_mh_hyd <- interaction_results_mh_hyd %>%
  mutate(Parameter = case_when(
    grepl("hydro.mort_NZMS", Parameter) ~ "S3_NZMS",
    grepl("hydro.mort_GAMM", Parameter) ~ "S3_GAMM",
    grepl("hydro.mort_CHIR", Parameter) ~ "S3_CHIR",
    grepl("hydro.mort_HYOS", Parameter) ~ "S3_HYOS",
    grepl("hydro.mort_BAET", Parameter) ~ "S3_BAET",
    TRUE ~ Parameter
  ))

# interaction_results_mf_hfe <- interaction_results_mf_hfe %>%
#   mutate(StageGroup = case_when(
#     grepl("S3_NZMS", StageGroup) ~ "NZMS",
#     grepl("S3_GAMM", StageGroup) ~ "GAMM",
#     grepl("S3_CHIR", StageGroup) ~ "CHIR",
#     grepl("S3_HYOS", StageGroup) ~ "HYOS",
#     grepl("S3_BAET", StageGroup) ~ "BAET",
#     grepl("S2_NZMS", StageGroup) ~ "NZMS",
#     grepl("S2_GAMM", StageGroup) ~ "GAMM",
#     grepl("S2_CHIR", StageGroup) ~ "CHIR",
#     grepl("S2_HYOS", StageGroup) ~ "HYOS",
#     grepl("S2_BAET", StageGroup) ~ "BAET",
#     grepl("S1_NZMS", StageGroup) ~ "NZMS",
#     grepl("S1_GAMM", StageGroup) ~ "GAMM",
#     grepl("S1_CHIR", StageGroup) ~ "CHIR",
#     grepl("S1_HYOS", StageGroup) ~ "HYOS",
#     grepl("S1_BAET", StageGroup) ~ "BAET",
#     TRUE ~ StageGroup
#   ))


# Extract the taxon type from parameter names
interaction_results_mh_hyd <- interaction_results_mh_hyd %>%
  mutate(Taxon = case_when(
    grepl("HYOS", Parameter) ~ "HYOS",
    grepl("GAMM", Parameter) ~ "GAMM",
    grepl("NZMS", Parameter) ~ "NZMS",
    grepl("CHIR", Parameter) ~ "CHIR",
    grepl("BAET", Parameter) ~ "BAET",
    TRUE ~ "Other"
  ))


# Create an edge list
edges_mh_hyd <- interaction_results_mh_hyd %>%
  select(Parameter, StageGroup, Slope, meanBio, TemperatureFactor, Taxon)

colnames(edges_mh_hyd) <- c("from", "to", "weight", "size", "TemperatureFactor", "taxon")  # Rename for igraph format

# Create igraph object
network_graph_mh_hyd <- graph_from_data_frame(d = edges_mh_hyd, vertices = nodes_all, directed = TRUE)
network_graph_mh_hyd <- add_edges(network_graph_mh_hyd, c("S2_BAET", "S2_CHIR", "S2_BAET", "S2_BAET", "S2_CHIR", "S2_BAET"))
E(network_graph_mh_hyd)$TemperatureFactor[6:8] <- c(1.1, 1.2, 1.5)
E(network_graph_mh_hyd)$size[which(is.na(E(network_graph_mh_hyd)$size)==T)] <- 0
# vertex size is scaled by the mean abundance of the vertex spp found in 'mean_ab_l[[i]]' 
# to do this i'm using the match function
# Assign vertex size based on meanBio from the edge list
V(network_graph_mh_hyd)$size <- edges_mh_hyd$size[match(V(network_graph_mh_hyd)$name, edges_mh_hyd$to)] * 1.5
V(network_graph_mh_hyd)$size[which(is.na(V(network_graph_mh_hyd)$size)==T)] <- min(V(network_graph_mh_hyd)$size, na.rm = T)/2
# Assign a layout manually for better visualization
network_layout_mh_hyd <- create_layout(network_graph_mh_hyd, layout = "circle")

# Define taxon colors
taxon_colors <- c(
  "HYOS" = "#AA3377",
  "GAMM" = "#CCBB44",
  "NZMS" = "#4477AA",
  "CHIR" = "#228833",
  "BAET" = "#66CCEE"
)
# Create a vector of labels for the nodes
node_labels <- rep(c("S1", "S2", "S3"), length.out = 60)

# Assign these labels to the vertices of the graph
V(network_graph_mh_hyd)$label <- node_labels * 3

# Plot the network
hyd_mort_net_hyd <- ggraph(network_layout_mh_hyd, layout = "circle") +
  geom_edge_arc(aes(edge_width = abs(weight), color = taxon),
                arrow = arrow(length = unit(4, "mm")),
                end_cap = circle(4, "mm"),
                start_cap = circle(4, "mm"),
                alpha = 0.8,  
                curvature = 0.25, 
                show.legend = FALSE) +
  geom_node_point(aes(color = Taxon, size = size)) +
  geom_node_text(aes(label = node_labels,  x = 1.2 * x,
                     y = 1.2 * y)) +  # Labels
  scale_edge_width(range = c(0.5, 3)) +
  scale_edge_color_manual(values = taxon_colors) +
  scale_color_manual(values = taxon_colors, 
                     labels=c(expression(paste(italic("Baetidae")," spp.")), 
                              expression(paste(italic("Chironomidae"), " spp.")), 
                              expression(italic("G. lacustris")), 
                              expression(paste(italic("Hydropsyche"), " spp.")), 
                              expression(italic("P. antipodarum")))) +
  scale_size_continuous(range = c(2, 10)) +  # Scale node size
  theme(
    panel.background = element_rect(fill = "white", color = "black"),  
    panel.grid = element_blank(), 
    strip.background = element_rect(fill = "white", color = "black"), 
    strip.text = element_text(color = "black"),  
    legend.key = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "Hydropeaking Sensitivity Analysis with HFE", 
       color = "Taxon", 
       edge_width = "Slope") +
  guides(size = "none") +  # Remove size legend
  facet_wrap(~TemperatureFactor, labeller = labeller(TemperatureFactor = c(
    "1" = "Baseline", 
    "1.1" = "+1°C", 
    "1.2" = "+2°C", 
    "1.5" = "+5°C"
  )))

# Print the plot
print(hyd_mort_net_hyd)

ggsave(filename = "hyd_net.png", hyd_mort_net_hyd, device = "png", 
       width = 6.25 , height = 5.5, units = c("in"), dpi= "retina")


Multispp_mh_sens_hyd %>% filter(Parameter == "hydro.mort_BAET", TemperatureFactor == 1, SensitivityIncrement == -0.001)
