##############################
#Figure Generation for Manuscript 2
#############################
library(ggpubr)
library(tidyverse)
library(patchwork)
library(svglite)
library(dataRetrieval)
library(cowplot)
library(png)
library(grid)
library(lm.beta)
library(effectsize)
library(flextable)
library(officer)
# Set global options for digits
options(digits = 3)


source("1spFunctions.R")
temp <- readNWISdv("09380000", "00010", "2015-10-01", "2022-05-01")

temperature <-ggplot(temp, aes(x = as.Date(Date), y = X_00010_00003)) +
  geom_line(aes(color = "Colorado \n River at \n Lee's \n Ferry"), linewidth = 1) +  # Correct legend mapping
  labs(x = " ", y = "Temperature", title = "Water Temperature", color = " ") +  # Legend title
  theme_bw() +
  scale_color_manual(values = c("Colorado \n River at \n Lee's \n Ferry" = "black")) +  # Assign legend color
  scale_x_date(date_labels = "%Y") +
  theme(
    text = element_text(size = 13),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12.5), 
    axis.text.y = element_text(size = 13)
  )

  
# correlation plots
x11()
Fig1 <- ggarrange(NZMSts, BAETts, GAMMts, HYOSts, CHIRts, temperature, 
          labels = c("a", "b", "c", "d", "e", "f"),
          ncol = 2, nrow = 3, common.legend =F)
ggsave(filename = "fig3.1.png", Fig1, device = "png", dpi = "retina", height = 15, width = 10)
# # hydropeaking intensity 
# source("BAET_Hydropeaking.R")
# source("HYOS_Hydropeaking.R")
# source("NZMS_Hydropeaking.R")
# source("GAMM_Hydropeaking.R")
# source("CHIR_Hydropeaking.R")
# 
# hyos_hydropeaking <- read_csv("hyos_hydropeaking_results.csv")
# hyos_hydropeaking$Taxa <- rep("HYOS", times = length(hyos_hydropeaking$Hydropeak))
# chir_hydropeaking <- read_csv("chir_hydropeaking_results.csv")
# chir_hydropeaking$Taxa <- rep("CHIR", times = length(chir_hydropeaking$Hydropeak))
# baet_hydropeaking_results <- read_csv("baet_hydropeaking_results.csv")
# baet_hydropeaking_results$Taxa <- rep("BAET", times = length(baet_hydropeaking_results$Hydropeak))
# nzms_hydropeaking <- read_csv("nzms_hydropeaking_results.csv")
# nzms_hydropeaking$Taxa <- rep("NZMS", times = length(nzms_hydropeaking$Hydropeak))
# gamm_hydropeaking <- read_csv("gamm_hydropeaking_results.csv")
# gamm_hydropeaking$Taxa <- rep("GAMM", times = length(gamm_hydropeaking$Hydropeak))
# #Hydro_Abund <- as.data.frame(rbind(BAET_hyd_means, HYOS_hyd_means, NZMS_hyd_means, CHIR_hyd_means, GAMM_hyd_means))
# Hydro_Abund <- as.data.frame(rbind(baet_hydropeaking_results, hyos_hydropeaking, nzms_hydropeaking, chir_hydropeaking, gamm_hydropeaking))
# # combine all abundance data
# Hydro_Abund$Hydropeak <- as.numeric(Hydro_Abund$Hydropeak)
# 
# # make sure in the correct format
# Hydro_Abund$MeanAbund <- as.numeric(Hydro_Abund$MeanAbund)
# Hydro_Abund$SdAbund <- as.numeric(Hydro_Abund$SdAbund)
# Hydro_Abund$Taxa <- as.factor(Hydro_Abund$Taxa)
# #plot
# hyd_abund <- ggplot(data = Hydro_Abund, aes(Hydropeak, log(MeanAbund), group = Taxa, color = Taxa))+ 
#   geom_line(linewidth = 1, alpha = 0.8)+
#   geom_ribbon(data = Hydro_Abund, aes(ymin = log(MeanAbund - SdAbund),
#                   ymax = log(MeanAbund + SdAbund), fill = Taxa),
#                          alpha = .15,
#                          color= "transparent",
#                          show.legend = F) +
#   scale_color_manual(name = "Taxon", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")), expression(paste(italic("Hydropsyche"), " spp.")), expression(italic("P. antipodarum"))), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
#   scale_fill_manual(name = "Taxon", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
#   geom_vline(aes(xintercept = 0.01), linetype = "dotted", 
#               linewidth=1)+
#   geom_vline(aes(xintercept = 0.17 ), linetype="dashed", 
#               linewidth=1)+
#   geom_vline(aes(xintercept = 0.55 ), linetype = "dotdash", 
#             linewidth=1)+
#   scale_y_continuous(breaks = scales::pretty_breaks(n = 4), limit = c(0, 11))+
#   theme_bw()+
#   xlab("Hydropeaking Intensity")+
#   ylab("Log Abundance per Reach")+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# 
# # read in mean timestep biomass data
# # Hydro_TS_Biomass <- rbind(BAET_hyd_size, HYOS_hyd_size, NZMS_hyd_size, CHIR_hyd_size, GAMM_hyd_size)
# # # make sure in correct format
# # Hydro_TS_Biomass$hydropeak <- as.numeric(Hydro_TS_Biomass$hydropeak)
# # Hydro_TS_Biomass$sizemeans <- as.numeric(Hydro_TS_Biomass$sizemeans)
# # Hydro_TS_Biomass$sizesd <- as.numeric(Hydro_TS_Biomass$sizesd)
# # Hydro_TS_Biomass$V4 <- as.factor(Hydro_TS_Biomass$V4)
# #plot
# hyd_ts <- ggplot(data = Hydro_Abund, aes(Hydropeak, log(SizeMean), group = Taxa, color = Taxa)) + 
#   geom_ribbon(aes(ymin = log(SizeMean - SizeSd),
#                     ymax = log(SizeMean + SizeSd), fill = Taxa),
#                 color = "transparent",
#                 alpha = .15,
#                 show.legend = F) +
#   geom_line(linewidth = 1, alpha = 0.8)+
#   scale_color_manual(name = "Taxon", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")), expression(paste(italic("Hydropsyche"), " spp.")), expression(italic("P. antipodarum"))), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
#   scale_fill_manual(name = "Taxon", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
#   geom_vline(aes(xintercept = 0.01), linewidth = 1, linetype = "dotted")+
#   geom_vline(aes(xintercept = 0.17 ), linewidth = 1, linetype="dashed")+ 
#   geom_vline(aes(xintercept = 0.55 ), linewidth = 1, linetype = "dotdash")+
#      theme_bw()+
#     xlab("Hydropeaking Intensity")+
#    ylab("Log Average Timestep Biomass (mg) per Reach")+
#   scale_y_continuous(breaks = scales::pretty_breaks(n = 4), limit = c(0, 12))+
#    theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#     axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
#                                    
# # read in mean annual productivity/biomass data
# #Hydro_Yr_Prod <- as.data.frame(rbind(BAET_hyd_yrprod, HYOS_hyd_yrprod, CHIR_hyd_yrprod, NZMS_hyd_yrprod, GAMM_hyd_yrprod))                                 
# # make sure in correct format
# 
# # Hydro_Yr_Prod$hydropeak <- as.numeric(Hydro_Yr_Prod$hydropeak)
# # Hydro_Yr_Prod$S3Yrprod <- as.numeric(Hydro_Yr_Prod$S3Yrprod)
# # Hydro_Yr_Prod$V3 <- as.factor(Hydro_Yr_Prod$V3)
# 
# Hyd_yr <- ggplot(data = Hydro_Abund, aes(Hydropeak, log(S3Yrprod), group = Taxa, color = Taxa)) + 
#   geom_ribbon(aes(ymin = log(S3Yrprod - S3Yrprodsd),
#                     ymax = log(S3Yrprod + S3Yrprodsd), fill = Taxa),
#                 color = "transparent",
#                 alpha = .15,
#                 show.legend = F) +
#   geom_line(linewidth = 1, alpha = 0.8)+
#   scale_color_manual(name = "Taxon", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")), expression(paste(italic("Hydropsyche"), " spp.")), expression(italic("P. antipodarum"))), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
#   scale_fill_manual(name = "Taxon", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
#   geom_vline(aes(xintercept = 0.01), linewidth = 1, linetype = "dotted")+
#   geom_vline(aes(xintercept = 0.17 ), linewidth = 1, linetype="dashed")+ 
#   geom_vline(aes(xintercept = 0.55 ), linewidth = 1, linetype = "dotdash")+
#   theme_bw()+
#   scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
#   xlab("Hydropeaking Intensity")+
#   ylab("Log Annual Emergent Biomass (mg) per Reach")+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# 
# Fig3.4 <- ggarrange(hyd_abund, hyd_ts, Hyd_yr,
#           labels = c("a", "b", "c"),
#           ncol = 3, nrow = 1, common.legend = T)

#------------------------------------------------
# read in Lees Ferry temp and discharge data from 2007 to 2023
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
temps <- average.yearly.temp(tempdata = temp, temp.column_name = "X_00010_00003", date.column_name = "Date")
# create a timeseries of average temperatures 100 years long
temps <- rep.avg.year(temps, n = 100, change.in.temp = 0, years.at.temp = 0)


# Define the temperature multipliers and labels
scenarios <- data.frame(
  multiplier = c(1, 1.1, 1.2, 1.5),
  label = c("base", "+10%", "+20%", "+50%")
)

# Apply modifications to create the different temperature regimes
temperature.regime <- map2_dfr(scenarios$multiplier, scenarios$label, ~ {
  temps[1:25, ] %>%
    mutate(Temperature = Temperature * .x, se = .y)
})

regime <- ggplot(data = temperature.regime, aes(x = as.Date(dts), y = Temperature, color = se))+
  geom_line(size = 1)+
  theme_bw()+
  scale_color_manual(name = "Warming Regime", labels=c("+10%", "+20%", "+50%", "Baseline"), values=c("#abd9e9", "#fdae61","#d7191c","#2c7bb6"))+
  scale_x_date(date_labels="%b", date_breaks  ="4 month")+
  xlab("")+
  theme(text = element_text(size = 11), axis.text.x = element_text(hjust = 1, size = 10, angle = 0), 
        axis.text.y = element_text(size = 10), legend.key = element_rect(fill = "transparent"))

hyos_abund <- read_csv("HYOS_temp_abund.csv")
baet_abund <- read_csv("BAET_temp_abund.csv")
chir_abund <- read_csv("CHIR_temp_abund.csv")
gamm_abund <- read_csv("GAMM_temp_abund.csv")
nzms_abund <- read_csv("NZMS_temp_abund.csv")
# 
#----------------
# with spike

# read in LF temp and discharge data from 2007 to 2023
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
# calculate average yearly temperatures
temps <- average.yearly.temp(tempdata = temp, temp.column_name = "X_00010_00003", date.column_name = "Date")
# create summertime spike (up to 21 C, then scale from there)
temps$Temperature[16:22] <- c(14, 16, 18, 21, 21, 18, 16, 14)

# create a timeseries of average temperatures 100 years long
temps <- rep.avg.year(temps, n = 100, change.in.temp = 0, years.at.temp = 0)

# Define the temperature multipliers and labels
scenarios <- data.frame(
  multiplier = c(1, 1.1, 1.2, 1.5),
  label = c("base", "+10%", "+20%", "+50%")
)

# Apply modifications to create the different temperature regimes
temperature.regime <- map2_dfr(scenarios$multiplier, scenarios$label, ~ {
  temps[1:25, ] %>%
    mutate(Temperature = Temperature * .x, se = .y)
})



#add temperature spike in September
spikeregime <- ggplot(data = temperature.regime, aes(x = as.Date(dts), y = Temperature, color = se))+
  geom_line(size = 1)+
  theme_bw()+
  scale_color_manual(name = "Warming Regime", labels=c("+10%", "+20%", "+50%", "Baseline"), values=c("#abd9e9", "#fdae61","#d7191c","#2c7bb6"))+
  scale_x_date(date_labels="%b", date_breaks  ="4 month")+
  xlab("")+
  theme(text = element_text(size = 11), axis.text.x = element_text(hjust = 1, size = 10, angle = 0), 
        axis.text.y = element_text(size = 10), legend.key = element_rect(fill = "transparent"))

########################
source("GAMMSurvivorship.R")
source("HYOSSurvivorship.R")
source("BAETSurvivorship.R")
source("NZMSSurvivorship.R")
source("CHIRSurvivorship.R")

gamm_surv <- ggplot(data = GAMMSurvRates, aes(x = Temperature, y = Survival, color = "Taxa"))+
  geom_point(aes(color = "data"), size = 1.5, alpha = 1, show.legend = T)+
  geom_line(data = tempGamm, aes(x = tem, y = V2), linewidth = 1, show.legend = T, alpha = 0.8)+
  theme_bw()+
  ylim(c(0,1))+
  ylab("G. lacustris Survival")+
  scale_color_manual(name = " ", labels=c("Empirical Data", "Fitted"), values=c("#CCBB44","black"))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13),legend.position = "bottom", legend.key = element_rect(fill = "transparent"))

hyos_surv <- ggplot(data = HYOSSurvRates, aes(x = Temperature, y = Survival, color = "Taxa"))+
  geom_point(aes(color = "data"), size = 1.5, alpha = 1, show.legend = T)+
  geom_line(data = tempHyos, aes(x = tem, y = V2), linewidth = 1, show.legend = T, alpha = 0.8)+
  theme_bw()+
  ylab("Hydropsyche spp. Survival")+
  scale_color_manual(name = " ", labels=c("Empirical Data", "Fitted"), values=c("#AA3377","black"))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13),legend.position = "bottom", legend.key = element_rect(fill = "transparent"))

baet_surv <- ggplot(data = BAETSurvRate, aes(x = Temperature, y = Survival, color = "Taxa"))+
  geom_point(aes(color = "data"), size = 1.5, alpha = 1, show.legend = T)+
  geom_line(data = tempBaet, aes(x = tem, y = V2), linewidth = 1, show.legend = T, alpha = 0.8)+
  theme_bw()+
  ylab("Baetidae spp. Survival")+
  scale_color_manual(name = " ", labels=c("Empirical Data", "Fitted"), values=c("#66CCEE","black"))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13),legend.position = "bottom", legend.key = element_rect(fill = "transparent"))

nzms_surv <- ggplot(data = NZMSSurvRates, aes(x = Temperature, y = Survival, color = "Taxa"))+
  geom_point(aes(color = "data"), size = 1.5, alpha = 1, show.legend = T)+
  geom_line(data = tempHyos, aes(x = tem, y = V2), linewidth = 1, show.legend = T, alpha = 0.8)+
  theme_bw()+
  ylab("P. antipordarum Survival")+
  scale_color_manual(name = " ", labels=c("Empirical Data", "Fitted"), values=c("#4477AA","black"))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13),legend.position = "bottom", legend.key = element_rect(fill = "transparent"))

chir_surv <- ggplot(data = CHIRSurvRate, aes(x = Temp, y = Survival, color = "Taxa"))+
  geom_point(aes(color = "data"), size = 1.5, alpha = 1, show.legend = T)+
  geom_line(data = tempHyos, aes(x = tem, y = V2), linewidth = 1, show.legend = T, alpha = 0.8)+
  theme_bw()+
  ylim(c(0,1))+
  ylab("Chironomidae spp. Survival")+
  xlab("Temperature")+
  scale_color_manual(name = " ", labels=c("Empirical Data", "Fitted"), values=c("#228833","black"))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13),legend.position = "bottom", legend.key = element_rect(fill = "transparent"))

tempsurv_plot <- ggarrange(nzms_surv, baet_surv, gamm_surv, hyos_surv, chir_surv,  
          labels = c("a", "b", "c", "d", "e"),
          ncol = 3, nrow = 2, common.legend =F)

ggsave("TempSurvPlot.png", plot = tempsurv_plot, device = "png", width = 8, height = 8, dpi = "retina")
# plot(GAMMSurvRates$Temperature, GAMMSurvRates$Survival, col = "red", pch = 16, xlab = "Temperature", ylab = "Survival", xlim = c(0,40), ylim = c(0, 1))
# lines(tem, TempSurv_GAMM(tem))
#http://127.0.0.1:35901/graphics/63814419-21c3-4b1a-a8b6-ed037422e6cd.png

######
# taxon by taxon 
# Combine them with a new column to differentiate

# hydropysche spp
hyos_abund <- read_csv("HYOS_temp_abund.csv")
hyos_pc <- read_csv("HYOS_temp_percapita.csv")
hyos_biomass <- read_csv("HYOS_temp_biomass.csv")
hyos_abund_sp <- read_csv("HYOS_temp_abund_spike.csv")
hyos_pc_sp <- read_csv("HYOS_temp_percapita_spike.csv")
hyos_biomass_sp <- read_csv("HYOS_temp_biomass_spike.csv")
hyos_abund_t.h <- read_csv("HYOS_temp_hyd_abund.csv")
hyos_biomass_t.h <- read_csv("HYOS_temp_hyd_biomass.csv")
hyos_abund_t.h_spike <- read_csv("HYOS_temp_hyd_abund_spike.csv")
hyos_biomass_t.h_spike <- read_csv("HYOS_temp_hyd_biomass_spike.csv")
hyos_abund_HFE <- read_csv("HYOS_temp_abund_HFE.csv")
hyos_biomass_HFE <- read_csv("HYOS_temp_biomass_HFE.csv")
hyos_biomass_HFE_sp <- read_csv("HYOS_temp_biomass_HFE_spike.csv")
hyos_abund_HFE_sp <- read_csv("HYOS_temp_abund_HFE_spike.csv")
hyos_abund_hyd_HFE <- read_csv("HYOS_temp_hyd_abund_HFE.csv")
hyos_biomass_hyd_HFE <- read_csv("HYOS_temp_hyd_biomass_HFE.csv")
hyos_abund_hyd_HFE.sp <- read_csv("HYOS_temp_hyd_abund_HFE_spike.csv")
hyos_biomass_hyd_HFE.sp <- read_csv("HYOS_temp_hyd_biomass_HFE_spike.csv")
hyos_pc_hfe <- read_csv("HYOS_temp_percapita_HFE.csv")
hyos_pc_hfe_sp <- read_csv("HYOS_temp_percapita_HFE_spike.csv")
hyos_pa <- read_csv("HYOS_temp_propadult.csv")
hyos_pa_hfe <- read_csv("HYOS_temp_propadult_HFE.csv")
hyos_pa_sp <- read_csv("HYOS_temp_propadult_spike.csv")
hyos_pa_hfe_sp <- read_csv("HYOS_temp_propadult_HFE_spike.csv")




# abundance
hyos_abund_combo <- bind_rows(
  hyos_abund %>% mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  hyos_abund_sp %>% mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  hyos_abund_t.h %>% mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  hyos_abund_t.h_spike %>% mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  hyos_abund_HFE %>% mutate(source = "Temperature & HFE",  temp_spike = 0, hydropeaking = 0, HFE = 1),
  hyos_abund_HFE_sp %>% mutate(source = "Temperature & spike & HFE",  temp_spike = 1, hydropeaking = 0, HFE = 1),
  hyos_abund_hyd_HFE %>% mutate(source = "Temperature & hydropeaking & HFE",  temp_spike = 0, hydropeaking = 1, HFE = 1),
  hyos_abund_hyd_HFE.sp %>% mutate(source = "Temperature & spike & hydropeaking & HFE",  temp_spike = 1, hydropeaking = 1, HFE = 1)
)

# biomass
hyos_biomass_combo <- bind_rows(
  hyos_biomass %>% mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  hyos_biomass_sp %>% mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  hyos_biomass_t.h %>% mutate(source = "Temperature & spike & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  hyos_biomass_t.h_spike %>% mutate(source = "Temperature & HFE",  temp_spike = 1, hydropeaking = 1, HFE = 0),
  hyos_biomass_HFE %>% mutate(source = "Temperature & spike & HFE",  temp_spike = 0, hydropeaking = 0, HFE = 1),
  hyos_biomass_HFE_sp %>% mutate(source = "Temperature & hydropeaking & HFE",  temp_spike =1, hydropeaking = 0, HFE = 1),
  hyos_biomass_hyd_HFE %>% mutate(source = "Temperature & hydropeaking & HFE",  temp_spike = 0, hydropeaking = 1, HFE = 1),
  hyos_biomass_hyd_HFE.sp %>% mutate(source = "Temperature & spike & hydropeaking & HFE",  temp_spike = 1, hydropeaking = 1, HFE = 1)
)

hyos_pc_combo <- bind_rows(
  hyos_pc %>% mutate(source = "Temperature"),
  hyos_pc_sp %>% mutate(source = "Temperature & spike"),
)


paired.colors <- c("#56B4E9", "#365C8DFF",
                   "#FDE725FF","#E69F00", 
                   "#9FDA3AFF","#009E73", 
                   "#CC79A7","#46337EFF")

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

hyos_abund_combo$abundance <- scale(hyos_abund_combo$abundance)
hyos_abund_combo$source <- ordered(hyos_abund_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE", 
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"))  
h.abund <- ggplot(data = hyos_abund_combo, aes(y = abundance, x = as.factor(temperature), color = source, group = source)) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, 
               position = position_dodge(0.5), show.legend = F) +
  stat_summary(fun = mean, geom = "point", size = 3, position = position_dodge(0.5), show.legend = F) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels=label.names, values=paired.colors)+
  xlab("Temperature")+
  labs(y=expression(paste(italic("Hydropsyche"), "  Scaled Abundance")))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# hyos biomass
hyos_biomass_combo$biomass <- scale(hyos_biomass_combo$biomass)
hyos_biomass_combo$source <- ordered(hyos_biomass_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE", 
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"))  
h.biomass <- ggplot(data = hyos_biomass_combo, aes(y = biomass, x = as.factor(temperature), color = source, group = source)) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, 
               position = position_dodge(0.5)) +
  stat_summary(fun = mean, geom = "point", size = 3, position = position_dodge(0.5)) +
  stat_summary(fun = mean, geom = "line", size = 0.75, position = position_dodge(0.5)) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels=label.names, values=paired.colors)+
  xlab("Temperature")+
  labs(y=expression(paste(italic("Hydropsyche"), " spp. Scaled Biomass")))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 11), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

h_plot <- wrap_plots(h.abund, h.biomass) +
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides = 'collect') &
  theme(
    legend.text = element_text(size = 11, lineheight = 0.9),
    legend.spacing.y = unit(0.5, "cm"),
    legend.key.height = unit(1.25, "cm")
  )
# Read the PNG image
caddis <- readPNG("caddisfly.png")
caddis_grob <- rasterGrob(caddis, x = 0.78, y = 0.92, width = 0.135, height = 0.135)  # Adjust x, y, width, and height

# Overlay the image on the plot
h.fig <- ggdraw() +
  draw_plot(h_plot) +
  draw_grob(caddis_grob)

ggsave("hfig.png", plot = h.fig, device = "png", width = 8.44, height = 6, units = c("in"), dpi = "retina")

h.pc <- ggplot(data = hyos_pc_combo, aes(x = as.factor(temperature), y = percapita, group = source, color = source))+
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 3)+
  stat_summary(fun = mean, geom = "line", size = 0.75) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels= c("Temperature", 
                                           "Temperature +\nSpike"), values= c("#56B4E9", "#365C8DFF"))+
  xlab("Temperature")+
  guides(color=guide_legend(title="Scenario"))+
  labs(y="Adult Biomass (mg)", title = expression(paste(italic("Hydropsyche"), " spp.")))+
  theme(text = element_text(size = 11), axis.text.x = element_text(hjust = 1, angle=45, size = 10), 
        axis.text.y = element_text(size = 11), legend.key = element_rect(fill = "transparent"))


# Fit linear models with ALL interactions (full factorial design)
hyosmod_biomass <- lm(biomass ~ (temperature + temp_spike + hydropeaking + HFE)^4, data = hyos_biomass_combo)
hyosmod_abund <- lm(abundance ~ (temperature + temp_spike + hydropeaking + HFE)^4, data = hyos_abund_combo)
sum_biomass_lm <- lm.beta(hyosmod_biomass)
sum_abund_lm <- lm.beta(hyosmod_abund)
# Extract the standardized coefficients
#lm.beta stores the standardized coefficients in the model object as 'standardized.coefficients'
std_coefs.b <- sum_biomass_lm$standardized.coefficients
std_coefs.a <- sum_abund_lm$standardized.coefficients
# Remove the intercept (if present) so we only have the effects and interactions
std_coefs.b <- std_coefs.b[-1]
std_coefs.a <- std_coefs.a[-1]

# # Create a data frame with effect names and their beta values
coef_df.b <- data.frame(
  Stressor = names(std_coefs.b),
  Beta = std_coefs.b,
  AbsBeta = abs(std_coefs.b)  # for sorting purposes
)


coef_df.a <- data.frame(
  Stressor = names(std_coefs.a),
  Beta = std_coefs.a,
  AbsBeta = abs(std_coefs.a)  # for sorting purposes
)
# Sort the data frame by absolute beta values in decreasing order
coef_df_sorted.hyos.b <- coef_df.b %>% arrange(desc(AbsBeta))
coef_df_sorted.hyos.a <- coef_df.a %>% arrange(desc(AbsBeta))

# Define effect size categories
coef_df_sorted.hyos.b$Effect_Size <- cut(abs(coef_df_sorted.hyos.b$Beta),
                                  breaks = c(0, 0.1, 0.3, 0.5, Inf),
                                  labels = c("Small", "Medium", "Large", "Very Large"),
                                  right = FALSE)

# Define effect size categories
coef_df_sorted.hyos.a$Effect_Size <- cut(abs(coef_df_sorted.hyos.a$Beta),
                                    breaks = c(0, 0.1, 0.3, 0.5, Inf),
                                    labels = c("Small", "Medium", "Large", "Very Large"),
                                    right = FALSE)

# Now plot the sorted coefficients using ggplot2
# ggplot(coef_df_sorted.hyos.b, aes(x = reorder(Effect, Beta), y = Beta, fill = Effect_Size)) +
#   geom_bar(stat = "identity") +
#   coord_flip() +  # Flip for better readability
#   labs(title = "Standardized Beta Coefficients",
#        x = "Predictor",
#        y = "Standardized Beta",
#        fill = "Effect Size") +
#   scale_fill_manual(values = c("Small" = "blue", "Medium" = "orange", "Large" = "red", "Very Large" = "darkred")) +
#   theme_minimal()

# Now plot the sorted coefficients using ggplot2
# ggplot(coef_df_sorted.hyos.a, aes(x = reorder(Effect, Beta), y = Beta, fill = Effect_Size)) +
#   geom_bar(stat = "identity") +
#   coord_flip() +  # flips the axes so the labels are easier to read
#   labs(
#     title = "Standardized Beta Coefficients for Abundance",
#     x = "Effect",
#     y = "Standardized Beta Coefficient"
#   ) +
#   scale_fill_manual(values = c("Small" = "blue", "Medium" = "orange", "Large" = "red", "Very Large" = "darkred")) +
#   theme_bw()


# because if categorical, fit manova
hyosmod_biomass <- aov(biomass ~ (temperature + temp_spike + hydropeaking + HFE)^4, data = hyos_biomass_combo)
hyosmod_abund <- aov(abundance ~ (temperature + temp_spike + hydropeaking + HFE)^4, data = hyos_abund_combo)

# Extract partial eta squared values
eta_sq_df.hyos.b <- eta_squared(hyosmod_biomass, partial = TRUE) %>%
  as.data.frame() %>%
  arrange(desc(Eta2_partial))  # Order from largest to smallest effect

eta_sq_df.hyos.a <- eta_squared(hyosmod_abund, partial = TRUE) %>%
  as.data.frame() %>%
  arrange(desc(Eta2_partial))  # Order from largest to smallest effect
# Create the bar plot
ggplot(eta_sq_df.hyos.b, aes(x = reorder(Parameter, Eta2_partial), y = Eta2_partial, fill = Eta2_partial)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +  # Flip for better readability
  scale_fill_viridis_c() +  # Color scale
  labs(title = "Hydropsyche spp.",
       x = "Stressor(s)",
       y = "Partial Eta Squared") +
  theme_bw()

ggplot(eta_sq_df.hyos.a, aes(x = reorder(Parameter, Eta2_partial), y = Eta2_partial, fill = Eta2_partial)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +  # Flip for better readability
  scale_fill_viridis_c() +  # Color scale
  labs(title = "Hydropsyche spp.",
       x = "Stressor(s)",
       y = "Partial Eta Squared") +
  theme_bw()


# Merge Dataframes
results_df.hyos.b <- left_join(eta_sq_df.hyos.b, coef_df_sorted.hyos.b, by = c("Parameter" = "Stressor")) %>%
  mutate(Direction = ifelse(Beta > 0, "Positive", "Negative"))

# Sort by Partial Eta Squared (Descending)
results_df.hyos.b <- results_df.hyos.b %>% arrange(desc(Eta2_partial))
write.table(results_df.hyos.b, file = "partialetasquare_hyos_biomass.csv", sep=",")

# Format Table with `flextable`
results_hyos.b <- flextable(results_df.hyos.b) %>%
  theme_vanilla() %>%  # Clean formatting
  set_header_labels(
    Parameter = "Stressor(s)",
    Eta2_partial = "Partial Eta²",
    Beta = "β Coefficient",
    Direction = "Effect Direction"
  ) %>%
  colformat_num(j = c("Eta2_partial", "Beta")) %>%  
  autofit()  


# Merge Dataframes
results_df.hyos.a <- left_join(eta_sq_df.hyos.a, coef_df_sorted.hyos.a, by = c("Parameter" = "Stressor")) %>%
  mutate(Direction = ifelse(Beta > 0, "Positive", "Negative"))

# Sort by Partial Eta Squared (Descending)
results_df.hyos.a <- results_df.hyos.a %>% arrange(desc(Eta2_partial))
write.table(results_df.hyos.a, file = "partialetasquare_hyos_abund.csv", sep= " , ")

# Format Table with `flextable`
results_hyos.a <- flextable(results_df.hyos.a) %>%
  theme_vanilla() %>%  # Clean formatting
  set_header_labels(
    Parameter = "Stressor(s)",
    Eta2_partial = "Partial Eta²",
    Beta = "β Coefficient",
    Direction = "Effect Direction"
  ) %>%
  colformat_num(j = c("Eta2_partial", "Beta")) %>%  
  autofit()  # Adjust column widths


## baetidae
baet_abund <- read_csv("BAET_temp_abund.csv")
baet_pc <- read_csv("BAET_temp_percapita.csv")
baet_biomass <- read_csv("BAET_temp_biomass.csv")
baet_abund_sp <- read_csv("BAET_temp_abund_spike.csv")
baet_pc_sp <- read_csv("BAET_temp_percapita_spike.csv")
baet_abund_t.h <- read_csv("BAET_temp_hyd_abund.csv")
baet_biomass_sp <- read_csv("BAET_temp_biomass_spike.csv")
baet_biomass_t.h <- read_csv("BAET_temp_hyd_biomass.csv")
baet_abund_t.h_spike <- read_csv("BAET_temp_hyd_abund_spike.csv")
baet_biomass_t.h_spike <- read_csv("BAET_temp_hyd_biomass_spike.csv")
baet_abund_HFE <- read_csv("BAET_temp_abund_HFE.csv")
baet_biomass_HFE <- read_csv("BAET_temp_biomass_HFE.csv")
baet_abund_HFE_sp <- read_csv("BAET_temp_abund_HFE_spike.csv")
baet_biomass_HFE.sp <- read_csv("BAET_temp_biomass_HFE_spike.csv")
baet_abund_HFE_sp <- read_csv("BAET_temp_abund_HFE_spike.csv")
baet_abund_hyd_HFE <- read_csv("BAET_temp_hyd_abund_HFE.csv")
baet_biomass_hyd_HFE <- read_csv("BAET_temp_hyd_biomass_HFE.csv")
baet_abund_hyd_HFE.sp <- read_csv("BAET_temp_hyd_abund_HFE_spike.csv")
baet_biomass_hyd_HFE.sp <- read_csv("BAET_temp_hyd_biomass_HFE_spike.csv")


# abundance
baet_abund_combo <- bind_rows(
  baet_abund %>% mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  baet_abund_sp %>% mutate(source = "Temperature", temp_spike = 1, hydropeaking = 0, HFE = 0),
  baet_abund_t.h %>% mutate(source = "Temperature", temp_spike = 0, hydropeaking = 1, HFE = 0),
  baet_abund_t.h_spike %>% mutate(source = "Temperature", temp_spike =1, hydropeaking = 1, HFE = 0),
  baet_abund_HFE %>% mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 1),
  baet_abund_HFE_sp %>% mutate(source = "Temperature", temp_spike = 1, hydropeaking = 0, HFE = 1),
  baet_abund_hyd_HFE %>% mutate(source = "Temperature", temp_spike = 0, hydropeaking = 1, HFE = 1),
  baet_abund_hyd_HFE.sp %>% mutate(source = "Temperature", temp_spike = 1, hydropeaking = 1, HFE = 1),
)

# biomass
baet_biomass_combo <- bind_rows(
  baet_biomass %>% mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  baet_biomass_sp %>% mutate(source = "Temperature", temp_spike = 1, hydropeaking = 0, HFE = 0),
  baet_biomass_t.h %>% mutate(source = "Temperature", temp_spike = 0, hydropeaking = 1, HFE = 0),
  baet_biomass_t.h_spike %>%mutate(source = "Temperature", temp_spike =1, hydropeaking = 1, HFE = 0),
  baet_biomass_HFE %>%mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 1),
  baet_biomass_HFE.sp %>%mutate(source = "Temperature", temp_spike = 1, hydropeaking = 0, HFE = 1),
  baet_biomass_hyd_HFE %>% mutate(source = "Temperature", temp_spike = 0, hydropeaking = 1, HFE = 1),
  baet_biomass_hyd_HFE.sp %>% mutate(source = "Temperature", temp_spike = 1, hydropeaking = 1, HFE = 1),
)

baet_pc_combo <- bind_rows(
  baet_pc %>% mutate(source = "Temperature"),
  baet_pc_sp %>% mutate(source = "Temperature & spike"),
)

baet_abund_combo$abundance <- scale(baet_abund_combo$abundance) 
baet_abund_combo$source <- ordered(baet_abund_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE", 
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"))  
b.abund <- ggplot(data = baet_abund_combo, aes(y = abundance, x = as.factor(temperature), color = source, group = source)) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, 
               position = position_dodge(0.5), show.legend = F) +
  stat_summary(fun = mean, geom = "point", size = 3, position = position_dodge(0.5), show.legend = F) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels=label.names, values=paired.colors)+
  xlab("Temperature")+
  labs(y=expression(paste(italic("Baetidae"), "  Scaled Abundance")))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# hyos biomass
baet_biomass_combo$biomass <- scale(baet_biomass_combo$biomass)
baet_biomass_combo$source <- ordered(baet_biomass_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE", 
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"))
b.biomass <- ggplot(data = baet_biomass_combo, aes(y = biomass, x = as.factor(temperature), color = source, group = source)) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, 
               position = position_dodge(0.5)) +
  stat_summary(fun = mean, geom = "point", size = 3, position = position_dodge(0.5)) +
  stat_summary(fun = mean, geom = "line", size = 0.75, position = position_dodge(0.5)) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels=label.names, values=paired.colors)+
  xlab("Temperature")+
  labs(y=expression(paste(italic("Baetidae"), " spp. Scaled Biomass")))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 11), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

b_plot <- wrap_plots(b.abund, b.biomass) +
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides = 'collect') &
  theme(
    legend.text = element_text(size = 11, lineheight = 0.9),
    legend.spacing.y = unit(0.5, "cm"),
    legend.key.height = unit(1.25, "cm")
  )
# Read the PNG image
mayfly <- readPNG("mayfly.png")
mayfly_grob <- rasterGrob(mayfly, x = 0.825, y = 0.92, width = 0.135, height = 0.135)  # Adjust x, y, width, and height

# Overlay the image on the plot
b.fig <- ggdraw() +
  draw_plot(b_plot) +
  draw_grob(mayfly_grob)

ggsave("bfig.png", plot = b.fig, device = "png", width = 8.44, height = 6, units = c("in"), dpi = "retina")


baet_pc_combo <- bind_rows(
  baet_pc %>% mutate(source = "Temperature"),
  baet_pc_sp %>% mutate(source = "Temperature & spike"),
)
b.pc <- ggplot(data = baet_pc_combo, aes(x = as.factor(temperature), y = percapita, group = source, color = source))+
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, show.legend = F) +
  stat_summary(fun = mean, geom = "point", size = 3, show.legend = F)+
  stat_summary(fun = mean, geom = "line", size = 0.75, show.legend = F) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels= c("Temperature", 
                                           "Temperature +\nSpike"), values= c("#56B4E9", "#365C8DFF"))+
  xlab("Temperature")+
  labs(y="Adult Biomass (mg)", title = expression(paste(italic("Baetidae"), " spp.")))+
  theme(text = element_text(size = 11), axis.text.x = element_text(hjust = 1, angle=45, size = 10), 
        axis.text.y = element_text(size = 11), legend.key = element_rect(fill = "transparent"))

# Fit linear models with ALL interactions (full factorial design)
baetmod_biomass <- lm(biomass ~ (temperature + temp_spike + hydropeaking + HFE)^4, data = baet_biomass_combo)
baetmod_abund <- lm(abundance ~ (temperature + temp_spike + hydropeaking + HFE)^4, data = baet_abund_combo)

sum_biomass_lm <- lm.beta(baetmod_biomass)
sum_abund_lm <- lm.beta(baetmod_abund)
# Extract the standardized coefficients
# lm.beta stores the standardized coefficients in the model object as 'standardized.coefficients'
std_coefs.b <- sum_biomass_lm$standardized.coefficients
std_coefs.a <- sum_abund_lm$standardized.coefficients
# Remove the intercept (if present) so we only have the effects and interactions
std_coefs.b <- std_coefs.b[-1]
std_coefs.a <- std_coefs.a[-1]

# Create a data frame with effect names and their beta values
coef_df.b <- data.frame(
  Stressor = names(std_coefs.b),
  Beta = std_coefs.b,
  AbsBeta = abs(std_coefs.b)  # for sorting purposes
)


coef_df.a <- data.frame(
  Stressor = names(std_coefs.a),
  Beta = std_coefs.a,
  AbsBeta = abs(std_coefs.a)  # for sorting purposes
)
# Sort the data frame by absolute beta values in decreasing order
coef_df_sorted.baet.b <- coef_df.b %>% arrange(desc(AbsBeta))
coef_df_sorted.baet.a <- coef_df.a %>% arrange(desc(AbsBeta))

# Define effect size categories
coef_df_sorted.baet.b$Effect_Size <- cut(abs(coef_df_sorted.baet.b$Beta),
                                         breaks = c(0, 0.1, 0.3, 0.5, Inf),
                                         labels = c("Small", "Medium", "Large", "Very Large"),
                                         right = FALSE)
# Define effect size categories
coef_df_sorted.baet.a$Effect_Size <- cut(abs(coef_df_sorted.baet.a$Beta),
                                         breaks = c(0, 0.1, 0.3, 0.5, Inf),
                                         labels = c("Small", "Medium", "Large", "Very Large"),
                                         right = FALSE)
# Now plot the sorted coefficients using ggplot2
# ggplot(coef_df_sorted.baet.b, aes(x = reorder(Effect, Beta), y = Beta, fill = Effect_Size)) +
#   geom_bar(stat = "identity") +
#   coord_flip() +  # flips the axes so the labels are easier to read
#   labs(
#     title = "Standardized Beta Coefficients for Biomass",
#     x = "Effect",
#     y = "Standardized Beta Coefficient"
#   ) +
#   scale_fill_manual(values = c("Small" = "blue", "Medium" = "orange", "Large" = "red", "Very Large" = "darkred")) +
#   theme_bw()
# 
# 
# 
# # Now plot the sorted coefficients using ggplot2
# ggplot(coef_df_sorted.baet.a, aes(x = reorder(Effect, Beta), y = Beta, fill = Effect_Size)) +
#   geom_bar(stat = "identity") +
#   coord_flip() +  # flips the axes so the labels are easier to read
#   labs(
#     title = "Standardized Beta Coefficients for Abundance",
#     x = "Effect",
#     y = "Standardized Beta Coefficient"
#   ) +
#   scale_fill_manual(values = c("Small" = "blue", "Medium" = "orange", "Large" = "red", "Very Large" = "darkred")) +
#   theme_bw()

# because if categorical, fit manova
baetmod_biomass <- aov(biomass ~ (temperature + temp_spike + hydropeaking + HFE)^4, data = baet_biomass_combo)
baetmod_abund <- aov(abundance ~ (temperature + temp_spike + hydropeaking + HFE)^4, data = baet_abund_combo)

# Extract partial eta squared values
eta_sq_df.baet.b <- eta_squared(baetmod_biomass, partial = TRUE) %>%
  as.data.frame() %>%
  arrange(desc(Eta2_partial))  # Order from largest to smallest effect

eta_sq_df.baet.a <- eta_squared(baetmod_abund, partial = TRUE) %>%
  as.data.frame() %>%
  arrange(desc(Eta2_partial))  # Order from largest to smallest effect
# Create the bar plot
ggplot(eta_sq_df.baet.b, aes(x = reorder(Parameter, Eta2_partial), y = Eta2_partial, fill = Eta2_partial)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +  # Flip for better readability
  scale_fill_viridis_c() +  # Color scale
  labs(title = "Hydropsyche spp.",
       x = "Stressor(s)",
       y = "Partial Eta Squared") +
  theme_bw()

ggplot(eta_sq_df.baet.a, aes(x = reorder(Parameter, Eta2_partial), y = Eta2_partial, fill = Eta2_partial)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +  # Flip for better readability
  scale_fill_viridis_c() +  # Color scale
  labs(title = "Hydropsyche spp.",
       x = "Stressor(s)",
       y = "Partial Eta Squared") +
  theme_bw()


# Merge Dataframes
results_df.baet.b <- left_join(eta_sq_df.baet.b, coef_df_sorted.baet.b, by = c("Parameter" = "Stressor")) %>%
  mutate(Direction = ifelse(Beta > 0, "Positive", "Negative"))

# Sort by Partial Eta Squared (Descending)
results_df.baet.b <- results_df.baet.b %>% arrange(desc(Eta2_partial))
write.table(results_df.baet.b, file = "partialetasquare_baet_biomass.csv", sep = " , ")

# Format Table with `flextable`
results_baet.b <- flextable(results_df.baet.b) %>%
  theme_vanilla() %>%  # Clean formatting
  set_header_labels(
    Parameter = "Stressor(s)",
    Eta2_partial = "Partial Eta²",
    Beta = "β Coefficient",
    Direction = "Effect Direction"
  ) %>%
  colformat_num(j = c("Eta2_partial", "Beta")) %>%  
  autofit()  # Adjust column widths

# Merge Dataframes
results_df.baet.a <- left_join(eta_sq_df.baet.a, coef_df_sorted.baet.a, by = c("Parameter" = "Stressor")) %>%
  mutate(Direction = ifelse(Beta > 0, "Positive", "Negative"))

# Sort by Partial Eta Squared (Descending)
results_df.baet.a <- results_df.baet.a %>% arrange(desc(Eta2_partial))

write.table(results_df.baet.a, file = "partialetasquare_baet_abund.csv", sep = " , ")
# Format Table with `flextable`
results_baet.a <- flextable(results_df.baet.a) %>%
  theme_vanilla() %>%  # Clean formatting
  set_header_labels(
    Parameter = "Stressor(s)",
    Eta2_partial = "Partial Eta²",
    Beta = "β Coefficient",
    Direction = "Effect Direction"
  ) %>%
  colformat_num(j = c("Eta2_partial", "Beta")) %>%  
  autofit()  # Adjust column widths

## Chironomidae spp
chir_abund <- read_csv("CHIR_temp_abund.csv")
chir_pc <- read_csv("CHIR_temp_percapita.csv")
chir_biomass <- read_csv("CHIR_temp_biomass.csv")
chir_abund_sp <- read_csv("CHIR_temp_abund_spike.csv")
chir_pc_sp <- read_csv("CHIR_temp_percapita_spike.csv")
chir_biomass_sp <- read_csv("CHIR_temp_biomass_spike.csv")
chir_abund_t.h <- read_csv("CHIR_temp_hyd_abund.csv")
chir_biomass_t.h <- read_csv("CHIR_temp_hyd_biomass.csv")
chir_abund_t.h_spike <- read_csv("CHIR_temp_hyd_abund_spike.csv")
chir_biomass_t.h_spike <- read_csv("CHIR_temp_hyd_biomass_spike.csv")
chir_abund_HFE <- read_csv("CHIR_temp_abund_HFE.csv")
chir_biomass_HFE <- read_csv("CHIR_temp_biomass_HFE.csv")
chir_abund_HFE_sp <- read_csv("CHIR_temp_abund_HFE_spike.csv")
chir_biomass_HFE.sp <- read_csv("CHIR_temp_biomass_HFE_spike.csv")
chir_abund_HFE_sp <- read_csv("CHIR_temp_abund_HFE_spike.csv")
chir_abund_hyd_HFE <- read_csv("CHIR_temp_hyd_abund_HFE.csv")
chir_biomass_hyd_HFE <- read_csv("CHIR_temp_hyd_biomass_HFE.csv")
chir_abund_hyd_HFE.sp <- read_csv("CHIR_temp_hyd_abund_HFE_spike.csv")
chir_biomass_hyd_HFE.sp <- read_csv("CHIR_temp_hyd_biomass_HFE_spike.csv")

# abundance
chir_abund_combo <- bind_rows(
  chir_abund %>% mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  chir_abund_sp %>% mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  chir_abund_t.h %>% mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  chir_abund_t.h_spike %>% mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  chir_abund_HFE %>% mutate(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
  chir_abund_HFE_sp %>% mutate(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
  chir_abund_hyd_HFE %>% mutate(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
  chir_abund_hyd_HFE.sp %>% mutate(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1),
)

# biomass
chir_biomass_combo <- bind_rows(
  chir_biomass %>% mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  chir_biomass_sp %>% mutate(source = "Temperature & spike",temp_spike = 1, hydropeaking = 0, HFE = 0),
  chir_biomass_t.h %>% mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  chir_biomass_t.h_spike %>% mutate(source =  "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  chir_biomass_HFE %>% mutate(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
  chir_biomass_HFE.sp %>% mutate(source =  "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
  chir_biomass_hyd_HFE %>% mutate(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
  chir_biomass_hyd_HFE.sp %>% mutate(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1),
)

chir_pc_combo <- bind_rows(
  chir_pc %>% mutate(source = "Temperature"),
  chir_pc_sp %>% mutate(source = "Temperature & spike"),
  #chir_pc_hfe %>% mutate(source = "Temperature & HFE"),
  #chir_pc_hfe_sp %>% mutate(source = "Temperature & spike & HFE")
)

chir_abund_combo$abundance <- scale(chir_abund_combo$abundance) 
chir_abund_combo$source <- ordered(chir_abund_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE", 
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"))
c.abund <- ggplot(data = chir_abund_combo, aes(y = abundance, x = as.factor(temperature), color = source, group = source)) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, 
               position = position_dodge(0.5), show.legend = F) +
  stat_summary(fun = mean, geom = "point", size = 3, position = position_dodge(0.5), show.legend = F) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels=label.names, values=paired.colors)+
  xlab("Temperature")+
  labs(y=expression(paste(italic("Chironomidae"), "  Scaled Abundance")))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# baet biomass
chir_biomass_combo$biomass <- scale(chir_biomass_combo$biomass)
chir_biomass_combo$source <- ordered(chir_biomass_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE", 
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"))
c.biomass <- ggplot(data = chir_biomass_combo, aes(y = biomass, x = as.factor(temperature), color = source, group = source)) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, 
               position = position_dodge(0.5)) +
  stat_summary(fun = mean, geom = "point", size = 3, position = position_dodge(0.5)) +
  stat_summary(fun = mean, geom = "line", size = 0.75, position = position_dodge(0.5)) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels=label.names, values=paired.colors)+
  xlab("Temperature")+
  labs(y=expression(paste(italic("Chironomidae"), " spp. Scaled Biomass")))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 11), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

c_plot <- wrap_plots(c.abund, c.biomass) +
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides = 'collect') &
  theme(
    legend.text = element_text(size = 11, lineheight = 0.9),
    legend.spacing.y = unit(0.5, "cm"),
    legend.key.height = unit(1.25, "cm")
  )
# Read the PNG image
chiro <- readPNG("chironomidae.png")
chiro_grob <- rasterGrob(chiro, x = 0.825, y = 0.92, width = 0.135, height = 0.135)  # Adjust x, y, width, and height

# Overlay the image on the plot
c.fig <- ggdraw() +
  draw_plot(c_plot) +
  draw_grob(chiro_grob)
ggsave("cfig.png", plot = c.fig, device = "png", width = 8.44, height = 6, units = c("in"), dpi = "retina")


chir_pc_combo <- bind_rows(
  baet_pc %>% mutate(source = "Temperature"),
  baet_pc_sp %>% mutate(source = "Temperature & spike"),
)

c.pc <- ggplot(data = chir_pc_combo, aes(x = as.factor(temperature), y = percapita, group = source, color = source))+
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, show.legend = F) +
  stat_summary(fun = mean, geom = "point", size = 3, show.legend = F)+
  stat_summary(fun = mean, geom = "line", size = 0.75, show.legend = F) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels= c("Temperature", 
                                           "Temperature +\nSpike"), values= c("#56B4E9", "#365C8DFF"))+
  xlab("Temperature")+
  labs(y="Adult Biomass (mg)", title = expression(paste(italic("Chironomidae"), " spp.")))+
  theme(text = element_text(size = 11), axis.text.x = element_text(hjust = 1, angle=45, size = 10), 
        axis.text.y = element_text(size = 11), legend.key = element_rect(fill = "transparent"))

# Fit linear models with ALL interactions (full factorial design)
chirmod_biomass <- lm(biomass ~ (temperature + temp_spike + hydropeaking + HFE)^4, data = chir_biomass_combo)
chirmod_abund <- lm(abundance ~ (temperature + temp_spike + hydropeaking + HFE)^4, data = chir_abund_combo)

sum_biomass_lm <- lm.beta(chirmod_biomass)
sum_abund_lm <- lm.beta(chirmod_abund)
# Extract the standardized coefficients
# lm.beta stores the standardized coefficients in the model object as 'standardized.coefficients'
std_coefs.b <- sum_biomass_lm$standardized.coefficients
std_coefs.a <- sum_abund_lm$standardized.coefficients
# Remove the intercept (if present) so we only have the effects and interactions
std_coefs.b <- std_coefs.b[-1]
std_coefs.a <- std_coefs.a[-1]

# Create a data frame with effect names and their beta values
coef_df.b <- data.frame(
  Stressor = names(std_coefs.b),
  Beta = std_coefs.b,
  AbsBeta = abs(std_coefs.b)  # for sorting purposes
)


coef_df.a <- data.frame(
  Stressor = names(std_coefs.a),
  Beta = std_coefs.a,
  AbsBeta = abs(std_coefs.a)  # for sorting purposes
)
# Sort the data frame by absolute beta values in decreasing order
coef_df_sorted.chir.b <- coef_df.b %>% arrange(desc(AbsBeta))
coef_df_sorted.chir.a <- coef_df.a %>% arrange(desc(AbsBeta))

# Define effect size categories
coef_df_sorted.chir.b$Effect_Size <- cut(abs(coef_df_sorted.chir.b$Beta),
                                         breaks = c(0, 0.1, 0.3, 0.5, Inf),
                                         labels = c("Small", "Medium", "Large", "Very Large"),
                                         right = FALSE)

# Define effect size categories
coef_df_sorted.chir.a$Effect_Size <- cut(abs(coef_df_sorted.chir.a$Beta),
                                         breaks = c(0, 0.1, 0.3, 0.5, Inf),
                                         labels = c("Small", "Medium", "Large", "Very Large"),
                                         right = FALSE)

# # Now plot the sorted coefficients using ggplot2
# ggplot(coef_df_sorted.chir.b, aes(x = reorder(Effect, Beta), y = Beta, fill = Effect_Size)) +
#   geom_bar(stat = "identity") +
#   coord_flip() +  # flips the axes so the labels are easier to read
#   labs(
#     title = "Standardized Beta Coefficients for Biomass",
#     x = "Effect",
#     y = "Standardized Beta Coefficient"
#   ) +
#   scale_fill_manual(values = c("Small" = "blue", "Medium" = "orange", "Large" = "red", "Very Large" = "darkred")) +
#   theme_bw()
# 
# # Now plot the sorted coefficients using ggplot2
# ggplot(coef_df_sorted.chir.a, aes(x = reorder(Effect, Beta), y = Beta, fill = Effect_Size)) +
#   geom_bar(stat = "identity") +
#   coord_flip() +  # flips the axes so the labels are easier to read
#   labs(
#     title = "Standardized Beta Coefficients for Abundance",
#     x = "Effect",
#     y = "Standardized Beta Coefficient"
#   ) +
#   scale_fill_manual(values = c("Small" = "blue", "Medium" = "orange", "Large" = "red", "Very Large" = "darkred")) +
#   theme_bw()
chirmod_biomass <- aov(biomass ~ (temperature + temp_spike + hydropeaking + HFE)^4, data = chir_biomass_combo)
chirmod_abund <- aov(abundance ~ (temperature + temp_spike + hydropeaking + HFE)^4, data = chir_abund_combo)

# Extract partial eta squared values
eta_sq_df.chir.b <- eta_squared(chirmod_biomass, partial = TRUE) %>%
  as.data.frame() %>%
  arrange(desc(Eta2_partial))  # Order from largest to smallest effect

eta_sq_df.chir.a <- eta_squared(chirmod_abund, partial = TRUE) %>%
  as.data.frame() %>%
  arrange(desc(Eta2_partial))  # Order from largest to smallest effect
# Create the bar plot
ggplot(eta_sq_df.chir.b, aes(x = reorder(Parameter, Eta2_partial), y = Eta2_partial, fill = Eta2_partial)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +  # Flip for better readability
  scale_fill_viridis_c() +  # Color scale
  labs(title = "Hydropsyche spp.",
       x = "Stressor(s)",
       y = "Partial Eta Squared") +
  theme_bw()

ggplot(eta_sq_df.chir.a, aes(x = reorder(Parameter, Eta2_partial), y = Eta2_partial, fill = Eta2_partial)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +  # Flip for better readability
  scale_fill_viridis_c() +  # Color scale
  labs(title = "Hydropsyche spp.",
       x = "Stressor(s)",
       y = "Partial Eta Squared") +
  theme_bw()


# Merge Dataframes
results_df.chir.b <- left_join(eta_sq_df.chir.b, coef_df_sorted.chir.b, by = c("Parameter" = "Stressor")) %>%
  mutate(Direction = ifelse(Beta > 0, "Positive", "Negative"))

# Sort by Partial Eta Squared (Descending)
results_df.chir.b <- results_df.chir.b %>% arrange(desc(Eta2_partial))
write.table(results_df.chir.b, file = "partialetasquare_chir_biomass.csv", sep = " , ")

# Format Table with `flextable`
results_chir.b <- flextable(results_df.chir.b) %>%
  theme_vanilla() %>%  # Clean formatting
  set_header_labels(
    Parameter = "Stressor(s)",
    Eta2_partial = "Partial Eta²",
    Beta = "β Coefficient",
    Direction = "Effect Direction"
  ) %>%
  colformat_num(j = c("Eta2_partial", "Beta")) %>%  
  autofit()  # Adjust column widths

# Merge Dataframes
results_df.chir.a <- left_join(eta_sq_df.chir.a, coef_df_sorted.chir.a, by = c("Parameter" = "Stressor")) %>%
  mutate(Direction = ifelse(Beta > 0, "Positive", "Negative"))

# Sort by Partial Eta Squared (Descending)
results_df.chir.a <- results_df.chir.a %>% arrange(desc(Eta2_partial))
write.table(results_df.chir.a, file = "partialetasquare_chir_abund.csv", sep = " , ")

# Format Table with `flextable`
results_chir.a <- flextable(results_df.chir.a) %>%
  theme_vanilla() %>%  # Clean formatting
  set_header_labels(
    Parameter = "Stressor(s)",
    Eta2_partial = "Partial Eta²",
    Beta = "β Coefficient",
    Direction = "Effect Direction"
  ) %>%
  colformat_num(j = c("Eta2_partial", "Beta")) %>%  
  autofit()  # Adjust column widths

# G. lacustris
gamm_abund <- read_csv("GAMM_temp_abund.csv")
gamm_pc <- read_csv("GAMM_temp_percapita.csv")
gamm_biomass <- read_csv("GAMM_temp_biomass.csv")
gamm_abund_sp <- read_csv("GAMM_temp_abund_spike.csv")
gamm_pc_sp <- read_csv("GAMM_temp_percapita_spike.csv")
gamm_biomass_sp <- read_csv("GAMM_temp_biomass_spike.csv")
gamm_abund_t.h <- read_csv("GAMM_temp_hyd_abund.csv")
gamm_biomass_t.h <- read_csv("GAMM_temp_hyd_biomass.csv")
gamm_abund_t.h_spike <- read_csv("GAMM_temp_hyd_abund_spike.csv")
gamm_biomass_t.h_spike <- read_csv("GAMM_temp_hyd_biomass_spike.csv")
gamm_abund_HFE <- read_csv("GAMM_temp_abund_HFE.csv")
gamm_biomass_HFE <- read_csv("GAMM_temp_biomass_HFE.csv")
gamm_abund_HFE_sp <- read_csv("GAMM_temp_abund_HFE_spike.csv")
gamm_biomass_HFE.sp <- read_csv("GAMM_temp_biomass_HFE_spike.csv")
gamm_abund_HFE_sp <- read_csv("GAMM_temp_abund_HFE_spike.csv")
gamm_abund_hyd_HFE <- read_csv("GAMM_temp_hyd_abund_HFE.csv")
gamm_biomass_hyd_HFE <- read_csv("GAMM_temp_hyd_biomass_HFE.csv")
gamm_abund_hyd_HFE.sp <- read_csv("GAMM_temp_hyd_abund_HFE_spike.csv")
gamm_biomass_hyd_HFE.sp <- read_csv("GAMM_temp_hyd_biomass_HFE_spike.csv")

# abundance
gamm_abund_combo <- bind_rows(
  gamm_abund %>% mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  gamm_abund_sp %>%  mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  gamm_abund_t.h %>% mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  gamm_abund_t.h_spike %>% mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  gamm_abund_HFE %>% mutate(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
  gamm_abund_HFE_sp %>% mutate(source = "Temperature & spike & HFE",  temp_spike = 1, hydropeaking = 0, HFE = 1),
  gamm_abund_hyd_HFE %>% mutate(source = "Temperature & hydropeaking & HFE",  temp_spike = 0, hydropeaking = 1, HFE = 1),
  gamm_abund_hyd_HFE.sp %>% mutate(source = "Temperature & spike & hydropeaking & HFE",  temp_spike = 1, hydropeaking = 1, HFE = 1)
)

# biomass
gamm_biomass_combo <- bind_rows(
  gamm_biomass %>%  mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  gamm_biomass_sp %>%  mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  gamm_biomass_t.h %>% mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  gamm_biomass_t.h_spike %>% mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  gamm_biomass_HFE %>% mutate(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
  gamm_biomass_HFE.sp %>% mutate(source = "Temperature & spike & HFE",  temp_spike = 1, hydropeaking = 0, HFE = 1),
  gamm_biomass_hyd_HFE %>% mutate(source = "Temperature & hydropeaking & HFE",  temp_spike = 0, hydropeaking = 1, HFE = 1),
  gamm_biomass_hyd_HFE.sp %>% mutate(source = "Temperature & spike & hydropeaking & HFE",  temp_spike = 1, hydropeaking = 1, HFE = 1)
)

gamm_pc_combo <- bind_rows(
  gamm_pc %>% mutate(source = "Temperature"),
  gamm_pc_sp %>% mutate(source = "Temperature & spike"),
  #gamm_pc_hfe %>% mutate(source = "Temperature & HFE"),
  #gamm_pc_hfe_sp %>% mutate(source = "Temperature & spike & HFE")
)

gamm_abund_combo$abundance <- scale(gamm_abund_combo$abundance)
gamm_abund_combo$source <- ordered(gamm_abund_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE", 
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"))  

g.abund <- ggplot(data = gamm_abund_combo, aes(y = abundance, x = as.factor(temperature), color = source, group = source)) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, 
               position = position_dodge(0.5), show.legend = F) +
  stat_summary(fun = mean, geom = "point", size = 3, position = position_dodge(0.5), show.legend = F) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels=label.names, values=paired.colors)+
  xlab("Temperature")+
  labs(y=expression(paste(italic("G. lacustris"), " Scaled Abundance")))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# gamm biomass
gamm_biomass_combo$biomass <- scale(gamm_biomass_combo$biomass)
gamm_biomass_combo$source <- ordered(gamm_biomass_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE", 
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"))  
g.biomass <- ggplot(data = gamm_biomass_combo, aes(y = biomass, x = as.factor(temperature), color = source, group = source)) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, 
               position = position_dodge(0.5)) +
  stat_summary(fun = mean, geom = "point", size = 3, position = position_dodge(0.5)) +
  stat_summary(fun = mean, geom = "line", size = 0.75, position = position_dodge(0.5)) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels=label.names, values=paired.colors)+
  xlab("Temperature")+
  labs(y=expression(paste(italic("G. lacustris"), " Scaled Biomass")))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 11), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

g_plot <- wrap_plots(g.abund, g.biomass) +
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides = 'collect') &
  theme(
    legend.text = element_text(size = 11, lineheight = 0.9),
    legend.spacing.y = unit(0.5, "cm"),
    legend.key.height = unit(1.25, "cm")
  )
# Read the PNG image
gamm <- readPNG("Gammarus.png")
gamm_grob <- rasterGrob(gamm, x = 0.825, y = 0.92, width = 0.135, height = 0.135)  # Adjust x, y, width, and height

# Overlay the image on the plot
g.fig <- ggdraw() +
  draw_plot(g_plot) +
  draw_grob(gamm_grob)

ggsave("gfig.png", plot = g.fig, device = "png", width = 8.44, height = 6, units = c("in"), dpi = "retina")


g.pc <- ggplot(data = gamm_pc_combo, aes(x = as.factor(temperature), y = percapita, group = source, color = source))+
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, show.legend = F) +
  stat_summary(fun = mean, geom = "point", size = 3, show.legend = F)+
  stat_summary(fun = mean, geom = "line", size = 0.75, show.legend = F) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels= c("Temperature", 
                                           "Temperature +\nSpike"), values= c("#56B4E9", "#365C8DFF"))+
  xlab("Temperature")+
  labs(y="Proportion Adult", title = expression(paste(italic("G. lacustris"))))+
  theme(text = element_text(size = 11), axis.text.x = element_text(hjust = 1, angle=45, size = 10), 
        axis.text.y = element_text(size = 11), legend.key = element_rect(fill = "transparent"))

# Fit linear models with ALL interactions (full factorial design)
gammmod_biomass <- lm(biomass ~ (temperature + temp_spike + hydropeaking + HFE)^4, data = gamm_biomass_combo)
gammmod_abund <- lm(abundance ~ (temperature + temp_spike + hydropeaking + HFE)^4, data = gamm_abund_combo)

sum_biomass_lm <- lm.beta(gammmod_biomass)
sum_abund_lm <- lm.beta(gammmod_abund)
# Extract the standardized coefficients
# lm.beta stores the standardized coefficients in the model object as 'standardized.coefficients'
std_coefs.b <- sum_biomass_lm$standardized.coefficients
std_coefs.a <- sum_abund_lm$standardized.coefficients
# Remove the intercept (if present) so we only have the effects and interactions
std_coefs.b <- std_coefs.b[-1]
std_coefs.a <- std_coefs.a[-1]

# Create a data frame with effect names and their beta values
coef_df.b <- data.frame(
  Stressor = names(std_coefs.b),
  Beta = std_coefs.b,
  AbsBeta = abs(std_coefs.b)  # for sorting purposes
)


coef_df.a <- data.frame(
  Stressor = names(std_coefs.a),
  Beta = std_coefs.a,
  AbsBeta = abs(std_coefs.a)  # for sorting purposes
)
# Sort the data frame by absolute beta values in decreasing order
coef_df_sorted.gamm.b <- coef_df.b %>% arrange(desc(AbsBeta))
coef_df_sorted.gamm.a <- coef_df.a %>% arrange(desc(AbsBeta))

# Define effect size categories
coef_df_sorted.gamm.b$Effect_Size <- cut(abs(coef_df_sorted.gamm.b$Beta),
                                         breaks = c(0, 0.1, 0.3, 0.5, Inf),
                                         labels = c("Small", "Medium", "Large", "Very Large"),
                                         right = FALSE)

# Define effect size categories
coef_df_sorted.gamm.a$Effect_Size <- cut(abs(coef_df_sorted.gamm.a$Beta),
                                         breaks = c(0, 0.1, 0.3, 0.5, Inf),
                                         labels = c("Small", "Medium", "Large", "Very Large"),
                                         right = FALSE)



# Now plot the sorted coefficients using ggplot2
# ggplot(coef_df_sorted.gamm.b, aes(x = reorder(Effect, Beta), y = Beta, fill = Effect_Size)) +
#   geom_bar(stat = "identity") +
#   coord_flip() +  # flips the axes so the labels are easier to read
#   labs(
#     title = "Standardized Beta Coefficients for Biomass",
#     x = "Effect",
#     y = "Standardized Beta Coefficient"
#   ) +
#   scale_fill_manual(values = c("Small" = "blue", "Medium" = "orange", "Large" = "red", "Very Large" = "darkred")) +
#   theme_bw()
# 
# 
# # Now plot the sorted coefficients using ggplot2
# ggplot(coef_df_sorted.gamm.a, aes(x = reorder(Effect, Beta), y = Beta, fill = Effect_Size)) +
#   geom_bar(stat = "identity") +
#   coord_flip() +  # flips the axes so the labels are easier to read
#   labs(
#     title = "Standardized Beta Coefficients for Abundance",
#     x = "Effect",
#     y = "Standardized Beta Coefficient"
#   ) +
#   scale_fill_manual(values = c("Small" = "blue", "Medium" = "orange", "Large" = "red", "Very Large" = "darkred")) +
#   theme_bw()

gammmod_biomass <- aov(biomass ~ (temperature + temp_spike + hydropeaking + HFE)^4, data = gamm_biomass_combo)
gammmod_abund <- aov(abundance ~ (temperature + temp_spike + hydropeaking + HFE)^4, data = gamm_abund_combo)

# Extract partial eta squared values
eta_sq_df.gamm.b <- eta_squared(gammmod_biomass, partial = TRUE) %>%
  as.data.frame() %>%
  arrange(desc(Eta2_partial))  # Order from largest to smallest effect

eta_sq_df.gamm.a <- eta_squared(gammmod_abund, partial = TRUE) %>%
  as.data.frame() %>%
  arrange(desc(Eta2_partial))  # Order from largest to smallest effect
# Create the bar plot
ggplot(eta_sq_df.gamm.b, aes(x = reorder(Parameter, Eta2_partial), y = Eta2_partial, fill = Eta2_partial)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +  # Flip for better readability
  scale_fill_viridis_c() +  # Color scale
  labs(title = "Hydropsyche spp.",
       x = "Stressor(s)",
       y = "Partial Eta Squared") +
  theme_bw()

ggplot(eta_sq_df.gamm.a, aes(x = reorder(Parameter, Eta2_partial), y = Eta2_partial, fill = Eta2_partial)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +  # Flip for better readability
  scale_fill_viridis_c() +  # Color scale
  labs(title = "Hydropsyche spp.",
       x = "Stressor(s)",
       y = "Partial Eta Squared") +
  theme_bw()


# Merge Dataframes
results_df.gamm.b <- left_join(eta_sq_df.gamm.b, coef_df_sorted.gamm.b, by = c("Parameter" = "Stressor")) %>%
  mutate(Direction = ifelse(Beta > 0, "Positive", "Negative"))

# Sort by Partial Eta Squared (Descending)
results_df.gamm.b <- results_df.gamm.b %>% arrange(desc(Eta2_partial))
write.table(results_df.gamm.b, file = "partialetasquare_gamm_biomass.csv", sep = " , ")

# Format Table with `flextable`
results_gamm.b <- flextable(results_df.gamm.b) %>%
  theme_vanilla() %>%  # Clean formatting
  set_header_labels(
    Parameter = "Stressor(s)",
    Eta2_partial = "Partial Eta²",
    Beta = "β Coefficient",
    Direction = "Effect Direction"
  ) %>%
  colformat_num(j = c("Eta2_partial", "Beta")) %>%  
  autofit()  # Adjust column widths

# Merge Dataframes
results_df.gamm.a <- left_join(eta_sq_df.gamm.a, coef_df_sorted.gamm.a, by = c("Parameter" = "Stressor")) %>%
  mutate(Direction = ifelse(Beta > 0, "Positive", "Negative"))

# Sort by Partial Eta Squared (Descending)
results_df.gamm.a <- results_df.gamm.a %>% arrange(desc(Eta2_partial))
write.table(results_df.gamm.a, file = "partialetasquare_gamm_abund.csv", sep = " , ")

# Format Table with `flextable`
results_gamm.a <- flextable(results_df.gamm.a) %>%
  theme_vanilla() %>%  # Clean formatting
  set_header_labels(
    Parameter = "Stressor(s)",
    Eta2_partial = "Partial Eta²",
    Beta = "β Coefficient",
    Direction = "Effect Direction"
  ) %>%
  colformat_num(j = c("Eta2_partial", "Beta")) %>%  
  autofit()  # Adjust column widths


nzms_abund <- read_csv("NZMS_temp_abund.csv")
nzms_pc <- read_csv("NZMS_temp_percapita.csv")
nzms_biomass <- read_csv("NZMS_temp_biomass.csv")
nzms_abund_sp <- read_csv("NZMS_temp_abund_spike.csv")
nzms_pc_sp <- read_csv("NZMS_temp_percapita_spike.csv")
nzms_biomass_sp <- read_csv("NZMS_temp_biomass_spike.csv")
nzms_abund_t.h <- read_csv("NZMS_temp_hyd_abund.csv")
nzms_biomass_t.h <- read_csv("NZMS_temp_hyd_biomass.csv")
nzms_abund_t.h_spike <- read_csv("NZMS_temp_hyd_abund_spike.csv")
nzms_biomass_t.h_spike <- read_csv("NZMS_temp_hyd_biomass_spike.csv")
nzms_abund_HFE <- read_csv("NZMS_temp_abund_HFE.csv")
nzms_biomass_HFE <- read_csv("NZMS_temp_biomass_HFE.csv")
nzms_abund_HFE_sp<- read_csv("NZMS_temp_abund_HFE_spike.csv")
nzms_biomass_HFE.sp <- read_csv("NZMS_temp_biomass_HFE_spike.csv")
nzms_abund_HFE_sp<- read_csv("NZMS_temp_abund_HFE_spike.csv")
nzms_abund_hyd_HFE <- read_csv("NZMS_temp_hyd_abund_HFE.csv")
nzms_biomass_hyd_HFE <- read_csv("NZMS_temp_hyd_biomass_HFE.csv")
nzms_abund_hyd_HFE.sp <- read_csv("NZMS_temp_hyd_abund_HFE_spike.csv")
nzms_biomass_hyd_HFE.sp <- read_csv("NZMS_temp_hyd_biomass_HFE_spike.csv")

nzms_abund <- read_csv("NZMS_temp_abund_2000.csv")
nzms_pc <- read_csv("NZMS_temp_percapita_2000.csv")
nzms_biomass <- read_csv("NZMS_temp_biomass_2000.csv")
nzms_abund_sp <- read_csv("NZMS_temp_abund_spike_2000.csv")
nzms_pc_sp <- read_csv("NZMS_temp_percapita_spike_2000.csv")
nzms_biomass_sp <- read_csv("NZMS_temp_biomass_spike_2000.csv")
nzms_abund_t.h <- read_csv("NZMS_temp_hyd_abund_2000.csv")
nzms_biomass_t.h <- read_csv("NZMS_temp_hyd_biomass_2000.csv")
nzms_abund_t.h_spike <- read_csv("NZMS_temp_hyd_abund_spike_2000.csv")
nzms_biomass_t.h_spike <- read_csv("NZMS_temp_hyd_biomass_spike_2000.csv")
nzms_abund_HFE <- read_csv("NZMS_temp_abund_HFE_2000.csv")
nzms_biomass_HFE <- read_csv("NZMS_temp_biomass_HFE_2000.csv")
nzms_abund_HFE_sp<- read_csv("NZMS_temp_abund_HFE_spike_2000.csv")
nzms_biomass_HFE.sp <- read_csv("NZMS_temp_biomass_HFE_spike_2000.csv")
nzms_abund_HFE_sp<- read_csv("NZMS_temp_abund_HFE_spike_2000.csv")
nzms_abund_hyd_HFE <- read_csv("NZMS_temp_hyd_abund_HFE_2000.csv")
nzms_biomass_hyd_HFE <- read_csv("NZMS_temp_hyd_biomass_HFE_2000.csv")
nzms_abund_hyd_HFE.sp <- read_csv("NZMS_temp_hyd_abund_HFE_spike_2000.csv")
nzms_biomass_hyd_HFE.sp <- read_csv("NZMS_temp_hyd_biomass_HFE_spike_2000.csv")

#abundance
nzms_abund_combo <- bind_rows(
  nzms_abund %>% mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  nzms_abund_sp %>% mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  nzms_abund_t.h %>% mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  nzms_abund_t.h_spike %>% mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  nzms_abund_HFE %>% mutate(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
  nzms_abund_HFE_sp %>% mutate(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
  nzms_abund_hyd_HFE %>% mutate(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
  nzms_abund_hyd_HFE.sp %>% mutate(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
)

# biomass
nzms_biomass_combo <- bind_rows(
  nzms_biomass %>% mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
  nzms_biomass_sp %>% mutate(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
  nzms_biomass_t.h %>% mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
  nzms_biomass_t.h_spike %>% mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
  nzms_biomass_HFE %>% mutate(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
  nzms_biomass_HFE.sp %>% mutate(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
  nzms_biomass_hyd_HFE %>% mutate(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
  nzms_biomass_hyd_HFE.sp %>% mutate(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
)

nzms_pc_combo <- bind_rows(
  nzms_pc %>% mutate(source = "Temperature"),
  nzms_pc_sp %>% mutate(source = "Temperature & spike"),
  #nzms_pc_hfe %>% mutate(source = "Temperature & HFE"),
  #nzms_pc_hfe_sp %>% mutate(source = "Temperature & spike & HFE")
)

nzms_abund_combo$abundance <- scale(nzms_abund_combo$abundance)
nzms_abund_combo$source <- ordered(nzms_abund_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE", 
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"))  

n.abund <- ggplot(data = nzms_abund_combo, aes(y = abundance, x = as.factor(temperature), color = source, group = source)) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, 
               position = position_dodge(0.5), show.legend = F) +
  stat_summary(fun = mean, geom = "point", size = 3, position = position_dodge(0.5), show.legend = F) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels=label.names, values=paired.colors)+
  xlab("Temperature")+
  labs(y=expression(paste(italic("P. antipodarum"), " Scaled Abundance")))+
  theme(text = element_text(size = 13), axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# nzms biomass
nzms_biomass_combo$biomass <- scale(nzms_biomass_combo$biomass)
nzms_biomass_combo$source <- ordered(nzms_biomass_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE", 
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"))  
n.biomass <- ggplot(data = nzms_biomass_combo, aes(y = biomass, x = as.factor(temperature), color = source, group = source)) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2, 
               position = position_dodge(0.5)) +
  stat_summary(fun = mean, geom = "point", size = 3, position = position_dodge(0.5)) +
  stat_summary(fun = mean, geom = "line", size = 0.75, position = position_dodge(0.5)) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels=label.names, values=paired.colors)+
  xlab("Temperature")+
  labs(y=expression(paste(italic("P. antipodarum"), " Scaled Biomass")))+
  theme(text = element_text(size = 11), axis.text.x = element_text(hjust = 1, angle=45, size = 11), 
        axis.text.y = element_text(size = 11), legend.key = element_rect(fill = "transparent"))

n_plot <- wrap_plots(n.abund, n.biomass) +
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides = 'collect') &
  theme(
    legend.text = element_text(size = 11, lineheight = 0.9),
    legend.spacing.y = unit(0.5, "cm"),
    legend.key.height = unit(1.25, "cm")
  )
# Read the PNG image
nzms <- readPNG("nzms.png")
nzms_grob <- rasterGrob(nzms, x = 0.825, y = 0.92, width = 0.135, height = 0.135)  # Adjust x, y, width, and height

# Overlay the image on the plot
n.fig <- ggdraw() +
  draw_plot(n_plot) +
  draw_grob(nzms_grob)

ggsave("nfig.png", plot = n.fig, device = "png", width = 8.44, height = 6, units = c("in"), dpi = "retina")


n.pc <- ggplot(data = nzms_pc_combo, aes(x = as.factor(temperature), y = percapita, group = source, color = source))+
  #stat_summary(fun.data = mean_sdl, geom = "errorbar", linewidth = 0.75, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 3, show.legend = F)+
  stat_summary(fun = mean, geom = "line", size = 0.75, show.legend = F) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+
  theme_bw()+
  scale_color_manual(name = " ", labels= c("Temperature", 
                                           "Temperature +\nSpike"), values= c("#56B4E9", "#365C8DFF"))+
  xlab("Temperature")+
  labs(y="Proportion Adult", title = expression(paste(italic("P. antipodarum"))))+
  theme(text = element_text(size = 11), axis.text.x = element_text(hjust = 1, angle=45, size = 10), 
        axis.text.y = element_text(size = 11), legend.key = element_rect(color = "transparent"))


# Fit linear models with ALL interactions (full factorial design)
nzmsmod_biomass <- lm(biomass ~ (temperature + temp_spike + hydropeaking + HFE)^4, data = nzms_biomass_combo)
nzmsmod_abund <- lm(abundance ~ (temperature + temp_spike + hydropeaking + HFE)^4, data = nzms_abund_combo)

sum_biomass_lm <- lm.beta(nzmsmod_biomass)
sum_abund_lm <- lm.beta(nzmsmod_abund)
# Extract the standardized coefficients
# lm.beta stores the standardized coefficients in the model object as 'standardized.coefficients'
std_coefs.b <- sum_biomass_lm$standardized.coefficients
std_coefs.a <- sum_abund_lm$standardized.coefficients
# Remove the intercept (if present) so we only have the effects and interactions
std_coefs.b <- std_coefs.b[-1]
std_coefs.a <- std_coefs.a[-1]

# Create a data frame with effect names and their beta values
coef_df.b <- data.frame(
  Stressor = names(std_coefs.b),
  Beta = std_coefs.b,
  AbsBeta = abs(std_coefs.b)  # for sorting purposes
)


coef_df.a <- data.frame(
  Stressor = names(std_coefs.a),
  Beta = std_coefs.a,
  AbsBeta = abs(std_coefs.a)  # for sorting purposes
)
# Sort the data frame by absolute beta values in decreasing order
coef_df_sorted.nzms.b <- coef_df.b %>% arrange(desc(AbsBeta))
coef_df_sorted.nzms.a <- coef_df.a %>% arrange(desc(AbsBeta))

# Define effect size categories
coef_df_sorted.nzms.b$Effect_Size <- cut(abs(coef_df_sorted.nzms.b$Beta),
                                         breaks = c(0, 0.1, 0.3, 0.5, Inf),
                                         labels = c("Small", "Medium", "Large", "Very Large"),
                                         right = FALSE)

# Define effect size categories
coef_df_sorted.nzms.a$Effect_Size <- cut(abs(coef_df_sorted.nzms.a$Beta),
                                         breaks = c(0, 0.1, 0.3, 0.5, Inf),
                                         labels = c("Small", "Medium", "Large", "Very Large"),
                                         right = FALSE)

# # Now plot the sorted coefficients using ggplot2
# ggplot(coef_df_sorted.nzms.b, aes(x = reorder(Effect, Beta), y = Beta, fill = Effect_Size)) +
#   geom_bar(stat = "identity") +
#   coord_flip() +  # flips the axes so the labels are easier to read
#   labs(
#     title = "Standardized Beta Coefficients for Biomass",
#     x = "Effect",
#     y = "Standardized Beta Coefficient"
#   ) +
#   scale_fill_manual(values = c("Small" = "blue", "Medium" = "orange", "Large" = "red", "Very Large" = "darkred")) +
#   theme_bw()
# 
# print(coef_df_sorted.nzms.b)
# # Now plot the sorted coefficients using ggplot2
# ggplot(coef_df_sorted.nzms.a, aes(x = reorder(Effect, Beta), y = Beta, fill = Effect_Size)) +
#   geom_bar(stat = "identity") +
#   coord_flip() +  # flips the axes so the labels are easier to read
#   labs(
#     title = "Standardized Beta Coefficients for Abundance",
#     x = "Effect",
#     y = "Standardized Beta Coefficient"
#   ) +
#   scale_fill_manual(values = c("Small" = "blue", "Medium" = "orange", "Large" = "red", "Very Large" = "darkred")) +
#   theme_bw()

nzmsmod_biomass <- aov(biomass ~ (temperature + temp_spike + hydropeaking + HFE)^4, data = nzms_biomass_combo)
nzmsmod_abund <- aov(abundance ~ (temperature + temp_spike + hydropeaking + HFE)^4, data = nzms_abund_combo)

# Extract partial eta squared values
eta_sq_df.nzms.b <- eta_squared(nzmsmod_biomass, partial = TRUE) %>%
  as.data.frame() %>%
  arrange(desc(Eta2_partial))  # Order from largest to smallest effect

eta_sq_df.nzms.a <- eta_squared(nzmsmod_abund, partial = TRUE) %>%
  as.data.frame() %>%
  arrange(desc(Eta2_partial))  # Order from largest to smallest effect
# Create the bar plot
ggplot(eta_sq_df.nzms.b, aes(x = reorder(Parameter, Eta2_partial), y = Eta2_partial, fill = Eta2_partial)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +  # Flip for better readability
  scale_fill_viridis_c() +  # Color scale
  labs(title = "Hydropsyche spp.",
       x = "Stressor(s)",
       y = "Partial Eta Squared") +
  theme_bw()

ggplot(eta_sq_df.nzms.a, aes(x = reorder(Parameter, Eta2_partial), y = Eta2_partial, fill = Eta2_partial)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +  # Flip for better readability
  scale_fill_viridis_c() +  # Color scale
  labs(title = "Hydropsyche spp.",
       x = "Stressor(s)",
       y = "Partial Eta Squared") +
  theme_bw()


# Merge Dataframes
results_df.nzms.b <- left_join(eta_sq_df.nzms.b, coef_df_sorted.nzms.b, by = c("Parameter" = "Stressor")) %>%
  mutate(Direction = ifelse(Beta > 0, "Positive", "Negative"))
write.table(results_df.nzms.b, file = "partialetasquare_nzms_biomass.csv", sep = " , ")

# Sort by Partial Eta Squared (Descending)
results_df.nzms.b <- results_df.nzms.b %>% arrange(desc(Eta2_partial))

# Format Table with `flextable`
results_nzms.b <- flextable(results_df.nzms.b) %>%
  theme_vanilla() %>%  # Clean formatting
  set_header_labels(
    Parameter = "Stressor(s)",
    Eta2_partial = "Partial Eta²",
    Beta = "β Coefficient",
    Direction = "Effect Direction"
  ) %>%
  colformat_num(j = c("Eta2_partial", "Beta")) %>%  
  autofit()  # Adjust column widths

# Merge Dataframes
results_df.nzms.a <- left_join(eta_sq_df.nzms.a, coef_df_sorted.nzms.a, by = c("Parameter" = "Stressor")) %>%
  mutate(Direction = ifelse(Beta > 0, "Positive", "Negative"))

# Sort by Partial Eta Squared (Descending)
results_df.nzms.a <- results_df.nzms.a %>% arrange(desc(Eta2_partial))
write.table(results_df.nzms.a, file = "partialetasquare_nzms_abund.csv", sep = " , ")

# Format Table with `flextable`
results_nzms.a <- flextable(results_df.nzms.a) %>%
  theme_vanilla() %>%  # Clean formatting
  set_header_labels(
    Parameter = "Stressor(s)",
    Eta2_partial = "Partial Eta²",
    Beta = "β Coefficient",
    Direction = "Effect Direction"
  ) %>%
  colformat_num(j = c("Eta2_partial", "Beta")) %>%  
  autofit()  # Adjust column widths



# List of all data frames
df_list <- list(results_df.hyos.a, results_df.hyos.b, results_df.baet.a, results_df.baet.b, 
                results_df.chir.a, results_df.chir.b, results_df.gamm.a, results_df.gamm.b, 
                results_df.nzms.a, results_df.nzms.b)

# Combine all data frames into one
combined_df <- bind_rows(df_list, .id = "Taxon_Response")  # Add identifier

# Extract top stressor for each taxon-response based on highest eta squared
top_stressor_df <- combined_df %>%
  group_by(Taxon_Response) %>%
  filter(!grepl(":", Parameter)) %>%  # Exclude interaction terms
  slice_max(Eta2_partial, n = 1) %>%
  rename(Top_Stressor = Parameter, Top_Stressor_Eta2 = Eta2_partial, Top_Stressor_Direction = Direction)

# Extract strongest interaction for each taxon-response
top_interaction_df <- combined_df %>%
  group_by(Taxon_Response) %>%
  filter(grepl(":", Parameter)) %>%  # Only keep interaction terms
  slice_max(Eta2_partial, n = 1) %>%
  rename(Strongest_Interaction = Parameter, Interaction_Eta2 = Eta2_partial, Interaction_Direction = Direction)

# Merge top stressor and strongest interaction data
summary_df <- (left_join(top_stressor_df, top_interaction_df, by = "Taxon_Response"))
summary_df <- summary_df %>% arrange((as.numeric(Taxon_Response)))
# Manually assign Taxon & Response (Ensure the length matches your data)
summary_df$Taxon_Response <- c(
  "Hydropsyche spp. Abundance", "Hydropsyche spp. Biomass",
  "Baetidae spp. Abundance", "Baetidae spp. Biomass",
  "Chironomidae spp. Abundance", "Chironomidae spp. Biomass",
  "G. lacustris Abundance", "G. lacustris Biomass",
  "P. antipodarum Abundance", "P. antipodarum Biomass"
)
library(flextable)

# Select and rename relevant columns for clarity
summary_table <- summary_df %>%
  select(
    Taxon_Response,
    Top_Stressor, Top_Stressor_Eta2, Top_Stressor_Direction,
    Strongest_Interaction, Interaction_Eta2, Interaction_Direction
  ) %>%
  rename(
    `Taxon & Response` = Taxon_Response,
    `Strongest Predictor` = Top_Stressor,
    `Partial Eta² (Predictor)` = Top_Stressor_Eta2,
    `Effect Direction (Predictor)` = Top_Stressor_Direction,
    `Key Interaction` = Strongest_Interaction,
    `Partial Eta² (Interaction)` = Interaction_Eta2,
    `Effect Direction (Interaction)` = Interaction_Direction
  )

# Create a formatted flextable
results_flextable <- flextable(summary_table) %>%
  theme_vanilla() %>%
  colformat_num(j = c("Partial Eta² (Predictor)", "Partial Eta² (Interaction)")) %>%
  autofit()

# Print the table
results_flextable
save_as_docx(results_flextable, path = "eta_sq_results_table.docx")

# read in Lees Ferry temp and discharge data from 2007 to 2023
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
temps <- average.yearly.temp(tempdata = temp, temp.column_name = "X_00010_00003", date.column_name = "Date")
# create a timeseries of average temperatures 100 years long
temps <- rep.avg.year(temps, n = 100, change.in.temp = 0, years.at.temp = 0)

# Define the temperature multipliers and labels
scenarios <- data.frame(
  multiplier = c(1, 1.1, 1.2, 1.5),
  label = c("base", "+10%", "+20%", "+50%")
)


# Apply modifications to create the different temperature regimes
temperature.regime <- map2_dfr(scenarios$multiplier, scenarios$label, ~ {
  temps[1:25, ] %>%
    mutate(Temperature = Temperature * .x, se = .y)
})

regime <- ggplot(data = temperature.regime, aes(x = as.Date(dts), y = Temperature, color = se))+
  geom_line(size = 1)+
  theme_bw()+
  scale_color_manual(name = "Warming Regime", labels=c("+10%", "+20%", "+50%", "Baseline"), values=c("#abd9e9", "#fdae61","#d7191c","#2c7bb6"))+
  scale_x_date(date_labels="%b", date_breaks  ="4 month")+
  xlab("")+
  theme(text = element_text(size = 11), axis.text.x = element_text(hjust = 1, size = 10, angle = 0), 
        axis.text.y = element_text(size = 10), legend.key = element_rect(fill = "transparent"))


temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
# calculate average yearly temperatures
temps <- average.yearly.temp(tempdata = temp, temp.column_name = "X_00010_00003", date.column_name = "Date")
# create summertime spike (up to 21 C, then scale from there)
temps$Temperature[16:23] <- c(14, 16, 18, 21, 21, 18, 16, 14)

# create a timeseries of average temperatures 100 years long
temps <- rep.avg.year(temps, n = 100, change.in.temp = 0, years.at.temp = 0)

# Define the temperature multipliers and labels
scenarios <- data.frame(
  multiplier = c(1, 1.1, 1.2, 1.5),
  label = c("base", "+10%", "+20%", "+50%")
)

# Apply modifications to create the different temperature regimes
temperature.regime2 <- map2_dfr(scenarios$multiplier, scenarios$label, ~ {
  temps[1:25, ] %>%
    mutate(Temperature = Temperature * .x, se = .y)
})


#add temperature spike in September
spikeregime <- ggplot(data = temperature.regime2, aes(x = as.Date(dts), y = Temperature, color = se))+
  geom_line(size = 1)+
  theme_bw()+
  scale_color_manual(name = "Warming Regime", labels=c("+10%", "+20%", "+50%", "Baseline"), values=c("#abd9e9", "#fdae61","#d7191c","#2c7bb6"))+
  scale_x_date(date_labels="%b", date_breaks  ="4 month")+
  xlab("")+
  theme(text = element_text(size = 11), axis.text.x = element_text(hjust = 1, size = 10, angle = 0), 
        axis.text.y = element_text(size = 10), legend.key = element_rect(fill = "transparent"))

temp_combo <- bind_rows(
  temperature.regime %>% mutate(source = "Temperature"),
  temperature.regime2 %>% mutate(source = "Temperature & spike"),
)



#add temperature spike in September
regimes <- ggplot(data = temp_combo, aes(x = as.Date(dts), y = Temperature, color = se, linetype = source))+
  geom_line(size = 1)+
  theme_bw()+
  scale_color_manual(name = "Warming Regime", labels=c("+10%", "+20%", "+50%", "Baseline"), values=c("#abd9e9", "#fdae61","#d7191c","#2c7bb6"))+
  scale_x_date(date_labels="%b", date_breaks  ="4 month")+
  xlab("Date")+
  labs(title = "Temperature")+
  guides(linetype=guide_legend(title="Scenario"))+
  theme(text = element_text(size = 11), axis.text.x = element_text(hjust = 1, size = 10, angle = 0), 
        axis.text.y = element_text(size = 11), legend.key = element_rect(fill = "transparent"))


Fig2.2 <- wrap_plots(b.pc, h.pc, c.pc, g.pc, n.pc, regimes)+
   plot_layout(guides = 'collect', nrow = 3)+
  plot_annotation(tag_levels = 'a')

ggsave("fig2.2.png", Fig2.2, device = "png", height = 11, width = 8.5, unit = c("in"), dpi = "retina")

