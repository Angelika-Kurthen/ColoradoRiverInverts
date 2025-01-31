##############################
#Figure Generation for Manuscript 2
#############################
library(ggpubr)
library(patchwork)
#install.packages("svglite")
library(svglite)

# source("CHIR_Validation.R")
# soruce("HYOS_Validation.R") etc

# correlation plots
x11()
Fig1 <- ggarrange(NZMSts, BAETts, GAMMts, HYOSts, CHIRts,  
          labels = c("a", "b", "c", "d", "e"),hjust = 0, vjust = 0.5,
          ncol = 2, nrow = 3, common.legend = T)

# hydropeaking intensity 
source("BAET_Hydropeaking.R")
source("HYOS_Hydropeaking.R")
source("NZMS_Hydropeaking.R")
source("GAMM_Hydropeaking.R")
source("CHIR_Hydropeaking.R")

hyos_hydropeaking <- read_csv("hyos_hydropeaking_results.csv")
hyos_hydropeaking$Taxa <- rep("HYOS", times = length(hyos_hydropeaking$Hydropeak))
chir_hydropeaking <- read_csv("chir_hydropeaking_results.csv")
chir_hydropeaking$Taxa <- rep("CHIR", times = length(chir_hydropeaking$Hydropeak))
baet_hydropeaking_results <- read_csv("baet_hydropeaking_results.csv")
baet_hydropeaking_results$Taxa <- rep("BAET", times = length(baet_hydropeaking_results$Hydropeak))
nzms_hydropeaking <- read_csv("nzms_hydropeaking_results.csv")
nzms_hydropeaking$Taxa <- rep("NZMS", times = length(nzms_hydropeaking$Hydropeak))
gamm_hydropeaking <- read_csv("gamm_hydropeaking_results.csv")
gamm_hydropeaking$Taxa <- rep("GAMM", times = length(gamm_hydropeaking$Hydropeak))
#Hydro_Abund <- as.data.frame(rbind(BAET_hyd_means, HYOS_hyd_means, NZMS_hyd_means, CHIR_hyd_means, GAMM_hyd_means))
Hydro_Abund <- as.data.frame(rbind(baet_hydropeaking_results, hyos_hydropeaking, nzms_hydropeaking, chir_hydropeaking, gamm_hydropeaking))
# combine all abundance data
Hydro_Abund$Hydropeak <- as.numeric(Hydro_Abund$Hydropeak)

# make sure in the correct format
Hydro_Abund$MeanAbund <- as.numeric(Hydro_Abund$MeanAbund)
Hydro_Abund$SdAbund <- as.numeric(Hydro_Abund$SdAbund)
Hydro_Abund$Taxa <- as.factor(Hydro_Abund$Taxa)
#plot
hyd_abund <- ggplot(data = Hydro_Abund, aes(Hydropeak, log(MeanAbund), group = Taxa, color = Taxa))+ 
  geom_line(linewidth = 1, alpha = 0.8)+
  geom_ribbon(data = Hydro_Abund, aes(ymin = log(MeanAbund - SdAbund),
                  ymax = log(MeanAbund + SdAbund), fill = Taxa),
                         alpha = .15,
                         color= "transparent",
                         show.legend = F) +
  scale_color_manual(name = "Taxon", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")), expression(paste(italic("Hydropsyche"), " spp.")), expression(italic("P. antipodarum"))), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  scale_fill_manual(name = "Taxon", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  geom_vline(aes(xintercept = 0.01), linetype = "dotted", 
              linewidth=1)+
  geom_vline(aes(xintercept = 0.17 ), linetype="dashed", 
              linewidth=1)+
  geom_vline(aes(xintercept = 0.55 ), linetype = "dotdash", 
            linewidth=1)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), limit = c(0, 11))+
  theme_bw()+
  xlab("Hydropeaking Intensity")+
  ylab("Log Abundance per Reach")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

# read in mean timestep biomass data
Hydro_TS_Biomass <- rbind(BAET_hyd_size, HYOS_hyd_size, NZMS_hyd_size, CHIR_hyd_size, GAMM_hyd_size)
# make sure in correct format
Hydro_TS_Biomass$hydropeak <- as.numeric(Hydro_TS_Biomass$hydropeak)
Hydro_TS_Biomass$sizemeans <- as.numeric(Hydro_TS_Biomass$sizemeans)
Hydro_TS_Biomass$sizesd <- as.numeric(Hydro_TS_Biomass$sizesd)
Hydro_TS_Biomass$V4 <- as.factor(Hydro_TS_Biomass$V4)
#plot
hyd_ts <- ggplot(data = Hydro_Abund, aes(Hydropeak, log(SizeMean), group = Taxa, color = Taxa)) + 
  geom_ribbon(aes(ymin = log(SizeMean - SizeSd),
                    ymax = log(SizeMean + SizeSd), fill = Taxa),
                color = "transparent",
                alpha = .15,
                show.legend = F) +
  geom_line(linewidth = 1, alpha = 0.8)+
  scale_color_manual(name = "Taxon", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")), expression(paste(italic("Hydropsyche"), " spp.")), expression(italic("P. antipodarum"))), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  scale_fill_manual(name = "Taxon", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  geom_vline(aes(xintercept = 0.01), linewidth = 1, linetype = "dotted")+
  geom_vline(aes(xintercept = 0.17 ), linewidth = 1, linetype="dashed")+ 
  geom_vline(aes(xintercept = 0.55 ), linewidth = 1, linetype = "dotdash")+
     theme_bw()+
    xlab("Hydropeaking Intensity")+
   ylab("Log Average Timestep Biomass (mg) per Reach")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), limit = c(0, 12))+
   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
    axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
                                   
# read in mean annual productivity/biomass data
#Hydro_Yr_Prod <- as.data.frame(rbind(BAET_hyd_yrprod, HYOS_hyd_yrprod, CHIR_hyd_yrprod, NZMS_hyd_yrprod, GAMM_hyd_yrprod))                                 
# make sure in correct format

# Hydro_Yr_Prod$hydropeak <- as.numeric(Hydro_Yr_Prod$hydropeak)
# Hydro_Yr_Prod$S3Yrprod <- as.numeric(Hydro_Yr_Prod$S3Yrprod)
# Hydro_Yr_Prod$V3 <- as.factor(Hydro_Yr_Prod$V3)

Hyd_yr <- ggplot(data = Hydro_Abund, aes(Hydropeak, log(S3Yrprod), group = Taxa, color = Taxa)) + 
  # geom_ribbon(aes(ymin = log(sizemeans - sizesd),
  #                   ymax = log(sizemeans + sizesd), fill = V3),
  #               color = "transparent",
  #               alpha = .15,
  #               show.legend = F) +
  geom_line(linewidth = 1, alpha = 0.8)+
  scale_color_manual(name = "Taxon", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")), expression(paste(italic("Hydropsyche"), " spp.")), expression(italic("P. antipodarum"))), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  scale_fill_manual(name = "Taxon", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  geom_vline(aes(xintercept = 0.01), linewidth = 1, linetype = "dotted")+
  geom_vline(aes(xintercept = 0.17 ), linewidth = 1, linetype="dashed")+ 
  geom_vline(aes(xintercept = 0.55 ), linewidth = 1, linetype = "dotdash")+
  theme_bw()+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  xlab("Hydropeaking Intensity")+
  ylab("Log Annual Emergent Biomass (mg) per Reach")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

Fig3.4 <- ggarrange(hyd_abund, hyd_ts, Hyd_yr,
          labels = c("a", "b", "c"),
          ncol = 3, nrow = 1, common.legend = T)


#------------------------------------------------
# temperature ramp

source("NZMS_TempToggle.R")
source("HYOS_TempToggle.R")
source("CHIR_TempToggle.R")
source("BAET_TempToggle.R")
source("GAMM_TempToggle.R")

# combine all abundance data
Temp_Abund <- as.data.frame(rbind(BAET_te_means, HYOS_te_means, NZMS_temp_means, CHIR_te_means, GAMM_te_means))

Temp_Abund$temp_seq <- as.numeric(Temp_Abund$temp_seq)

# make sure in the correct format
Temp_Abund$means <- as.numeric(Temp_Abund$means)
Temp_Abund$sd <- as.numeric(Temp_Abund$sd)
Temp_Abund$V4 <- as.factor(Temp_Abund$V4)
#plot 
temp_abund <- ggplot(data = Temp_Abund, aes(temp_seq, (means), color = V4, group = V4))+ 
  geom_line(linewidth = 1, alpha = 0.8)+
  geom_ribbon(data = Temp_Abund, aes(ymin = (means - sd),
                                      ymax = (means + sd), fill = V4),
              alpha = .15,
              color= "transparent",
              show.legend = F) +
  scale_color_manual(name = "Taxon", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")), expression(paste(italic("Hydropsyche"), " spp.")), expression(italic("P. antipodarum"))), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  scale_fill_manual(name = "Taxon", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  geom_vline(aes(xintercept = 1), linetype = "dotted", 
             linewidth=1)+
  geom_vline(aes(xintercept = 2 ), linetype="dashed",
             linewidth=1)+
  # geom_vline(aes(xintercept = 0.55 ), linetype = "dotdash", 
  #            linewidth=1)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  theme_bw()+
  xlab("Temperature Increase")+
  ylab("Abundance per Reach")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

# read in mean timestep biomass da
Temp_TS_Biomass <- rbind(BAET_te_size, HYOS_te_size, NZMS_temp_size, CHIR_te_size, GAMM_te_size)
# make sure in correct format
Temp_TS_Biomass$temp_seq <- as.numeric(Temp_TS_Biomass$temp_seq)
Temp_TS_Biomass$sizemeans <- as.numeric(Temp_TS_Biomass$sizemeans)
Temp_TS_Biomass$sizesd <- as.numeric(Temp_TS_Biomass$sizesd)
Temp_TS_Biomass$V4 <- as.factor(Temp_TS_Biomass$V4)
#plot
te_ts <- ggplot(data = Temp_TS_Biomass, aes(temp_seq, (sizemeans), group = V4, color = V4)) + 
  geom_ribbon(aes(ymin = (sizemeans - sizesd),
                  ymax = (sizemeans + sizesd), fill = V4),
              color = "transparent",
              alpha = .15,
              show.legend = F) +
  geom_line(linewidth = 1, alpha = 0.8)+
  scale_color_manual(name = "Taxon", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")), expression(paste(italic("Hydropsyche"), " spp.")), expression(italic("P. antipodarum"))), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  scale_fill_manual(name = "Taxon", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  geom_vline(aes(xintercept = 1), linewidth = 1, linetype = "dotted")+
  geom_vline(aes(xintercept = 2 ), linewidth = 1, linetype="dashed")+ 
  #geom_vline(aes(xintercept = 0.55 ), linewidth = 1, linetype = "dotdash")+
  theme_bw()+
  xlab("Temperature Increase")+
  ylab("Average Timestep Biomass (mg) per Reach")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

# read in mean annual productivity/biomass data
Temp_Yr_Prod <- as.data.frame(rbind(BAET_te_yrprod, HYOS_te_yrprod, CHIR_te_yrprod, NZMS_temp_yrprod, GAMM_te_yrprod))                                 
# make sure in correct format

Temp_Yr_Prod$temp_seq <- as.numeric(Temp_Yr_Prod$temp_seq)
Temp_Yr_Prod$S3Yrprod <- as.numeric(Temp_Yr_Prod$S3Yrprod)
Temp_Yr_Prod$S3Yrprod_sd <- as.numeric(Temp_Yr_Prod$S3Yrprod_sd)
Temp_Yr_Prod$V4 <- as.factor(Temp_Yr_Prod$V4)

temp_yr <- ggplot(data = Temp_Yr_Prod, aes(temp_seq, (S3Yrprod), group = V4, color = V4)) + 
  geom_ribbon(aes(ymin = (S3Yrprod - S3Yrprod_sd),
                  ymax = (S3Yrprod + S3Yrprod_sd), fill = V4),
              color = "transparent",
              alpha = .15,
              show.legend = F) +
  geom_line(linewidth = 1, alpha = 0.8)+
  scale_color_manual(name = "Taxon", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")), expression(paste(italic("Hydropsyche"), " spp.")), expression(italic("P. antipodarum"))), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  scale_fill_manual(name = "Taxon", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  geom_vline(aes(xintercept = 1), linewidth = 1, linetype = "dotted")+
  geom_vline(aes(xintercept = 2), linewidth = 1, linetype="dashed")+
  #geom_vline(aes(xintercept = 0.55 ), linewidth = 1, linetype = "dotdash")+
  theme_bw()+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  xlab("Temperature Increase")+
  ylab("Annual Emergent Biomass (mg) per Reach")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

Fig3 <- ggarrange(temp_abund, te_ts, temp_yr,
          labels = c("a", "b", "c"),
          ncol = 3, nrow = 1, common.legend = T)
ggsave(filename = "Fig3.3.png", plot = Fig3, device = "png", width = 12, height = 5.5, units = "in", dpi = "retina")



#####################


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

regime <- ggplot(data = temperature.regime, aes(x = dts, y = Temperature, color = se))+
  geom_line(size = 1)+
  theme_bw()+
  scale_color_manual(name = "Warming Regime", labels=c("+10%", "+20%", "+50%", "Baseline"), values=c("#abd9e9", "#fdae61","#d7191c","#2c7bb6"))+
  scale_x_date(date_labels="%B", date_breaks  ="2 month")+
  xlab("")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5, angle = 45), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

hyos_abund <- read_csv("HYOS_temp_abund.csv")
baet_abund <- read_csv("BAET_temp_abund.csv")
chir_abund <- read_csv("CHIR_temp_abund.csv")
gamm_abund <- read_csv("GAMM_temp_abund.csv")
nzms_abund <- read_csv("NZMS_temp_abund.csv")
# 
# temp_abund <- as.data.frame(rbind(baet_abund, hyos_abund, chir_abund, gamm_abund, nzms_abund))
# 
# temp_ab <- ggplot(data = temp_abund, aes(x = as.factor(temperature), y = (as.numeric(abundance)), fill = taxa))+
#   geom_boxplot()+
#   facet_wrap(.~taxa, scale= "free_y")
# 
h_abund <- ggplot(data = hyos_abund, aes(x = as.factor(temperature), y = (as.numeric(abundance)), fill = taxa))+
  geom_boxplot()+
  scale_fill_manual(name = " ", labels=c("Hydropsyche spp.","Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#AA3377","#66CCEE", "#228833", "#CCBB44", "#4477AA"))+
  theme_bw()+
  xlab("Temperature")+
  ylab("Abundance")
b_abund <- ggplot(data = baet_abund, aes(x = as.factor(temperature), y = (as.numeric(abundance)), fill = taxa))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(name = " ", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  xlab("Temperature")+
  ylab("Abundance")
c_abund <- ggplot(data = chir_abund, aes(x = as.factor(temperature), y = (as.numeric(abundance)), fill = taxa))+
  geom_boxplot()+
  scale_fill_manual(name = " ", labels=c( "Chironomidae spp.","Baetidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#228833","#66CCEE", "#CCBB44", "#AA3377", "#4477AA"))+
  theme_bw()+
  xlab("Temperature")+
  ylab("Abundance")
g_abund <- ggplot(data = gamm_abund, aes(x = as.factor(temperature), y = (as.numeric(abundance)), fill = taxa))+
  geom_boxplot()+
  theme_bw()+
scale_fill_manual(name = " ", labels=c( "G. lacustris", "Chironomidae spp.","Baetidae spp.", "Hydropsyche spp.", "P. antipodarum"), values=c("#CCBB44", "#228833","#66CCEE", "#AA3377", "#4477AA"))+
  xlab("Temperature")+
  ylab("Abundance")
n_abund <-ggplot(data = nzms_abund, aes(x = as.factor(temperature), y = (as.numeric(abundance)), fill = taxa))+
  geom_boxplot()+
  xlab("Temperature")+
  ylab("Abundance")+
  theme_bw()+
scale_fill_manual(name = " ", labels=c("P. antipodarum", "G. lacustris", "Chironomidae spp.","Baetidae spp.", "Hydropsyche spp.", "P. antipodarum"), values=c("#4477AA","#CCBB44", "#228833","#66CCEE", "#AA3377"))

wrap_plots(b_abund, h_abund, c_abund, g_abund, n_abund, regime)+
plot_layout(guides = 'collect')


hyos_biomass <- read_csv("HYOS_temp_biomass.csv")
baet_biomass <- read_csv("BAET_temp_biomass.csv")
chir_biomass <- read_csv("CHIR_temp_biomass.csv")
gamm_biomass <- read_csv("GAMM_temp_biomass.csv")
nzms_biomass <- read_csv("NZMS_temp_biomass.csv")
# 
# temp_biomass <- as.data.frame(rbind(baet_biomass, hyos_biomass, chir_biomass, gamm_biomass, nzms_biomass))
# 
# tempbio <- ggplot(data = temp_biomass, aes(x = as.factor(temperature), y = (as.numeric(biomass)), fill = taxa))+
#   geom_boxplot()+
#   facet_wrap(.~taxa, scale = "free_y")
 
h_biomass <- ggplot(data = hyos_biomass, aes(x = as.factor(temperature), y = (as.numeric(biomass)), color = taxa, fill = taxa))+
  geom_boxplot()+
  scale_fill_manual(name = " ", labels=c("Hydropsyche spp.","Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#AA3377","#66CCEE", "#228833", "#CCBB44", "#4477AA"))+
  scale_color_manual(name = " ", labels=c("Hydropsyche spp.","Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#AA3377","#66CCEE", "#228833", "#CCBB44", "#4477AA"))+
  theme_bw()+
  xlab("Temperature")+
  ylab("Biomass (mg)")
b_biomass <- ggplot(data = baet_biomass, aes(x = as.factor(temperature), y = (as.numeric(biomass)), color = taxa, fill = taxa))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(name = " ", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  scale_color_manual(name = " ", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  xlab("Temperature")+
  ylab("Biomass (mg)")
c_biomass <- ggplot(data = chir_biomass, aes(x = as.factor(temperature), y = (as.numeric(biomass)), color  = taxa, fill = taxa))+
  geom_boxplot()+
  scale_fill_manual(name = " ", labels=c( "Chironomidae spp.","Baetidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#228833","#66CCEE", "#CCBB44", "#AA3377", "#4477AA"))+
  scale_color_manual(name = " ", labels=c( "Chironomidae spp.","Baetidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#228833","#66CCEE", "#CCBB44", "#AA3377", "#4477AA"))+
  theme_bw()+
  xlab("Temperature")+
  ylab("Biomass (mg)")
g_biomass <- ggplot(data = gamm_biomass, aes(x = as.factor(temperature), y = (as.numeric(biomass)), color = taxa, fill = taxa))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(name = " ", labels=c( "G. lacustris", "Chironomidae spp.","Baetidae spp.", "Hydropsyche spp.", "P. antipodarum"), values=c("#CCBB44", "#228833","#66CCEE", "#AA3377", "#4477AA"))+
  scale_color_manual(name = " ", labels=c( "G. lacustris", "Chironomidae spp.","Baetidae spp.", "Hydropsyche spp.", "P. antipodarum"), values=c("#CCBB44", "#228833","#66CCEE", "#AA3377", "#4477AA"))+
  xlab("Temperature")+
  ylab("Biomass (mg)")
n_biomass <-ggplot(data = nzms_biomass, aes(x = as.factor(temperature), y = (as.numeric(biomass)), color = taxa, fill = taxa))+
  geom_boxplot()+
  xlab("Temperature")+
  ylab("Biomass (mg)")+
  theme_bw()+
  scale_color_manual(name = " ", labels=c("P. antipodarum", "G. lacustris", "Chironomidae spp.","Baetidae spp.", "Hydropsyche spp.", "P. antipodarum"), values=c("#4477AA","#CCBB44", "#228833","#66CCEE", "#AA3377"))+
  scale_fill_manual(name = " ", labels=c("P. antipodarum", "G. lacustris", "Chironomidae spp.","Baetidae spp.", "Hydropsyche spp.", "P. antipodarum"), values=c("#4477AA","#CCBB44", "#228833","#66CCEE", "#AA3377"))

wrap_plots(b_biomass, h_biomass, c_biomass, g_biomass, n_biomass, regime)+
  plot_layout(guides = 'collect')
