###########################
##Code to produce figures for Manuscript 1
##########################



# code for temperature regime shift 
source("A_sp_Temp_Toggle.R")
source("Bsp_Temp_Toggle.R")
source("CspTempToggle.R")
source("D_sp_TempToggle.R")
temp_df <- rbind(a_temp_adjust_df, b_temp_adjust_df ,c_temp_adjust_df ,d_temp_adjust_df)

ggplot(data = temp_df, aes(temp_regime, temp_means/10000, color = V3))+
  geom_line(size = 1, alpha = 0.8)+
  scale_color_manual(name = "Strategy", labels=c("Stonefly", "Mayfly", "Caddisfly", "Beetle"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377"))+
  theme_bw()+
  geom_vline(xintercept = mean(a_temp_adjust_df$temp_regime), linetype="dotted", 
             size=1)+
  xlab("Mean Annual Water Temperature in C")+
  ylab("Relativized Abundance")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

# code for fecundity sensitivity analysis
source("A_sp_Fecundity_Toggle.R")
source("B_sp_fecundity_Toggle.R")
source("C_sp_FecundityToggle.R")
source("D_sp_Fecundity_Toggle.R")

fec_df <- rbind(a_fec_df, b_fec_df, c_fec_df, d_fec_df)
ggplot(data = fec_df, aes(fec_seq, y= fec_means/10000, color = V3))+
  geom_point(size = 1, alpha = 0.5)+
  stat_smooth(method = "lm",
              position = "identity", 
              formula = y ~ x, se = F) +
  scale_color_manual(name = "Strategy", labels=c("Stonefly", "Mayfly", "Caddisfly", "Beetle"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377"))+
  geom_vline(aes(xintercept = mean(a_fec_df$fec_seq), color = "A"), linetype = "dotdash", 
           size=1)+
  geom_vline(aes(xintercept = mean(b_fec_df$fec_seq),color = "B" ), linetype="dotted", 
             size=1)+
  geom_vline(aes(xintercept = mean(c_fec_df$fec_seq),color = "C" ), linetype = "dotted", 
             size=1)+
  geom_vline(aes(xintercept = mean(d_fec_df$fec_seq), color = "D"), linetype="dotted", 
             size=1)+
  theme_bw()+
  xlab("Fecundity (# of eggs)")+
  ylab("Relativized Abundance")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

# code for degree day sensitivity analysis
source("A_sp_DD_toggle.R")
source("B_sp_DD_toggle.R")
source("CspDDToggle.R")
source("D_sp_DD_Toggle.R")

dd_df <- rbind(add_df, bdd_df, cdd_df, ddd_df)
ggplot(data = dd_df, aes(dd_seq, dd_means/10000, color = V3)) + 
  geom_point(size = 1, alpha = 0.5)+
  stat_smooth(method = "lm", 
              position = "identity",
              formula = y~x, se = F)+
  scale_color_manual(name = "Strategy", labels=c("Stonefly", "Mayfly", "Caddisfly", "Beetle"), values=c("#66CCEE", "#228833", "#EE6677", "#AA3377"))+
  geom_vline(aes(xintercept = mean(add_df$dd_seq), color = "A"), linetype = "dotdash", 
             size=1)+
  geom_vline(aes(xintercept = mean(bdd_df$dd_seq),color = "B" ), linetype="dotted", 
             size=1)+
  geom_vline(aes(xintercept = mean(cdd_df$dd_seq),color = "C" ), linetype = "dotted", 
             size=1)+
  geom_vline(aes(xintercept = mean(ddd_df$dd_seq), color = "D"), linetype="dotted", 
             size=1)+
  theme_bw()+
  xlab("Degree Days to Emergence")+
  ylab("Relativized Abundance")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

# code for julian date based timing
source("SpA_JulianPulse.R")
source("SpB_JulianPulse.R")
source("SpC_JulianPulse.R")
source("SpD_JulianPulse.R")

julianshort <-rbind(ashort_df, bshort_df, cshort_df, dshort_df) 
julianresil <- rbind(aresil_df, bresil_df, cresil_df, dresil_df)
julianlong <- rbind(along_df, blong_df, clong_df, dlong_df)

julianshort$ma <- rollmean(julianshort)

ggplot(data = julianshort, aes(x = all.dates, y = (short_mean), color = V3))+
  geom_point(size = 1, alpha = 0.5)+
  stat_smooth(geom = "smooth", formula =y ~x)+
  scale_color_manual(name = "Taxa", values=c("#66CCEE", "#228833", "#EE6677", "#AA3377"))+
  theme_bw()+
  #stat_poly_eq(parse=T, aes(label = ..eq.label..), formula=y ~ poly(x,3))+
  xlab("Julian Date")+
  ylab("Log Abundance")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5),
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

ggplot(data = julianresil, aes(x = all.dates, y = resil_mean/10000, color = V3))+
  geom_point(size = 1, alpha = 0.5)+
  stat_smooth(geom = "smooth", 
              position = "identity", se = F)+
  scale_color_manual(name = "Strategy", labels=c("Stonefly", "Mayfly", "Caddisfly", "Beetle"), values=c("#66CCEE", "#228833", "#EE6677", "#AA3377"))+
  theme_bw()+
  xlab("Julian Date")+
  ylab("Relativized Abundance")+
  scale_x_date(date_labels="%B", date_breaks  ="1 month")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle=45, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

ggplot(data = julianlong, aes(x = all.dates, y = long_mean/10000, color = V3))+
  geom_point(size = 1, alpha = 0.5)+
  stat_smooth(geom = "smooth",
              position = "identity",
              formula = y~x, se = F)+
  scale_color_manual(name = "Strategy", labels=c("Stonefly", "Mayfly", "Caddisfly", "Beetle"), values=c("#66CCEE", "#228833", "#EE6677", "#AA3377"))+
  theme_bw()+
  xlab("Julian Date")+
  ylab("Relativized Abundance")+
  scale_x_date(date_labels="%B", date_breaks  ="1 month")+
  theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

# code for pulse magnitude figure
source("SpA_PulseMagnitude.R")
source("SpB_PulseMagnitude.R")
source("SpC_PulseMagnitude.R")
source("SpD_PulseMagnitude.R")

# magnitude_df <- rbind(a_magnitude_df, b_magnitude_df, c_magnitude_df, d_magnitude_df)
# magnitude_df$magnitudes <- as.numeric(magnitude_df$magnitudes)
# magnitude_df$immediate_response <- as.numeric(magnitude_df$immediate_response)


max_df <- rbind(a_short_df, b_short_df, c_short_df, d_short_df)
max_df$magnitudes <- as.numeric(max_df$magnitudes)
max_df$short_response <- as.numeric(max_df$short_response)
ggplot(data = max_df, aes(x = magnitudes, y = log(short_response), color = V3))+
  geom_point(size = 1, alpha = 0.5)+
  stat_smooth(geom = "smooth",
position = "identity",
formula = y~x, se = F)+
  scale_color_manual(name = "Strategy", labels=c("Stonefly", "Mayfly", "Caddisfly", "Beetle"), values=c("#66CCEE", "#228833", "#EE6677", "#AA3377"))+
  xlab("Pulse Disturbance Magnitude (proportion bankfull discharge)")+
  ylab("Log Post-Pulse Abundance")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))



# code for pulse frequency figure
source("SpA_Frequency.R")
source("SpB_Frequency.R")
source("SpC_Frequency.R")
source("SpD_Frequency.R")


immediate_df <- rbind(a_immediate_df, b_immediate_df, c_immediate_df, d_immediate_df)
immediate_df$immediate <- as.numeric(immediate_df$immediate)
immediate_df$V2 <- as.numeric(immediate_df$V2)
short_df <- rbind(a_short_df, b_short_df, c_short_df, d_short_df)
short_df$short <- as.numeric(short_df$short)
short_df$V2 <- as.numeric(short_df$V2)
long_df <- rbind(a_long_df, b_long_df, c_long_df, d_long_df)
long_df$long <- as.numeric(long_df$long)
long_df$V2 <- as.numeric(long_df$V2)

ggplot(data = immediate_df, aes(x = V2, y = immediate/10000, color = V3))+
  scale_color_manual(name = "Strategy", labels=c("Stonefly", "Mayfly", "Caddisfly", "Beetle"), values=c("#66CCEE", "#228833", "#EE6677", "#AA3377"))+
  #geom_point(size = 1, alpha = 0.5)+
  stat_smooth(se= F)+
  xlab("Annual Frequency of Pulse Disturbance")+
  ylab("Log Abundance")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

ggplot(data = short_df, aes(x = V2, y = short/10000, color = V3))+
  scale_color_manual(name = "Strategy", labels=c("Stonefly", "Mayfly", "Caddisfly", "Beetle"), values=c("#66CCEE", "#228833", "#EE6677", "#AA3377"))+
  geom_point(size = 1, alpha = 0.5)+
  stat_smooth(se= F)+
  xlab("Annual Frequency of Pulse Disturbance")+
  ylab("Log Abundance")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

ggplot(data = long_df, aes(x = V2, y = long/10000, color = V3))+
  scale_color_manual(name = "Strategy", labels=c("Stonefly", "Mayfly", "Caddisfly", "Beetle"), values=c("#66CCEE", "#228833", "#EE6677", "#AA3377"))+
  geom_point(size = 1, alpha = 0.5)+
  stat_smooth(se= F)+
  xlab("Annual Frequency of Pulse Disturbance")+
  ylab("Relativized Abundance")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

# code for temperature regime used in most runs where temp isn't adjusted
source("SpA_Frequency.R")
temp$dts <- as.Date(temp$dts, origin = "1970-01-01")
ggplot(data = temp[1:27,], aes(dts, Temperature))+
  geom_line(size = 1)+
  xlab("Month")+
  ylab("Temperature C")+
  theme_bw()+
  scale_x_date(date_labels="%B", date_breaks  ="2 months")+
  theme(text = element_text(size = 14), axis.text.x = element_text(size = 12.5), 
  axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

# code for Temperature-Mortality relationship  
source("NegExpSurv.R")
ggplot(data = tempsurvdf, aes(x = tem, y = temSurv))+
  geom_line(size = 1)+
  xlab("Temperature C")+
  ylab("Survival")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

# code for flood pulse mortalty curves
df <- as.data.frame(rbind(med.df, low.df))

ggplot(df, aes(x = Q, y = Survival, col = Response))+
  geom_line(linewidth = 1)+
  coord_cartesian(ylim = c(0,1)) +
  ylab('Immediate Post-Disturbance Survival') +
  scale_color_manual(name = " ", labels=c("High Survival", "Low Survival"), values=c("#feb078","#2c115f"))+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
  xlab('Max Event Discharge/Bankfull Discharge')


# annual 
source("Annual.R")
ggplot(data = annual, aes(x = Date, y  = 
                            log(Abundance), color = Taxa))+
  #geom_point(size = 1, alpha = 0.5)+
  geom_line(size = 1, alpha = 0.8)+
  #stat_smooth(size= 1, span = 0.4, se =F)+
  xlab("Month")+
  ylab("Log Abundance")+
  ylim(c(5, 18))+
  #ylim(c(0, 5000000))+
  scale_color_manual(name = "Strategy", labels=c("Stonefly", "Mayfly", "Caddisfly", "Beetle"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377"))+
  scale_x_date(date_labels="%B", date_breaks  ="1 month")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle = 45, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

# disturbance
source("Annual.R")
threeyear
ggplot(data = threeyear, aes(x = Date, y  =
                           log(Abundance), color = Taxa))+
  #geom_point(size = 1, alpha = 0.5)+
  geom_line(size = 0.8, alpha = 0.7)+
 # stat_smooth(size= 1, span = 0.4, se =F)+
  ylim(c(5, 18))+
  #ylim(c(0, 5000000))+
  xlab("Month")+
  ylab("Log Abundance")+
  geom_vline(xintercept = as.numeric(as.Date("2035-05-08")), 
             color = "black", 
             lwd = 1,
             linetype = "dotted") +
  scale_color_manual(name = "Strategy", labels=c("Stonefly", "Mayfly", "Caddisfly", "Beetle"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377"))+
  scale_x_date(date_labels="%B", date_breaks  ="3 month")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle = 45, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

threeyear[which(threeyear$Date == "2034-05-09"),]
threeyear[which(threeyear$Date == "2036-05-07"),]
threeyear[which(threeyear$Date == "2037-05-06"),]

source("Annual.R")

ggplot(data = pulse, aes(x = Date, y  =
                           log(Abundance), color = Taxa))+
  #geom_point(size = 1, alpha = 0.5)+
  geom_line(size = 1, alpha = 0.8)+
  # stat_smooth(size= 1, span = 0.4, se =F)+
  ylim(c(5, 18))+
  #ylim(c(0, 5000000))+
  xlab("Month")+
  ylab("Log Abundance")+
  geom_vline(xintercept = as.numeric(as.Date("2035-05-08")), 
             color = "black", 
             lwd = 1,
             linetype = "dotted") +
  scale_color_manual(name = "Strategy", labels=c("Stonefly", "Mayfly", "Caddisfly", "Beetle"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377"))+
  scale_x_date(date_labels="%B", date_breaks  ="1 month")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle = 45, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))



# press disturbance magnitude
source("A_sp_press_mag.R")
source("B_sp_press_mag.R")
source("C_sp_press_mag.R")
source("D_sp_press_mag.R")

press_mag_df <- rbind(a_magnitude_df, b_magnitude_df, c_magnitude_df, d_magnitude_df)
press_mag_df$magnitudes <- as.numeric(press_mag_df$magnitudes)
press_mag_df$mag_response <- as.numeric(press_mag_df$mag_response)
ggplot(data = press_mag_df, aes(x = magnitudes, y = mag_response/10000, color = V3))+
  geom_point(size = 1, alpha = 0.5)+
  #geom_line()+
  stat_smooth(size = 1, span = 0.4, se = F)+
  scale_color_manual(name = "Strategy", labels=c("Stonefly", "Mayfly", "Caddisfly", "Beetle"), values=c("#66CCEE", "#228833", "#EE6677", "#AA3377"))+
  xlab("Press Magnitude (Hydropeaking Index)")+
  ylab("Relatived Abundance")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))


# Press v Pulse Magnitude

source("SpA_PressVPulse.R")
source("SpB_PressVPulse.R")
 source("SpC_PressVPulse.R")
source("SPD_PressVPulse.R")


press_pulse <- rbind(a_immediate_df, b_immediate_df, c_immediate_df, d_immediate_df)

TaxaNames <- list(
  'A'="Stonefly",
  'B'="Mayfly",
  'C'="Caddisfly",
  'D'="Beetle"
)

Taxa_labeller <- function(variable,value){
  return(TaxaNames[value])
}

ggplot(data = press_pulse, aes(x = Press_mag, y = Pulse_mag))+
  geom_raster(aes(fill = log(abundance)), interpolate = F)+
  scale_fill_viridis_c(option = "magma") +
  geom_point(data = subset(press_pulse, abundance == 0), aes(x= Press_mag, y = Pulse_mag, shape = "Extirpated"))+
  scale_color_grey()+
  labs(shape = "") +
  facet_wrap(~Taxa, nrow =2, labeller = Taxa_labeller)+
  theme_bw()+
  xlab("Press Magnitude")+
  scale_y_discrete(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))+
  ylab("Pulse Magnitude")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
  axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
  guides(fill=guide_legend(title="Log Abundance"))+
  theme(strip.text.x = element_text(size = 14), 
    strip.background = element_rect(
      color="black", fill="white", linetype="solid"))+
  theme(legend.margin = margin(-1,0,0,0, unit="cm"))




press_pulse_short <- rbind(a_short_df, b_short_df, c_short_df, d_short_df)
ggplot(data = press_pulse_short, aes(x = Press_mag, y = Pulse_mag))+
  geom_raster(aes(fill = log(abundance)), interpolate = F)+
  scale_fill_viridis_c(option = "magma") +
  geom_point(data = subset(press_pulse, abundance == 0), aes(x= Press_mag, y = Pulse_mag, shape = "Extirpated"))+
  scale_color_grey()+
  labs(shape = "") +
  facet_wrap(~Taxa, ncol = 2, labeller = Taxa_labeller)+
  theme_bw()+
  xlab("Press Magnitude")+
  scale_y_discrete(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))+
  ylab("Pulse Magnitude")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
  guides(fill=guide_legend(title="Log Abundance"))+
  theme(strip.text.x = element_text(size = 14), 
        strip.background = element_rect(
          color="black", fill="white", linetype="solid"))+
  theme(legend.margin = margin(-1,0,0,0, unit="cm"))

press_pulse_max <- rbind(a_max_df, b_max_df, c_max_df, d_max_df)
ggplot(data = press_pulse_max, aes(x = Press_mag, y = Pulse_mag))+
  geom_raster(aes(fill = log(abundance)), interpolate = F)+
  scale_fill_viridis_c(option = "magma") +
  geom_point(data = subset(press_pulse_max, abundance == 0), aes(x= Press_mag, y = Pulse_mag, shape = "Extirpated"))+
  scale_color_grey()+
  labs(shape = "") +
  facet_wrap(~Taxa, ncol = 2, labeller = Taxa_labeller)+
  theme_bw()+
  xlab("Press Magnitude")+
  scale_y_discrete(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))+
  ylab("Pulse Magnitude")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
  guides(fill=guide_legend(title="Log Abundance"))+
  theme(strip.text.x = element_text(size = 14), 
        strip.background = element_rect(
          color="black", fill="white", linetype="solid"))+
  theme(legend.margin = margin(-1,0,0,0, unit="cm"))

# Pulse Disturbance Frequency vs Magnitude

# either run these (take a long time)
source("SpA_Freq_V_Mag.R")
source("SpB_Freq_V_Mag.R")
source("SpC_Freq_V_Mag.R")
source("SpD_Freq_V_Mag.R")

#or run them on the HPC and read in csv(still take a long time)


#code to make heatmap for K in response to Disturbance and time post disturbance
source("Kwireplot.R")
ggplot(data = KQT, aes(x = t , y = Q))+
  geom_raster(aes(fill = K), interpolate = F)+
  scale_fill_viridis_c(option = "magma") +
  scale_color_grey()+
  labs(shape = "") +
  theme_bw()+
  xlab("Timesteps Post Disturbance")+
  ylab("Disturbance Magnitude")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
      axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
  guides(fill=guide_legend(title="K (carrying capacity)"))+
  theme(strip.text.x = element_text(size = 14), 
        strip.background = element_rect(
          color="black", fill="white", linetype="solid"))+
  theme(legend.margin = margin(-1,0,0,0, unit="cm"))


spA <- read.csv(file = "SpA_FreqVMag_short.csv")
ggplot(data = spA, aes(x = frequency, y = magnitude))+
  geom_raster(aes(fill = log(abundance)), interpolate = F)+
  scale_fill_viridis_c(option = "magma")

spAim <- read.csv(file = "SpA_FreqVMag_immediate.csv")
ggplot(data = spAim, aes(x = frequency, y = magnitude))+
  geom_raster(aes(fill = log(abundance)), interpolate = F)+
  scale_fill_viridis_c(option = "magma")

spB <- read.csv(file = "SpB_FreqVMag_short.csv")
ggplot(data = spB, aes(x = frequency, y = magnitude))+
  geom_raster(aes(fill = log(abundance)), interpolate = F)+
  scale_fill_viridis_c(option = "magma")

spBim <- read.csv(file = "SpB_FreqVMag_immediate.csv")
ggplot(data = spBim, aes(x = frequency, y = magnitude))+
  geom_raster(aes(fill= log(abundance)), interpolate  =F)+
  scale_fill_viridis_c(option = "magma")
