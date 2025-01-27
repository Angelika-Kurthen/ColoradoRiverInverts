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

Hydro_Abund <- as.data.frame(rbind(BAET_hyd_means, HYOS_hyd_means, NZMS_hyd_means, CHIR_hyd_means, GAMM_hyd_means))

# combine all abundance data
Hydro_Abund$hydropeak <- as.numeric(Hydro_Abund$hydropeak)

# make sure in the correct format
Hydro_Abund$means <- as.numeric(Hydro_Abund$means)
Hydro_Abund$sd <- as.numeric(Hydro_Abund$sd)
Hydro_Abund$V4 <- as.factor(Hydro_Abund$V4)
#plot
hyd_abund <- ggplot(data = Hydro_Abund, aes(hydropeak, log(means), group = V4, color = V4))+ 
  geom_line(linewidth = 1, alpha = 0.8)+
  geom_ribbon(data = Hydro_Abund, aes(ymin = log(means - sd),
                  ymax = log(means + sd), fill = V4),
                         alpha = .15,
                         color= "transparent",
                         show.legend = F) +
  scale_color_manual(name = "Taxon", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
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
  ylab("Log Abundance")+
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
hyd_ts <- ggplot(data = Hydro_TS_Biomass, aes(hydropeak, log(sizemeans), group = V4, color = V4)) + 
  geom_ribbon(aes(ymin = log(sizemeans - sizesd),
                    ymax = log(sizemeans + sizesd), fill = V4),
                color = "transparent",
                alpha = .15,
                show.legend = F) +
  geom_line(linewidth = 1, alpha = 0.8)+
    scale_color_manual(name = "Taxon", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  scale_fill_manual(name = "Taxon", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  geom_vline(aes(xintercept = 0.01), linewidth = 1, linetype = "dotted")+
  geom_vline(aes(xintercept = 0.17 ), linewidth = 1, linetype="dashed")+ 
  geom_vline(aes(xintercept = 0.55 ), linewidth = 1, linetype = "dotdash")+
     theme_bw()+
    xlab("Hydropeaking Intensity")+
   ylab("Log Average Timestep Biomass (mg)")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), limit = c(0, 12))+
   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
    axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
                                   
# read in mean annual productivity/biomass data
Hydro_Yr_Prod <- as.data.frame(rbind(BAET_hyd_yrprod, HYOS_hyd_yrprod, CHIR_hyd_yrprod, NZMS_hyd_yrprod, GAMM_hyd_yrprod))                                 
# make sure in correct format

Hydro_Yr_Prod$hydropeak <- as.numeric(Hydro_Yr_Prod$hydropeak)
Hydro_Yr_Prod$S3Yrprod <- as.numeric(Hydro_Yr_Prod$S3Yrprod)
Hydro_Yr_Prod$V3 <- as.factor(Hydro_Yr_Prod$V3)

Hyd_yr <- ggplot(data = Hydro_Yr_Prod, aes(hydropeak, log(S3Yrprod), group = V3, color = V3)) + 
  geom_ribbon(aes(ymin = log(sizemeans - sizesd),
                    ymax = log(sizemeans + sizesd), fill = V3),
                color = "transparent",
                alpha = .15,
                show.legend = F) +
  geom_line(linewidth = 1, alpha = 0.8)+
  scale_color_manual(name = "Taxon", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  scale_fill_manual(name = "Taxon", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  geom_vline(aes(xintercept = 0.01), linewidth = 1, linetype = "dotted")+
  geom_vline(aes(xintercept = 0.17 ), linewidth = 1, linetype="dashed")+ 
  geom_vline(aes(xintercept = 0.55 ), linewidth = 1, linetype = "dotdash")+
  theme_bw()+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  xlab("Hydropeaking Intensity")+
  ylab("Log Annual Emergent Biomass (mg)")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

ggarrange(hyd_abund, hyd_ts, Hyd_yr,
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
temp_abund <- ggplot(data = Temp_Abund, aes(temp_seq, log(means), color = V4, group = V4))+ 
  geom_line(linewidth = 1, alpha = 0.8)+
  geom_ribbon(data = Temp_Abund, aes(ymin = log(means - sd),
                                      ymax = log(means + sd), fill = V4),
              alpha = .15,
              color= "transparent",
              show.legend = F) +
  scale_color_manual(name = "Taxon", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  scale_fill_manual(name = "Taxon", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  geom_vline(aes(xintercept = 1), linetype = "dotted", 
             linewidth=1)+
  geom_vline(aes(xintercept = 2 ), linetype="dashed",
             linewidth=1)+
  # geom_vline(aes(xintercept = 0.55 ), linetype = "dotdash", 
  #            linewidth=1)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), limit = c(0, 12))+
  theme_bw()+
  xlab("Temperature Increase")+
  ylab("Log Abundance")+
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
te_ts <- ggplot(data = Temp_TS_Biomass, aes(temp_seq, log(sizemeans), group = V4, color = V4)) + 
  geom_ribbon(aes(ymin = log(sizemeans - sizesd),
                  ymax = log(sizemeans + sizesd), fill = V4),
              color = "transparent",
              alpha = .15,
              show.legend = F) +
  geom_line(linewidth = 1, alpha = 0.8)+
  scale_color_manual(name = "Taxon", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  scale_fill_manual(name = "Taxon", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  geom_vline(aes(xintercept = 1), linewidth = 1, linetype = "dotted")+
  geom_vline(aes(xintercept = 2 ), linewidth = 1, linetype="dashed")+ 
  #geom_vline(aes(xintercept = 0.55 ), linewidth = 1, linetype = "dotdash")+
  theme_bw()+
  xlab("Temperature Increase")+
  ylab("Log Average Timestep Biomass (mg)")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), limit = c(0, 12))+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

# read in mean annual productivity/biomass data
Temp_Yr_Prod <- as.data.frame(rbind(BAET_te_yrprod, HYOS_te_yrprod, CHIR_te_yrprod, NZMS_temp_yrprod, GAMM_te_yrprod))                                 
# make sure in correct format

Temp_Yr_Prod$temp_seq <- as.numeric(Temp_Yr_Prod$temp_seq)
Temp_Yr_Prod$S3Yrprod <- as.numeric(Temp_Yr_Prod$S3Yrprod)
Temp_Yr_Prod$S3Yrprod_sd <- as.numeric(Temp_Yr_Prod$S3Yrprod_sd)
Temp_Yr_Prod$V4 <- as.factor(Temp_Yr_Prod$V4)

temp_yr <- ggplot(data = Temp_Yr_Prod, aes(temp_seq, log(S3Yrprod), group = V4, color = V4)) + 
  geom_ribbon(aes(ymin = log(S3Yrprod - S3Yrprod_sd),
                  ymax = log(S3Yrprod + S3Yrprod_sd), fill = V4),
              color = "transparent",
              alpha = .15,
              show.legend = F) +
  geom_line(linewidth = 1, alpha = 0.8)+
  scale_color_manual(name = "Taxon", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  scale_fill_manual(name = "Taxon", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
  geom_vline(aes(xintercept = 1), linewidth = 1, linetype = "dotted")+
  geom_vline(aes(xintercept = 2), linewidth = 1, linetype="dashed")+
  #geom_vline(aes(xintercept = 0.55 ), linewidth = 1, linetype = "dotdash")+
  theme_bw()+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), limit = c(0,12))+
  xlab("Temperature Increase")+
  ylab("Log Annual Emergent Biomass (mg)")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

Fig3 <- ggarrange(temp_abund, te_ts, temp_yr,
          labels = c("a", "b", "c"),
          ncol = 3, nrow = 1, common.legend = T)
ggsave(filename = "Fig3.3.png", plot = Fig3, device = "png", width = 12, height = 5, units = "in", dpi = "retina")

