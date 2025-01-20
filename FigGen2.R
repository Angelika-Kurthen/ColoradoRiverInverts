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
  # geom_ribbon(data = Hydro_Abund, aes(ymin = scale(means - sd),
  #                 ymax = scale(means + sd), fill = V4),
  #                        alpha = .15,
  #                        show.legend = T) +
  scale_color_manual(name = "Taxon", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
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
hyd_ts <- ggplot(data = Hydro_TS_Biomass, aes(hydropeak, log(sizemeans), color = V4)) + 
  # geom_ribbon(aes(ymin = sizemeans - sizesd,
  #                   ymax = sizemeans + sizesd),
  #               colour = V3,
  #               alpha = .15,
  #               show.legend = T) +
  geom_line(linewidth = 1, alpha = 0.8)+
    scale_color_manual(name = "Taxon", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
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
  # geom_ribbon(aes(ymin = sizemeans - sizesd,
  #                   ymax = sizemeans + sizesd),
  #               colour = V3,
  #               alpha = .15,
  #               show.legend = T) +
  geom_line(linewidth = 1, alpha = 0.8)+
  scale_color_manual(name = "Taxon", labels=c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA"))+
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

