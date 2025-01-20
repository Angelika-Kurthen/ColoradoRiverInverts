###########
## Figure generation for 3rd manuscript
##########
library(ggpubr)
library(patchwork)
#install.packages("svglite")
library(svglite)

#source("Multispp.R")
source("MultisppHydropeak.R")

# make sure in the correct format
hydropeak.abunds$`rep(hydropeak[hyd], times = 5)` <- as.numeric(hydropeak.abunds$`rep(hydropeak[hyd], times = 5)`)
hydropeak.abunds$avg.abund<- as.numeric(hydropeak.abunds$avg.abund)
hydropeak.abunds$taxa <- as.factor(hydropeak.abunds$taxa)

#plot
hyd_abund_multi <- ggplot(data = hydropeak.abunds, aes(`rep(hydropeak[hyd], times = 5)`, avg.abund, group = taxa, color = taxa))+ 
  geom_line(linewidth = 1, alpha = 0.8)+
  # geom_ribbon(data = Hydro_Abund, aes(ymin = scale(means - sd),
  #                 ymax = scale(means + sd), fill = V4),
  #                        alpha = .15,
  #                        show.legend = T) +
  scale_color_manual(name = "Taxon", labels=c("Hydropsyche spp.", "Baetidae spp.", "P. anitpodarum", "Chironomidae spp.", "G. lacustris"), values=c("#AA3377","#66CCEE","#4477AA","#228833", "#CCBB44"))+
  geom_vline(aes(xintercept = 0.01), linetype = "dotted", 
             linewidth=1)+
  geom_vline(aes(xintercept = 0.17 ), linetype="dashed", 
             linewidth=1)+
  geom_vline(aes(xintercept = 0.55 ), linetype = "dotdash", 
             linewidth=1)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  theme_bw()+
  xlab("Hydropeaking Intensity")+
  ylab("Relative Abundance")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

# make sure in correct format
hydropeak.biomass$`rep(hydropeak[hyd], times = 5)` <- as.numeric(hydropeak.biomass$`rep(hydropeak[hyd]`)
hydropeak.biomass$avg.biomass <- as.numeric(hydropeak.biomass$avg.biomass)
hydropeak.biomass$taxa <- as.factor(hydropeak.biomass$taxa) 
#plot
hyd_ts_multi <- ggplot(data = hydropeak.biomass, aes(`rep(hydropeak[hyd], times = 5)`,avg.biomass , color = taxa)) + 
  # geom_ribbon(aes(ymin = sizemeans - sizesd,
  #                   ymax = sizemeans + sizesd),
  #               colour = V3,
  #               alpha = .15,
  #               show.legend = T) +
  geom_line(linewidth = 1, alpha = 0.8)+
  scale_color_manual(name = "Taxon", labels=c("Hydropsyche spp.", "Baetidae spp.", "P. anitpodarum", "Chironomidae spp.", "G. lacustris"), values=c("#AA3377","#66CCEE","#4477AA","#228833", "#CCBB44"))+
  geom_vline(aes(xintercept = 0.01), linewidth = 1, linetype = "dotted")+
  geom_vline(aes(xintercept = 0.17 ), linewidth = 1, linetype="dashed")+ 
  geom_vline(aes(xintercept = 0.55 ), linewidth = 1, linetype = "dotdash")+
  theme_bw()+
  xlab("Hydropeaking Intensity")+
  ylab("Relative Biomass (mg)")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

#make sure in correct format
# can also make this emergent by only looking at HYOS, CHIR, and BAET
hydropeak.avg.annual$`rep(hydropeak[hyd], times = 5)` <- as.numeric(hydropeak.avg.annual$`rep(hydropeak[hyd], times = 5)`)
hydropeak.avg.annual$taxa <- as.factor(hydropeak.avg.annual$taxa)
hydropeak.avg.annual$`means.s3.biomass$mean.S3.biomass` <- as.numeric(hydropeak.avg.annual$`means.s3.biomass$mean.S3.biomass`)
hyd_yr_multi <- ggplot(data = hydropeak.avg.annual, aes(`rep(hydropeak[hyd], times = 5)`, log(`means.s3.biomass$mean.S3.biomass`), group = taxa, color = taxa)) + 
  # geom_ribbon(aes(ymin = sizemeans - sizesd,
  #                   ymax = sizemeans + sizesd),
  #               colour = V3,
  #               alpha = .15,
  #               show.legend = T) +
  geom_line(linewidth = 1, alpha = 0.8)+
  scale_color_manual(name = "Taxon", labels=c("Hydropsyche spp.", "Baetidae spp.", "P. anitpodarum", "Chironomidae spp.", "G. lacustris"), values=c("#AA3377","#66CCEE","#4477AA","#228833", "#CCBB44"))+
  geom_vline(aes(xintercept = 0.01), linewidth = 1, linetype = "dotted")+
  geom_vline(aes(xintercept = 0.17 ), linewidth = 1, linetype="dashed")+ 
  geom_vline(aes(xintercept = 0.55 ), linewidth = 1, linetype = "dotdash")+
  theme_bw()+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  xlab("Hydropeaking Intensity")+
  ylab("Log Annual S3 Biomass (mg)")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

ggarrange(hyd_abund_multi, hyd_ts_multi , hyd_yr_multi,
          labels = c("a", "b", "c"),
          ncol = 3, nrow = 1, common.legend = T)


