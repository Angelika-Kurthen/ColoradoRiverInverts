###########################################
##Code to produce figures for Manuscript 1
###########################################

library(ggpubr)
library(patchwork)
library(svglite)
library(renv)

renv::init()
renv::snapshot()


strategy_colors <- c(
  Boom = "#228833",
  Fast = "#CCBB44",
  Moderate = "#66CCEE",
  Slow = "#AA3377"
)

strategy_scale <- scale_color_manual(
  name = "Strategy",
  values = strategy_colors
)

#----------------------------------------
# Produce Figure 1
#----------------------------------------
source("LifeHistoriesMatrixModels/Scripts/Annual.R")

yr4temp <- ggplot(data = subset(temp, dts >= "2031-07-21" & dts <= "2033-12-20"), aes(as.Date(dts), Temperature))+
  geom_line(size = 0.8)+
  xlab("Month")+
  ylab("Temperature C")+
  theme_bw()+
  scale_x_date(date_labels="%B", date_breaks  ="3 months")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5, angle = 45), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))



threeyear$Taxa <- factor(threeyear$Taxa, levels = c("B", "C", "A", "D"))
logthreeyearplot <- ggplot(data = threeyear, aes(x = Date, y  =
                                                   (Abundance), color = Taxa))+
  geom_line(size = 0.8, alpha = 0.7)+
  
  xlab("Month")+
  ylab("Abundance")+
  scale_color_manual(name = "Strategy", labels=c("Boom", "Fast", "Moderate", "Slow"), values=c("#228833", "#CCBB44","#66CCEE", "#AA3377"))+
  scale_x_date(date_labels="%B", date_breaks  ="3 month")+
  theme_bw()+
  theme(text = element_text(size = 13.5), axis.text.x = element_text(hjust = 1, angle = 45, size = 12.5),
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))


source("LifeHistoriesMatrixModels/Scripts/DensIndependence.R")
dens.ind$Taxa <- factor(dens.ind$Taxa, levels = c("B", "C", "A", "D"))
DensInd <- ggplot(data = dens.ind, aes(x = Date, y  =log(Abundance), color = Taxa))+
  geom_point(size = 1, alpha = 0.5)+
  geom_line(size = 1)+
  xlab("Month")+
  ylab("Log Abundance")+
  scale_color_manual(name = "Strategy", labels=c("Boom", "Fast", "Moderate", "Slow"), values=c("#228833", "#CCBB44","#66CCEE", "#AA3377"))+
  scale_x_date(date_labels="%B", date_breaks  ="3 month")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle = 45, size = 12.5),
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

source("LifeHistoriesMatrixModels/Scripts/DensityDependence.R")
dens.dep$Taxa <- factor(dens.dep$Taxa, levels = c("B","C", "A", "D"))
DensDep <- ggplot(data = dens.dep, aes(x = Date, y  =(Abundance), color = Taxa))+
  geom_point(size = 1, alpha = 0.5)+
  geom_line(size = 1, alpha = 0.7)+
  xlab("Month")+
  ylab("Abundance")+
  scale_color_manual(name = "Strategy", labels=c("Boom", "Fast", "Moderate", "Slow"), values=c("#228833", "#CCBB44","#66CCEE", "#AA3377"))+
  scale_y_continuous( breaks = c(0, 1e5, 2e5, 3e5), labels = scales::scientific)+
  scale_x_date(date_labels="%B", date_breaks  ="3 month")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle = 45, size = 12.5),
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

Fig1 <- ggarrange(DensDep, logthreeyearplot,
                  labels = c("a", "b"), hjust = 0, vjust = 0.5,
                  ncol = 1, nrow = 2, common.legend = T)
ggsave(filename = "Fig1.png", plot= Fig1, width = 6.5, height= 8.5, device = "png", dpi = "retina" )


#-----------------------------------------
# Produce Figure 2
#-----------------------------------------

source("LifeHistoriesMatrixModels/Scripts/SpA_Multivolt.R")
source("LifeHistoriesMatrixModels/Scripts/SpB_Multivolt.R")
source("LifeHistoriesMatrixModels/Scripts/SpC_Multivolt.R")
source("LifeHistoriesMatrixModels/Scripts/SPD_Multivolt.R")

oneyear <- rbind(B.oneyear, C.oneyear, A.oneyear, D.oneyear)


# Define arrow positions
arrow_data <- data.frame(
  Strategy = c("Boom", "Boom", "Boom", "Boom", "Boom", "Boom", "Boom", "Boom",
               "Fast", "Fast", "Fast", "Fast", "Fast", "Fast",
               "Moderate","Moderate","Moderate", "Slow","Slow", "Slow"),  # Adjust these to match your facet variable levels
  MeanTemp = c(1,1,1, 2, 2,2,2,2,
               1,1,2,2,2,2, 
               1,2,2,1,2,2),
  x_start = as.Date(c("2035-07-17", "2035-08-28", "2035-10-23", "2035-02-27", "2035-05-22", "2035-07-31", "2035-09-25", "2035-12-28",
                      "2035-08-14", "2035-06-19", "2035-01-16", "2035-05-08", "2035-07-17", "2035-09-11",
                      "2035-09-11", "2035-04-24", "2035-10-09", "2035-08-28", "2035-05-08","2035-10-09")),  
  y_end = c(9.5, 9.5, 9.5, 9.8, 9.8, 9.8, 9.8, 9.8,
            8.75, 8.75, 10.05, 10.5, 10.5, 10.5, 
            7, 7, 7, 7, 7, 7),  # Adjust y positions
  x_end =as.Date(c("2035-07-17", "2035-08-28", "2035-10-23", "2035-02-27", "2035-05-22", "2035-07-31", "2035-09-25", "2035-12-28",
                   "2035-08-14", "2035-06-19", "2035-01-16", "2035-05-08", "2035-07-17", "2035-09-11",
                   "2035-09-11", "2035-04-24", "2035-10-09", "2035-08-28", "2035-05-08","2035-10-09")),  
  y_start = c(11.5, 11.5, 11.5, 11.8, 11.8, 11.8, 11.8, 11.8,
              10.75, 10.75, 12.05, 12.05,12.05, 12.05, 
              9, 9, 9, 9, 9, 9)  # Adjust arrow end points
)

Fig2 <- ggplot(data = oneyear, aes(x = Date, y = log(Abund), group = as.factor(MeanTemp), color = as.factor(MeanTemp)))+
  geom_line(size = 1, alpha = 0.8)+
  scale_color_manual(name = "Mean Temperature (C)", labels = c("12", "20"), values = c("#4477AA", "#EE6677"))+
  ylab("Log Adult Abundance") + 
  theme_bw()+
  scale_x_date(date_labels="%B", date_breaks  ="2 month")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle=45, size = 12.5), 
        axis.text.y = element_text(size = 13),legend.position = "bottom", legend.key = element_rect(fill = "transparent"))+
  facet_wrap(~Strategy)+
  geom_segment(data = arrow_data, aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = as.factor(MeanTemp)), 
               arrow = arrow(type = "closed", length = unit(0.1, "inches")), 
               inherit.aes = FALSE, show.legend = F)
ggsave(filename = "Fig2.png", Fig2, height = 5, width = 5, device = "png", dpi = "retina")


#----------------------------------
# Produce Figure 3
#----------------------------------
# code for temperature regime shift 
source("LifeHistoriesMatrixModels/Scripts/A_sp_Temp_Toggle.R")
source("LifeHistoriesMatrixModels/Scripts/Bsp_Temp_Toggle.R")
source("LifeHistoriesMatrixModels/Scripts/CspTempToggle.R")
source("LifeHistoriesMatrixModels/Scripts/D_sp_TempToggle.R")

temp_df <- rbind(a_temp_adjust_df, b_temp_adjust_df ,c_temp_adjust_df ,d_temp_adjust_df)
# Reorder factor levels for correct legend order
temp_df$V3 <- factor(temp_df$V3, levels = c("B", "C", "A", "D"))

# create the plot
abund <- ggplot(data = temp_df, aes(temp_regime, log(temp_means), color = V3))+
  geom_line(size = 1, alpha = 0.8)+
  scale_color_manual(name = "Strategy", labels=c("Boom", "Fast", "Moderate", "Slow"), values=c("#228833", "#CCBB44","#66CCEE", "#AA3377"))+
  theme_bw()+
  geom_vline(xintercept = mean(a_temp_adjust_df$temp_regime), linetype="dotted", 
             size=1)+
  xlab(" ")+
  scale_y_continuous(labels = scales::number_format(accuracy = 1))+
  ylab("Log Abundance")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"), plot.margin = margin(5,5,5,20))


size_df <- rbind(a_size_df, b_size_df, c_size_df, d_size_df)
# Reorder factor levels for correct legend order
size_df$V4 <- factor(size_df$V4, levels = c("B", "C", "A", "D"))

biomass <- ggplot(data = size_df, aes(temp_regime, stage3s_means, color = V4))+
  geom_line(size = 1, alpha = 0.8)+
  scale_color_manual(name = "Strategy", labels=c("Boom", "Fast", "Moderate", "Slow"), values=c("#228833", "#CCBB44","#66CCEE", "#AA3377"))+
  theme_bw()+
  geom_vline(xintercept = mean(a_temp_adjust_df$temp_regime), linetype="dotted", 
             size=1)+
  xlab(" ")+
  ylab("Individual Biomass (mg)")+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"), plot.margin = margin(5,5,5,20))

size_df$totbiomass <- temp_df$temp_means * size_df$size_means

totbiomass <- ggplot(data = size_df, aes(temp_regime, totbiomass/1000, color = V4))+
  geom_line(size = 1, alpha = 0.8)+
  scale_color_manual(name = "Strategy", labels=c("Boom", "Fast", "Moderate", "Slow"), values=c("#228833", "#CCBB44","#66CCEE", "#AA3377"))+
  theme_bw()+
  geom_vline(xintercept = mean(a_temp_adjust_df$temp_regime), linetype="dotted", 
             size=1)+
  xlab("Mean Annual Water Temperature in C")+
  ylab("Standing Biomass (g)")+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), n.breaks = 3)+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"), plot.margin = margin(5,5,5,20))

Fig3<- ggarrange(abund, biomass, totbiomass, 
          labels = c("a", "b", "c"),
          ncol = 1, nrow = 3, common.legend = T)

ggsave(filename = "Fig3.png", plot = Fig3, device = "png", width = 6.5, height = 8.5, dpi = "retina")

#----------------------------------
## Produce Figure 4
#----------------------------------
temp_dist <- bind_rows(temp_dist_a, temp_dist_b, temp_dist_c, temp_dist_d, .id = "taxa")

# deltatemp <- subset(temp_dist, season == 3)
# temp_dist <- subset(temp_dist, season < 3)

supp.labs <- c("Winter Disturbance", "Summer Disturbance", "\u0394 Abundance")
names(supp.labs) <- c("1", "2", "3")
#temp_dist[which(is.na(temp_dist$V3)), ] <- -Inf

temp_dist$V3[which(is.na(temp_dist$V3))] <- -Inf
temp_dist$taxa <- factor(temp_dist$taxa, levels = c("2", "3", "1", "4"))

d <- ggplot(data = temp_dist, aes(x = temp_regime, y = V3, color = taxa))+
  geom_line(linewidth = 1, alpha = 0.8)+
  scale_color_manual(name = "Strategy", labels=c("Boom", "Fast", "Moderate", "Slow"), values=c("#228833", "#CCBB44","#66CCEE", "#AA3377"))+
  geom_point()+
  theme_bw()+
  xlab("Mean Annual Water Temperature in C")+
  ylab("Log Abundance")+
  scale_y_continuous(labels = scales::number_format(accuracy = 1))+
  facet_grid(.~season, scales = "free_y", labeller = labeller(season = supp.labs))+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"), plot.margin = margin(5,5,5,20))

temp_size <- bind_rows(temp_size_a, temp_size_b, temp_size_c, temp_size_d, .id = "taxa")
temp_dist$taxa <- factor(temp_dist$taxa, levels = c("2", "3", "1", "4"))

supp.labs <- c("Winter Disturbance", "Summer Disturbance", "\u0394 Biomass")
names(supp.labs) <- c("1", "2", "3")

es <- ggplot(data = temp_size, aes(x = temp_regime, y = size_means/1000, color = taxa))+
  geom_line(linewidth = 1, alpha = 0.8)+
  scale_color_manual(name = "Strategy", labels=c("Boom", "Fast", "Moderate", "Slow"), values=c("#228833", "#CCBB44","#66CCEE", "#AA3377"))+
  geom_point()+
  theme_bw()+
  xlab("Mean Annual Water Temperature in C")+
  ylab("Standing Biomass (g)")+
  facet_grid(.~season, scales = "free_y", labeller = labeller(season = supp.labs))+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"), plot.margin = margin(5,5,5,20))

Fig4 <- ggarrange(d, es, nrow  = 2, labels = c("a","b"))

ggsave(filename = "Fig4.png", plot = Fig4, device = "png", width = 7.5, height = 8.5, dpi = "retina")

#--------------------------------
## Produce Figure S1
#--------------------------------
source("LifeHistoriesMatrixModels/Scripts/HilbertMetric.R")

chaostestplot <- ggplot(data = chaostestdf, aes(x = fecs, y = chaos1))+
  geom_point(alpha = 0.8)+
  xlab("Stage 3 Boom Fecundity")+
  ylab("Chaos 0 - 1 Test Index")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5),
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
  theme_bw()

chaosts<- ggplot(data = outchaos2[250:500,], aes(x = as.numeric(timesteps), y = log(mean.abund), group = 1))+
  geom_line(linewidth = 1, col = "#4477AA")+
  xlab("Timestep")+
  ylab("Log Abundance")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5),
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
  scale_x_continuous(breaks = seq(250, 2600, by = 250))+
  theme_bw()

chaosn <- ggplot(data = chaosdf, aes(x = log(V1), y = log(V2)))+
  geom_point(col = "#4477AA", alpha = 0.8)+
  geom_line(linewidth = 1, col = "#4477AA")+
  geom_path(col = "#4477AA")+
  xlab(bquote(log(N[boom](t))))+
  ylab(bquote(log(N[boom](t+1))))+
  annotate("text", x= 12, y= 14.05485, label = "F3 = 5200")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5),
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
  theme_bw()

stablets <- ggplot(data = outstable2[250:500,], aes(x = as.numeric(timesteps), y = log(mean.abund), group = 1))+
  geom_line(linewidth = 1, col = "#EE6677")+
  xlab("Timestep")+
  ylab("Log Abundance")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5),
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
  scale_x_continuous(breaks = seq(250, 2600, by = 250))+
  theme_bw()

stablen <- ggplot(data = stabledf, aes(x = log(V1), y = log(V2)))+
  geom_point(alpha = 0.8, col = "#EE6677")+
  geom_line(linewidth = 1,col = "#EE6677")+
  geom_path(linewidth = 1, col = "#EE6677")+
  xlab(bquote(log(N[boom](t)))  )+
  ylab(bquote(log(N[boom](t+1)))  )+
  annotate("text", x= 11, y= 12.40901, label = "F3 = 1200")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5),
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
  theme_bw()

FigS3 <- ggarrange(chaostestplot,
          ggarrange(chaosn, chaosts, stablen, stablets, ncol = 2, nrow = 2, vjust = 0.5, labels = c("b", "c", "d", "e")),
          labels = "a",
          nrow = 2, common.legend = T)

ggsave(filename = "FigS3.png", FigS3, height = 8.5, width = 6.5, device = "png", dpi = "retina")

source("LifeHistoriesMatrixModels/Scripts/A_sp_Fecundity_Toggle.R")
source("LifeHistoriesMatrixModels/Scripts/B_sp_fecundity_Toggle.R")
source("LifeHistoriesMatrixModels/Scripts/C_sp_FecundityToggle.R")
source("LifeHistoriesMatrixModels/Scripts/D_sp_Fecundity_Toggle.R")

fec_df <- rbind(a_fec_df, b_fec_df, c_fec_df, d_fec_df)
fec_df$V3 <- factor(fec_df$V3, levels = c("B", "C", "A", "D"))
FigS3 <- ggplot(data = fec_df, aes(fec_seq, y= fec_means/10000, color = V3))+
  geom_point(size = 1, alpha = 0.5)+
  stat_smooth(method = "lm",
              position = "identity", 
              formula = y ~ x, se = F) +
  scale_color_manual(name = "Strategy", labels=c("Boom", "Fast", "Moderate", "Slow"), values=c("#228833", "#CCBB44","#66CCEE", "#AA3377"))+
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

ggsave(filename = "FigS3.png", plot = FigS3, device = "png", width = 6, height = 5, dpi = "retina")



# code for degree day sensitivity analysis
source("LifeHistoriesMatrixModels/Scripts/A_sp_DD_toggle.R")
source("LifeHistoriesMatrixModels/Scripts/B_sp_DD_toggle.R")
source("LifeHistoriesMatrixModels/Scripts/CspDDToggle.R")
source("LifeHistoriesMatrixModels/Scripts/D_sp_DD_Toggle.R")

dd_df <- rbind(add_df, bdd_df, cdd_df, ddd_df)
dd_df$V3 <- factor(dd_df$V3, levels = c("B", "C", "A", "D"))
FigS4 <- ggplot(data = dd_df, aes(dd_seq, dd_means/10000, color = V3)) + 
  geom_point(size = 1, alpha = 0.5)+
  stat_smooth(method = "lm", 
              position = "identity",
              formula = y~x, se = F)+
  scale_color_manual(name = "Strategy", labels=c("Boom", "Fast", "Moderate", "Slow"), values=c("#228833", "#CCBB44","#66CCEE", "#AA3377"))+
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

ggsave(filename = "FigS4.png", plot = FigS4, device = "png", width = 6, height = 5, dpi = "retina")

#code to make heatmap for K in response to Disturbance and time post disturbance
source("LifeHistoriesMatrixModels/Scripts/Kwireplot.R")
FigS5 <- ggplot(data = KQT, aes(x = t , y = Q))+
  geom_raster(aes(fill = K), interpolate = T)+
  scale_fill_viridis_c(option = "magma") +
  scale_color_grey()+
  labs(shape = "") +
  theme_bw()+
  xlab("LifeHistoriesMatrixModels/Scripts/Timesteps Post Disturbance")+
  ylab("LifeHistoriesMatrixModels/Scripts/Disturbance Magnitude")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
  guides(fill=guide_legend(title="K (carrying capacity)"))+
  theme(strip.text.x = element_text(size = 14), 
        strip.background = element_rect(
          color="black", fill="white", linetype="solid"))+
  theme(legend.margin = margin(-1,0,0,0, unit="cm"))

ggsave(filename = "FigS5.png", FigS5, height = 5, width = 6, device = "png", dpi = "retina")


# code for temperature regime used in most runs where temp isn't adjusted
source("LifeHistoriesMatrixModels/Scripts/SpA_PulseMagnitude.R")
temp$dts <- as.Date(temp$dts, origin = "1970-01-01")

# code for Temperature-Mortality relationship  
source("LifeHistoriesMatrixModels/Scripts/NegExpSurv.R")
FigS6 <- ggplot(data = tempsurvdf, aes(x = tem, y = temSurv))+
  geom_line(size = 1)+
  xlab("Temperature C")+
  ylab("Survival")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

ggsave(filename = "FigS6.png", plot= FigS6, width = 7, height = 5, device= "png", dpi = "retina")
