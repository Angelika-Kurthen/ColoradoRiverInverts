###########################
##Code to produce figures for Manuscript 1
##########################

library(ggpubr)
library(patchwork)
#install.packages("svglite")
library(svglite)
# code for temperature regime shift 
source("A_sp_Temp_Toggle.R")
source("Bsp_Temp_Toggle.R")
source("CspTempToggle.R")
source("D_sp_TempToggle.R")
temp_df <- rbind(a_temp_adjust_df, b_temp_adjust_df ,c_temp_adjust_df ,d_temp_adjust_df)
# Reorder factor levels for correct legend order
temp_df$V3 <- factor(temp_df$V3, levels = c("B", "C", "A", "D"))

# create the plot
abund <- ggplot(data = temp_df, aes(temp_regime, temp_means/10000, color = V3))+
  geom_line(size = 1, alpha = 0.8)+
  scale_color_manual(name = "Strategy", labels=c("Boom", "Fast", "Moderate", "Slow"), values=c("#228833", "#CCBB44","#66CCEE", "#AA3377"))+
  theme_bw()+
  geom_vline(xintercept = mean(a_temp_adjust_df$temp_regime), linetype="dotted", 
             size=1)+
  xlab(" ")+
  scale_y_continuous(labels = scales::number_format(accuracy = 1))+
  ylab("Relativized Abundance")+
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
  ylab("Final Body Mass (mg)")+
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
  ylab("Total Biomass (g)")+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), n.breaks = 3)+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"), plot.margin = margin(5,5,5,20))
x11()

Fig2<- ggarrange(abund, biomass, totbiomass, 
          labels = c("a", "b", "c"),
          ncol = 1, nrow = 3, common.legend = T)

ggsave(filename = "Fig2stage3.png", plot = Fig2, device = "png", width = 6.5, height = 8.5, dpi = "retina")
ggsave(filename = "Fig2stagebiomass.png", plot = Fig2, device = "png", width = 6.5, height = 8.5, dpi = "retina")

ggsave(filename = "Fig2.png", plot = Fig2, device = "png", width = 6.5, height = 8.5, dpi = "retina")


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

# supp.labs <- c("\u0394 Abundance")
# names(supp.labs) <- c("3")
# 
# es <- ggplot(data = deltatemp, aes(x = temp_regime, y = short, color = taxa))+
#   geom_line(linewidth = 1, alpha = 0.8)+
#   scale_color_manual(name = "Strategy", labels=c("Stonefly", "Mayfly", "Caddisfly", "Beetle"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377"))+
#   geom_point()+
#   theme_bw()+
#   xlab("Mean Annual Water Temperature in C")+
#   ylab("Abundance")+
#   facet_grid(.~season, scales = "free_y", labeller = labeller(season = supp.labs))+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"), plot.margin = margin(5,5,5,20))

temp_size <- bind_rows(temp_size_a, temp_size_b, temp_size_c, temp_size_d, .id = "taxa")
temp_dist$taxa <- factor(temp_dist$taxa, levels = c("2", "3", "1", "4"))

# deltasize <- subset(temp_size, season == 3)
# temp_size <- subset(temp_size, season < 3)
supp.labs <- c("Winter Disturbance", "Summer Disturbance", "\u0394 Biomass")
names(supp.labs) <- c("1", "2", "3")

es <- ggplot(data = temp_size, aes(x = temp_regime, y = size_means/1000, color = taxa))+
  geom_line(linewidth = 1, alpha = 0.8)+
  scale_color_manual(name = "Strategy", labels=c("Boom", "Fast", "Moderate", "Slow"), values=c("#228833", "#CCBB44","#66CCEE", "#AA3377"))+
  geom_point()+
  theme_bw()+
  xlab("Mean Annual Water Temperature in C")+
  ylab("Total Biomass (g)")+
  facet_grid(.~season, scales = "free_y", labeller = labeller(season = supp.labs))+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"), plot.margin = margin(5,5,5,20))
x11()
Fig4 <- ggarrange(d, es, nrow  = 2, labels = c("a","b"))
ggsave(filename = "Fig4.png", plot = Fig4, device = "png", width = 7.5, height = 8.5, dpi = "retina")
ggsave(filename = "Fig4stage3.png", plot = Fig4, device = "png", width = 7.5, height = 8.5, dpi = "retina")

# code for fecundity sensitivity analysis
source("A_sp_Fecundity_Toggle.R")
source("B_sp_fecundity_Toggle.R")
source("C_sp_FecundityToggle.R")
source("D_sp_Fecundity_Toggle.R")

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
source("A_sp_DD_toggle.R")
source("B_sp_DD_toggle.R")
source("CspDDToggle.R")
source("D_sp_DD_Toggle.R")

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

# code for julian date based timing
# source("SpA_JulianPulse.R")
# source("SpB_JulianPulse.R")
# source("SpC_JulianPulse.R")
# source("SpD_JulianPulse.R")
# library(zoo)
# julianshort <-rbind(ashort_df, bshort_df, cshort_df, dshort_df) 
# julianresil <- rbind(aresil_df, bresil_df, cresil_df, dresil_df)
# julianlong <- rbind(along_df, blong_df, clong_df, dlong_df)

# ggplot(data = julianshort, aes(x = all.dates, y = (short_mean), color = V3))+
#   geom_point(size = 1, alpha = 0.5)+
#   stat_smooth(geom = "smooth", formula =y ~x)+
#   scale_color_manual(name = "Taxa", values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377"))+
#   theme_bw()+
#   #stat_poly_eq(parse=T, aes(label = ..eq.label..), formula=y ~ poly(x,3))+
#   xlab("Julian Date")+
#   ylab("Log Abundance")+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5),
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

# ggplot(data = julianresil, aes(x = all.dates, y = resil_mean/10000, color = V3))+
#   geom_point(size = 1, alpha = 0.5)+
#   stat_smooth(geom = "smooth", 
#               position = "identity", se = F)+
#   scale_color_manual(name = "Strategy", labels=c("Stonefly", "Mayfly", "Caddisfly", "Beetle"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377"))+
#   theme_bw()+
#   xlab("Julian Date")+
#   ylab("Relativized Abundance")+
#   scale_x_date(date_labels="%B", date_breaks  ="1 month")+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle=45, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

# ggplot(data = julianlong, aes(x = all.dates, y = long_mean/10000, color = V3))+
#   geom_point(size = 1, alpha = 0.5)+
#   stat_smooth(geom = "smooth",
#               position = "identity",
#               formula = y~x, se = F)+
#   scale_color_manual(name = "Strategy", labels=c("Stonefly", "Mayfly", "Caddisfly", "Beetle"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377"))+
#   theme_bw()+
#   xlab("Julian Date")+
#   ylab("Relativized Abundance")+
#   scale_x_date(date_labels="%B", date_breaks  ="1 month")+
#   theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

# code for pulse magnitude figure
# source("SpA_PulseMagnitude.R")
# source("SpB_PulseMagnitude.R")
# source("SpC_PulseMagnitude.R")
# source("SpD_PulseMagnitude.R")

# magnitude_df <- rbind(a_magnitude_df, b_magnitude_df, c_magnitude_df, d_magnitude_df)
# magnitude_df$magnitudes <- as.numeric(magnitude_df$magnitudes)
# magnitude_df$immediate_response <- as.numeric(magnitude_df$immediate_response)


# max_df <- rbind(a_short_df, b_short_df, c_short_df, d_short_df)
# max_df$magnitudes <- as.numeric(max_df$magnitudes)
# max_df$short_response <- as.numeric(max_df$short_response)
# ggplot(data = max_df, aes(x = magnitudes, y = log(short_response), color = V3))+
#   geom_point(size = 1, alpha = 0.5)+
#   stat_smooth(geom = "smooth",
# position = "identity",
# formula = y~x, se = F)+
#   scale_color_manual(name = "Strategy", labels=c("Stonefly", "Mayfly", "Caddisfly", "Beetle"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377"))+
#   xlab("Pulse Disturbance Magnitude (proportion bankfull discharge)")+
#   ylab("Log Post-Pulse Abundance")+
#   theme_bw()+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))



# code for pulse frequency figure
# source("SpA_Frequency.R")
# source("SpB_Frequency.R")
# source("SpC_Frequency.R")
# source("SpD_Frequency.R")
# 
# 
# immediate_df <- rbind(a_immediate_df, b_immediate_df, c_immediate_df, d_immediate_df)
# immediate_df$immediate <- as.numeric(immediate_df$immediate)
# immediate_df$V2 <- as.numeric(immediate_df$V2)
# short_df <- rbind(a_short_df, b_short_df, c_short_df, d_short_df)
# short_df$short <- as.numeric(short_df$short)
# short_df$short_response <- as.numeric(short_df$short_response)
# long_df <- rbind(a_long_df, b_long_df, c_long_df, d_long_df)
# long_df$long <- as.numeric(long_df$long)
# long_df$V2 <- as.numeric(long_df$V2)
# 
# ggplot(data = immediate_df, aes(x = V2, y = immediate/10000, color = V3))+
#   scale_color_manual(name = "Strategy", labels=c("Stonefly", "Mayfly", "Caddisfly", "Beetle"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377"))+
#   #geom_point(size = 1, alpha = 0.5)+
#   stat_smooth(se= F)+
#   xlab("Annual Frequency of Pulse Disturbance")+
#   ylab("Log Abundance")+
#   theme_bw()+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# 
# ggplot(data = short_df, aes(x = all.dates, y = short/10000, color = V3))+
#   scale_color_manual(name = "Strategy", labels=c("Stonefly", "Mayfly", "Caddisfly", "Beetle"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377"))+
#   geom_point(size = 1, alpha = 0.5)+
#   stat_smooth(se= F)+
#   xlab("Annual Frequency of Pulse Disturbance")+
#   ylab("Log Abundance")+
#   theme_bw()+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# 
# ggplot(data = long_df, aes(x = V2, y = long/10000, color = V3))+
#   scale_color_manual(name = "Strategy", labels=c("Stonefly", "Mayfly", "Caddisfly", "Beetle"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377"))+
#   geom_point(size = 1, alpha = 0.5)+
#   stat_smooth(se= F)+
#   xlab("Annual Frequency of Pulse Disturbance")+
#   ylab("Relativized Abundance")+
#   theme_bw()+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

# code for temperature regime used in most runs where temp isn't adjusted
source("SpA_PulseMagnitude.R")
temp$dts <- as.Date(temp$dts, origin = "1970-01-01")
# tempgraph <- ggplot(data = temp[1:27,], aes(dts, Temperature))+
#   geom_line(size = 1)+
#   xlab("Month")+
#   ylab("Temperature C")+
#   theme_bw()+
#   scale_x_date(date_labels="%B", date_breaks  ="2 months")+
#   theme(text = element_text(size = 14), axis.text.x = element_text(size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# 


# code for Temperature-Mortality relationship  
source("NegExpSurv.R")
FigS6 <- ggplot(data = tempsurvdf, aes(x = tem, y = temSurv))+
  geom_line(size = 1)+
  xlab("Temperature C")+
  ylab("Survival")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

ggsave(filename = "FigS6.png", plot= FigS6, width = 5, height = 5, device= "png", dpi = "retina")
# code for flood pulse mortalty curves
# df <- as.data.frame(rbind(med.df, low.df))
# 
# ggplot(df, aes(x = Q, y = Survival, col = Response))+
#   geom_line(linewidth = 1)+
#   coord_cartesian(ylim = c(0,1)) +
#   ylab('Immediate Post-Disturbance Survival') +
#   scale_color_manual(name = " ", labels=c("High Survival", "Low Survival"), values=c("#feb078","#2c115f"))+
#   theme_bw()+
#   theme(text = element_text(size = 14), axis.text.x = element_text(size = 12.5, angle = 45), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
#   xlab('Max Event Discharge/Bankfull Discharge')
# 

# annual 
# source("Annual.R")
# ggplot(data = annual, aes(x = Date, y  = 
#                             log(Abundance), color = Taxa))+
#   #geom_point(size = 1, alpha = 0.5)+
#   geom_line(size = 1, alpha = 0.8)+
#   #stat_smooth(size= 1, span = 0.4, se =F)+
#   xlab("Month")+
#   ylab("Log Abundance")+
#   ylim(c(5, 18))+
#   #ylim(c(0, 5000000))+
#   scale_color_manual(name = "Strategy", labels=c("Stonefly", "Mayfly", "Caddisfly", "Beetle"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377"))+
#   scale_x_date(date_labels="%B", date_breaks  ="1 month")+
#   theme_bw()+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle = 45, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
# 
# disturbance
source("Annual.R")
# 
# threeyearplot <- ggplot(data = threeyear, aes(x = Date, y  =
#                            (Abundance)/100000, color = Taxa))+
#   #geom_point(size = 1, alpha = 0.5)+
#   geom_line(size = 0.8, alpha = 0.7)+
#  # stat_smooth(size= 1, span = 0.4, se =F)+
#   #ylim(c(5, 18))+
#   ylim(c(0, 10))+
#   xlab("Month")+
#   ylab("Relative Abundance")+
#   geom_vline(xintercept = as.numeric(as.Date("2035-05-08")), 
#              color = "black", 
#              lwd = 1,
#              linetype = "dotted") +
#   geom_vline(xintercept = as.numeric(as.Date("2035-05-22")), 
#              color = "black", 
#              lwd = 1,
#              linetype = "dotted") +
#   geom_vline(xintercept = as.numeric(as.Date("2035-06-05")), 
#              color = "black", 
#              lwd = 1,
#              linetype = "dotted") +
#   geom_vline(xintercept = as.numeric(as.Date("2035-06-19")), 
#              color = "black", 
#              lwd = 1,
#              linetype = "dotted") +  
#   scale_color_manual(name = "Strategy", labels=c("Stonefly", "Mayfly", "Caddisfly", "Beetle"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377"))+
#   scale_x_date(date_labels="%B", date_breaks  ="3 month")+
#   theme_bw()+
#   theme(text = element_text(size = 13.54), axis.text.x = element_text(hjust = 1, angle = 45, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
#   #inset_element(yr4temp, 0.45, 0.55, 1, 1)+
#   #plot_annotation(tag_levels = 'a')
# 

yr4temp <- ggplot(data = subset(temp, dts >= "2031-07-21" & dts <= "2033-12-20"), aes(as.Date(dts), Temperature))+
  geom_line(size = 0.8)+
  xlab("Month")+
  ylab("Temperature C")+
  theme_bw()+
  scale_x_date(date_labels="%B", date_breaks  ="3 months")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5, angle = 45), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
#,  plot.background = element_rect(colour = "black", fill="white", size=0.5))


threeyear$Taxa <- factor(threeyear$Taxa, levels = c("B", "C", "A", "D"))
logthreeyearplot <- ggplot(data = threeyear, aes(x = Date, y  =
                                                (Abundance), color = Taxa))+
  #geom_point(size = 1, alpha = 0.5)+
  geom_line(size = 0.8, alpha = 0.7)+
  # stat_smooth(size= 1, span = 0.4, se =F)+
  #ylim(c(5, 18))+
  #ylim(c(0, 10))+
  xlab("Month")+
  ylab("Abundance")+
  scale_color_manual(name = "Strategy", labels=c("Boom", "Fast", "Moderate", "Slow"), values=c("#228833", "#CCBB44","#66CCEE", "#AA3377"))+
  scale_x_date(date_labels="%B", date_breaks  ="3 month")+
  theme_bw()+
  theme(text = element_text(size = 13.5), axis.text.x = element_text(hjust = 1, angle = 45, size = 12.5),
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))


source("DensIndependence.R")
dens.ind$Taxa <- factor(dens.ind$Taxa, levels = c("B", "C", "A", "D"))
DensInd <- ggplot(data = dens.ind, aes(x = Date, y  =log(Abundance), color = Taxa))+
  geom_point(size = 1, alpha = 0.5)+
  geom_line(size = 1)+
  #stat_smooth(size= 1, span = 0.4, se =F)+
  xlab("Month")+
  ylab("Log Abundance")+
  scale_color_manual(name = "Strategy", labels=c("Boom", "Fast", "Moderate", "Slow"), values=c("#228833", "#CCBB44","#66CCEE", "#AA3377"))+
  scale_x_date(date_labels="%B", date_breaks  ="3 month")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle = 45, size = 12.5),
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

source("DensityDependence.R")
dens.dep$Taxa <- factor(dens.dep$Taxa, levels = c("B","C", "A", "D"))
DensDep <- ggplot(data = dens.dep, aes(x = Date, y  =(Abundance), color = Taxa))+
  geom_point(size = 1, alpha = 0.5)+
  geom_line(size = 1)+
  #stat_smooth(size= 1, span = 0.4, se =F)+
  xlab("Month")+
  ylab("Abundance")+
  scale_color_manual(name = "Strategy", labels=c("Boom", "Fast", "Moderate", "Slow"), values=c("#228833", "#CCBB44","#66CCEE", "#AA3377"))+
  scale_x_date(date_labels="%B", date_breaks  ="3 month")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle = 45, size = 12.5),
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

Fig1 <- ggarrange(DensInd, DensDep, logthreeyearplot, yr4temp,
          labels = c("a", "b", "c","d"), hjust = 0, vjust = 0.5,
          ncol = 2, nrow = 2, common.legend = T)

ggsave(filename = "Fig1.png", plot= Fig1, width = 6.5, height= 8.5, device = "png", dpi = "retina" )
 
# threeyear[which(threeyear$Date == "2034-05-09"),]
# threeyear[which(threeyear$Date == "2036-05-07"),]
# threeyear[which(threeyear$Date == "2037-05-06"),]

# source("Annual.R")
# 
# ggplot(data = pulse, aes(x = Date, y  =
#                            log(Abundance), color = Taxa))+
#   #geom_point(size = 1, alpha = 0.5)+
#   geom_line(size = 1, alpha = 0.8)+
#   # stat_smooth(size= 1, span = 0.4, se =F)+
#   ylim(c(5, 18))+
#   #ylim(c(0, 5000000))+
#   xlab("Month")+
#   ylab("Log Abundance")+
#   geom_vline(xintercept = as.numeric(as.Date("2035-05-08")), 
#              color = "black", 
#              lwd = 1,
#              linetype = "dotted") +
#   scale_color_manual(name = "Strategy", labels=c("Stonefly", "Mayfly", "Caddisfly", "Beetle"), values=c("#66CCEE", "#228833", "#CCBB44", "#AA3377"))+
#   scale_x_date(date_labels="%B", date_breaks  ="1 month")+
#   theme_bw()+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle = 45, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))



# press disturbance magnitude
source("A_sp_press_mag.R")
source("B_sp_press_mag.R")
source("C_sp_press_mag.R")
source("D_sp_press_mag.R")

press_mag_df <- rbind(a_magnitude_df, b_magnitude_df, c_magnitude_df, d_magnitude_df)
press_mag_df$magnitudes <- as.numeric(press_mag_df$magnitudes)
press_mag_df$mag_response <- as.numeric(press_mag_df$mag_response)
press_mag_df$V3 <- factor(press_mag_df$V3, levels = c("B", "C", "A", "D"))
FigS2 <- ggplot(data = press_mag_df, aes(x = magnitudes, y = mag_response/10000, color = V3))+
  geom_point(size = 1, alpha = 0.5)+
  geom_line(linewidth = 1, alpha = 0.8)+
  #stat_smooth(size = 1, span = 0.3, se = F)+
  scale_color_manual(name = "Strategy", labels=c("Boom", "Fast", "Moderate", "Slow"), values=c("#228833", "#CCBB44","#66CCEE", "#AA3377"))+
  xlab("Press Magnitude (Hydropeaking Index)")+
  ylab("Relatived Abundance")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

ggsave("FigS2.png", plot = FigS2, width = 6, height = 5, device = "png", dpi = "retina")

# Press v Pulse Magnitude
# 
# source("SpA_PressVPulse.R")
# source("SpB_PressVPulse.R")
#  source("SpC_PressVPulse.R")
# source("SPD_PressVPulse.R")
# 

# press_pulse <- rbind(a_immediate_df, b_immediate_df, c_immediate_df, d_immediate_df)
# 
# TaxaNames <- list(
#   'A'="Stonefly",
#   'B'="Mayfly",
#   'C'="Caddisfly",
#   'D'="Beetle"
# )
# 
# Taxa_labeller <- function(variable,value){
#   return(TaxaNames[value])
# }
# 
# ggplot(data = press_pulse, aes(x = Press_mag, y = Pulse_mag))+
#   geom_raster(aes(fill = log(abundance)), interpolate = F)+
#   scale_fill_viridis_c(option = "magma") +
#   geom_point(data = subset(press_pulse, abundance == 0), aes(x= Press_mag, y = Pulse_mag, shape = "Extirpated"))+
#   scale_color_grey()+
#   labs(shape = "") +
#   facet_wrap(~Taxa, nrow =2, labeller = Taxa_labeller)+
#   theme_bw()+
#   xlab("Press Magnitude")+
#   scale_y_discrete(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))+
#   ylab("Pulse Magnitude")+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#   axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
#   guides(fill=guide_legend(title="Log Abundance"))+
#   theme(strip.text.x = element_text(size = 14), 
#     strip.background = element_rect(
#       color="black", fill="white", linetype="solid"))+
#   theme(legend.margin = margin(-1,0,0,0, unit="cm"))
# 
# 
# 
# 
# press_pulse_short <- rbind(a_short_df, b_short_df, c_short_df, d_short_df)
# ggplot(data = press_pulse_short, aes(x = Press_mag, y = Pulse_mag))+
#   geom_raster(aes(fill = log(abundance)), interpolate = F)+
#   scale_fill_viridis_c(option = "magma") +
#   geom_point(data = subset(press_pulse, abundance == 0), aes(x= Press_mag, y = Pulse_mag, shape = "Extirpated"))+
#   scale_color_grey()+
#   labs(shape = "") +
#   facet_wrap(~Taxa, ncol = 2, labeller = Taxa_labeller)+
#   theme_bw()+
#   xlab("Press Magnitude")+
#   scale_y_discrete(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))+
#   ylab("Pulse Magnitude")+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
#   guides(fill=guide_legend(title="Log Abundance"))+
#   theme(strip.text.x = element_text(size = 14), 
#         strip.background = element_rect(
#           color="black", fill="white", linetype="solid"))+
#   theme(legend.margin = margin(-1,0,0,0, unit="cm"))
# 
# press_pulse_max <- rbind(a_max_df, b_max_df, c_max_df, d_max_df)
# ggplot(data = press_pulse_max, aes(x = Press_mag, y = Pulse_mag))+
#   geom_raster(aes(fill = log(abundance)), interpolate = F)+
#   scale_fill_viridis_c(option = "magma") +
#   geom_point(data = subset(press_pulse_max, abundance == 0), aes(x= Press_mag, y = Pulse_mag, shape = "Extirpated"))+
#   scale_color_grey()+
#   labs(shape = "") +
#   facet_wrap(~Taxa, ncol = 2, labeller = Taxa_labeller)+
#   theme_bw()+
#   xlab("Press Magnitude")+
#   scale_y_discrete(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))+
#   ylab("Pulse Magnitude")+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
#   guides(fill=guide_legend(title="Log Abundance"))+
#   theme(strip.text.x = element_text(size = 14), 
#         strip.background = element_rect(
#           color="black", fill="white", linetype="solid"))+
#   theme(legend.margin = margin(-1,0,0,0, unit="cm"))
# 
# Pulse Disturbance Frequency vs Magnitude

# either run these (take a long time)
# source("SpA_Freq_V_Mag.R")
# source("SpB_Freq_V_Mag.R")
# source("SpC_Freq_V_Mag.R")
# source("SpD_Freq_V_Mag.R")

#or run them on the HPC and read in csv(still take a long time)


#code to make heatmap for K in response to Disturbance and time post disturbance
source("Kwireplot.R")
FigS5 <- ggplot(data = KQT, aes(x = t , y = Q))+
  geom_raster(aes(fill = K), interpolate = T)+
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

ggsave(filename = "FigS5.png", FigS5, height = 5, width = 6, device = "png", dpi = "retina")

# code to make heatmaps for pulse Freq v Mag
# 
# spA <- read.csv(file = "SpA_FreqVMag_short.csv")
# ggplot(data = spA, aes(x = frequency, y = magnitude))+
#   geom_raster(aes(fill = log(abundance)), interpolate = F)+
#   scale_fill_viridis_c(option = "magma")
# 
# spB <- read.csv(file = "SpB_FreqVMag_short.csv")
# ggplot(data = spB, aes(x = frequency, y = magnitude))+
#   geom_raster(aes(fill = log(abundance)), interpolate = F)+
#   scale_fill_viridis_c(option = "magma")
# 
# spC <- read.csv(file = "SpC_FreqVMag_short.csv")
# ggplot(data = spC, aes(x = frequency, y = magnitude))+
#   geom_raster(aes(fill = log(abundance)), interpolate = F)+
#   scale_fill_viridis_c(option = "magma")
# 
# spD <- read.csv(file = "SpD_FreqVMag_short.csv")
# ggplot(data = spD, aes(x = frequency, y = magnitude))+
#   geom_raster(aes(fill= log(abundance)), interpolate  =F)+
#   scale_fill_viridis_c(option = "magma")
# 

# AppendMe <- function(dfNames) {
#   do.call(rbind, lapply(dfNames, function(x) {
#     cbind(get(x), source = x)
#   }))
# }
# 
# freq_mag <- AppendMe(c("spA", "spB", "spC", "spD"))
# 
# TaxaNames <- list(
#   'spA'="Stonefly",
#   'spB'="Mayfly",
#   'spC'="Caddisfly",
#   'spD'="Beetle"
# )
# 
# Taxa_labeller <- function(variable,value){
#   return(TaxaNames[value])
# }
# ggplot(data = freq_mag, aes(x = frequency, y = magnitude))+
#   facet_wrap(~source, ncol = 2, labeller = Taxa_labeller)+
#   geom_raster(aes(fill = log(abundance)), interpolate = F)+
#   scale_fill_viridis_c(option = "magma") +
#   labs(shape = "") +
#   theme_bw()+
#   xlab("Pulse Frequency")+
#   #scale_y_discrete(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))+
#   ylab("Pulse Magnitude")+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
#   guides(fill=guide_legend(title="Log Abundance"))+
#   theme(strip.text.x = element_text(size = 14), 
#         strip.background = element_rect(
#           color="black", fill="white", linetype="solid"))+
#   theme(legend.margin = margin(-1,0,0,0, unit="cm"))

# multivoltinism
source("SpA_Multivolt.R")
source("SpB_Multivolt.R")
source("SpC_Multivolt.R")
source("SPD_Multivolt.R")

oneyear <- rbind(B.oneyear, C.oneyear, A.oneyear, D.oneyear)

Fig3 <- ggplot(data = oneyear, aes(x = Date, y = log(Abund), group = as.factor(MeanTemp), color = as.factor(MeanTemp)))+
  geom_line(size = 1, alpha = 0.8)+
  scale_color_manual(name = "Mean Temperature (C)", labels = c("12", "20"), values = c("#4477AA", "#EE6677"))+
  ylab("Log Adult Abundance") + 
  theme_bw()+
  scale_x_date(date_labels="%B", date_breaks  ="2 month")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle=45, size = 12.5), 
        axis.text.y = element_text(size = 13),legend.position = "bottom", legend.key = element_rect(fill = "transparent"))+
  facet_wrap(~Strategy)
ggsave(filename = "Fig3.png", Fig3, height = 5, width = 5, device = "png", dpi = "retina")



# # chaos plots
source("HilbertMetric.R")

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

x11()
FigS1 <- ggarrange(chaostestplot,
          ggarrange(chaosn, chaosts, stablen, stablets, ncol = 2, nrow = 2, vjust = 0.5, labels = c("b", "c", "d", "e")),
          labels = "a",
          nrow = 2, common.legend = T)

ggsave(filename = "FigS1.png", FigS1, height = 8.5, width = 6.5, device = "png", dpi = "retina")

# 
#####
#Table 1
source("NegExpSurv.R")
#make some copies 
high.df1<-high.df
high.df2<-high.df
high.df1$Response <- rep("A", times = 2000 )
high.df2$Response <- rep("B", times = 2000)

#combine into one large dataset
flow.response <- rbind(high.df, high.df1, high.df2, low.df)

flow <- ggplot(data = flow.response, aes(x = Q, y = Survival, color = Response))+
  geom_line(linewidth = 1, aes(linetype = Response), alpha = 0.8)+
  scale_color_manual(name = "Strategy", labels=c( "Boom", "Fast","Moderate", "Slow"), values=c("#228833","#CCBB44","#66CCEE",  "#AA3377"))+
  scale_linetype_manual(name = "Strategy", labels=c( "Boom", "Fast","Moderate", "Slow"), values=c(2,4,3,1))+
  xlab("Flow Magnitude")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5),
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
  theme_bw()
####
Fecundity <- c(1200, 500, 500, 300)
DegreeDays <- c(500, 900, 1200, 1500)
Strategy <- c("Boom", "Fast", "Moderate", "Slow")
Y1 <- c(0, 0.8, 0, 0.8)
Y2 <- c(0.2, 1, 0.2, 1)
traits <- as.data.frame(cbind(as.factor(Strategy), 
                              as.numeric(DegreeDays), 
                              as.numeric(Fecundity), 
                              as.numeric(Y1), 
                              as.numeric(Y2)))

fec <- ggplot(data = traits, aes(x = Strategy, y = Fecundity, color = Strategy))+
  geom_point(size = 3)+
  ylab("Average Fecundity")+
  xlab("Life History Strategy")+
  scale_color_manual(name = "Strategy", labels=c( "Boom", "Fast","Moderate", "Slow"), values=c("#228833","#CCBB44","#66CCEE",  "#AA3377"))+
  geom_text(
    aes(label = paste(Fecundity)),
    parse = TRUE,
    color = "black",
    nudge_y = 100,
  )+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5),
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
  theme_bw()


  
dd <- ggplot(data= traits, aes(x = Strategy, y = DegreeDays, color = Strategy))+
  geom_point(size = 3)+
  ylab("Degree Days")+
  xlab("Life History Strategy")+
  geom_text(
    aes(label = paste(DegreeDays)),
    parse = TRUE,
    color = "black",
    nudge_y = 100,
  )+
  scale_color_manual(name = "Strategy", labels=c( "Boom", "Fast","Moderate", "Slow"), values=c("#228833","#CCBB44","#66CCEE",  "#AA3377"))+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5),
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
  theme_bw()

ovi <- ggplot(data = traits, aes(x = Strategy, color = Strategy))+
  geom_segment(aes(x = Strategy, xend = Strategy, y = Y1, yend = Y2), linewidth = 3)+
  ylab("Oviposition Location")+
  xlab("Life History Strategy")+
  scale_color_manual(name = "Strategy", labels=c( "Boom", "Fast","Moderate", "Slow"), values=c("#228833","#CCBB44","#66CCEE",  "#AA3377"))+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5),
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
  theme_bw()

ggarrange(fec, dd, ovi, flow, ncol = 2, nrow = 2, common.legend = T, labels = c("a", "b", "c", "d"))

  
