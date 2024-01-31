
install.packages("patchwork")
library(patchwork)

## check to make sure temp regime is the same
source("AspTempToggle.R")
source("Bsp_Temp_Toggle.R")
source("CspTempToggle.R")
source("D_sp_TempToggle.R")


(atemp + btemp)/(ctemp + dtemp) + plot_annotation(tag_levels = 'A')


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
  scale_color_manual(name = "Taxa", values=c("#66CCEE", "#228833", "#EE6677", "#AA3377"))+
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

source("A_sp_DD_Toggle.R")
source("B_sp_DD_toggle.R")
source("CspDDToggle.R")
source("D_sp_DD_Toggle.R")

(add + bdd)/(cdd + ddd) + plot_annotation(tag_levels = "A")
summary(add_lm)
summary(bdd_lm)
summary(cdd_lm)
summary(ddd_lm)


source("SpA_JulianPulse.R")
source("SpB_JulianPulse.R")
source("SpC_JulianPulse.R")
source("SpD_JulianPulse.R")

(ashort + bshort)/(cshort + dshort)+ plot_annotation(tag_levels = "A")
(aresil + bresil)/(cresil + dresil)+ plot_annotation(tag_levels = "A")
(along + blong)/(clong + dlong)+ plot_annotation(tag_levels = "A")


source("SpA_PulseMagnitude.R")
source("SpB_PulseMagnitude.R")
source("SpC_PulseMagnitude.R")
source("SpD_PulseMagnitude.R")
(amag + bmag)/(cmag + dmag)+plot_annotation(tag_levels = "A")


source("SpA_Frequency.R")
source("SpB_Frequency.R")
source("SpC_Frequency.R")
source("SpD_Frequency.R")
(a_imm + b_imm)/(c_imm + d_imm) + plot_annotation(tag_levels = "A")
(a_short + b_short)/(c_short + d_short) + plot_annotation(tag_levels = "A")
(a_long + b_long)/(c_long + d_long) + plot_annotation(tag_levels = "A")

temp$dts <- as.Date(temp$dts, origin = "1970-01-01")
ggplot(data = temp[1:27,], aes(dts, Temperature))+
  geom_line(size = 1)+
  xlab("Month")+
  ylab("Temperature C")+
  theme_bw()+
  scale_x_date(date_labels="%B", date_breaks  ="2 months")+
  theme(text = element_text(size = 14), axis.text.x = element_text(size = 12.5), 
  axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

  
ggplot(data = tempsurvdf, aes(x = tem, y = temSurv))+
  geom_line(size = 1)+
  xlab("Temperature C")+
  ylab("Survival")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
       