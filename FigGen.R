
install.packages("patchwork")
library(patchwork)

## check to make sure temp regime is the same
source("AspTempToggle.R")
source("Bsp_Temp_Toggle.R")
source("CspTempToggle.R")
source("D_sp_TempToggle.R")


(atemp + btemp)/(ctemp + dtemp) + plot_annotation(tag_levels = 'A')


source("A_sp_FecundityToggle.R")
source("B_sp_fecundity_Toggle.R")
source("C_sp_FecundityToggle.R")
source("D_sp_Fecundity_Toggle.R")

(afec + bfec)/(cfec + dfec) + plot_annotation(tag_levels = "A")
summary(afec_lm)
summary(bfec_lm)
summary(cfec_lm)
summary(dfec_lm)


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
       