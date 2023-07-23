##############################
# Code to calculate negative exponential fit of survivorships
##############################


# Code for HPC - tidyverse has some issues on our HPC because one of the packages is deprecated

library(minpack.lm, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")


flow.surv.fit <- function(magnitude, mortality, Qmin){
  x <- magnitude # rename x and y
  y <- 1 - mortality # convert to survival
  df <- as.data.frame(cbind(x,y))
  df <- rbind(df, c(Qmin, 1))  # make sure we specify that we have 100% survival at Qmin (0.25)
  nls.fit <- nlsLM(formula = (y ~ k* exp(-h*x)), data = df, start = c(k = 1, h = 1)) # fit to negative exponential function
  return(nls.fit)
}


flow.surv.rate <- function(h, k, max, min, interval, Qmin) {
  Q <- seq(min, max, by = interval)
  surv <- k*(exp(-h*Q))
  surv.df <- as.data.frame(cbind(Q, surv))
  surv.df$surv[which(surv.df$Q <= Qmin)] <- 1
  return(surv.df)
}


mag <- c(0.5, 1, 3)
mort <- c(0.85, 0.999, 0.9999999)
high <- flow.surv.fit(mag, mort, 0.25)
high.df <- flow.surv.rate(high$m$getPars()[2], high$m$getPars()[1],2, 0.001, 0.001, 0.25 )
high.df <- as.data.frame(cbind(high.df, rep("High Response", times = length(high.df$Q))))
colnames(high.df) <- c("Q", "Survival", "Response")


mag <- c(0.75, 1, 3)
mort <- c(0.6, 0.65, 0.75)
med <- flow.surv.fit(mag, mort, 0.25)
med.df <- flow.surv.rate(med$m$getPars()[2], med$m$getPars()[1],2, 0.001, 0.001, 0.25 )
med.df <- as.data.frame(cbind(med.df, rep("Medium Response", times = length(med.df$Q))))
colnames(med.df) <- c("Q", "Survival", "Response")

mag <- c(0.75, 1, 2, 3)
mort <- c(0.4, 0.45, 0.55, 0.6)
low <- flow.surv.fit(mag, mort, 0.25)
low.df <- flow.surv.rate(low$m$getPars()[2], low$m$getPars()[1],2, 0.001, 0.001, 0.25 )
low.df <- as.data.frame(cbind(low.df, rep("Low Response", times = length(low.df$Q))))
colnames(low.df) <- c("Q", "Survival", "Response")

df <- as.data.frame(rbind(high.df, med.df, low.df))
ggplot(df, aes(x = Q, y = Survival, col = Response))+
  geom_line()+
  #geom_point(data = HYOSVitalRates, aes(x = `Max Event Discharge/Bankfull Discharge` , y = 1-(Mortality), color = Citation))+
  #coord_cartesian(ylim = c(0,1)) +
  ylab('Immediate Post-Disturbance Survival') +
  theme_bw()+
  xlab('`Max Event Discharge/Bankfull Discharge`')
