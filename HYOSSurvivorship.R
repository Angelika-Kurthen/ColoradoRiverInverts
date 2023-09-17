#####################################
## Code to fit curve to Hydropsyche spp mortality rates to data
####################################
library(readxl)
library(minpack.lm)
library(tidyverse)
library(car)
library(boot)
library(data.table)
HYOSVitalRates <- read_excel("VitalRates.xlsx", sheet = "Hydropsyche Mortality Rates")
HYOSVitalRates <- as.data.frame(HYOSVitalRates)

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

surv.fit.HYOS <- flow.surv.fit(HYOSVitalRates$`Max Event Discharge/Bankfull Discharge`, HYOSVitalRates$Mortality, 0.25)
surv.df.HYOS <- flow.surv.rate(surv.fit.HYOS$m$getPars()[2] , surv.fit.HYOS$m$getPars()[1], 2, 0.001, 0.001, 0.25)

ggplot(surv.df.HYOS, aes(x = Q, y = surv))+
  geom_line()+
  geom_point(data = HYOSVitalRates, aes(x = `Max Event Discharge/Bankfull Discharge` , y = 1-(Mortality), color = Citation))+
  #coord_cartesian(ylim = c(0,1)) +
  ylab('Immediate Post-Disturbance Survival') +
  theme_bw()+
  xlab('`Max Event Discharge/Bankfull Discharge`')

# Calculate temperature dependent development time
HYOSDevRates <- read_excel("VitalRates.xlsx", sheet = "Hydropsyche Development Rates")
HYOSDevRates <- as.data.frame(HYOSDevRates)

polyfit <- nlsLM(logit(MaturationRate) ~ a*Temperature^2 + b*Temperature + c, data = HYOSDevRates, start = c(a = 1, b = 1, c = 1))

devtime <- function(x){
  y = -0.01385  *x^2+ 0.30973*x -5.72982 
  return(inv.logit(y))
}

HYOSSurvRates <- read_excel("VitalRates.xlsx", sheet = "Hydropsyche Survival Rates ")
HYOSSurvRates <- as.data.frame(HYOSSurvRates)

fit <- nlsLM(logit(Survival)~ a*Temperature^2 + b*Temperature + c, data = HYOSSurvRates, start = c(a= 1, b=1, c = 1))

TempSurv <- function(x){
  a <- -0.09934 *x^2 +  3.44127  *x -15.47038 
  return(inv.logit(a))
}
