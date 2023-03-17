############################################################
## Code to fit curve to P. anitpodarum mortality rates to data
#############################################################

library(readxl)
library(minpack.lm)
library(tidyverse)
NZMSVitalRates <- read_excel("VitalRates.xlsx", sheet = "NZMS Mortality Rates")
NZMSVitalRates <- as.data.frame(NZMSVitalRates)

flow.surv.fit <- function(magnitude, mortality, Qmin){
  x <- magnitude # rename x and y
  y <- 1 - mortality # convert to survival
  df <- as.data.frame(cbind(x,y))
  df <- rbind(df, c(Qmin, 1.375))  # make sure we specify that we have 100% survival at Qmin (0.25)
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

surv.fit.NZMS <- flow.surv.fit(NZMSVitalRates$`Max Event Discharge/Bankfull Discharge`, NZMSVitalRates$Mortality, 0.25)
surv.df.NZMS <- flow.surv.rate(surv.fit.NZMS$m$getPars()[2] , surv.fit.NZMS$m$getPars()[1], 2, 0.001, 0.001, 0.25)

ggplot(surv.df.NZMS, aes(x = Q, y = surv))+
  geom_line()+
  geom_point(data = NZMSVitalRates, aes(x = `Max Event Discharge/Bankfull Discharge` , y = 1-(Mortality), color = Citation))+
  coord_cartesian(ylim = c(0,1)) +
  ylab('Immediate Post-Disturbance Survival') +
  theme_bw()+
  xlab('`Max Event Discharge/Bankfull Discharge`')

# calculating temperature fecundity relationship for NZMS using a power function
# log(y)=log(b*x^z)=log(b)+z*log(x)
# In this case lm(log(y)~log(x) ) solve your problem

# max recruitment between 16 and 19 C (Dybahl and Kane)
# recruitment stops below 9C (Bennett )
# above 27 C everything stops working well (Dybahl and Kane)
x <- c(9, 16, 17.5, 19, 27)
y <- c(0.001, 1, 1, 1, 0.001)
data <- as.data.frame(cbind(x,y))
nlsLM(y~ -b*(x - 17.5)^4 + 1,start = list(b = 0),data=data)

# eq = y = -0.0001427(x - 17.5)^4 + 1 * F 

