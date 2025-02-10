#########################
## Code to Validate NZMS model
###########################
library(purrr)
library(tidyverse)
library(lubridate)
library(plyr)
library(dplyr)
library(ggplot2)
# data retrieval tool from USGS
library(dataRetrieval)
#install.packages("devtools")
library(devtools)
#install_github(repo = "jmuehlbauer-usgs/R-packages", subdir = "foodbase")
library(foodbase)


source("NZMS_1sp_Model.R")

discharge <- readNWISdv("09380000", "00060", "2007-10-01", "2023-05-01")
flow.magnitude <- TimestepDischarge(discharge, 85000)
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
temps <- TimestepTemperature(temp)

set.seed(333)
out <- NZMSmodel(flow.data = flow.magnitude$Discharge, temp.data = temps, disturbanceK = 9000, baselineK = 5000, Qmin = 0.3, extinct = 50, iteration = 1000, peaklist = 0.17, peakeach = length(temps$Temperature))


# adults<-as.data.frame(cbind(as.Date(temps$dts), out[1:length(temps$dts),2:3,1]))
# colnames(adults) <- c("Time","Adult")
# adults$Time <- as.Date(adults$Time, origin = "1970-01-01")

means.list.NZMS <- mean.data.frame(out,burnin = 200, iteration= 1000)
means.list.NZMS <- cbind(means.list.NZMS, temps$dts[199:length(temps$dts)])
means.list.NZMS$`temps$dts` <- as.Date(means.list.NZMS$`temps$dts`)
# plot abundance over time

drift.data.total <- readDB(gear = "Drift", type = "Sample", updater = F)

drift.LF <- drift.data.total[which(drift.data.total$RiverMile >= -6 & drift.data.total$RiverMile <= 0),]

NZMS.samp.LF <- sampspec(samp = drift.LF, species = "NZMS", stats = T)
NZMS.samp <- merge(NZMS.samp.LF$Statistics, NZMS.samp.LF$Samples, by = "BarcodeID", all = T)
NZMS.samp <- NZMS.samp[which(NZMS.samp$GearID == 4),] 
NZMS.samp$Density <- NZMS.samp$CountTotal/NZMS.samp$Volume

NZMS.samp <- merge(NZMS.samp, discharge[, 3:4], by = "Date")
NZMS.samp$Density <- (NZMS.samp$Density/(9e-15 *(NZMS.samp$X_00060_00003 * 0.02831683)^4.1))
NZMS.samp <- aggregate(NZMS.samp$Density, list(NZMS.samp$Date), FUN = mean)

means <- vector()
for (i in 1:length(temps$dts)){
  d <- NZMS.samp[which(NZMS.samp$Group.1 %within% interval(temps$dts[i], temps$dts[i+1]) == T),]
  if (is.nan(mean(d$x)) == T) {
    s = NA
  } else {
    s<- mean(d$x)}
  means <- append(means, s)
}

# test for temporal autocorrelation
#acf(na.omit(means)) #we have some

NZMS.samp.sum <- na.omit(as.data.frame(cbind(as.Date(means.list.NZMS$`temps$dts`, format = "%Y-%m-%d"), means[201:406])))
NZMS.samp.sum$V1 <- as.Date(NZMS.samp.sum$V1, origin = "1970-01-01")

# culling data by lags - Temporal autocorrelation: a neglected factor in the study of behavioral repeatability and plasticity  

cor.df <- left_join(NZMS.samp.sum, means.list.NZMS, by=c('V1'="temps$dts"), copy = T)
cor.lm <- lm((cor.df$mean.abund) ~ (cor.df$V2))
cor.test((cor.df$V2+1), (cor.df$mean.abund+1), method = "spearman")


# culling data by lags - Temporal autocorrelation: a neglected factor in the study of behavioral repeatability and plasticity  
NZMS.samp.sum1 <- NZMS.samp.sum %>% slice(which(row_number() %% 3 == 0))
NZMS.samp.sum2 <- NZMS.samp.sum %>%  slice(which(row_number() %% 3 == 1))
NZMS.samp.sum3 <- NZMS.samp.sum %>%  slice(which(row_number() %% 3 == 2))

cor.df1 <- left_join(NZMS.samp.sum1, means.list.NZMS, by=c('V1'="temps$dts"), copy = T)
cor.lm1 <- lm(cor.df1$mean.abund ~ cor.df1$V2)
summary(cor.lm1)
plot(cor.df1$V2, cor.df1$mean.abund)
rho1 <- cor.test((cor.df1$V2+1), (cor.df1$mean.abund+1), method = "spearman")

cor.df2 <- left_join(NZMS.samp.sum2, means.list.NZMS, by=c('V1'="temps$dts"), copy = T)
cor.lm2 <- lm(cor.df2$mean.abund ~ cor.df2$V2)
rho2 <- cor.test((cor.df2$V2+1), (cor.df2$mean.abund+1), method = "spearman")

cor.df3 <- left_join(NZMS.samp.sum3, means.list.NZMS, by=c('V1'="temps$dts"), copy = T)
cor.lm3 <- lm(cor.df3$mean.abund ~ cor.df3$V2)
summary(cor.lm3)
rho3 <- cor.test((cor.df3$V2+1), (cor.df3$mean.abund+1), method = "spearman")

rho <- mean(c(rho1$estimate, rho2$estimate, rho3$estimate))
sd <- sd(c(rho1$estimate, rho2$estimate, rho3$estimate))

# summary(cor.df)
# ggplot(data = cor.df, aes(x = (V2) , y = (mean.abund)))+
#   geom_point()+
#   stat_smooth(method = "lm",
#               formula = y ~ x,
#               geom = "smooth")+
#   geom_text(x = 3e+05, y = 6100, label = "")+
#   labs(y = "NZMS Model Output", x = "NZMS Emprical Data")


rmse.nzms<- sqrt(mean((cor.df$V2 - cor.df$mean.abund)^2))

rmse.nzms.scale <- sqrt(mean((scale(cor.df$V2) - scale(cor.df$mean.abund))^2))
# coverage 
coverage <- mean(scale(cor.df$V2) >= (scale(cor.df$mean.abund) - (1.96*rmse.nzms.scale)) & scale(cor.df$V2) <= (scale(cor.df$mean.abund) + (1.96*rmse.nzms.scale)))
colors <- c("black","#4477AA")
linetypes <- c("solid", "solid")
NZMSts <- ggplot(data = cor.df, aes(x = V1, y = scale(mean.abund), group = 1, color = "Model")) +
  geom_ribbon(aes(ymin = scale(mean.abund) - 1.96 * rmse.nzms.scale,
                  ymax = scale(mean.abund) + 1.96 * rmse.nzms.scale),
              colour = 'transparent',
              alpha = .1,
              fill = "black",
              show.legend = F)+
  geom_line(show.legend = T, linewidth = 1, alpha = 0.8) +
  geom_line(data = cor.df, aes(x =V1, y = scale(V2), color = "P. antipodarum"), linewidth = 1, alpha = 0.8, show.legend = T)+
  #geom_point(data = NZMS.samp.sum[125,], aes(x = V1, y = scale(means), color = "Empirical"), show.legend = T)+
  labs(y=expression(paste(italic("P. antipodarum"), " Abund.")))+
  ylim(c(-4,7))+
  geom_text(mapping = aes(x = as.Date("2019-06-01"), y =5, label = paste('rho', "==", 0.63)), parse = T, color = "black", size = 4.5)+
  geom_text(mapping = aes(x = as.Date("2019-06-01"), y =5.75, label = paste('C = 94%')), color = "black", size = 4.5)+
  geom_text(mapping = aes(x = as.Date("2019-06-01"), y =6.5, label = paste('Scaled RMSE = 1.19')), color = "black", size = 4.5)+
  xlab("")+
  labs(colour=" ")+
  theme_bw()+
  scale_color_manual(values = colors)+
  theme(text = element_text(size = 13), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), )+
  scale_x_date(date_labels="%Y")
##################
# N mix models
##################
library(rjags)
library(jagsUI)
library(MCMCvis)

drift.data.total <- readDB(gear = "Drift", type = "Sample", updater = F)
# get drift data from between Lees Ferry and RM -6
drift.LF <- drift.data.total[which(drift.data.total$RiverMile >= -6 & drift.data.total$RiverMile <= 0),]
#specify
NZMS.samp.LF <- sampspec(samp = drift.LF, species = c("NZMS"), stats = T)
# pull stats and merge with sample info
NZMS.samp <- merge(NZMS.samp.LF$Statistics, NZMS.samp.LF$Samples, by = "BarcodeID", all = T)
# make sure we are using the same gear
NZMS.samp <- NZMS.samp[which(NZMS.samp$GearID == 4),] 
NZMS.samp <- NZMS.samp[which(NZMS.samp$FlagStrange == 0), ]

vals <- vector()
for (i in 1:length(temps$dts)){
  d <- NZMS.samp[which(NZMS.samp$Date %within% interval(temps$dts[i], temps$dts[i+1]-1) == T),]
  if (length(d$CountTotal) > 0) {
    s<- rep(i, times = length(d$CountTotal))
    vals <- append(vals, s)}
}

# add to data frame
NZMS.samp <- cbind(NZMS.samp, vals)

# now we need into include mean water temperature and discharge for each 
NZMS.samp <- cbind(NZMS.samp, temps$Temperature[NZMS.samp$vals])
NZMS.samp <- cbind(NZMS.samp, flow.magnitude$Discharge[NZMS.samp$vals])
max_visits <- vector() # vector to put max obs per site
means <- vector()
for (i in 1:length(temps$dts)){  # cycle through all the 14 day timesteps that we have model output for
  # pull abundances between each each timestep
  d <- NZMS.samp[which(NZMS.samp$Date >= temps$dts[i] & NZMS.samp$Date < temps$dts[i+1]), ]
  max_visits[i] <- length(d$Date) # number of observations
  means[i] <- mean(d$CountTotal) # means - this is just for checking NAs
}
# make data frame with the timestep, the mean
df <- cbind.data.frame(temps$dts, means)
# remove anythig where there are NA values of Density (means that volume or count is missing)
df <- df[!is.na(df$means), ]
# phenology may also play into dynamics, so include month column as well
#month <- month(NZMS.samp$Date)
#install.packages("aspace")
library(aspace)
df$circdate <- sin(as_radians((lubridate::yday(df$`temps$dts`)/365)*360))

# define our RxJ matrix
R <- length(temps$dts)
J <- max(max_visits)
site_mat <- matrix(data = NA, nrow = R, ncol = J)
dens_mat <- matrix(data = NA, nrow = R, ncol = J)
obs_intercept <- matrix(data = 1, nrow = R, ncol = J)
# make vector for flows at each timestep
# make RxJ matrix full of densities
# make RxJ matrix full of raw counts
# make RxJ matrix full of volumes sampled for each abundance
flows <- vector()
temperature <- vector()
volumes <- matrix(data = NA, nrow = R, ncol = J)
for (i in 1:length(temps$dts)){
  d <- NZMS.samp[which(NZMS.samp$Date >= temps$dts[i] & NZMS.samp$Date < temps$dts[i+1]), ]
  flows[i] <- mean(d$`flow.magnitude$Discharge[NZMS.samp$vals]`)
  #date <- df[]
  temperature[i] <- mean(d$`temps$Temperature[NZMS.samp$vals]`)
  #windspeed[i, ] <- c(d$WindSpeed, rep(NA, times = (J- length(d$CountTotal))))
  site_mat[i, ] <- c(d$CountTotal, rep(NA, times = (J- length(d$CountTotal))))
  #habitat[i, ] <- c(d$Habitat, rep(NA, times = (J- length(d$CountTotal))))
  #time[i, ] <- c(d$TimeElapsed,rep(NA, times = (J- length(d$CountTotal))))
  #weather[i, ] <- c(d$Weather, rep(NA, times = (J- length(d$CountTotal))))
  volumes[i, ] <- c((d$Volume),rep(NA, times = (J- length(d$CountTotal))))
}

# we need to remove all timesteps that are just NAs
nodata <- which(is.na(site_mat[,1]))
# first identify all the timsteps that don't have data (so we can match them up later)
site_mat <- as.matrix(site_mat[-nodata,]) # count data
#dens_mat <- as.matrix(dens_mat[-nodata, ]) # density data
obs_intercept <- as.matrix(obs_intercept[-nodata,]) # intercept for obs cov

flows <- as.data.frame(scale(flows[-nodata])) # site cov flow 
temperature <- as.data.frame(scale(temperature[-nodata])) # site cov temp
circdate <- as.data.frame(df$circdate)
#windspeed <- as.matrix((scale(windspeed[-nodata,])))
##windspeed[is.na(windspeed)] <- mean(windspeed, na.rm = TRUE) # replace NAs with mean duration time since NAs not allowed in predictors or offsets
#habitat <- as.matrix(habitat[-nodata,])
#Mode <- function(x) {
#  ux <- unique(x)
#  ux[which.max(tabulate(match(x, ux)))]
#}
#habitat[is.na(habitat)] <- Mode(na.omit(habitat)) # replace NAs with mean duration time since NAs not allowed in predictors or offsets
#weather <- as.matrix(weather[-nodata,])
#weather[is.na(weather)] <- Mode(na.omit(weather))

site_intercept <- rep(1, times = length(flows$V1)) 
site_covs<- as.matrix(cbind(site_intercept, flows, circdate)) #flows,temperature, circdate)
obs_covs <- array(data= NA, dim = c(length(flows$V1),J,1))
obs_covs[,,1] <- obs_intercept                                  

#offset
offset <- as.matrix(scale(log(volumes[-nodata, ])))
offset[is.na(offset)] <- mean(offset, na.rm = TRUE) # replace NAs with mean duration time since NAs not allowed in predictors or offsets
# 

sink("N-mixturePoisNZMS.jags")
cat("
model{
    # Priors
    for(i in 1:nAlpha){ # nAlpha is the number of site predictor variables
      alpha[i] ~ dnorm(0, 0.1) # alphas are the site covariates
    }

    for(i in 1:nBeta){ # nBeta is number of site x observation predictors
      beta[i] ~ dnorm(0, 0.1) # betas are the observation covariates
    }

    # Likelihood
    for(r in 1:R){
     N[r] ~ dpois(lambda[r]) #start with pulling from Poisson
     log(lambda[r]) <- sum(alpha * XN[r, ]) #XN is matrix of site covariates

     for(j in 1:J){
      y[r, j] ~ dbinom(p[r, j], N[r]) #binomial for observed counts
      logit(p[r, j]) <- sum(off[r, j] + beta * Xp[r,j,]) #XP(obs covariates)

        ## Expected count at site r, sample j
        exp[r,j] <- N[r] * p[r, j]

        ## Discrepancy
        ## (note small value added to denominator to avoid potential divide by zero)
        ## This is the X2 + descrepancy
        E[r, j] <- pow((y[r, j] - exp[r, j]), 2) / (exp[r, j] + 0.5)

        ## Simulate new count from model
        y.rep[r, j] ~ dbinom(p[r, j], N[r])

        ## X2
        E.rep[r, j] <- pow((y.rep[r, j] - exp[r, j]), 2) / (exp[r, j] + 0.5)
     }
    }
     # chi-squared test statistics
    fit <- sum(E[,])
    fit.rep <- sum(E.rep[ , ])
    } # End model
", fill = TRUE)
sink()

jags_data <- list(y = (site_mat),
                  XN = site_covs,
                  Xp = (obs_covs),
                  J = dim(site_mat)[2], #visits
                  R = dim(site_mat)[1], #sites
                  off = (offset), #obs offset
                  nAlpha = dim(site_covs)[2],
                  nBeta = dim(obs_covs)[3])

nAlpha <- dim(site_covs)[2]
nBeta <- dim(obs_covs)[3]
jags_inits <- function(){
  list(
    N = apply(jags_data$y, 1, max, na.rm=TRUE),
    alpha=runif(nAlpha,-1,1),
    beta=runif(nBeta,-1,1))}

parameters <- c("alpha", "beta", "lambda", "p", "N")


nc <- 3
ni <- 10000
nb <- 2500
nt <- 1

Nmix_fit <- jags.model("N-mixturePoisNZMS.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)

Nmix_fit_UI1 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixturePoisNZMS.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)
#
print(Nmix_fit_UI1)

zm = coda.samples(Nmix_fit, variable.names = c("lambda"), n.iter = ni, n.thin = nt)


lam <- MCMCpstr(zm, "lambda")
# N <- MCMCpstr(zm, "N")
# N <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(N)))
# N$V1 <- as.Date(N$V1, origin = "1970-01-01")


lam <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(lam)))
lam$V1 <- as.Date(lam$V1, origin = "1970-01-01")

# y.rep <- MCMCpstr(zm, "y.rep")
# exp <- MCMCpstr(zm, "exp")
#
# plot(unlist(y.rep), unlist(site_mat))
#
# plot(unlist(y.rep), unlist(site_mat))
# abline(0, 1)
#
# fit <- MCMCchains(zm, "fit")
# fit.rep <- MCMCchains(zm, "fit.rep")
# # mean(fit > fit.rep) # close to 1 so bad fit?
# plot(fit.rep ~ fit)
# abline(0, 1) # 1 to 1 line not even there


# fit_df <- data.frame(y = c(c(unlist(site_mat)), c(unlist(y.rep))),
#                      data = rep(c("Observed", "Simulated"), each = length(site_mat)))
# library(ggplot2)
# ggplot(fit_df, aes(x = y, fill = data)) + geom_histogram() + facet_grid(.~data)# not getting all the 0s and missing the really high #s

#cor.df <- left_join(N, means.list.NZMS, by=c('V1'="Date"), copy = T)
cor.df <- na.omit(left_join(lam, means.list.NZMS, by=c('V1'="temps$dts"), copy = T))

cor.df1 <- cor.df %>% slice(which(row_number() %% 2 == 0))
rho1 <- cor.test((cor.df1$V2), (cor.df1$mean.abund), method = "spearman")

cor.df2 <- cor.df %>%  slice(which(row_number() %% 2 == 1))
rho2 <- cor.test(cor.df2$V2, cor.df2$mean.abund, method = "spearman")

cor.df3 <- cor.df %>%  slice(which(row_number() %% 3 == 2))
rho3 <- cor.test(cor.df3$V2, cor.df3$mean.abund, method = "spearman")

Poislam  <- mean(c(rho1$estimate, rho2$estimate, rho3$estimate))
#cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2)
#PoisN <- cor.test((cor.df$V2.x), (cor.df$mean.abund), method = "spearman")

#cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2)
#Poislam <- cor.test((cor.df$V2.x), (cor.df$mean.abund), method = "spearman")

rm(zm)
#rm(Nmix_fit_UI1)
rm(Nmix_fit)

sink("N-mixtureZIPNZMS.jags")
cat("
model{
    # Priors
    omega ~ dbeta(1,1)

    for(i in 1:nAlpha){ # nAlpha is the number of site predictor         variables
      alpha[i] ~ dnorm(0, 0.1) # alphas are the site covariates
    }

    for(i in 1:nBeta){ # nBeta is number of site x observation predictors
      beta[i] ~ dnorm(0, 0.1) # betas are the observation covariates
    }

    # Likelihood
    for(r in 1:R){
     z[r] ~ dbern(omega) # either there or not
     N[r] ~ dpois(lambda[r] * z[r]) #start with pulling from Poisson with z variable
     log(lambda[r]) <- sum(alpha * XN[r, ]) #XN is matrix of site covariates

     for(j in 1:J){
      y[r, j] ~ dbinom(p[r, j], N[r]) #binomial for observed counts
      logit(p[r, j]) <- sum(off[r, j] + beta * Xp[r,j,]) #XP(obs covariates)

        ## Expected count at site r, sample j
        exp[r,j] <- N[r] * p[r, j]

        ## Discrepancy
        ## (note small value added to denominator to avoid potential divide by zero)
        ## This is the X2 + descrepancy
        E[r, j] <- pow((y[r, j] - exp[r, j]), 2) / (exp[r, j] + 0.5)

        ## Simulate new count from model
        y.rep[r, j] ~ dbinom(p[r, j], N[r])

        ## X2
        E.rep[r, j] <- pow((y.rep[r, j] - exp[r, j]), 2) / (exp[r, j] + 0.5)
     }
    }
     # chi-squared test statistics
    fit <- sum(E[,])
    fit.rep <- sum(E.rep[,])
    } # End model
", fill = TRUE)
sink()

jags_data <- list(y = (site_mat),
                  XN = site_covs,
                  Xp = (obs_covs),
                  J = dim(site_mat)[2], #visits
                  R = dim(site_mat)[1], #sites
                  off = (offset), #obs offset
                  nAlpha = dim(site_covs)[2],
                  nBeta = dim(obs_covs)[3])

nAlpha <- dim(site_covs)[2]
nBeta <- dim(obs_covs)[3]
jags_inits <- function(){
  list(
    omega = runif(1, 0.5, 0.7),
    N = apply(jags_data$y, 1, max, na.rm=TRUE),
    alpha=runif(nAlpha,-1,1),
    beta=runif(nBeta,-1,1))}

parameters <- c("omega", "alpha", "beta", "lambda", "p", "N")

nc <- 3
ni <- 10000
nb <- 2500
nt <- 1

Nmix_fit <- jags.model("N-mixtureZIPNZMS.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)

Nmix_fit_UI2 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixtureZIPNZMS.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)
#
print(Nmix_fit_UI2)

zm = coda.samples(Nmix_fit, variable.names = c("lambda"), n.iter = ni, n.thin = nt)


lam <- MCMCpstr(zm, "lambda")
# N <- MCMCpstr(zm, "N")
# N <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(N)))
# N$V1 <- as.Date(N$V1, origin = "1970-01-01")


lam <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(lam)))
lam$V1 <- as.Date(lam$V1, origin = "1970-01-01")

# cor.df <- left_join(N, means.list.NZMS, by=c('V1'="Date"), copy = T)
# cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2.x)
# ZipN <-  cor.test((cor.df$V2.x), (cor.df$mean.abund), method = "spearman")

cor.df <- na.omit(left_join(lam, means.list.NZMS, by=c('V1'="temps$dts"), copy = T))


cor.df1 <- cor.df %>% slice(which(row_number() %% 2 == 0))
rho1 <- cor.test((cor.df1$V2), (cor.df1$mean.abund), method = "spearman")

cor.df2 <- cor.df %>%  slice(which(row_number() %% 2 == 1))
rho2 <- cor.test(cor.df2$V2, cor.df2$mean.abund, method = "spearman")

cor.df3 <- cor.df %>%  slice(which(row_number() %% 3 == 2))
rho3 <- cor.test(cor.df3$V2, cor.df3$mean.abund, method = "spearman")

Ziplam  <- mean(c(rho1$estimate, rho2$estimate, rho3$estimate))
#cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2)
rm(zm)
rm(Nmix_fit)
#rm(Nmix_fit_UI2)


sink("N-mixtureZIPoverdispNZMS.jags")
cat("
model{
    # Priors
    omega ~ dbeta(1,1)

    tau.p <- pow(sd.p, -2)
    sd.p ~ dunif(0,3)

    for(i in 1:nAlpha){ # nAlpha is the number of site predictor         variables
      alpha[i] ~ dnorm(0, 0.1) # alphas are the site covariates
    }

    for(i in 1:nBeta){ # nBeta is number of site x observation predictors
      beta[i] ~ dnorm(0, 0.1) # betas are the observation covariates
    }

    # Likelihood
    for(r in 1:R){
     z[r] ~ dbern(omega) # either there or not
     N[r] ~ dpois(lambda[r] * z[r]) #start with pulling from Poisson with z variable
     log(lambda[r]) <- sum(alpha * XN[r, ]) #XN is matrix of site covariates

     for(j in 1:J){
      y[r, j] ~ dbinom(p[r, j], N[r]) #binomial for observed counts
      logit(p[r, j]) <- lp[r,j]
      mu.lp[r, j] <- sum(off[r, j] + beta * Xp[r,j,]) #XP(obs covariates)
      lp[r,j] ~ dnorm(mu.lp[r,j], tau.p) #sample effect based on mean p

        ## Expected count at site r, sample j
        exp[r,j] <- N[r] * p[r, j]

        ## Discrepancy
        ## (note small value added to denominator to avoid potential divide by zero)
        ## This is the X2 + descrepancy
        E[r, j] <- pow((y[r, j] - exp[r, j]), 2) / (exp[r, j] + 0.5)

        ## Simulate new count from model
        y.rep[r, j] ~ dbinom(p[r, j], N[r])

        ## X2
        E.rep[r, j] <- pow((y.rep[r, j] - exp[r, j]), 2) / (exp[r, j] + 0.5)
     }
    }
     # chi-squared test statistics
    fit <- sum(E[,])
    fit.rep <- sum(E.rep[,])
    } # End model
", fill = TRUE)
sink()

jags_data <- list(y = (site_mat),
                  XN = site_covs,
                  Xp = (obs_covs),
                  J = dim(site_mat)[2], #visits
                  R = dim(site_mat)[1], #sites
                  off = (offset), #obs offset
                  nAlpha = dim(site_covs)[2],
                  nBeta = dim(obs_covs)[3])

nAlpha <- dim(site_covs)[2]
nBeta <- dim(obs_covs)[3]
jags_inits <- function(){
  list(
    omega = runif(1, 0.5, 0.7),
    sd.p = runif(1, 0.3, 0.7),
    N = apply(jags_data$y, 1, max, na.rm=TRUE),
    alpha=runif(nAlpha,-1,1),
    beta=runif(nBeta,-1,1))}

parameters <- c("omega", "alpha", "beta", "lambda", "p", "N")

nc <- 3
ni <- 10000
nb <- 2500
nt <- 1

Nmix_fit <- jags.model("N-mixtureZIPoverdispNZMS.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)

Nmix_fit_UI3 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixtureZIPoverdispNZMS.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)

print(Nmix_fit_UI3)

zm = coda.samples(Nmix_fit, variable.names = c("lambda"), n.iter = ni, n.thin = nt)


lam <- MCMCpstr(zm, "lambda")
# p <- MCMCpstr(zm, "p")
# N <- MCMCpstr(zm, "N")
# N <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(N)))
# N$V1 <- as.Date(N$V1, origin = "1970-01-01")
# 

lam <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(lam)))
lam$V1 <- as.Date(lam$V1, origin = "1970-01-01")

#plot(unlist(y.rep), unlist(site_mat))
#abline(0, 1)

#
# mean(fit > fit.rep) # close to 1 so bad fit?
# plot(fit.rep ~ fit)
# abline(0, 1) # 1 to 1 line not even there


# fit_df <- data.frame(y = c(c(unlist(site_mat)), c(unlist(y.rep))),
#                     data = rep(c("Observed", "Simulated"), each = length(site_mat)))
# library(ggplot2)
# ggplot(fit_df, aes(x = y, fill = data)) + geom_histogram() + facet_grid(.~data) #still not getting all the 0s and missing the really high #s

cor.df <- na.omit(left_join(lam, means.list.NZMS, by=c('V1'="temps$dts"), copy = T))

cor.df1 <- cor.df %>% slice(which(row_number() %% 2 == 0))
rho1 <- cor.test((cor.df1$V2), (cor.df1$mean.abund), method = "spearman")

cor.df2 <- cor.df %>%  slice(which(row_number() %% 2 == 1))
rho2 <- cor.test(cor.df2$V2, cor.df2$mean.abund, method = "spearman")

cor.df3 <- cor.df %>%  slice(which(row_number() %% 3 == 2))
rho3 <- cor.test(cor.df3$V2, cor.df3$mean.abund, method = "spearman")

Zip_ovdlam  <- mean(c(rho1$estimate, rho2$estimate, rho3$estimate))
#cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2)
# cor.df <- left_join(N, means.list.NZMS, by=c('V1'="V2"), copy = T)
# #cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2.x)
# Zip_ovdN <- cor.test((cor.df$V2), (cor.df$mean.abund), method = "spearman")
# 
# cor.df <- left_join(lam, means.list.NZMS, by=c('V1'="V2"), copy = T)
# #cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2.x)
# Zip_ovdlam <- cor.test((cor.df$V2), (cor.df$mean.abund), method = "spearman")

rm(zm)
#rm(Nmix_fit_UI3)
rm(Nmix_fit)

sink("N-mixturePoisoverdispNZMS.jags")
cat("
model{
    # Priors
    for(i in 1:nAlpha){ # nAlpha is the number of site predictor variables
      alpha[i] ~ dnorm(0, 0.1) # alphas are the site covariates
    }

    for(i in 1:nBeta){ # nBeta is number of site x observation predictors
      beta[i] ~ dnorm(0, 0.1) # betas are the observation covariates
    }

    tau.p <- pow(sd.p, -2)
    sd.p ~ dunif(0,3)

    # Likelihood
    for(r in 1:R){
     N[r] ~ dpois(lambda[r]) #start with pulling from Poisson
     log(lambda[r]) <- sum(alpha * XN[r, ]) #XN is matrix of site covariates

     for(j in 1:J){
      y[r, j] ~ dbinom(p[r, j], N[r]) #binomial for observed counts
      logit(p[r, j]) <- lp[r,j]
      mu.lp[r, j] <- sum(off[r, j] + beta * Xp[r,j,]) #XP(obs covariates)
      lp[r,j] ~ dnorm(mu.lp[r,j], tau.p) #sample effect based on mean p

        ## Expected count at site r, sample j
        exp[r,j] <- N[r] * p[r, j]

        ## Discrepancy
        ## (note small value added to denominator to avoid potential divide by zero)
        ## This is the X2 + descrepancy
        E[r, j] <- pow((y[r, j] - exp[r, j]), 2) / (exp[r, j] + 0.5)

        ## Simulate new count from model
        y.rep[r, j] ~ dbinom(p[r, j], N[r])

        ## X2
        E.rep[r, j] <- pow((y.rep[r, j] - exp[r, j]), 2) / (exp[r, j] + 0.5)
     }
    }
     # chi-squared test statistics
    fit <- sum(E[,])
    fit.rep <- sum(E.rep[,])
    } # End model
", fill = TRUE)
sink()

jags_data <- list(y = (site_mat),
                  XN = site_covs,
                  Xp = (obs_covs),
                  J = dim(site_mat)[2], #visits
                  R = dim(site_mat)[1], #sites
                  off = (offset), #obs offset
                  nAlpha = dim(site_covs)[2],
                  nBeta = dim(obs_covs)[3])

nAlpha <- dim(site_covs)[2]
nBeta <- dim(obs_covs)[3]
jags_inits <- function(){
  list(
    sd.p = runif(1, 0.5, 0.8),
    N = apply(jags_data$y, 1, max, na.rm=TRUE),
    alpha=runif(nAlpha,-1,1),
    beta=runif(nBeta,-1,1))}

parameters <- c("omega", "alpha", "beta", "lambda", "p", "N")

nc <- 3
ni <- 10000
nb <- 2500
nt <- 1

Nmix_fit <- jags.model("N-mixturePoisoverdispNZMS.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)

Nmix_fit_UI4 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixturePoisoverdispNZMS.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)

print(Nmix_fit_UI4)
#rm(Nmix_fit_UI4)
zm = coda.samples(Nmix_fit, variable.names = c("lambda"), n.iter = ni, n.thin = nt)

lam <- MCMCpstr(zm, "lambda")
# p <- MCMCpstr(zm, "p")
# N <- MCMCpstr(zm, "N")
# N <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(N)))
# N$V1 <- as.Date(N$V1, origin = "1970-01-01")
# 

lam <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(lam)))
lam$V1 <- as.Date(lam$V1, origin = "1970-01-01")

# close to 1 so bad fit?
# plot(fit.rep ~ fit)
# abline(0, 1) # 1 to 1 line not even there

#cor.df <- left_join(N, means.list.NZMS, by=c('V1'="V2"), copy = T)
#cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2.x)
#Pois_ovdN <- cor.test((cor.df$V2), (cor.df$mean.abund), method = "spearman")


#cor.df <- left_join(lam, means.list.NZMS, by=c('V1'="V2"), copy = T)
#cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2.x)
#Pois_ovdlam <- cor.test((cor.df$V2), (cor.df$mean.abund), method = "spearman")


cor.df <- na.omit(left_join(lam, means.list.NZMS, by=c('V1'="temps$dts"), copy = T))

cor.df1 <- cor.df %>% slice(which(row_number() %% 2 == 0))
rho1 <- cor.test((cor.df1$V2), (cor.df1$mean.abund), method = "spearman")

cor.df2 <- cor.df %>%  slice(which(row_number() %% 2 == 1))
rho2 <- cor.test(cor.df2$V2, cor.df2$mean.abund, method = "spearman")


cor.df3 <- cor.df %>%  slice(which(row_number() %% 3 == 2))
rho3 <- cor.test(cor.df3$V2, cor.df3$mean.abund, method = "spearman")

Pois_ovdlam  <- mean(c(rho1$estimate, rho2$estimate, rho3$estimate))
#cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2)


rm(zm)
rm(Nmix_fit)
#rm(Nmix_fit_UI4)

sink("N-mixtureNBNZMS.jags")
cat("
model{

# State model
for (r in 1:R){
 N[r] <- n[r] * z[r]
  n[r] ~ dnegbin(s[r], th)
  s[r] <- th / (th + lambda[r])
  log(lambda[r]) <- sum(alpha * XN[r, ]) #XN is matrix of site covariates
   z[r] ~ dbern(omega)
}

omega ~ dbeta(1,1)
th ~ dgamma(0.01, 0.01)
phi <- 1/th
theta <- th

  for(i in 1:nAlpha){ # nAlpha is the number of site predictor variables
      alpha[i] ~ dnorm(0, 0.01) # alphas are the site covariates
    }


# Detection model
for (r in 1:R){
  for (j in 1:J){
    logit(p[r,j]) <- max(1e-5, min (0.999999, sum(off[r, j] + (beta * Xp[r,j,]))))
    y[r,j] ~ dbinom(p[r,j], N[r])
  }
}

 for(i in 1:nBeta){ # nBeta is number of site x observation predictors
      beta[i] ~ dnorm(0, 0.01) # betas are the observation covariates
    }


# Fit statistic for real data
for (r in 1:R){
  for (j in 1:J){
    yhat[r,j] <- N[r] * p[r,j] + 0.001 # add small value to avoid divide by zero
    chi2[r,j] <- (y[r,j] - yhat[r,j])^2 / yhat[r,j]
  }
}
fit <- sum(chi2)

# Fit statistic for simulated data
for (r in 1:R){
  for (j in 1:J){
    y_new[r,j] ~ dbinom(p[r,j], N[r]) # simulate new datapoint
    chi2_new[r,j] <- (y_new[r,j] - yhat[r,j])^2 / yhat[r,j]
  }
}
fit_new <- sum(chi2_new)

sumN <- sum(N[])

}
", fill = TRUE)
sink()

jags_data <- list(y = (site_mat),
                  XN = site_covs,
                  Xp = (obs_covs),
                  J = dim(site_mat)[2], #visits
                  R = dim(site_mat)[1], #sites
                  off = (offset), #obs offset
                  nAlpha = dim(site_covs)[2],
                  nBeta = dim(obs_covs)[3])



nAlpha <- dim(site_covs)[2]
nBeta <- dim(obs_covs)[3]
jags_inits <- function(){
  list(
    omega = runif(1, 0.5, 0.7),
    n = apply(jags_data$y, 1, max, na.rm=TRUE),
    alpha=runif(nAlpha,-1,1),
    beta=runif(nBeta,-1,1))}

parameters <- c("alpha", "beta", "lambda", "p", "N", "theta", "phi", "fit", "fit_new")

nc <- 3
ni <- 10000
nb <- 2500
nt <- 1

Nmix_fit <- jags.model("N-mixtureNBNZMS.jags",data = jags_data, inits = jags_inits, n.chains = nc, n.adapt = 1000)

update(Nmix_fit, n.iter = 1000)



zm = coda.samples(Nmix_fit, variable.names = c("lambda", "N"), n.iter = ni, n.thin = nt)


lam <- MCMCpstr(zm, "lambda")
# N <- MCMCpstr(zm, "N")
# N <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(N)))
# N$V1 <- as.Date(N$V1, origin = "1970-01-01")
# 

lam <- as.data.frame(cbind(as.Date(temps$dts[-nodata]), unlist(lam)))
lam$V1 <- as.Date(lam$V1, origin = "1970-01-01")


Nmix_fit_UI5 <- jagsUI::jags(data = jags_data, inits = jags_inits, parameters.to.save = parameters, model.file = "N-mixtureNBNZMS.jags",  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)

# cor.df <- left_join(lam, means.list.NZMS, by=c('V1'="V2"), copy = T)
# #cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2.x)
# nblam <- cor.test((cor.df$V2), (cor.df$mean.abund), method = "spearman")
# 
# 
# cor.df <- left_join(N, means.list.NZMS, by=c('V1'="V2"), copy = T)
# #cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2.x)
# nbN <- cor.test((cor.df$V2), (cor.df$mean.abund), method = "spearman")

cor.df <- na.omit(left_join(lam, means.list.NZMS, by=c('V1'="temps$dts"), copy = T))

cor.df1 <- cor.df %>% slice(which(row_number() %% 2 == 0))
rho1 <- cor.test((cor.df1$V2), (cor.df1$mean.abund), method = "spearman")

cor.df2 <- cor.df %>%  slice(which(row_number() %% 2 == 1))
rho2 <- cor.test(cor.df2$V2, cor.df2$mean.abund, method = "spearman")


cor.df3 <- cor.df %>%  slice(which(row_number() %% 3 == 2))
rho3 <- cor.test(cor.df3$V2, cor.df3$mean.abund, method = "spearman")

nblam  <- mean(c(rho1$estimate, rho2$estimate, rho3$estimate))
#cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2)
print(Nmix_fit_UI5)


rm(zm)
rm(Nmix_fit)
#rm(Nmix_fit_UI5)

