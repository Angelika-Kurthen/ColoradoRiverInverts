##################
# NZMS N Mixture
#################
#library(dataRetrieval, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
library(unmarked, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")


# HPC and dataRetrieval don't like each other

#discharge <- readNWISdv("09380000", "00060", "2007-10-01", "2023-05-01")
#flow.magnitude <- TimestepDischarge(discharge, 85000)

write.csv(flow.magnitude, file = "flow.magnitude.csv")
flow.magnitude <- read.csv(file = "flow.magnitude.csv", header = T)

# temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
# temps <- TimestepTemperature(temp)

write.csv(temps, file = "temps.csv")
temps <- read.csv(file = "temps.csv", header = T)


NZMS.samp <- read.csv("NZMSsamp.csv", header = T)

# need to fit our data to an unmarked frame

# observations need to be in MxJ matrix, where M is # of sites and J is max number of obs per site

# we measured over temps$dts (those are our timesteps)

max_visits <- vector()
means <- vector()
se <- vector()
for (i in 1:length(temps$dts)){
  d <- NZMS.samp[which(NZMS.samp$Date >= temps$dts[i] & NZMS.samp$Date < temps$dts[i+1]), ]
  max_visits[i] <- length(d$Date)
  means[i] <- mean(d$Density)
  se[i] <- sd(d$Density)/(length(d$Density))
}

df <- cbind.data.frame(temps$dts, means, se)
df <- df[!is.na(df$means), ]

R <- length(temps$dts)
J <- max(max_visits)

site_mat <- matrix(data = NA, nrow = R, ncol = J)
dens_mat <- matrix(data = NA, nrow = R, ncol = J)
flows <- vector()
volumes <- matrix(data = NA, nrow = R, ncol = J)
for (i in 1:length(temps$dts)){
  d <- NZMS.samp[which(NZMS.samp$Date >= temps$dts[i] & NZMS.samp$Date < temps$dts[i+1]), ]
  flows[i] <- mean(d$X_00060_00003)
  site_mat[i, ] <- c(d$CountTotal, rep(NA, times = (J- length(d$CountTotal))))
  dens_mat[i, ] <- c(as.integer(d$Density), rep (NA, times = (J - length(d$Density))))
  volumes[i, ] <- c(scale(d$Volume),rep(NA, times = (J- length(d$CountTotal))))
}

# we need to remove all timesteps that are just NAs
# first identify all the timsteps that don't have data (so we can match them up later)
nodata <- which(is.na(site_mat[,1]))
length(temps$dts) - length(nodata)
site_mat <- as.matrix(site_mat[-nodata,])
dens_mat <- as.matrix(dens_mat[-nodata, ])
flows <- as.data.frame(flows[-nodata])
volumes <- as.matrix(volumes[-nodata, ])
dimnames(volumes) <- list(temps$dts[-nodata], seq(1:48))


volumes <- list(volumes)
names(volumes) <- c("vol")

dat <- unmarkedFramePCount(y = site_mat, siteCovs = flows, obsCovs = volumes)
dat_dens <- unmarkedFramePCount(y = dens_mat, siteCovs = flows, obsCovs = volumes)

# Model with no covariates or random effect on abundance
# mod_null <- stan_pcount(~ 1 ~1, dat, K=10453,
#                         chains=3, iter=2000, seed=123, prior_coef_det = logistic(0.69, 0.97), log_lik = F)
# mod_vol <- stan_pcount(~scale(vol)~1, dat, K = 10453, chains = 3, iter = 2000, seed= 123,prior_coef_det = logistic(0, 1), log_lik= F)

mod_um_vol <- pcount(~ vol ~1 , dat, K = 10453)
mod_off_vol <- pcount(~ offset(vol) ~ 1, dat, K = 10453)

mod_um_vol_NB <- pcount(~ 1 + offset(vol) ~1 , dat, K = 10453, mixture = "NB")
mod_um_vol_ZIP <- pcount(~ 1 + offset(vol) ~1 , dat, K = 10453, mixture = "ZIP")
mod_um_dens <- pcount(~1 ~1, dat_dens, K = 1176244)
mod_um_dens_NB <- pcount(~1 ~1, dat_dens, K = 1176244, mixture = "NB")
mod_um_dens_ZIP <- pcount(~1 ~1, dat_dens, K = 1176244, mixture = "ZIP")


#dens_mod_mull <- stan_pcount(~1 ~1, dat_dens, K = 1176244, chains = 3, iter = 2000, seed = 123, log_lik = T)
fm <- fitList(mod_um_dens, mod_um_dens_NB)
modSel(fm)
# Model selectionfit
# Function returning three fit-statistics.
fitstats <- function(fm, na.rm=TRUE) {
  observed <- getY(fm@data)
  expected <- fitted(fm)
  resids <- residuals(fm)
  sse <- sum(resids^2, na.rm=na.rm)
  chisq <- sum((observed - expected)^2 / expected, na.rm=na.rm)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm=na.rm)
  out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
  return(out)
}
pb <- parboot(mod_um_vol, fitstats, nsim=25, report=1)
pb_dens <- parboot(mod_um_dens, fitstats, nsim = 25, report = 1)
pb_plot_mod_vol <- plot(pb, main="")
pb_plot_mod_dens <- plot(pb_dens, main = "")
png(filename = "pb_plot_mod_vol.png")
plot(pb_plot_mod_vol)
dev.off()
png(filename = "pb_plot_mod_dens.png")
plot(pb_plot_mod_dens)
dev.off()
# Finite-sample inference for a derived parameter.
# Population size in sampled area
Nhat <- function(fm) {
  sum(bup(ranef(fm, K=50)))
}
set.seed(345)



pb_dens.N <- parboot(mod_um_vol_NB, fitstats, nsim=25, report=5)



pb_vol.N <- parboot(mod_um_NB, Nhat, nsim = 25, report = 5)
# Compare to empirical Bayes confidence intervals
colSums(confint(ranef(mod_um_vol_NB, K=10453)))
colSums(confint(ranef(mod_um_dens, K = 50)))

# Estimates of conditional abundance distribution at each site
re <- ranef(mod_um_dens_NB)
re_dens <- ranef(mod_um_vol_NB)
# Best Unbiased Predictors
bup(re, stat="mean") # Posterior mean
bup(re, stat="mode") # Posterior mode
est <- bup(re) # estimated population size

est <- est[!is.nan(est)]
obs <- df$means[!is.nan(est)]
cor(est, obs)

mutate(df$means, dens =  )


confint(re, level=0.95) # 90% CI

# Plots
plot(re)

# Best Unbiased Predictors
bup(re_dens, stat="mean") # Posterior mean
bup(re_dens, stat="mode") # Posterior mode
bup(re_dens) # estimated population size
confint(re_dens, level=0.9) # 90% CI
# Plots
plot(re, subset=site %in% 1:25)
plot(re, subset=site %in% 26:50)
plot(re, subset= site %in% 51:100)
plot(re, subset = site %in% 101:142)
     
#dens_mod_mull <- stan_pcount(~1 ~1, dat_dens, K = 1176244, chains = 3, iter = 2000, seed = 123, log_lik = T)

# occupancy to get global p 
# site_mat[site_mat >0] <- 1
# occu_dat <- unmarkedFramePCount(y = site_mat, siteCovs = flows, obsCovs = volumes)
# 
# occu_mod_null <- stan_occu(~1 ~1, occu_dat, chains = 3, iter = 2000)
# occu_mod_vol <- stan_occu(~scale(vol)~1, occu_dat, chains = 3, iter = 2000)
# 
# #posterior
# names(mod_um)
# ## [1] "beta_state[(Intercept)]" "beta_det[(Intercept)]"
# occ_intercept <- extract(mod_um, "det")[[1]]
# hist(occ_intercept, freq=FALSE)
# lines(density(occ_intercept), col='red', lwd=2)
# 
# #Compare the models
# #First we combine the models into a fitList:
# mods <- fitList(mod_um, mod_um_vol)
# #Then we generate a model selection table:
# round(modSel(mods), 3)
# #the model with the largest elpd performed best
# plot_residuals(mod_um_vol, submodel="det")
# fit_top_gof <- gof(mod_um_vol, raws=100, quiet=TRUE)
# 
# prob_det <- predict(mod_um_vol, submodel="det")
# # from predicted detection probabilities, get range
# range(prob_det$Predicted, na.rm = T)
# 
# 
