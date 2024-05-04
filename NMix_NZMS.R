##################
# NZMS N Mixture
#################
# data retrieval tool from USGS
#install.packages("dataRetrieval")
library(dataRetrieval)
#install.packages("devtools")
#install_github(repo = "jmuehlbauer-usgs/R-packages", subdir = "foodbase")
library(devtools)
library(foodbase)
library(lubridate)
library(AICcmodavg)
library(unmarked)
library(ubms)

source("1spFunctions.R")
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
temps <- TimestepTemperature(temp)
discharge <- readNWISdv("09380000", "00060", "2007-10-01", "2023-05-01")

# pull drift data from DB
drift.data.total <- readDB(gear = "Drift", type = "Sample", updater = TRUE)
# get drift data from between Lees Ferry and RM -6
drift.LF <- drift.data.total[which(drift.data.total$RiverMile >= -6 & drift.data.total$RiverMile <= 0),]
#specify NZMS
NZMS.samp.LF <- sampspec(samp = drift.LF, species = "NZMS", stats = T)
# pull stats and merge with sample info
NZMS.samp <- merge(NZMS.samp.LF$Statistics, NZMS.samp.LF$Samples, by = "BarcodeID", all = T)
# make sure we are using the same gear
NZMS.samp <- NZMS.samp[which(NZMS.samp$GearID == 4),] 
# calculate density
NZMS.samp$Density <- NZMS.samp$CountTotal/NZMS.samp$Volume
NZMS.samp <- merge(NZMS.samp, discharge[, 3:4], by = "Date")
# for densities, apply equation to adjust for flow
#NZMS.samp$Density <- (NZMS.samp$Density/(9e-15 *(NZMS.samp$X_00060_00003 * 0.02831683)^4.1))
#write.csv(NZMS.samp, "NZMSsamp.csv")

# need to fit our data to an unmarked frame

# observations need to be in RxJ matrix, where R is # of sites and J is max number of obs per site
# we measured over temps$dts (those are our timesteps)

max_visits <- vector() # vector to put max obs per site
means <- vector()
for (i in 1:length(temps$dts)){  # cycle through all the 14 day timesteps that we have model output for
  # pull abundances between each each timestep
  d <- NZMS.samp[which(NZMS.samp$Date >= temps$dts[i] & NZMS.samp$Date < temps$dts[i+1]), ]
  max_visits[i] <- length(d$Date) # number of observations
  means[i] <- mean(d$CountTotal) # means
}

# make data frame with the timestep, the mean
df <- cbind.data.frame(temps$dts, means)
# remove anythig where there are NA values of Density (means that volume or count is missing)
df <- df[!is.na(df$means), ]

# define our RxJ matrix
R <- length(temps$dts)
J <- max(max_visits)
site_mat <- matrix(data = NA, nrow = R, ncol = J)
dens_mat <- matrix(data = NA, nrow = R, ncol = J)

# make vector for flows at each timestep
# make RxJ matrix full of densities
# make RxJ matrix full of raw counts
# make RxJ matrix full of volumes sampled for each abundance
flows <- vector()
volumes <- matrix(data = NA, nrow = R, ncol = J)
for (i in 1:length(temps$dts)){
  d <- NZMS.samp[which(NZMS.samp$Date >= temps$dts[i] & NZMS.samp$Date < temps$dts[i+1]), ]
  flows[i] <- mean(d$X_00060_00003)
  site_mat[i, ] <- c(d$CountTotal, rep(NA, times = (J- length(d$CountTotal))))
  dens_mat[i, ] <- c(as.integer(d$Density), rep (NA, times = (J - length(d$Density))))
  volumes[i, ] <- c((d$Volume),rep(NA, times = (J- length(d$CountTotal))))
}

# we need to remove all timesteps that are just NAs
nodata <- which(is.na(site_mat[,1]))
# first identify all the timsteps that don't have data (so we can match them up later)
site_mat <- as.matrix(site_mat[-nodata,])
dens_mat <- as.matrix(dens_mat[-nodata, ])
flows <- as.data.frame(flows[-nodata])
volumes <- as.matrix(volumes[-nodata, ])
dimnames(volumes) <- list(temps$dts[-nodata], seq(1:48))
volumes <- list(volumes)
names(volumes) <- c("vol")

# convert organized data into unmarked dataframe format, with flow as site covariate and volume of sample as observation covariate
dat <- unmarkedFramePCount(y = site_mat, siteCovs = data.frame(flows), obsCovs = volumes)

# pcount is N-mixture model to estimate counts at each site
# we are doing space-time approximation

# we want adjust our count data by volume (since volumes vary by each sample)
# these are basically our three different hypotheses
#Poisson mixture
mod_um_vol_P <- pcount(~ 1 + (vol) ~1 , dat, K = 11000, mixture = "P")
# Neg Binom mixture
mod_um_vol_NB <- pcount(~ 1 + (vol) ~1 , dat, K = 11000, mixture = "NB")
# Zero Inflated Poisson mixture
mod_um_vol_ZIP <- pcount(~ 1 + (vol) ~1 , dat, K =11000, mixture = "ZIP")

# now we want to compare the models
fm <- fitList(mod_um_vol_P, mod_um_vol_NB, mod_um_vol_ZIP)
modSel(fm)

# Model selectionfit
# suggests that NB is best model

Nmix.gof.test(mod_um_vol_NB, nsim = 10)
# Compare to empirical Bayes confidence intervals
colSums(confint(ranef(mod_um_vol_NB, K=11000)), na.rm = T)
# Estimates of conditional abundance distribution at each site
re <- ranef(mod_um_vol_NB)
# Best Unbiased Predictors
est <- bup(re) # estimated population size at each timestep
conf <- as.data.frame(confint(re, level=0.95)) # 90% CI

###############################
# occupancy to get global p with ubms
###############################
# need 0s and 1s
site_mat[site_mat >0] <- 1

# turn into unmarked data frame
occu_dat <- unmarkedFramePCount(y = site_mat, siteCovs = flows, obsCovs = volumes)
# 
occu_mod_null <- stan_occu(~1 ~1, occu_dat, chains = 3, iter = 2000)
occu_mod_vol <- stan_occu(~scale(vol)~1, occu_dat, chains = 3, iter = 2000)
# 
# #posterior
names(occu_mod_vol)
# ## [1] "beta_state[(Intercept)]" "beta_det[(Intercept)]"
occ_intercept <- extract(occu_mod_vol, "beta_det")[[1]]
hist(occ_intercept, freq=FALSE)
lines(density(occ_intercept), col='red', lwd=2)
# 
# #Compare the models
# #First we combine the models into a fitList:
mods <- fitList(occu_mod_null, occu_mod_vol)
# #Then we generate a model selection table:
round(modSel(mods), 3)
# #the model with the largest elpd performed best
plot_residuals(occu_mod_vol, submodel="det")
fit_top_gof <- gof(occu_mod_vol, raws=100, quiet=TRUE)
# 
prob_det <- predict(occu_mod_vol, submodel="det")
# # from predicted detection probabilities, get range
range(prob_det$Predicted, na.rm = T)
# 
# Now what? 
