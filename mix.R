install.packages("unmarked")
library(unmarked)
# Real data
data(mallard)
mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site,
                                  obsCovs = mallard.obs)
(fm.mallard <- pcount(~ ivel+ date + I(date^2) ~ length + elev + forest, mallardUMF, K=30))
(fm.mallard.nb <- pcount(~ date + I(date^2) ~ length + elev, mixture = "NB", mallardUMF, K=30))

NZMS.samp

forts <- vector()
for (i in 1:length(temps$dts)){
  d <- NZMS.samp[which(NZMS.samp$Date >= temps$dts[i] & NZMS.samp$Date < temps$dts[i+1]),]
  if (length(d$Date) == 0) {
    s = NA
  } else {
    s<- rep(i, length(d$Date))
  }
  forts <- append(forts, s)
}

NZMS.samp <- as.data.frame(cbind(NZMS.samp, na.omit(forts)))

View(aggregate(NZMS.samp$Density ~ NZMS.samp$Date, FUN = length))

install.packages("fitdistrplus")
library(fitdistrplus)

NZMS.samp <- group_by(NZMS.samp, na.omit(forts))
NZMS.sp <- split(NZMS.samp, f = na.omit(forts))
for i in seq(1:142){
  b <- NZMS.sp[[i]]$Density
  mean(b)
  plot <- distplot(table(b), "nbin")
  fit <- glm.nb(Counts ~ ., data = plot)
  sfit <- summary(fit)
  
}
hist(NZMS.sp$`12`$Density)
mean(NZMS.sp$`12`$Density)
plot <- distplot(table(NZMS.sp$`12`$Density), "nbin")
fit <- glm.nb(Counts ~ ., data = plot)
summary(fit)

rnbinom(1, size = 0.667, mu = 15093.47)


ggplot(data = NZMS.samp, aes(Date, fill = as.factor(RiverMile)))+
  geom_histogram(stackgroups = TRUE, binwidth = 14, dotsize = 1.25)





NZMS.samp 
site <- unique(NZMS.samp$RiverMile)
uniquebarcode<- unique(NZMS.samp$BarcodeID)
mat <- matrix(nrow = length(unique(NZMS.samp$RiverMile)), ncol = length(unique(NZMS.samp$BarcodeID)))
obstemp <-  matrix(nrow = length(unique(NZMS.samp$RiverMile)), ncol = length(unique(NZMS.samp$BarcodeID)))
obsflow <-  matrix(nrow = length(unique(NZMS.samp$RiverMile)), ncol = length(unique(NZMS.samp$BarcodeID)))
for (i in 1:length(uniquebarcode)) {
  a <- NZMS.samp[which(NZMS.samp$BarcodeID == uniquebarcode[i]),]
  b <- which(site == a$RiverMile)
  c <- last(which(temps$dts <= a$Date))
  obstemp[ ,i] <- rep(temps$Temperature[c], 27)
  obsflow[ ,i] <- rep(flow.magnitude$Discharge[c], 27)
  mat[b, i] <- a$CountTotal
}

obs <- list(obstemp, obsflow)
names(obs) <- c('temp', 'flowmePCount(mat,obsCovs = obs)
fm <- pcount(~temp')
dat <- unmarkedFra + flow ~1, dat, K = 20353)
fm.nb <- pcount(~temp + flow ~1, dat,mixture = "NB", K = 20353) #
#AIC: 15566.42 
