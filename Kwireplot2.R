################################
## Wireplot for Q, K, and tau
################################
library(tidyverse)
library(ggplot2)
library(lattice)
library(RColorBrewer)

## checkpos makes sure that the K-occupied term is positive, assigns 0 if not
checkpos <- function(x) {
  ifelse(x < 0, 0, x)
}

## Plot setup
clrs <- colorRampPalette(brewer.pal(9, "YlOrRd"))
trellis.par.set("axis.line", list(col = NA, lty = 1, lwd = 1))
theme.novpadding <- list(
  layout.heights = list(top.padding = 0, bottom.padding = 0),
  layout.widths = list(left.padding = 0, right.padding = 0)
)

## Minimum threshold of what is considered a flood
Qmin <- 50000 
## Half saturation constant
a <- 75000
## Rate that K returns to pre-disturbance level
g <- 0.1

## Maximum flood size to run model to
Qmax = 100000

## Timesteps to run model out to - days
t = 75



dfK <- data.frame(Q = seq(0, Qmax, by = 200))
for(i in 1:length(dfK$Q)){
if (dfK$Q[i] < Qmin) {
  dfK$Q[i] <- 0
} else {
  dfK$Q[i] <- (dfK$Q[i] - Qmin)/(a + dfK$Q[i]- Qmin)
}}

KQT <- data.frame(Q = numeric(((Qmax/200) + 1) * (t + 1)))
## Filling in Q
KQT$Q <- rep(seq(0, Qmax, by =200), each = t + 1)

for(i in 1:length(KQT$Q)){
  if (KQT$Q[i] < Qmin) {
    KQT$Qf[i] <- 0
  } else {
    KQT$Qf[i] <- (KQT$Q[i] - Qmin)/(a + KQT$Q[i]- Qmin)
  }}


for (i in 1:length(KQT$Q)){
# Calculate K arrying capacity immediately following the disturbance
  KQT$K0[i]  <- 10000 + ((40000-10000)*KQT$Qf[i])
  }

KQT$t <- rep(seq(0, t), (Qmax/200) + 1)

for (i in 1:length(KQT$Q)){
KQT$K[i] <- 10000 + (KQT$K0[i] - 10000) * exp(-g * KQT$t[i])}

wireframe(K ~ Q + t, data = KQT,
          aspect = c(1, .4),
          drape = TRUE,
          shade = FALSE,
          colorkey = FALSE,
          col = alpha('#ffeda0', 0.08),
          scales = list(arrows = FALSE, col = 'black'),
          screen = list(z = -40, x = -70),
          par.settings = theme.novpadding,
          col.regions = clrs(1000),
          main = 'Baetis spp. K')
