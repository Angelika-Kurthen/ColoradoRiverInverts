
source("NZMS_1sp_Model.R")
library(plotly)
# 3d plot comparing different steady temps * Q for NZMS

flow.seq = seq(0, 1, by = 0.1)
temp.seq = seq(0, 30, by = 1)

outarray <- array(data = NA, dim = c(length(temp.seq), length(flow.seq)), dimnames = list(temp.seq, flow.seq))

for (flow in 1:length(flow.seq)){
  discharge <- rep(flow.seq[flow], times = 1000)
  for (temp in 1:length(temp.seq)){
    temps <- rep(temp.seq[temp], times = 1000)
    dts <- seq(1, 1000, by = 1)
    temps <- as.data.frame(cbind(dts, temps))
    colnames(temps) <- c("dts", "Temperature")
    Nout <- NZMSmodel(flow.data = discharge, temp.data = temps, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 500, iteration = 1, peaklist =  0, peakeach = length(temps$Temperature))
    outarray[temp, flow] <- mean(rowSums(Nout[200:300, ,]))
  }
}

View(as.list(outarray))
axx <- list(title = "Q")
axy <- list( title = "Temperature (C)")
axz <- list( title = "NZMS Abundace/Recruitment Limit")
p <- plot_ly(z = ~outarray/10000, x =flow.seq, type = "surface", contours = list(
  
  x = list(show = TRUE,start = 0, end = 1, size = 0.05, color = 'white', width = 1),
  
  z = list(show = TRUE,start = 0, end = 4, size = 0.25, color = 'white', width = 1)))
p <- p %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
hide_colorbar(p)
write.csv(outarray, "NZMSTvQ.csv", fileEncoding = "UTF-8")
