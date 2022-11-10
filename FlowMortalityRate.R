#####################################
## Flow Mortality Rate graph
###########################

# follows exp(-h*Q[t-1])
h =  0.0075
Q <- seq(from = 1, to = 400, by = 2)
mort <- 1- (exp(-h*Q))

df <- as.data.frame(cbind(Q, mort))
ggplot(df, aes(x = Q, y = mort))+
  geom_line()+
  coord_cartesian(ylim = c(0,1)) +
  ylab('Immediate Post-Disturbance Mortality (%)') +
  xlab('Sycamore Creek, AZ Discharge (cfs)')
