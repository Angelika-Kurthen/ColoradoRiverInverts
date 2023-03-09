





G1 = 0.9 * 1/14
G2 = 0.9 * 1/7
P1 = 0.9 * 6/7
P2 = 0.9 * 13/14
P3 = 0.9 * 6/7 



F2 = seq(0, 3, by = 1)
F3 <- seq(0, 100, by = 0.1)
time <- seq(1, 1000, by = 1)

for (y in 1:length(F2)){
  N <- array(data = 0, dim = c(length(time) + 1, 3))
  N[1, ] <- c(10, 10, 10)
  out.df <- matrix(NA, ncol = 2, nrow = 0)
for (i in 1:length(F3)){
  S1 <- c(P1, F2[y], F3[i])
  S2 <- c(G1, P2, 0)
  S3 <- c(0, G2, P3) 
  A <-  rbind( S1, S2, S3)
  Total.N <- vector()
  
  for (t in 2:length(time)){
    N[t,] <- N[t-1, ] %*% A
    Total.N <- append(Total.N, sum(N[t,]))
    uval <- unique(Total.N[940:1000])
  }
  out.df <- rbind(out.df, cbind(rep(F3[i], length(uval)), uval))
}
out.df <- as.data.frame(out.df)
colnames(out.df) <- c("F3", "N")
p <- ggplot(out.df, aes(x = F3, y = N)) + 
  geom_point(size = 0.5) + 
  ylab('log N')+
  ggtitle(paste0("F2=", F2[y]))  
png(filename = paste0(F2[y], "plot.png"))
print(p)
dev.off()
  }
