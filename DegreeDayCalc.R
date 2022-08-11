##################################################
# Calculate Degree Days based on Temp
#######################################
# calculate degree days
# that is basically  degree day = Sum( mean temp - threshold temp) over time period
# we can calculate threshold temp later (based on info regarding Cloean dipterum, 10 C)


degreeday <- aggregate(temp[sapply(temp,is.numeric)],
                       by=list(ID),
                       FUN=sum)
names(degreeday)[1:2] <-c("dts","DegreeDay")
# add the correct dates as the beginning of every period
# 
degreeday$dts <- as.POSIXct(temp$Date[(degreeday$dts*14)+1])
# order by date in chronological order
degreeday <- degreeday[order(degreeday$dts),]
degreeday <- degreeday[1:363,]

# there are less temperatures than discharge readings, so we will just repeat the same temp sequence for this exercise
DDs <- rbind(degreeday, degreeday, degreeday)
DDs <- DDs[1:length(out$Discharge), ]






emergetime <- vector()

for (t in timestep){
  if (t == 1) {print("t = 1")}
  else {
  degseq <- seq(t-1, 1, by = -1)
  #tdeg <- DDs$DegreeDay[degseq]
  vec <- 0
  for (s in degseq) {
    if(vec <= 266) { vec <- DDs$DegreeDay[s] + vec
    print(vec)}
    else {emerg <- t - s
    emergetime <- append(emergetime, emerg)
      break}
  
    }
    
  }
}



plot(timestep, temps$Temperature[1:26], type = "l", col = "red", xlab = "Time", ylim = c(5,15))
lines(timestep[7:26], emergetime, col = "blue")
legend(5, 14, legend=c("Time to Emerge", "Temperature (C)"),
       col=c("Blue", "REd"), lty=1, cex=0.8)
temps$Temperature
