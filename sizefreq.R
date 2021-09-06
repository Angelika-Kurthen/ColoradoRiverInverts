install.packages("devtools")
library(devtools)
install_github(repo = "jmuehlbauer-usgs/R-packages", subdir = "foodbase")
library(foodbase)
library(tidyr)
library(ggplot2)
drift <- readDB(gear = "Drift", type = "Sample", updater = T)
library(gam)
library(lubridate)

species <- c("NZMS","GAMM", "QUAG", "CHIA", "CHIP", "CHIL", "SIMA", "SIMP", "SIML", "HYDE", "BAEL")
SizeFreqPlot <- function(spp, loc){
   
t1 <- readDB(type = 'SpeciesList')

# if the location is the upper canyon, want df to be RM -14 through 75
  if (loc = "upper"){ 
    df <- drift[drift$RiverMile <= 75, ]
    df <- sampspec(samp = df)
  } 
# if the location is the middle canyon want df to be RM 75 through 150
  if (loc = "middle"){
    df <- drift[drift$RiverMile > 75 & drift$RiverMile < 150,]
    df <- sampspec(samp = df)
    }
# if the location is the lower canyon want the df to be all RMs after 150
  if (loc = "lower"){
    df <- drift[drift$RiverMile >= 150, ]
    df <- sampspec(samp = df)
  }

# pull species specific data
df <- sampspec(samp = df, species = spp, stats = T)

# get species data
df_sp <-  df$SpecDel

# convert Dates into ymd date format
df[[8]]$Date<-ymd(df[[8]]$Date)
# bind spp info to date info
df_sp <- cbind(df[[8]]$Date, df_sp)


# get rid of NAs
df_sp <- df_sp[-which((is.na(df_sp$CountTotal == "NA"))== T),]
# change name of Date column to "Date"
colnames(df_sp)[1] <- "Date"
# order dataframe by ascending Date
df_sp <- df_sp[order(df_sp$Date), ]

# if info about stats wanted, uncomment below
#dfSizeStats<- cbind(df[[8]]$Date, df$Statistics[,2:7])

# format plot view
par(oma=c(5,4,3,2.4),mar=c(0,0,0,0),yaxs='i')
p=2
# isolate all unique dates so we can combine info
Dates<-unique(df_sp$Date)
# get length of all the unique dates
nDates<-length(Dates) # number of dates
# create empty list
mxFreqsObs<-dim(0)
# absolutely no idea what this is for 
nf<-layout(matrix(seq(1,nDates),1,nDates),widths=1)

# now we cycle through each date and make a histogram of size frequency 
for(d in 1:length(Dates)){
  # get size and spp info for each date
  temp<-subset(df_sp,Date==ymd(Dates[d]))
  # make sure that the column names are in numerica form
  colnames(temp)[4:24] <- c(0:20)
  # only want the size info
  temp <- temp[,4:24]
  # remove 0s, make them NAs
  temp[temp==0]<-NA
  # sum the frequency of size classes for that unique date
  maxFreq<-max(apply(temp,2,sum,na.rm=TRUE))
  # 
  mxFreqsObs<-c(mxFreqsObs,maxFreq)
  # use this to set xlims
  mxFreq<-maxFreq+1
  # convert to matrix
  temp <- as.matrix(temp)
  # make barplot
  barplot(temp, axes=FALSE,xlim=c(0,mxFreq),ylim=c(0,10), horiz=TRUE,axisnames=FALSE,border=NA)
  # add Dates to bottom
  text(mxFreq/2, 0, labels = Dates[d], srt = 45, adj = c(1.1,1.1), xpd = NA, cex=1)
  # format upper axes so it is labled with the frequency of observations
  axis(3,at=c(mxFreq),label=c(mxFreq), cex.axis=.8,tcl=-0.2,mgp=c(0,0.3,0),las=2)
  # for the first graph, make y axis label
  if(d==1){
    axis(2,las=2,at=seq(0,10))
    mtext('Size (mm)',2,line=2.5)
  }
  # for last graph, make y axis tick marks
  if(d==length(Dates)){axis(4,las=2,at=seq(0,10))
  }
}
}
