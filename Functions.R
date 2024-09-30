######################################
## Colorado River Insects Functions
## Angelika Kurthen
######################################

# want to make a function that creates dataframes with river mile info, date, species sizes
TaxaLocationData <- function(spp, loc){
  
  # pull species specific data from the drift data
  df <- sampspec(samp = drift, species = spp, stats = F)
  df <- cbind(df[[1]][5], df[[1]][6], df$Specimen)
  
  # if the location is the upper canyon, want df to be RM -14 through 75
  if (loc == "upper"){ 
    df <- df[df$RiverMile < 75, ]
  } 
  # if the location is the middle canyon want df to be RM 75 through 150
  if (loc == "middle"){
    df <- df[df$RiverMile >= 75 & df$RiverMile < 150,]
  }
  # if the location is the lower canyon want the df to be all RMs after 150
  if (loc == "lower"){
    df <- df[df$RiverMile >= 150, ]
  }
  return(df)
}



SizeFreqPlotDaily <- function(df){
  # convert Dates into ymd date format
  df$Date<-ymd(df$Date)
  # all unique dates that samples were taken
  Dates<-unique(df$Date)
  nDates<-length(Dates)
  if (nDates > 200){ # if the number of dates is greater than 200, can't show data nicely 
    print("Too Many Dates - aggregate by week")
    df2 <- df %>% # aggregate by week
      group_by(year = year(Date), week = week(Date)) %>% 
      summarise_if(is.numeric, sum)
    df2$Date <- paste0(df2$year, "-", df2$week, "-", "1") # convert into year-week-weekday format
    df2$Date <- as.Date(df2$Date, format = "%Y-%U-%u") # convert to a Date (basically, first day of the week)
    
    Dates<-unique(df$Date)# get the unique dates
    nDates<-length(Dates)
    
    if (length(df2$week > 200)) {
      print("Too Many Dates - aggregate by month")
      df3 <- df %>% # aggregate by week
        group_by(year = year(Date), month = month(Date)) %>% 
        summarise_if(is.numeric, sum)
      df3$Date <- paste0(df3$year, "-", df3$month)
      df3$Date <- as.yearmon(df3$Date)
      
      Dates<-unique(df3$Date)
      nDates<-length(Dates)
      
      
      par(oma=c(5,4,3,2.4),mar=c(0,0,0,0),yaxs='i', mfrow = c(nDates, 1))
      p=2
      mxFreqs<-c(50,250,530,675,1000) # see hist of mxFreqsObs below to set
      mxFreqsObs<-dim(0)
      
      nf<-layout(matrix(seq(1,nDates),1,nDates),widths=1)
      for(d in 1:length(Dates)){
        temp<-subset(df3,Date == Dates[d])
        # make sure that the column names are in numerical form
        colnames(temp)[4:24] <- c(0:20)
        tab<-as.matrix(temp[4:24])
        tab[tab==0]<-NA
        # 		none<-apply(tab,1,sum,na.rm=TRUE)
        #  		tab[none==0,]<-0
        maxFreq<-max(apply(tab,2,sum,na.rm=TRUE))
        mxFreqsObs<-c(mxFreqsObs,maxFreq)
        # 		mxFreq<-mxFreqs[mxFreqs>maxFreq][1] # uncomment to scale to common xlim values
        mxFreq<-maxFreq+1 # uncomment to scale to individual xlim values
        barplot(tab,axes=FALSE,xlim=c(0,mxFreq),ylim=c(0,10), horiz=TRUE,axisnames=FALSE,border=NA)
        text(mxFreq/2, 0, labels = Dates[d], srt = 45, adj = c(1.1,1.1), xpd = NA, cex=0.5)
        axis(3,at=c(mxFreq),label=c(mxFreq), cex.axis=.8,tcl=-0.2,mgp=c(0,0.3,0),las=2)
        if(d==1){
          axis(2,las=2,at=seq(0,20))
          mtext('Size (mm)',2,line=2.5)		}
        if(d==length(Dates)){axis(4,las=2,at=seq(0,20))	}
      }
      box("inner", col="black")
      #dev.off()
      return()
    }
    
    par(oma=c(5,4,3,2.4),mar=c(0,0,0,0),yaxs='i', mfrow = c(nDates, 1))
    p=2
    mxFreqs<-c(50,250,530,675,1000) # see hist of mxFreqsObs below to set
    mxFreqsObs<-dim(0)
    
    nf<-layout(matrix(seq(1,nDates),1,nDates),widths=1)
    for(d in 1:length(Dates)){
      temp<-subset(df2,Date==ymd(Dates[d]))
      # make sure that the column names are in numerical form
      colnames(temp)[4:24] <- c(0:20)
      tab<-as.matrix(temp[4:24])
      tab[tab==0]<-NA
      # 		none<-apply(tab,1,sum,na.rm=TRUE)
      #  		tab[none==0,]<-0
      maxFreq<-max(apply(tab,2,sum,na.rm=TRUE))
      mxFreqsObs<-c(mxFreqsObs,maxFreq)
      # 		mxFreq<-mxFreqs[mxFreqs>maxFreq][1] # uncomment to scale to common xlim values
      mxFreq<-maxFreq+1 # uncomment to scale to individual xlim values
      barplot(tab,axes=FALSE,xlim=c(0,mxFreq),ylim=c(0,10), horiz=TRUE,axisnames=FALSE,border=NA)
      text(mxFreq/2, 0, labels = Dates[d], srt = 45, adj = c(1.1,1.1), xpd = NA, cex=0.5)
      axis(3,at=c(mxFreq),label=c(mxFreq), cex.axis=.8,tcl=-0.2,mgp=c(0,0.3,0),las=2)
      if(d==1){
        axis(2,las=2,at=seq(0,20))
        mtext('Size (mm)',2,line=2.5)		}
      if(d==length(Dates)){axis(4,las=2,at=seq(0,20))	}
    }
    box("inner", col="black")
    #dev.off()
    return()
  } else {
    # mxSize<-max(Obs$PredSize)+1
    # pdf('ExpPatch-Output/ExpPatch_Nucella_Size_Nostrina.pdf',height=4,width=8)
    # all unique dates that samples were take
    par(oma=c(5,4,3,2.4),mar=c(0,0,0,0),yaxs='i', mfrow = c(nDates, 1))
    p=2
    mxFreqs<-c(50,250,530,675,1000) # see hist of mxFreqsObs below to set
    mxFreqsObs<-dim(0)
    
    nf<-layout(matrix(seq(1,nDates),1,nDates),widths=1)
    for(d in 1:length(Dates)){
      temp<-subset(df,Date==ymd(Dates[d]))
      # make sure that the column names are in numerical form
      colnames(temp)[5:25] <- c(0:20)
      tab<-as.matrix(temp[5:25])
      tab[tab==0]<-NA
      # 		none<-apply(tab,1,sum,na.rm=TRUE)
      #  		tab[none==0,]<-0
      maxFreq<-max(apply(tab,2,sum,na.rm=TRUE))
      mxFreqsObs<-c(mxFreqsObs,maxFreq)
      # 		mxFreq<-mxFreqs[mxFreqs>maxFreq][1] # uncomment to scale to common xlim values
      mxFreq<-maxFreq+1 # uncomment to scale to individual xlim values
      barplot(tab,axes=FALSE,xlim=c(0,mxFreq),ylim=c(0,10), horiz=TRUE,axisnames=FALSE,border=NA)
      text(mxFreq/2, 0, labels = Dates[d], srt = 45, adj = c(1.1,1.1), xpd = NA, cex=0.5)
      axis(3,at=c(mxFreq),label=c(mxFreq), cex.axis=.8,tcl=-0.2,mgp=c(0,0.3,0),las=2)
      if(d==1){
        axis(2,las=2,at=seq(0,20))
        mtext('Size (mm)',2,line=2.5)		}
      if(d==length(Dates)){axis(4,las=2,at=seq(0,20))	}
    }
    box("inner", col="black")
    #dev.off()
  }
}
