#'  Plotting functions for Fst distributions after OutFLANK
#'  
#'This function takes the output of OutFLANK as
#'input with the OFoutput parameter.  It plots a histogram of the FST (by
#'default, the uncorrected FSTs used by OutFLANK) of loci and overlays the
#'inferred null histogram.
#'
#'
#'@param OFoutput The output of the function OutFLANK() 
#
#' @param withOutliers Determines whether the loci marked as outliers (with $OutlierFlag) are included in the histogram.
#' 
#' @param NoCorr Plots the distribution of FSTNoCorr when TRUE. Recommended, because this is the data used by OutFLANK to infer the distribution.
#' 
#' @param Hmin The minimum heterozygosity required before including  a locus in the plot.
#' 
#' @param binwidth The width of bins in the histogram.
#' 
#' @param Zoom If Zoom is set to TRUE, then the graph will zoom in on the right tail of the distirbution (based on argument RightZoomFraction)
#' 
#' @param RightZoomFraction Used when Zoom = TRUE. Defines the proportion of the distribution to plot.
#' 
#' @param titletext Allows a test string to be printed as a title on the graph
#' 
#' @return produces a historgram of the FST
#' @export


util.outflank.plotter <- function(OFoutput,withOutliers = TRUE, NoCorr= TRUE, Hmin=0.1, binwidth=0.005, Zoom = FALSE,RightZoomFraction = 0.05,titletext=NULL){
  data=OFoutput$results[which(OFoutput$results$He>Hmin),]
  if(NoCorr) {
    flist=data$FSTNoCorr
    fbar=sum(data$T1NoCorr)/sum(data$T2NoCorr)
    titletext= paste(c(titletext,"Fst without sample size correction"))
  }
  
  if(!NoCorr) {
    flist=data$FST
    fbar=OFoutput$FSTbar
    
    titletext= paste(c(titletext,"Fst with sample size correction"))
  }
  
  flist = flist[which(!is.na(flist))]
  keeperlist=which(!data$OutlierFlag)
  
  
  if(!withOutliers) flist = flist[keeperlist]
  
  if(Zoom) {FstDistPlotterZoom(df = OFoutput$dfInferred, FSTlist  = flist,  FSTbar = fbar, binwidth,titletext, RightZoomFraction)} else {
    FstDistPlotter(df = OFoutput$dfInferred, FSTlist = flist,  FSTbar = fbar, binwidth, titletext = titletext)}
  
}


################################

FstDistPlotter = function(df, FSTlist, FSTbar, binwidth=0.005,titletext=NULL){
  xPlotUpperBound=ceiling(max(FSTlist)*100)/100
  breakslist=seq(0,xPlotUpperBound+binwidth,by=binwidth)
  breaks = length(breakslist)
  
  x = breakslist
  y=rep(0,length(x))
  for(i in 1:breaks) y[i] = pchisq(((i-.5)*binwidth)/FSTbar*df , df=df) - pchisq((((i-1.5)*binwidth))/FSTbar*df , df=df)
  y=length(FSTlist)*y
  
  hist(FSTlist,col="darkgoldenrod1", breaks=breakslist, prob=F, xlab="Fst",  main=titletext)

  lines(x,y,col="darkblue", lwd=3)
}

###  FstDistPlotterZoom  #####
#This is a function that plots the right tail of the distribution.  

FstDistPlotterZoom = function(df, FSTlist,  FSTbar, binwidth = 0.005,titletext = NULL, RightZoomFraction = 0.1){
  
  FSTlistNoNA=FSTlist[which(!is.na(FSTlist))]
  
  xPlotUpperBound=ceiling(max(FSTlistNoNA)*100)/100
  xPlotLowerBound=floor(as.numeric(quantile(FSTlistNoNA, prob = 1 - RightZoomFraction, na.rm=TRUE)) * 100) / 100
  flist=FSTlistNoNA[which(FSTlistNoNA>xPlotLowerBound)]
  
  
  breakslist=seq(xPlotLowerBound,xPlotUpperBound,by=binwidth)
  breaks = length(breakslist)
  
  x = breakslist
  
  y=rep(0,length(x))
  for(i in 1:breaks) y[i] = pchisq((xPlotLowerBound + (i-.5)*binwidth)/FSTbar*df , df=df) - pchisq((xPlotLowerBound+ (i-1.5)*binwidth)/FSTbar*df , df=df)
  
  y=length(FSTlistNoNA)*y
  
  hist(flist,col="darkgoldenrod1", breaks=breakslist, prob=F, xlab="Fst",  main=titletext)
 
  lines(x,y,col="darkblue", lwd=3)
 
}

FstDistPlotterAddBadCurve = function(df, FSTlist,  FSTbar, binwidth = 0.005, RightZoomFraction = 0.99){
  
  FSTlistNoNA=FSTlist[which(!is.na(FSTlist))]
  
  xPlotUpperBound=ceiling(max(FSTlistNoNA)*100)/100
  xPlotLowerBound=floor(as.numeric(quantile(FSTlistNoNA, prob = 1 - RightZoomFraction, na.rm=TRUE)) * 100) / 100
  
  
  breakslist=seq(xPlotLowerBound,xPlotUpperBound,by=binwidth)
  breaks = length(breakslist)
  
  x = breakslist
  
  y=rep(0,length(x))
  for(i in 1:breaks) y[i] = pchisq((xPlotLowerBound + (i-.5)*binwidth)/FSTbar*df , df=df) - pchisq((xPlotLowerBound+ (i-1.5)*binwidth)/FSTbar*df , df=df)
  
  y=length(FSTlistNoNA)*y
  
  
  lines(x,y,col="red", lwd=3)
  
}


#  OutFLANKBadCurvePlotter draws a curve based on the same Fstbar but with soe differnet degrees of freedom
OutFLANKBadCurvePlotter = function(badDF,OFoutput,withOutliers = TRUE, NoCorr= TRUE, Hmin=0.1, binwidth=0.005, Zoom = FALSE,RightZoomFraction = 0.99,titletext=NULL){
  data=OFoutput$results[which(OFoutput$results$He>Hmin),]
  if(NoCorr) {
    flist=data$FSTNoCorr
    fbar=sum(data$T1NoCorr)/sum(data$T2NoCorr)
  }
  
  if(!NoCorr) {
    flist=data$FST
    fbar=OFoutput$FSTbar
    
  }
  
  flist = flist[which(!is.na(flist))]
  keeperlist=which(!data$OutlierFlag)
  
  
  if(!withOutliers) flist = flist[keeperlist]
  
  FstDistPlotterAddBadCurve(badDF, FSTlist  = flist,  FSTbar = fbar, binwidth,RightZoomFraction)
  
}

