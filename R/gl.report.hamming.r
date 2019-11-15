#' Calculates the pairwise Hamming distance between DArT trimmed DNA sequences
#'
#' Hamming distance is calculated as the number of base differences between two 
#' sequences which can be expressed as a count or a proportion. Typically, it is
#' calculated between two sequences of equal length. In the context of DArT
#' trimmed sequences, which differ in length but which are anchored to the left
#' by the restriction enzyme recognition sequence, it is sensible to compare the
#' two trimmed sequences starting from immediately after the common recognition
#' sequence and terminating at the last base of the shorter sequence. 
#' 
#' Hamming distance can be computed by exploiting the fact that the dot product 
#' of two binary vectors x and (1-y) counts the corresponding elements that are 
#' different between x and y. This approach can also be used for vectors that 
#' contain more than two possible values at each position (e.g. A, C, T or G).
#'
#' If a pair of DNA sequences are of differing length, the longer is truncated.
#'
#' The algorithm is that of Johann de Jong 
#' \url{https://johanndejong.wordpress.com/2015/10/02/faster-hamming-distance-in-r-2/}
#' as implimented in utils.hamming.r
#' 
#' A histogram and whiskerplot can be requested. Both display a user specified value
#' for the minumum acceptable Hamming distance.
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param rs -- number of bases in the restriction enzyme recognition sequence [default 5]
#' @param boxplot -- if 'standard', plots a standard box and whisker plot; 
#' if 'adjusted', plots a boxplot adjusted for skewed distributions [default 'adjusted']
#' @param range -- specifies the range for delimiting outliers [default = 1.5 interquartile ranges]
#' @param threshold minimum acceptable base pair difference for display on the whisker plot and histogram [default 3 bp]
#' @param taglength -- typical length of the sequence tags [default 69]
#' @param probar -- if TRUE, then a progress bar is desplayed on long loops [default TRUE]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return Tabulation of loc that will be lost on filtering, against values of the threshold
#' @importFrom adegenet glPlot
#' @importFrom graphics hist
#' @importFrom stats sd
#' @importFrom robustbase adjbox
#' @importFrom graphics par
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' out <- gl.report.hamming(testset.gl)

gl.report.hamming <- function(x, 
                              rs=5, 
                              boxplot="adjusted",
                              range=1.5,
                              threshold=3,
                              taglength=69,
                              probar=FALSE, 
                              verbose = 2) {
  
# TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  
# SET VERBOSITY
  
  if (is.null(verbose)){ 
    if(!is.null(x@other$verbose)){ 
      verbose <- x@other$verbose
    } else { 
      verbose <- 2
    }
  } 
  
  if (verbose < 0 | verbose > 5){
    cat(paste("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n"))
    verbose <- 2
  }
  
# FLAG SCRIPT START
  
  if (verbose >= 1){
    if(verbose==5){
      cat("Starting",funname,"[ Build =",build,"]\n")
    } else {
      cat("Starting",funname,"\n")
    }
  }
  
# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    stop("Fatal Error: genlight object required!\n")
  }
  
  if (all(x@ploidy == 1)){
    cat("  Processing Presence/Absence (SilicoDArT) data\n")
  } else if (all(x@ploidy == 2)){
    cat("  Processing a SNP dataset\n")
  } else {
    stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)!")
  }
  
# FUNCTION SPECIFIC ERROR CHECKING

  if(length(x@other$loc.metrics$TrimmedSequence) == 0) {
    stop("Fatal Error: Data must include Trimmed Sequences\n")
  }

  if (rs < 0 | rs > taglength){
    cat("Fatal Error: Length of restriction enzyme recognition sequence must be greater than zero, and less that the maximum length of a sequence tag; usually it is less than 9\n"); stop()
  }
    
  if (nLoc(x) == 1){
    stop("Fatal Error: Data must include more than one locus\n")
  }

# DO THE JOB

  s <- as.character(x@other$loc.metrics$TrimmedSequence)
  tld <- threshold/(taglength-rs)
  if( probar ) {
    pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
    getTxtProgressBar(pb)
  }

  if (verbose >= 2){
    cat("  Hamming distance ranges from zero (sequence identity) to 1 (no bases shared at any position)\n")
    cat("  Calculating pairwise Hamming distances between trimmed Reference sequence tags\n")
  }
  
  count=0
  nL <- nLoc(x)
  d <- rep(NA,(((nL-1)*nL)/2))
  
  for (i in 1:(nL-1)){
    for (j in ((i+1):nL)){
      count <- count + 1
      d[count] <- utils.hamming(s[i],s[j],r=rs)
    }
    if (probar) {setTxtProgressBar(pb, i/(nL-1))}
  }

  # Prepare for plotting
  if (all(x@ploidy==2)){
    title <- paste0("SNP data (DArTSeq)\nPairwise Hamming Distance between sequence tags")
  } else {
    title <- paste0("Fragment P/A data (SilicoDArT)\nPairwise Hamming Distance between sequence tags")
  }  
  if (verbose >= 2) {
    if (boxplot == "adjusted"){
      cat("  Plotting box and whisker plot (adjusted for skewness) and histogram of Hamming distance, showing a threshold of",threshold,"bp [HD",round(tld,2),"]\n")
    } else {
      cat("  Plotting box and whisker plot (standard) and histogram of Hamming distance, showing a threshold of",threshold,"bp [HD",round(tld,2),"]\n")
    }  
  }
  # Save the prior settings for mfrow, oma, mai and pty, and reassign
  op <- par(mfrow = c(2, 1), oma=c(1,1,1,1), mai=c(0.5,0.5,0.5,0.5),pty="m")
  # Set margins for first plot
  par(mai=c(1,0.5,0.5,0.5))
  # Plot Box-Whisker plot
  if (boxplot == "standard"){
    whisker <- boxplot(d, 
                       horizontal=TRUE, 
                       col='red', 
                       range=range,
                       ylim=c(0,1),
                       main = title)
    abline(v=tld,col="red")
   } else {
    whisker <- robustbase::adjbox(x=as.numeric(d),
                                  horizontal = TRUE,
                                  col='red',
                                  range=range,
                                  ylim=c(0,1),
                                  main = title)
    abline(v=tld,col="red")
  }  
  # Set margins for second plot
  par(mai=c(0.5,0.5,0,0.5))  
  # Plot Histogram
  hst <- hist(d,
       main="",
       xlab="",
       border="blue", 
       col="red",
       xlim=c(0,1),
       breaks=100)
  abline(v=tld,col="red")
  text(tld,max(hst$counts),pos=4,offset=1,paste(threshold,"bp"))
  
   cat("  No. of loci =", nLoc(x), "\n")
   cat("  No. of individuals =", nInd(x), "\n")
   cat("    Miniumum Hamming distance: ",round(min(d),2),"\n")
   cat("    Maximum Hamming distance: ",round(max(d),2),"\n")
   cat(paste0("    Mean Hamming Distance ",round(mean(d),2),"+/-",round(sd(d),3)," SD\n\n"))
   n.outliers <- sum(d[d<=(threshold/(taglength-rs))])
   cat("  No. of pairs with Hamming Distance less than or equal to",threshold,"base pairs:",n.outliers,"\n")

   # Determine the loss of loci for a given filter cut-off
   retained <- array(NA,21)
   pc.retained <- array(NA,21)
   filtered <- array(NA,21)
   pc.filtered <- array(NA,21)
   percentile <- array(NA,21)
   for (index in 1:21) {
     i <- (index-1)*5
     percentile[index] <- i/100
     retained[index] <- length(d[d >= percentile[index]])
     pc.retained[index] <- round(retained[index]*100/((((nL-1)*nL)/2)),1)
     filtered[index] <- ((((nL-1)*nL)/2)) - retained[index]
     pc.filtered[index] <- 100 - pc.retained[index]
   }
   df <- cbind(percentile,retained,pc.retained,filtered,pc.filtered)
   df <- data.frame(df)
   colnames(df) <- c("Threshold", "Retained", "Percent", "Filtered", "Percent")
   df <- df[order(-df$Threshold),]
   rownames(df) <- NULL
#   if (verbose >= 3){
#     cat("Note: The data below are for pairwise distances between",nL,"loci, for which there are",((((nL-1)*nL)/2)), "distances\n")
#     print(df)
#   } 
   
   # Reset the par options    
   par(op)   
   
# FLAG SCRIPT END

   if (verbose >= 1) {
     cat("Completed:",funname,"\n")
   }
   
   return(df)

}
