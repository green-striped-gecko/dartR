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
#' Hamming distance can be computed 
#' by exploiting the fact that the dot product of two binary vectors x and (1-y)
#' counts the corresponding elements that are different between x and y.
#' This approach can also be used for vectors that contain more than two possible 
#' values at each position (e.g. A, C, T or G).
#'
#' If a pair of DNA sequences are of differing length, the longer is truncated.
#'
#' The algorithm is that of Johann de Jong 
#' \url{https://johanndejong.wordpress.com/2015/10/02/faster-hamming-distance-in-r-2/}
#' as implimented in utils.hamming.r
#' 
#' A histogram and or a smearplot can be requested. Note that the smearplot is computationally intensive, and will take time to 
#' execute on large datasets.
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param rs -- number of bases in the restriction enzyme recognition sequence [default = 5]
#' @param plot specify if a histogram of Hamming distance is to be produced [default FALSE] 
#' @param smearplot if TRUE, will produce a smearplot of individuals against loci [default FALSE]
#' @param probar -- if TRUE, then a progress bar is desplayed on long loops [default = TRUE]
#' @return Tabulation of loc that will be lost on filtering, against values of the threshold
#' @importFrom adegenet glPlot
#' @importFrom graphics hist
#' @importFrom stats sd
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl.report.hamming(testset.gl)

gl.report.hamming <- function(x, rs=5, plot=FALSE, smearplot=FALSE, probar=TRUE) {
  
# TIDY UP FILE SPECS

  funname <- match.call()[[1]]

# FLAG SCRIPT START

    cat("Starting",funname,"\n")

# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }

  # Work around a bug in adegenet if genlight object is created by subsetting
    x@other$loc.metrics <- x@other$loc.metrics[1:nLoc(x),]

  # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nLoc(x))
      pop(x) <- as.factor(pop(x))
    }

  # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x, verbose=0)
    if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {cat("  Warning: genlight object contains monomorphic loci\n")}

# FUNCTION SPECIFIC ERROR CHECKING

  if(length(x@other$loc.metrics$TrimmedSequence) == 0) {
    cat("Fatal Error: Data must include Trimmed Sequences\n"); stop()
  }

  if (rs < 0 | rs > 69){
    cat("Fatal Error: Length of restriction enzyme recognition sequence must be greater than zero, and less that the maximum length of a sequence tag; usually it is less than 9\n"); stop()
  }

# DO THE JOB

  s <- as.character(x@other$loc.metrics$TrimmedSequence)

  if (nLoc(x)==1) {
    cat("Only one locus. No Hamming distances calculated.\n")
  } else {
    cat("Hamming distance ranges from zero (sequence identity) to 1 (no bases shared at any position)\n")
    cat("Calculating pairwise Hamming distances between trimmed reference sequence tags\n")
  count=0
  nL <- nLoc(x)
  
  # Calculate the number of iterations in loops below to set dimensions of d
  # niter = sum i=1 to nL-1 of (nL - i)
  # niter = sum i=1 to nL-1 of (nL) - sum i=1 to nL-1 of (i)
  # niter = nL(nL-1) - (nL-1)nL/2 [triangle number]
  # niter = nL(nL-1)/2 which seem intuitive
  d <- rep(NA,(((nL-1)*nL)/2))
  
  if( probar ) {
    pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
    getTxtProgressBar(pb)
  }
  for (i in 1:(nL-1)){
    for (j in ((i+1):nL)){
      count <- count + 1
      d[count] <- utils.hamming(s[i],s[j],r=rs)
    }
    if (probar) {setTxtProgressBar(pb, i/(nL-1))}
  }
  }

  # Plot a histogram of Hamming distance
  par(mfrow = c(2, 1),pty="m")
  
  if (plot) {
    hist(d, 
         main="Hamming distance between Loci taken pairwise", 
         xlab="Hamming distance", 
         border="blue", 
         col="red",
         xlim=c(0,1),
         breaks=100)
  }  
  if(smearplot){
    glPlot(x)
  } 
  
   xlimit <- min(d)
   
   cat("No. of loci =", nLoc(x), "\n")
   cat("No. of individuals =", nInd(x), "\n")
   cat("  Miniumum Hamming distance: ",round(min(d),2),"\n")
   cat("  Maximum Hamming distance: ",round(max(d),2),"\n")
   #cat("  Mean Hamming distance: ",round(mean(d),3),"\n\n")
   cat(paste0("  Mean Hamming Distance ",mn,"+/-",sdev," SD\n\n"))

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
   cat("Note: The data below are calculated from pairwise distances between",nL,"loci, for which there are",((((nL-1)*nL)/2)), "distances\n")
   print(df)
   
# FLAG SCRIPT END

    cat("Completed:",funname,"\n")

   return(df)
}
