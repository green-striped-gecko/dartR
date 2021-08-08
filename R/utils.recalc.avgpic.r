#' A utility script to recalculate the OneRatioRef, OneRatioSnp, PICRef, PICSnp, and AvgPIC by locus after some
#' populations have been deleted.
#'
#' The locus metadata supplied by DArT has OneRatioRef, OneRatioSnp, PICRef, PICSnp, and AvgPIC included,
#' but the allelec composition will change when some individuals are removed from the dataset and so the initial statistics will
#' no longer apply. This script recalculates these statistics and places the recalculated values in the appropriate place 
#' in the genlight object.
#' 
#' If the locus metadata OneRatioRef|Snp, PICRef|Snp and/or AvgPIC do not exist, the script
#' creates and populates them.
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return The modified genlight object
#' @author Custodian: Luis Mijangos (Post to \url{https://groups.google.com/d/forum/dartr})
#' @seealso \code{utils.recalc.metrics} for recalculating all metrics, \code{utils.recalc.callrate} for recalculating CallRate,
#' \code{utils.recalc.freqhomref} for recalculating frequency of homozygous reference, \code{utils.recalc.freqhomsnp} for recalculating frequency of homozygous alternate,
#' \code{utils.recalc.freqhet} for recalculating frequency of heterozygotes, \code{gl.recalc.maf} for recalculating minor allele frequency,
#' \code{gl.recalc.rdepth} for recalculating average read depth
#' @examples
#' #out <- utils.recalc.avgpic(testset.gl)

utils.recalc.avgpic <- function(x, verbose=NULL) {

# TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  hold <- x
  
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
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }
  
  if (all(x@ploidy == 1)){
      if (verbose >= 2){cat("  Processing  Presence/Absence (SilicoDArT) data, locus metric PIC\n")}
      data.type <- "SilicoDArT"
    } else if (all(x@ploidy == 2)){
      if (verbose >= 2){cat("  Processing a SNP dataset, locus metric AvgPIC\n")}
      data.type <- "SNP"
    } else {
      stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
  }

  # Check monomorphs have been removed up to date
  if (x@other$loc.metrics.flags$monomorphs == FALSE){
    if (verbose >= 2){
      cat("  Warning: Dataset contains monomorphic loci which will be included in the ",funname," calculations\n")
    }  
  }

# FUNCTION SPECIFIC ERROR CHECKING 

  if (data.type=="SNP"){
     if (is.null(x@other$loc.metrics$AvgPIC)) {
       x@other$loc.metrics$AvgPIC <- array(NA,nLoc(x))
       if (verbose >= 3){
         cat("  Locus metric AvgPIC does not exist, creating slot @other$loc.metrics$AvgPIC\n")
       }
     }
     if (is.null(x@other$loc.metrics$OneRatioRef)) {
       x@other$loc.metrics$OneRatioRef <- array(NA,nLoc(x))
       if (verbose >= 3){
         cat("  Locus metric OneRatioRef does not exist, creating slot @other$loc.metrics$OneRatioRef\n")
       }
     }
     if (is.null(x@other$loc.metrics$OneRatioSnp)) {
       x@other$loc.metrics$OneRatioSnp <- array(NA,nLoc(x))
       if (verbose >= 3){
         cat("  Locus metric OneRatioSnp does not exist, creating slot @other$loc.metrics$OneRatioSnp\n")
       }
     }
     if (is.null(x@other$loc.metrics$PICRef)) {
       x@other$loc.metrics$PICRef <- array(NA,nLoc(x))
       if (verbose >= 3){
         cat("  Locus metric PICRef does not exist, creating slot @other$loc.metrics$PICRef\n")
       }
     }
     if (is.null(x@other$loc.metrics$PICSnp)) {
       x@other$loc.metrics$PICSnp <- array(NA,nLoc(x))
       if (verbose >= 3){
         cat("  Locus metric PICSnp does not exist, creating slot @other$loc.metrics$PICSnp\n")
       }
     }
  }
  
  if (data.type=="SilicoDArT"){
    if (is.null(x@other$loc.metrics$PIC)) {
      x@other$loc.metrics$PIC <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric PIC does not exist, creating slot @other$loc.metrics$PIC\n")
      }
    }
    if (is.null(x@other$loc.metrics$OneRatio)) {
      x@other$loc.metrics$OneRatio <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric OneRatio does not exist, creating slot @other$loc.metrics$OneRatio\n")
      }
    }
  }

# DO THE JOB

  t <- as.matrix(x)
     
  if(data.type == "SNP"){
    
    if (verbose >= 2){
      cat("  Recalculating OneRatioRef, OneRatioSnp, PICRef, PICSnp, AvgPIC\n")
    }
     
     c0 <- colSums(t==0, na.rm=T)
     c1 <- colSums(t==1, na.rm=T)
     c2 <- colSums(t==2, na.rm=T)
     c <- (c0+c1+c2)
     x@other$loc.metrics$OneRatioRef <- (c0+c1)/c
     x@other$loc.metrics.flags$OneRatioRef <- TRUE
     
     x@other$loc.metrics$OneRatioSnp <- (c1+c2)/c
     x@other$loc.metrics.flags$OneRatioSnp <- TRUE
     
     OneRatioRef <- x@other$loc.metrics$OneRatioRef
     OneRatioSnp <- x@other$loc.metrics$OneRatioSnp
     ZeroRatioRef <- 1 - OneRatioRef
     ZeroRatioSnp <- 1 - OneRatioSnp
     x@other$loc.metrics$PICRef <- 1 - ((OneRatioRef*OneRatioRef) + (ZeroRatioRef*ZeroRatioRef))
     x@other$loc.metrics.flags$PICRef <- TRUE
     
     x@other$loc.metrics$PICSnp <- 1 - ((OneRatioSnp*OneRatioSnp) + (ZeroRatioSnp*ZeroRatioSnp))
     x@other$loc.metrics.flags$PICSnp <- TRUE
     
     x@other$loc.metrics$AvgPIC <- (x@other$loc.metrics$PICRef + x@other$loc.metrics$PICSnp)/2
     x@other$loc.metrics.flags$AvgPIC <- TRUE
     
  } 
  
  if(data.type == "SilicoDArT"){
    
    if (verbose >= 2){
      cat("  Recalculating OneRatio, PIC\n")
    }
    
    OneRatio <- colMeans(t==1, na.rm=T)
    x@other$loc.metrics$OneRatio <- OneRatio
    x@other$loc.metrics.flags$OneRatio <- TRUE
    
    ZeroRatio <- 1 - OneRatio

    x@other$loc.metrics$PIC <- 1 - ((OneRatio*OneRatio) + (ZeroRatio*ZeroRatio))
    x@other$loc.metrics.flags$PIC <- TRUE
    
  } 
  
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

   return(x)
}
