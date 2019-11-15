#' A utility script to reset to FALSE (or TRUE) the locus metric flags after some individuals or populations have been deleted.
#'
#' The locus metadata supplied by DArT has OneRatioRef, OneRatioSnp, PICRef, PICSnp, and AvgPIC included,
#' but the allelec composition will change when some individuals are removed from the dataset and so the initial statistics will
#' no longer apply. This applies also to some variable calculated by dartR (e.g. maf). This script resets the locus metrics
#' flags to FALSE to indicate that these statistics in the genlight object are no longer current. The verbosity default is also set, and
#' in the case of SilcoDArT, the flags PIC and OneRatio are also set.
#' 
#' If the locus metrics do not exist then they are added to the genlight object but not populated. If the locus metrics flags 
#' do not exist, then they are added to the genlight object and set to FALSE (or TRUE).
#'
#' @param x -- name of the genlight object containing the SNP data or tag presence/absence data (SilicoDArT) [required]
#' @param set -- set the flags to TRUE or FALSE [FALSE]
#' @param set.verbosity -- set the default verbosity for all functions, where verbosity is not specified
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default NULL]
#' @return The modified genlight object
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @seealso \code{utils.recalc.metrics} for recalculating all metrics, \code{utils.recalc.callrate} for recalculating CallRate,
#' \code{utils.recalc.freqhomref} for recalculating frequency of homozygous reference, \code{utils.recalc.freqhomsnp} for recalculating frequency of homozygous alternate,
#' \code{utils.recalc.freqhet} for recalculating frequency of heterozygotes, \code{gl.recalc.maf} for recalculating minor allele frequency,
#' \code{gl.recalc.rdepth} for recalculating average read depth
#' @examples
#' #result <- utils.reset.flags(testset.gl)

utils.reset.flags <- function(x, set=FALSE, set.verbosity=2, verbose=NULL) {
  
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
      cat("Starting",funname,"[Build =",build,"\n")
    } else {
      cat("Starting",funname,"\n")
    }
  }
  
# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    stop("  Fatal Error: genlight object required!\n")
  }
  
  if (all(x@ploidy == 1)){
    if (verbose >= 2){cat("  Processing  Presence/Absence (SilicoDArT) data\n")}
    data.type <- "SilicoDArT"
  } else if (all(x@ploidy == 2)){
    if (verbose >= 2){cat("  Processing a SNP dataset\n")}
    data.type <- "SNP"
  } else {
    stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
  }

# SCRIPT SPECIFIC ERROR TESTING

  if (set.verbosity < 0 | set.verbosity > 5){
    cat("  Warning: Parameter 'set.verbosity' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    set.verbosity <- 2
  }
  
# DO THE JOB
  
  if (data.type=="SNP"){
    if(verbose >= 2){
      cat("  Resetting flags for AvgPIC, OneRatioRef, OneRatioSnp, PICRef, PICSnp, CallRate, maf, FreqHets, FreqHomRef, FreqHomSnp, monomorphs, OneRatio, PIC to",set,"\n")
      cat("  Resetting SilicoDArT flags for OneRatio, PIC to FALSE\n")
    }
    
  #AvgPIC
    if (is.null(x@other$loc.metrics$AvgPIC)) {
      x@other$loc.metrics$AvgPIC <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric AvgPIC does not exist, creating slot @other$loc.metrics$AvgPIC\n")
      }
    }
    x@other$loc.metrics.flags$AvgPIC <- set
    
  #OneRatioRef  
    if (is.null(x@other$loc.metrics$OneRatioRef)) {
      x@other$loc.metrics$OneRatioRef <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric OneRatioRef does not exist, creating slot @other$loc.metrics$OneRatioRef\n")
      }
    }
    x@other$loc.metrics.flags$OneRatioRef <- set

  #OneRatioSnp        
    if (is.null(x@other$loc.metrics$OneRatioSnp)) {
      x@other$loc.metrics$OneRatioSnp <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric OneRatioSnp does not exist, creating slot @other$loc.metrics$OneRatioSnp\n")
      }
    }
    x@other$loc.metrics.flags$OneRatioSnp <- set

  #PICRef    
    if (is.null(x@other$loc.metrics$PICRef)) {
      x@other$loc.metrics$PICRef <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric PICRef does not exist, creating slot @other$loc.metrics$PICRef\n")
      }
    }
    x@other$loc.metrics.flags$PICRef <- set
    
  #PICSnp
    if (is.null(x@other$loc.metrics$PICSnp)) {
      x@other$loc.metrics$PICSnp <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric PICSnp does not exist, creating slot @other$loc.metrics$PICSnp\n")
      }
    }
    x@other$loc.metrics.flags$PICSnp <- set

  #CallRate
    if (is.null(x@other$loc.metrics$CallRate)) {
      x@other$loc.metrics$CallRate <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric CallRate does not exist, creating slot @other$loc.metrics$CallRate\n")
      }
    }
    x@other$loc.metrics.flags$CallRate <- set
    
  #FreqHomRef
    if (is.null(x@other$loc.metrics$FreqHomRef)) {
      x@other$loc.metrics$FreqHomRef <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric FreqHomRef does not exist, creating slot @other$loc.metrics$FreqHomRef\n")
      }
    }
    x@other$loc.metrics.flags$FreqHomRef <- set
    
  #FreqHomSnp
    if (is.null(x@other$loc.metrics$FreqHomSnp)) {
      x@other$loc.metrics$FreqHomSnp <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric FreqHomSnp does not exist, creating slot @other$loc.metrics$FreqHomSnp\n")
      }
    }
    x@other$loc.metrics.flags$FreqHomSnp <- set
    
  #FreqHets
    if (is.null(x@other$loc.metrics$FreqHets)) {
      x@other$loc.metrics$FreqHets <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric FreqHets does not exist, creating slot @other$loc.metrics$FreqHets\n")
      }
    }
    x@other$loc.metrics.flags$FreqHets <- set
    
  #monomorphs
    if (is.null(x@other$loc.metrics$monomorphs)) {
      x@other$loc.metrics$monomorphs <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric monomorphs does not exist, creating slot @other$loc.metrics$monomorphs\n")
      }
    }
    x@other$loc.metrics.flags$monomorphs <- set
    
    #maf
    if (is.null(x@other$loc.metrics$maf)) {
      x@other$loc.metrics$maf <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric maf does not exist, creating slot @other$loc.metrics$maf\n")
      }
    }
    x@other$loc.metrics.flags$maf <- set
    
  #OneRatio
    if (is.null(x@other$loc.metrics$OneRatio)) {
      x@other$loc.metrics$OneRatio <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric OneRatio does not exist, creating slot @other$loc.metrics$OneRatio\n")
      }
    }
    x@other$loc.metrics.flags$OneRatio <- FALSE
      
  #PIC
    if (is.null(x@other$loc.metrics$PIC)) {
      x@other$loc.metrics$PIC <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric PIC does not exist, creating slot @other$loc.metrics$PIC\n")
      }
    }
    x@other$loc.metrics.flags$PIC <- FALSE
    
  #verbosity
    if (is.null(x@other$loc.metrics$verbose)) {
      x@other$loc.metrics$verbose <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric 'PIC'verbose' does not exist, creating slot @other$loc.metrics$verbose\n")
      }
    }
    x@other$verbose <- set.verbosity
  }
  
  if (data.type=="SilicoDArT"){
    if(verbose >= 2){
      cat("  Resetting flags for CallRate, PIC, OneRatio, monomorphs to",set,"\n")
      cat("  Setting SNP flags for AvgPIC, OneRatioRef, OneRatioSnp, PICRef, PICSnp, maf, FreqHets, FreqHomRef, FreqHomSnp to FALSE\n")
    }
    
    #AvgPIC
    if (is.null(x@other$loc.metrics$AvgPIC)) {
      x@other$loc.metrics$AvgPIC <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric AvgPIC does not exist, creating slot @other$loc.metrics$AvgPIC\n")
      }
    }
    x@other$loc.metrics.flags$AvgPIC <- FALSE
    
    #OneRatioRef  
    if (is.null(x@other$loc.metrics$OneRatioRef)) {
      x@other$loc.metrics$OneRatioRef <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric OneRatioRef does not exist, creating slot @other$loc.metrics$OneRatioRef\n")
      }
    }
    x@other$loc.metrics.flags$OneRatioRef <- FALSE
    
    #OneRatioSnp        
    if (is.null(x@other$loc.metrics$OneRatioSnp)) {
      x@other$loc.metrics$OneRatioSnp <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric OneRatioSnp does not exist, creating slot @other$loc.metrics$OneRatioSnp\n")
      }
    }
    x@other$loc.metrics.flags$OneRatioSnp <- FALSE
    
    #PICRef    
    if (is.null(x@other$loc.metrics$PICRef)) {
      x@other$loc.metrics$PICRef <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric PICRef does not exist, creating slot @other$loc.metrics$PICRef\n")
      }
    }
    x@other$loc.metrics.flags$PICRef <- FALSE
    
    #PICSnp
    if (is.null(x@other$loc.metrics$PICSnp)) {
      x@other$loc.metrics$PICSnp <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric PICSnp does not exist, creating slot @other$loc.metrics$PICSnp\n")
      }
    }
    x@other$loc.metrics.flags$PICSnp <- FALSE
    
    #CallRate
    if (is.null(x@other$loc.metrics$CallRate)) {
      x@other$loc.metrics$CallRate <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric CallRate does not exist, creating slot @other$loc.metrics$CallRate\n")
      }
    }
    x@other$loc.metrics.flags$CallRate <- set
    
    #FreqHomRef
    if (is.null(x@other$loc.metrics$FreqHomRef)) {
      x@other$loc.metrics$FreqHomRef <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric FreqHomRef does not exist, creating slot @other$loc.metrics$FreqHomRef\n")
      }
    }
    x@other$loc.metrics.flags$FreqHomRef <- FALSE
    
    #FreqHomSnp
    if (is.null(x@other$loc.metrics$FreqHomSnp)) {
      x@other$loc.metrics$FreqHomSnp <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric FreqHomSnp does not exist, creating slot @other$loc.metrics$FreqHomSnp\n")
      }
    }
    x@other$loc.metrics.flags$FreqHomSnp <- FALSE
    
    #FreqHets
    if (is.null(x@other$loc.metrics$FreqHets)) {
      x@other$loc.metrics$FreqHets <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric FreqHets does not exist, creating slot @other$loc.metrics$FreqHets\n")
      }
    }
    x@other$loc.metrics.flags$FreqHets <- FALSE
    
    #monomorphs
    if (is.null(x@other$loc.metrics$monomorphs)) {
      x@other$loc.metrics$monomorphs <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric monomorphs does not exist, creating slot @other$loc.metrics$monomorphs\n")
      }
    }
    x@other$loc.metrics.flags$monomorphs <- set
    
    #maf
    if (is.null(x@other$loc.metrics$maf)) {
      x@other$loc.metrics$maf <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric maf does not exist, creating slot @other$loc.metrics$maf\n")
      }
    }
    x@other$loc.metrics.flags$maf <- FALSE
    
    #OneRatio
    if (is.null(x@other$loc.metrics$OneRatio)) {
      x@other$loc.metrics$OneRatio <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric OneRatio does not exist, creating slot @other$loc.metrics$OneRatio\n")
      }
    }
    x@other$loc.metrics.flags$OneRatio <- set
    
    #PIC
    if (is.null(x@other$loc.metrics$PIC)) {
      x@other$loc.metrics$PIC <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric PIC does not exist, creating slot @other$loc.metrics$PIC\n")
      }
    }
    x@other$loc.metrics.flags$PIC <- set
    
    #verbosity
    if (is.null(x@other$loc.metrics$verbose)) {
      x@other$loc.metrics$verbose <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric 'PIC'verbose' does not exist, creating slot @other$loc.metrics$verbose\n")
      }
    }
    x@other$verbose <- set.verbosity
    
  }

# ADD TO HISTORY not in utils functions
  
# FLAG SCRIPT END

  if (verbose >= 1) {
    cat("Completed:", funname, "\n")
  }
  
  return(x)
  
}
