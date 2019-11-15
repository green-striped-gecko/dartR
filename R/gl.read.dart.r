#' Import DarT data into R and conver it to a genlight object
#' 
#' This function is a wrapper function that allows you to convert you dart file into a genlight object in one step. In previous versions you had to use read.dart and then dart2genlight. In case you have individual metadata for each individual/sample you can specify as before in the dart2genlight command the file that combines the data.
#'
#' @param filename file containing the SNP data (csv file) [required]
#' @param ind.metafile file that contains additional information on individuals [required]
#' @param covfilename use ind.metafile parameter [deprectiated, NULL]
#' @param nas a character specifying NAs [nas = '-']
#' @param topskip a number specifying the number of rows to be skipped. If not provided the number of rows to be skipped are "guessed" by the number of rows with "*" at the beginning.
#' @param lastmetric specifies the last non genetic column (Default is "RepAvg"). Be sure to check if that is true, otherwise the number of individuals will not match. You can also specify the last column by a number.
#' @param recalc force the recalculation of locus metrics, in case individuals have been manually deleted from the input csv file [FALSE]
#' @param mono.rm force the removal of monomorphic loci (including all NAs), in case individuals have been manually deleted from the input csv file [FALSE]
#' @param probar show progress bar
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2, or as set by gl.set.verbose()]
#' @return a dart genlight object that contains individual metrics [if data were provided] and locus metrics [from a DArT report]. The dart genlight object can then be fed into a number of initial screening, export and export functions provided by the package. For some of the function it is necessary to have the metadata that was provided from DArT. Please check the vignette for more information. Additional information can also be found in the help documents for  \code{utils.read.dart}. 
#' @export
#' @examples{
#' dartfile <- system.file("extdata","testset_SNPs_2Row.csv", package="dartR")
#' metadata <- system.file("extdata","testset_metadata.csv", package="dartR")
#' gl <- gl.read.dart(dartfile, ind.metafile = metadata, probar=TRUE)
#' }


gl.read.dart <- function(filename, 
                         ind.metafile=NULL, 
                         recalc=FALSE, 
                         mono.rm=FALSE, 
                         nas = "-", 
                         topskip=NULL,  
                         lastmetric ="RepAvg", 
                         covfilename=NULL, 
                         probar=FALSE, 
                         verbose=NULL){
  
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
  
  if(verbose == 0 & probar){
    probar=FALSE
    #cat("  Note: Progress bar set to FALSE for verbosity==0\n")
  }

  # FLAG SCRIPT START
  
  if (verbose >= 1){
    if(verbose==5){
      cat("Starting",funname,"[ Build =",build,"]\n")
    } else {
      cat("Starting",funname,"\n")
    }
  }
  
# DO THE JOB
  
  # Deal with the redundant covfilename parameter
  if (is.null(ind.metafile)) {
    ind.metafile <- covfilename
  }
  
  dout <-utils.read.dart(filename = filename, nas=nas, topskip=topskip, lastmetric = lastmetric, verbose=verbose)
  glout <- utils.dart2genlight(dout, ind.metafile = ind.metafile, probar = probar, verbose=verbose)
  
  if (verbose>=2){
    cat(" Data read in. Please check carefully the output above\n")
  }  
  
  # Setting the recalc flags (TRUE=up-to-date, FALSE=no longer valid) for all locus metrics capable of being recalculated
  recalc.flags <-  c( "AvgPIC", "OneRatioRef","OneRatioSnp", "PICRef", "PICSnp", "CallRate",  "maf", "FreqHets" ,"FreqHomRef" , "FreqHomSnp", 
                      "monomorphs", "OneRatio", "PIC" )
  glout@other$loc.metrics.flags <-  data.frame(matrix(TRUE, nrow=1, ncol=length(recalc.flags)))
  names(glout@other$loc.metrics.flags) <- recalc.flags
  glout@other$verbose <- 2
  
  # Calculate locus metrics not provided by DArT
  
    # Calculate Read Depth
    if (is.null(glout@other$loc.metrics$rdepth)) {
      if (verbose >= 2){
        cat("  Read depth calculated and added to the locus metrics\n")
      }
      glout@other$loc.metrics$rdepth <- array(NA,nLoc(glout))
      for (i in 1:nLoc(glout)){
        called.ind <- round(nInd(glout)*glout@other$loc.metrics$CallRate[i],0)
        ref.count <- called.ind*glout@other$loc.metrics$OneRatioRef[i]
        alt.count <- called.ind*glout@other$loc.metrics$OneRatioSnp[i]
        sum.count.ref <- ref.count*glout@other$loc.metrics$AvgCountRef[i]
        sum.count.alt <- alt.count*glout@other$loc.metrics$AvgCountSnp[i]
        glout@other$loc.metrics$rdepth[i] <- round((sum.count.alt + sum.count.ref)/called.ind,1)
      }
    }
  
    # Calculate maf
    if (is.null(glout@other$loc.metrics$maf)) {
     utils.recalc.maf(glout, verbose=0)
     if (verbose >= 2){
       cat("  Minor Allele Frequency (maf) calculated and added to the locus metrics\n")
     }
    }
  
  # Calculate metrics provided by DArT, as a hedge against the user having deleted individuals from the input csv file
  
    if (recalc){
      if (verbose >= 2){
        cat("  Recalculating locus metrics provided by DArT (optionally specified)\n")
      }
      glout <- utils.recalc.avgpic(glout, verbose=0)
      glout <- utils.recalc.callrate(glout, verbose=0)
      glout <- utils.recalc.freqhets(glout, verbose=0)
      glout <- utils.recalc.freqhomref(glout, verbose=0)
      glout <- utils.recalc.freqhomsnp(glout, verbose=0)
    }
  
  # Remove monomorphs, which should not be present, but might have been introduced it the user deleted individuals from the input csv file
    
    glout@other$loc.metrics.flags$monomorphs <- FALSE
    if (mono.rm){
      if (verbose >= 2){
        cat("  Deleting monomorphic loci (optionally requested)\n")
      }
      glout <- gl.filter.monomorphs(glout,verbose=0)
    }
  
  # Set the SilicoDArT flags to FALSE
    glout@other$loc.metrics.flags$OneRatio <- FALSE
    glout@other$loc.metrics.flags$PIC <- FALSE

  # Provide a summary of the data
    if (verbose >= 3){
      cat("\nSummary of the SNP dataset\n")
      cat("  No. of loci:",nLoc(glout),"\n")  
      cat("  No. of individuals:",nInd(glout),"\n")
      cat("  No. of populations:",nPop(glout),"\n")
      if(!recalc){
        cat("  Locus metrics provided by DArT retained, not recalculated\n")
      }
      if(!mono.rm){
        cat("  Monomoporhic loci not deleted, assumed absent initially\n\n")
      }
    }
    
    # Create the history repository
    if (is.null(glout@other$history)) {
      glout@other$history <- list(match.call())
    }
    
# FLAG SCRIPT END
    
  if (verbose > 0) {
    cat(paste("Completed:",funname,"\n"))
  }
    
  return(glout)
}



