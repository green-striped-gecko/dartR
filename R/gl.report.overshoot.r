#' @name gl.report.overshoot
#' 
#' @title Reports loci for which the SNP has been trimmed from the sequence tag along with the adaptor 
#'
#' @description This function checks the position of the SNP within the trimmed sequence tag and identifies those for which the SNP position is outside
#' the trimmed sequence tag. This can happen, rarely, when the sequence containing the SNP resembles the adaptor.
#' 
#' @details The SNP genotype can still be used in most analyses, but functions like gl2fasta() will present challenges if the SNP has been trimmed from
#' the sequence tag.
#' 
#' @param x -- name of the genlight object [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return NULL
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl.report.overshoot(testset.gl)
#' 
#' #' @seealso \code{\link{gl.filter.overshoot}}, \code{\link{gl.access.report}}

gl.report.overshoot <- function(x, verbose=NULL) {

# TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  
  # GENERAL ERROR CHECKING
  
  datatype <- utils.check.gl(x)
  verbose <- utils.check.verbosity(verbose)
  
# FLAG SCRIPT START
  
  if (verbose >= 1){
    if(verbose==5){
      cat("Starting",funname,"[ Build =",build,"]\n")
    } else {
      cat("Starting",funname,"\n")
    }
  }
  
# STANDARD ERROR CHECKING
  
  if (datatype=="SilicoDArT"){
    cat("  Detected Presence/Absence (SilicoDArT) data\n")
    stop(error("Cannot identify overshoot arising from SNPS deleted with adaptors for fragment presence/absence data. 
               Please provide a SNP dataset.\n"))
  } 
  
# SCRIPT SPECIFIC ERROR CHECKING
  
  if(length(x@other$loc.metrics$TrimmedSequence) != nLoc(x)) {
    stop(error("Fatal Error: Data must include Trimmed Sequences for each loci in a column called 'TrimmedSequence' 
               in the @other$loc.metrics slot.\n"))
  }
  if(length(x@other$loc.metrics$SnpPosition) != nLoc(x)) {
    stop(error("Fatal Error: Data must include position information for each loci.\n"))
  }
  
# DO THE JOB

  if (verbose >=2) {cat(report("  Identifying loci for which the SNP has been trimmed with the adaptor\n"))}

  trimmed <- as.character(x@other$loc.metrics$TrimmedSequence)
  snpos <- x@other$loc.metrics$SnpPosition
  # Shift the index for snppos to start from 1 not zero
  snpos <- snpos +1
  # Pull those loci for which the SNP position is greater than the tag length
  xx <- x[,snpos > nchar(trimmed)]
  # Report the number of such loci
    cat("  No. of loci with SNP falling outside the trimmed sequence:",nLoc(xx),"\n")
    if(nLoc(xx) > 0){cat("\n",paste(locNames(xx),"\n"))}
    
    # creating temp file names
    temp_table <- tempfile(pattern = "table_")
    
    # saving to tempdir
    setwd(tempdir())
    saveRDS(df, file = temp_table)
    if(verbose>=2){cat(report("  Saving the list of recalitrant loci to the tempfile as",temp_table,"using saveRDS\n"))}
    if(verbose>=2){cat(report("  NOTE: Retrieve output files from tempdir using gl.access.report()\n"))}  

# FLAG SCRIPT END

  if (verbose >= 1) {
    cat("Completed:",funname,"\n")
  }
    
  return(NULL)

}

