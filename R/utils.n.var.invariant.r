#' An utility script to calculate the number of variant and invariant sites by
#' locus
#'
#' calculate the number of variant and invariant sites by locus and add them as
#' columns in \code{loc.metrics}. This can be useful to conduct further
#' filtering, for example where only loci with secondaries are wanted for
#' phylogenetic analyses.
#'
#' Invariant sites are the site (nucleotide) that are not polymorphic. When the
#' locus metadata supplied by DArT includes the sequence of the allele
#' (\code{TrimmedSequence}), it is used by this function to estimate the number
#' of sites that were sequenced in each tag (read). This script then subtracts
#' the number of polymorphic sites. The final number of variant and invariant
#' sites is the added to the table in \code{gl@others$loc.metrics}.
#'
#' \strong{NOTE}: It is important to realise that this function correctly
#' estimates the number of variant and invariant sites only when it is executed on
#' \code{genlight} objects before secondaries are removed.
#'
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'   progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return The modified genlight object
#' @author Carlo Pacioni (Post to \url{https://groups.google.com/d/forum/dartr})
#' @seealso
#' \code{\link{gl.filter.secondaries}},\code{\link{gl.report.heterozygosity}}
#' @export
#' @import data.table
#' @examples
#' #out <- utils.n.invariant(testset.gl)

utils.n.var.invariant <- function(x, verbose=NULL) {
  
  # SET VERBOSITY
  
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  
  funname <- match.call()[[1]]
  utils.flag.start(func=funname,build="Jody",v=verbose) # I'm not quite sure what the build here should be
  
  # CHECK DATATYPE 
  datatype <- utils.check.datatype(x,verbose=verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
  if (grepl(x@other$history, pattern = "gl.filter.secondaries")) {
    warning("This function should be run before removing secondaries. A gl.filter.secondaries call was found in the history. This may cause the results to be incorrect")
  }
  if(isFALSE("CloneID" %in% names(x$other$loc.metrics))) {
    # Extract the clone ID number
    a <- strsplit(as.character(x@other$loc.metrics$AlleleID),"\\|")
    b <- unlist(a)[ c(TRUE,FALSE,FALSE) ] # You may want to put these two lines in a FUN and use it here and in gl/report.scondaries
    x@other$loc.metrics$CloneID <- b
  }
  if(isFALSE("TrimmedSequence" %in% names(x$other$loc.metrics))) {
    stop("The column 'TrimmedSequence' is not present in x$other$loc.metrics, but it is needed for this function")
  }
  
  
  # DO THE JOB
  
  if (verbose >= 2) {cat("  Calculating n invariant sites\n")}
  proc.data <- data.table(x$other$loc.metrics) # using data.table
  proc.data[, lenTrimSeq := nchar(as.character(TrimmedSequence))]
  proc.data[, n.variant := .N, by=CloneID]
  proc.data[, n.invariant := lenTrimSeq - n.variant]
  x@other$loc.metrics <- setDF(proc.data)
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat(report("Completed:",funname,"\n"))
  }
  
  return(x)
}

