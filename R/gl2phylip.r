#' Create a Phylip input distance matrix from a genlight (SNP) \{adegenet\} object
#'
#' This function calculates and returns a matrix of Euclidean distances between populations
#' and produces an input file for the phylogenetic program Phylip (Joe Felsenstein).
#'
#' @param x Name of the genlight object containing the SNP data or a genind object containing presence absence data [required]
#' @param outfile Name of the file to become the input file for phylip [default phyinput.txt]
#' @param outpath -- path where to save the output file [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath="." when calling this function to direct output files to your working directory.
#' @param verbose -- specify the level of verbosity: 0, silent, fatal errors only; 1, flag function begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @param bstrap Number of bootstrap replicates [default 1]
#' @return Matrix of Euclidean distances between populations
#' @import utils
#' @importFrom stats dist
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' \donttest{
#' result <- gl2phylip(testset.gl, outfile="test.txt", bstrap=10)
#' }

gl2phylip <- function(x, outfile="phyinput.txt", outpath=tempdir(), bstrap=1, verbose=NULL) {

# TIDY UP FILE SPECS

  outfilespec <- file.path(outpath, outfile)
  funname <- match.call()[[1]]

# FLAG SCRIPT START
  # set verbosity
  if (is.null(verbose) & !is.null(x@other$verbose)) verbose=x@other$verbose
  if (is.null(verbose)) verbose=2
 

  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }

  if (verbose > 0) {
    cat("Starting",funname,"\n")
  }

# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }

  # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nInd(x))
      pop(x) <- as.factor(pop(x))
    }

  # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x, verbose=0)
    if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {cat("  Warning: genlight object contains monomorphic loci\n")}

# DO THE JOB

  # Convert gl object to a matrix of allele fequencies, locus by population
    if (verbose >= 2) {cat("Converting to a matrix of frequencies, locus by populations\n")}
    t=apply(as.matrix(x),2, tapply, pop(x), function(e) mean(e)/2)
  # Compute Euclidean distance
    if (verbose >= 2) {cat("Computing Euclidean distances\n")}
    d <- round(as.matrix(dist(t)),4)
    row.names(d) <- c(paste(row.names(d),"          "))
    row.names(d) <- substr(row.names(d),1,10)

  # Output phylip data file
    if (verbose >= 2) {
      cat("Writing the Phylip input file", outfilespec, "\n")
      if (bstrap > 1) {
        cat("Repeating calculations for", bstrap, "iterations\n")
      }
    }
    npops <- length(levels(factor(pop(x))))
    sink(outfilespec)
    cat(c("   ",npops,"\n"))
    for (i in 1:npops){
      cat(row.names(d)[i],d[i,],"\n")
    }

    # Check if bootstrap replicates are required
    if (bstrap > 1) {

    # Repeat for each bootstrap replicate
      for (j in (2:bstrap)) {

      # subsample the loci, with replication
        h <- seq(1:nLoc(x))
        newx <- x[,sample(h,size=nLoc(x),replace=TRUE)]

      # Convert gl object to a matrix of allele fequencies, locus by population
        t=apply(as.matrix(newx),2, tapply, pop(x), function(e) mean(e)/2)

      # Compute Euclidean distance
        d <- round(as.matrix(dist(t)),4)
        row.names(d) <- c(paste(row.names(d),"          "))
        row.names(d) <- substr(row.names(d),1,10)

      # Output phylip data file
        npops <- length(levels(factor(pop(x))))
        cat(c("   ",npops,"\n"))
        for (i in 1:npops){
          cat(row.names(d)[i],d[i,],"\n")
        }
      }
    }
    sink()
    if (verbose >= 2) {cat("Closing output file",outfile,"\n")}

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

    return(d)
}
