#' Output an nj tree to summarize genetic similarity among populations
#'
#' This function is a wrapper for the nj\{ape\} function applied to Euclidian
#' distances calculated from the genlight object.
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param outgroup -- Vector containing the population names that are the outgroups [Default NULL]
#' @param type -- Type of dendrogram phylogram|cladogram|fan|unrooted [Default Phylogram]
#' @param labelsize -- Size of the labels as a proportion of the graphics default [Default 0.7]
#' @param verbose -- specify the level of verbosity: 0, silent, fatal errors only; 1, flag function begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A tree file of class phylo
#' @import adegenet
#' @importFrom stringr str_pad
#' @importFrom ape nj root plot.phylo
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl.tree.nj(testset.gl,type="fan")

# Last amended 3-Feb-19

gl.tree.nj <- function(x, type="phylogram", outgroup=NULL, labelsize=0.7, verbose=2) {

# TIDY UP FILE SPECS

  funname <- match.call()[[1]]

# FLAG SCRIPT START

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

# DO THE JOB

  # Convert gl object to a matrix of allele fequencies, locus by population
    if (verbose >= 2) {cat("  Converting to a matrix of frequencies, locus by populations\n")}
    t=apply(as.matrix(x),2, tapply, pop(x), function(e) mean(e)/2)
  # Compute Euclidean distance
    if (verbose >= 2) {cat("  Computing Euclidean distances\n")}
    d <- round(as.matrix(dist(t)),4)
    row.names(d) <- c(paste(row.names(d),"          "))
    row.names(d) <- substr(row.names(d),1,10)
    
  # Plot the distances as an nj tree  
    tree <- nj(d)
    if (!is.null(outgroup)) {
      # Function plot.phylo{ape} has the labels all of the same length
      outgroup <- str_pad(outgroup, nchar(tree$tip.label[1]), side = c("right"), pad = " ")
      # Truncate to 10 characters
      outgroup <- substr(outgroup,1,10)
      # Root the tree
      rtree <- root(tree, outgroup)
      # Plot the tree
      plot.phylo(rtree, type=type, cex=labelsize)
      return(rtree)
    } else {
      # Just plot the tree unrooted
      plot.phylo(tree, type=type, cex=labelsize)
      return(tree)
    }

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
    
}
