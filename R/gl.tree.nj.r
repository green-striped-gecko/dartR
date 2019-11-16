#' Output an nj tree to summarize genetic similarity among populations
#'
#' This function is a wrapper for the nj\{ape\} function applied to Euclidian
#' distances calculated from the genlight object.
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param outgroup -- Vector containing the population names that are the outgroups [Default NULL]
#' @param type -- Type of dendrogram phylogram|cladogram|fan|unrooted [Default Phylogram]
#' @param labelsize -- Size of the labels as a proportion of the graphics default [Default 0.7]
#' @param treefile -- Name of the file for the tree topology using Newick format [Default NULL].
#' @param verbose -- specify the level of verbosity: 0, silent, fatal errors only; 1, flag function begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A tree file of class phylo
#' @importFrom stringr str_pad
#' @importFrom ape nj root plot.phylo write.tree
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' # SNP data
#'   gl.tree.nj(testset.gl,type="fan")
#' # Tag P/A data
#'   gl.tree.nj(testset.gs,type="fan")

gl.tree.nj <- function(x, 
                       type="phylogram", 
                       outgroup=NULL, 
                       labelsize=0.7, 
                       treefile=NULL, 
                       verbose=NULL) {

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
    if (verbose >= 2){cat("  Processing  Presence/Absence (SilicoDArT) data\n")}
  } else if (all(x@ploidy == 2)){
    if (verbose >= 2){cat("  Processing a SNP dataset\n")}
  } else {
    stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
  }

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
    tree <- ape::nj(d)
    if (!is.null(outgroup)) {
      # Function plot.phylo{ape} has the labels all of the same length
      outgroup <- stringr::str_pad(outgroup, nchar(tree$tip.label[1]), side = c("right"), pad = " ")
      # Truncate to 10 characters
      outgroup <- substr(outgroup,1,10)
      # Root the tree
      tree <- ape::root(tree, outgroup)
      # Plot the tree
      # Save the prior settings for mfrow, oma, mai and pty, and reassign
      op <- par(mfrow = c(1, 1), oma=c(1,1,1,1), mai=c(0,0,0,0),pty="m")
      ape::plot.phylo(tree, type=type, cex=labelsize)
    } else {
      # Just plot the tree unrooted
      op <- par(mfrow = c(1, 1), oma=c(1,1,1,1), mai=c(0,0,0,0),pty="m")
      ape::plot.phylo(tree, type=type, cex=labelsize)
    }
    
  # Output the tree file
    if(!is.null(treefile)){
      if(verbose>=2){cat("  Writing the tree topology to",treefile,"\n")}
      write.tree(tree,file=treefile)
    }  
    
  # Reset the par options    
    par(op)
    
# FLAG SCRIPT END
    
    if (verbose > 0) {
      cat("Completed:",funname,"\n")
    }
    
    return(tree)

}
