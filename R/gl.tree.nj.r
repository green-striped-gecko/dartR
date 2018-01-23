#' Output an nj tree to summarize genetic similarity among populations
#'
#' This function is a wrapper for the nj\{ape\} function applied to Euclidian
#' distances calculated from the genlight object.
#'
#' @param gl -- Name of the genlight object containing the SNP data or a genind object containing presence absence data [required]
#' @param outgroup -- Vector containing the population names that are the outgroups [Default NULL]
#' @param type -- Type of dendrogram phylogram|cladogram|fan|unrooted [Default Phylogram]
#' @param labelsize -- Size of the labels as a proportion of the graphics default [Default 0.7]
#' @return A tree file of type phylo
#' @import adegenet
#' @importFrom stringr str_pad
#' @importFrom ape nj root plot.phylo
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' gl.tree.nj(testset.gl,type="fan")

gl.tree.nj <- function(gl, type="phylogram",outgroup=NULL,labelsize=0.7) {
x <- gl

    if(class(x)!="genlight") {
      cat("Fatal Error: genlight object required!\n"); stop()
    }

  # Convert gl object to a matrix of allele fequencies, locus by population
    cat("Converting to a matrix of frequencies, locus by populations\n")
    t=apply(as.matrix(x),2, tapply, pop(x), function(e) mean(e)/2)
  # Compute Euclidean distance
    cat("Computing Euclidean distances\n")
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
}
