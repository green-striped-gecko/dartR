#' Create a Phylip input distance matrix from a genlight (SNP) or genind (SilicoDarT) \{adegenet\} object
#'
#' This function calculates and returns a matrix of Euclidean distances between populations
#' and produces an input file for the phylogenetic program Phylip (Joe Felsenstein).
#'
#' @param gl Name of the genlight object containing the SNP data or a genind object containing presence absence data [required]
#' @param outfile Name of the file to become the input file for phylip [default phyinput.txt]
#' @param outpath path where to save the output file (set to tempdir by default)
#' @param bstrap Number of bootstrap replicates [default 1]
#' @return Matrix of Euclidean distances between populations
#' @import adegenet utils
#' @importFrom stats dist
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' \donttest{
#' result <- gl2phylip(testset.gl, outfile="test.txt", bstrap=10)
#' }


gl2phylip <- function(gl, outfile="phyinput.txt", outpath=tempdir(), bstrap=1) {

    outfile <- file.path(outpath, outfile)
    x <- gl

    if(class(x)!="genlight") {
      cat("Fatal Error: genlight object required for gl2phylip!\n"); stop()
    }

  # Convert gl object to a matrix of allele fequencies, locus by population
    cat("Converting to a matrix of frequencies, locus by populations\n")
    t=apply(as.matrix(x),2, tapply, pop(x), function(e) mean(e)/2)
  # Compute Euclidean distance
    cat("Computing Euclidean distances\n")
    d <- round(as.matrix(dist(t)),4)
    row.names(d) <- c(paste(row.names(d),"          "))
    row.names(d) <- substr(row.names(d),1,10)

  # Output phylip data file
    cat("Writing the Phylip input file", outfile, "\n")
    npops <- length(levels(factor(pop(x))))
    sink(outfile)
    cat(c("   ",npops,"\n"))
    for (i in 1:npops){
      cat(row.names(d)[i],d[i,],"\n")
    }

    #cat("Repeating calculations for", bstrap, "iterations\n")
    # Check if bootstrap replicates are required
    if (bstrap > 1) {

      # Set up the progress counter
      #pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
      #getTxtProgressBar(pb)

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
        #sink("outfile", append=TRUE)
        cat(c("   ",npops,"\n"))
        for (i in 1:npops){
          cat(row.names(d)[i],d[i,],"\n")
        }
      }
      #setTxtProgressBar(pb, j)
    }
    sink()
    return(d)
}
