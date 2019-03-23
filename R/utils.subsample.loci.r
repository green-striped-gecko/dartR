#' Subsample n loci from a genlight object and return as a genlight object
#'
#' This is a support script, to subsample a genlight \{adegenet\} object based on loci. Two methods are used
#' to subsample, random and based on information content (avgPIC).
#'
#' @param x -- name of the genlight object containing the SNP genotypes by specimen and population [required]
#' @param k -- number of loci to include in the subsample [required]
#' @param method -- "random", in which case the loci are sampled at random; or avgPIC, in which case the top n loci
#' ranked on information content (AvgPIC) are chosen [default "random"]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A genlight object with k loci
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' result <- utils.subsample.loci(testset.gl, k=200, method="avgPIC")

utils.subsample.loci <- function(x, k, method="random", verbose=2) {

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
  
  if(method=="random") {
    if (verbose>=2){cat("    Subsampling at random,",k,"loci from",class(x),"object","\n")}
    #nblocks <- trunc((ncol(x)/n)+1)
    #blocks <- lapply(seploc(x, n.block=nblocks, random=TRUE, parallel=FALSE),as.matrix)
    blocks <- seploc(x, block.size=k, random=TRUE, parallel=FALSE)
    x.new <- blocks$block.1
    if (verbose>=3) {
      cat("      No. of loci retained =", nLoc(x.new),"\n")
      #cat("      Note: SNP metadata discarded\n")
    }
  } else if (method=="AvgPIC" | method=="avgpic" | method=='avgPIC' | method=='AvgPic'){
    x.new <- x[, order(-x@other$loc.metrics["AvgPIC"])]
    x.new <- x.new[,1:k]
    x.new@other$loc.metrics <- x.new@other$loc.metrics[1:k,]
    if (verbose>=3) {
      cat("      No. of loci retained =", nLoc(x.new),"\n")
      #cat("      Note: SNP metadata discarded\n")
    }
  } else {
    cat ("Fatal Error in gl.sample.loci.r: method must be random or avgpic\n"); stop()
  }

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

  return(x.new)

}

#test <- gl.subsample.loci(gl, 12, method="avgpic")
#as.matrix(test)[1:20,]

#as.matrix(x)[1:20,1:10]

#as.matrix(x.new)[1:20,1:10]

#x<-gl
#method<-"random"
#n=100
