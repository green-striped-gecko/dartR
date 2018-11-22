#' Subsample n loci from a genlight object and return as a genlight object
#'
#' This is a support script, to subsample a genlight \{adegenet\} object based on loci. Two methods are used
#' to subsample, random and based on information content (avgPIC)
#'
#' @param gl -- name of the genlight object containing the SNP genotypes by specimen and population [required]
#' @param n -- number of loci to include in the subsample [required]
#' @param method -- "random", in which case the loci are sampled at random; or avgPIC, in which case the top n loci
#' ranked on information content (AvgPIC) are chosen [default "random"]
#' @return A genlight object with n loci
#' @export
#' @author Arthur Georges (bugs? Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' result <- gl.subsample.loci(testset.gl, n=200, method="avgPIC")

gl.subsample.loci <- function(gl, n, method="random") {
x <- gl

  if(method=="random") {
    cat("Subsampling at random, approximately",n,"loci from",class(x),"object","\n")
    nblocks <- trunc((ncol(x)/n)+1)
    blocks <- lapply(seploc(x, n.block=nblocks, random=TRUE, parallel=FALSE),as.matrix)
    x.new <- blocks$block.1
    cat("   No. of loci retained =", ncol(x.new),"\n")
    cat("   Note: SNP metadata discarded\n")
  } else if (method=="AvgPIC" | method=="avgpic" | method=='avgPIC'){
    x.new <- x[, order(-x@other$loc.metrics["AvgPIC"])]
    x.new <- x.new[,1:n]
    cat("   No. of loci retained =", ncol(x.new),"\n")
    cat("   Note: SNP metadata discarded\n")
  } else {
    cat ("Fatal Error in gl.sample.loci.r: method must be random or repavg\n"); stop()
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
