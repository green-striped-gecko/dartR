#'  Calculate a statistic for each locus by group
#'  
#'  An internal function essentially to convey readability to rather contorted R code.
#'  It takes as input a genlight \{adegenet\} object with an index variable (say, population)
#'  and calculates the selected statistic for each locus, broken down by the groups defined
#'  by the index variable.
#'
#' @param gl -- name of the genlight object containing the SNP data [required]
#' @param method -- breakdown variable [default pop(x)]
#' @param stat -- statistic to calculate: mean [only mean(x)/2 currently implemented]
#' @return A matrix, populations (rows) by loci (columns), showing the statistic [mean/2]
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' #result <- gl.gene.freq(testset.gl, method=pop(gl), stat="mean")

gl.gene.freq <- function(gl, method=pop(gl), stat="mean") {
  
  # For each column of gl, apply the function tapply with index variable pop(gl)and function mean(x)/2
  if (stat=="mean"){
    t=apply(as.matrix(gl),2, tapply, pop(gl), function(x){mean(x, na.rm=TRUE)/2}) 
  } else {
    cat("Only stat=\"mean\" is implemented\n"); stop()
  }
  return(t)
} 
