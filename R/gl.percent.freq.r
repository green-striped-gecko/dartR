#' Generate percentage allele frequencies by locus and population
#'
#' This is a support script, to take SNP data or SilocoDArT presence/absence data grouped into populations in a genelight or genind object \{adegenet\}
#' and generate a table of allele frequencies for each population and locus
#'
#' @param gl -- name of the genlight containing the SNP genotypes or genind object containing the presence/absence data [required]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A matrix with allele frequencies (genelight) or presence/absence frequencies (genind) broken down by population and locus
#' @export
#' @import adegenet
#' @importFrom reshape2 melt
#' @importFrom plyr ddply
#' @author Arthur Georges (bugs? Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' m <-  gl.percent.freq(testset.gl)
#' m

gl.percent.freq<- function(gl, v=2) {

  # ERROR CHECKING
  
  if(class(gl)!="genlight") {
    cat("Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }
  
  if (v < 0 | v > 5){
    cat("    Warning: verbosity must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    v <- 2
  }  
  
  if (v > 0) {
    cat("Starting gl.percent.freq: Calculating allele frequencies for populations\n")
    cat("  This may take some time -- be patient\n")
  }
  
# Calculate the required statistics, to be backward compatible 
  nmissing <- apply(as.matrix(gl),2, tapply, pop(gl), function(x){sum(is.na(x))})
  nobs <- apply(as.matrix(gl),2, tapply, pop(gl), function(x){sum(!is.na(x))}) 
  n <- nmissing + nobs
  sum <- apply(as.matrix(gl),2, tapply, pop(gl), function(x){sum(x,na.rm=TRUE)})
  f <- apply(as.matrix(gl),2, tapply, pop(gl), function(x){mean(x, na.rm=TRUE)/2})
  f <- f*100

# Convert the matricies to long format pop, locus, value
  nobs <- melt(nobs, na.rm=FALSE)
  nmissing <- melt(nmissing, na.rm=FALSE)
  n <- melt(n, na.rm=FALSE)
  f <- melt(f, na.rm=FALSE)
  sum <- melt(sum, na.rm=FALSE)
  
  if(nPop(gl) == 1) {
    m <- cbind(levels(pop(gl)),rownames(sum),sum,nobs,nmissing,f,n)
  } else {
    m <- cbind(sum,nobs[,3],nmissing[,3],f[,3],n[,3])
  }  
  colnames(m) <- c("popn","locus","sum","nobs","nmissing","frequency","n")
  
  if(v > 0) {cat("Completed gl.percent.freq\n\n")}
  
  return(m)
  
}
