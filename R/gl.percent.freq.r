#' Generate percentage allele frequencies by locus and population
#'
#' This is a support script, to take SNP data or SilocoDArT presence/absence data grouped into populations in a genelight object \{adegenet\}
#' and generate a table of allele frequencies for each population and locus
#'
#' @param x -- name of the genlight object containing the SNP or Tag P/A (SilicoDArT) data [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return A matrix with allele (SNP data) or presence/absence frequencies (Tag P/A data) broken down by population and locus
#' @export
#' @importFrom plyr ddply
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' m <-  gl.percent.freq(testset.gl)
#' m

gl.percent.freq<- function(x, verbose=NULL) {

# CHECK IF PACKAGES ARE INSTALLED
  pkg <- "reshape2"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop("Package",pkg," needed for this function to work. Please   install it.") } else {
      
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
    data.type <- "SilicoDArT"
  } else if (all(x@ploidy == 2)){
    if (verbose >= 2){cat("  Processing a SNP dataset\n")}
    data.type <- "SNP"
  } else {
    stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
  }
    
# SCRIPT SPECIFIC ERROR CHECKING
    
  # Checking for and removing monomorphic loci
    if(!(x@other$loc.metrics.flags$monomorphs==TRUE)){
      if (verbose >= 1) {cat("Warning: Monomorphic loci retained, used in calculations\n")}
    } 

# DO THE JOB
  
  mat <- as.matrix(x)
  
  if(data.type=="SilicoDArT"){

    if (verbose >= 2) {
      cat("Starting gl.percent.freq: Calculating Tag P/A frequencies for populations\n")
    }
    # Treat SilicoDArT as biallelic, no heterozygotes
    mat[mat==1] <- 2

  } else {

    if (verbose >= 1) {
      cat("Starting gl.percent.freq: Calculating allele frequencies for populations\n")
    }
    if (verbose >= 3){
      cat("  This may take some time -- be patient\n")
    }
  }
  
# Calculate the required statistics, to be backward compatible 
  nmissing <- apply(mat,2, tapply, pop(x), function(x){sum(is.na(x))})
  nobs <- apply(mat,2, tapply, pop(x), function(x){sum(!is.na(x))})
  n <- nmissing + nobs
  sum <- apply(mat,2, tapply, pop(x), function(x){sum(x,na.rm=TRUE)})
  f <- apply(mat,2, tapply, pop(x), function(x){mean(x, na.rm=TRUE)/2})
  f <- f*100

# Convert the matricies to long format pop, locus, value
  nobs <- reshape2::melt(nobs, na.rm=FALSE)
  nmissing <- reshape2::melt(nmissing, na.rm=FALSE)
  n <- reshape2::melt(n, na.rm=FALSE)
  f <- reshape2::melt(f, na.rm=FALSE)
  sum <- reshape2::melt(sum, na.rm=FALSE)
  
  if(nPop(x) == 1) {
    m <- cbind(levels(pop(x)),rownames(sum),sum,nobs,nmissing,f,n)
  } else {
    m <- cbind(sum,nobs[,3],nmissing[,3],f[,3],n[,3])
  }  
  colnames(m) <- c("popn","locus","sum","nobs","nmissing","frequency","n")
  
# FLAG SCRIPT END

  if (verbose >= 1) {
    cat("Completed:",funname,"\n")
  }
  
  return(m)
  
  }
}
