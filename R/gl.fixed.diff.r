#' Generate a matrix of fixed differences from a genelight or genind object \{adegenet\}
#'
#' This script takes SNP data grouped into populations in a genlight object (DArTSeq)
#' and generates a matrix of fixed differences between populations taken pairwise
#'
#' A fixed difference at a locus occurs when two populations share no alleles. The challenge with this approach
#' is that when sample sizes are finite, fixed differences will occur through sampling error, compounded when
#' many loci are examined. Simulations suggest that sample sizes of n1=5 and n2=5 is adequate to reduce the
#' probability of [experiment-wide] type 1 error to negligible levels [ploidy=2]. A warning is issued if comparison
#' between two populations involves sample sizes less than 5, taking into account allele drop-out. The minimum sample 
#' size for scoring fixed differences between two populations can be set with the parameter nlimit.
#' An absolute fixed difference is as defined above. However, one might wish to score fixed differences at some lower
#' level of allele frequency difference, say where percent allele fequencies are 95,5 and 5,95 rather than 100:0 and 0:100.
#' This adjustment can be done with the tloc parameter. For example, tloc=0.05 means that SNP allele frequencies of 
#' 95,5 and 5,95 percent will be regarded as fixed when comparing two populations at a locus.
#'
#' @param gl -- name of the genlight object containing SNP genotypes or a genind object containing presence/absence data [required]
#' @param tloc -- threshold defining a fixed difference (e.g. 0.05 implies 95:5 vs 5:95 is fixed) [default 0]
#' @param nlimit -- number of individuals with non-missing SNP scores in the two populations combined, required for 
#' scoring of fixed differences [default 2]
#' @param pc -- logical value indicating whether to report fixed difference counts (pc=FALSE) or percentages (pc=TRUE) [default TRUE]
#' @param v -- verbosity = 0, silent; 1, brief; 2, verbose [default 1]
#' @return Matrix of percent fixed differences (lower matrix), number of loci (upper matrix)
#' @import adegenet utils
#' @export
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' #only used the first 20 individuals due to runtime reasons 
#' mat <- gl.fixed.diff(testset.gl[1:20,], tloc=0.05, nlimit=5)
#' @seealso \code{\link{is.fixed}}

gl.fixed.diff <- function(gl, tloc=0, nlimit=2, pc=FALSE, v=1) {
x <- gl

# Checking parameter values

  if(class(gl) == "genlight") {
    cat("Analysing a genlight object\n")
  } else {
    cat("Fatal Error: Specify a genlight object\n")
    stop()
  }
  if (tloc < 0 ) {
    cat("Error: Parameter tloc should be positive in the range 0 to 0.5, reset to zero\n")
    tloc <- 0
  }
  if (tloc > 0.5 ) {
    cat("Error: Parameter tloc should be positive in the range 0 to 0.5, reset to 0.5\n")
    tloc <- 0.5
  }
  if (nlimit < 2) {
    cat("Error: Minimum value for nlimit is 2 individuals, reset to 2\n")
    nlimit=2
  }
  if (v < 0 || v > 2) {
    cat("Error: Verbosity must be set to one of 0, silent; 1, brief; 2, verbose. Reset to default of 1\n")
    v=1
  }
  # Checking for and removing monomorphic loci
    if (v > 0 && pc==TRUE) {
      x2 <- gl.filter.monomorphs(x,v=0)
      if (nLoc(x2) < nLoc(x)) {
        cat("Warning: Monomorphic loci are present and will be used in percentage fixed difference calculations\n")
      }
      rm(x2)
    }

  # Calculate percent allele frequencies
    gl.mat.sum <- gl.percent.freq(x)

  # GENERATE A MATRIX OF PAIRWISE FIXED DIFFERENCES
    
  # Report samples sizes for each population
    if (v > 0){
      cat("Populations, aggregations and sample sizes")
      print(table(pop(x)))
      if (min(table(pop(x))) < 5 ){
        cat("Warning: fixed differences can arise through sampling error if sample sizes are small\n")
        cat("  Some sample sizes are small (N < 5)\n")
      }
      cat("\n")
      cat("Calculating pairwise fixed differences\n")
    }  
   
  # Establish an array to hold the fixed differences and sample sizes
    npops<-nlevels(gl.mat.sum$popn)
    nloci<-nlevels(gl.mat.sum$locus)
    fixed <- array(-1, c(npops, npops))
    loc.count <- array(-1, c(npops, npops))

    # Set up the progress counter

    if (v > 0){
      pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
      getTxtProgressBar(pb)
    }
    # Cycle through the data to sum the fixed differences into a square matrix
    flag <- 0
    for (i in 1:nloci) {                           # For each locus
      countj<-0
      startj <- (i-1)*(npops)+1
      endj <- i*npops
      for (j in startj:endj){    # For each population against that locus
        n1 <- gl.mat.sum$nobs[j]
        countj<-countj+1
        countk<-0
        startk <- (i-1)*(npops)+1
        endk <- i*npops
        for (k in startk:endk) { # For each population, compare pairwise
          n2 <- gl.mat.sum$nobs[k]
          countk<-countk+1
          if (!is.na(n1+n2)) {
          if ((n1+n2) >= nlimit) {
          # Compare and if not missing, increment
            cf <- is.fixed(gl.mat.sum$frequency[j],gl.mat.sum$frequency[k],t=tloc)
            if (!is.na(cf)) {
              if (fixed[countj,countk] == -1) {
                fixed[countj,countk] <- cf
                loc.count[countj,countk] <- 1
              } else {
                fixed[countj,countk] <- fixed[countj,countk] + cf
                loc.count[countj,countk] <- loc.count[countj,countk] + 1
              }
            }
          }
          }
        }
      }

      if (v > 0){setTxtProgressBar(pb, i/nloci)}
    }

  # Cycle through the populations to determine sample sizes
    
    ind.count <- array(-1, c(npops, npops))
    flag <- 0
    t  <- table(pop(x))
    for (i in 1:(length(t)-1)) {
      for (j in i:length(t)) {
        ind.count[i,j] <- t[i] + t[j]
        if (t[i] < 5 | t[j] < 5) { 
        flag <- 1
        }
      }
    }
    if (flag == 1 && v > 0) {

      cat("\n   Warning: Some comparisons involve sample sizes were less than 5.\n")
      cat("   Compounded Type I error rate may be high. Consider a priori amalgamation.\n")
    }

  # Convert missing values to NA
    fixed[fixed == -1] <- NA
    loc.count[loc.count == -1] <- NA
  # Convert to percentages if requested  
    if ( pc ) {fixed <- round(fixed*100/loc.count,4)}
  # Tidy up adding row and column names
    rownames(fixed)<-levels(gl.mat.sum$popn)
    colnames(fixed)<-levels(gl.mat.sum$popn)
    fixed <- fixed[order(rownames(fixed)), order(colnames(fixed))]
    rownames(loc.count)<-levels(gl.mat.sum$popn)
    colnames(loc.count)<-levels(gl.mat.sum$popn)
    loc.count <- loc.count[order(rownames(loc.count)), order(colnames(loc.count))]
    
    # Put the SNP locus counts in the upper matrix, percent fixed differences in the lower matrix

    if (npops == 1) {
      if (v > 0){cat("All populations amalgamated into one\n")}
      fixed = 0
    } else {
      for (i in 1:npops-1) {
        for (j in (i+1):npops) {
          fixed[i,j] <- round(loc.count[i,j],0)
        }
      }
    }
    
  # Return the matrix
    return(fixed)
}
