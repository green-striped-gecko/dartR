#' Generate a matrix of fixed differences from a genelight or genind object \{adegenet\}
#'
#' This script takes SNP data grouped into populations in a genlight object (DArTSeq)
#' or presence absence data (SilicoDArT) grouped into populations in a genind object
#' and generates a matrix of fixed differences between populations taken pairwise
#'
#' A fixed difference at a locus occurs when two populations share no alleles. The challenge with this approach
#' is that when sample sizes are finite, fixed differences will occur through sampling error, compounded when
#' many loci are examined. Simulations suggest that sample sizes of n1=5 and n2=5 is adequate to reduce the
#' probability of [experiment-wide] type 1 error to negligible levels [ploidy=2]. A warning is issued if comparison
#' between two populations involves sample sizes less than 5, taking into account allele drop-out.
#'
#' Tollerance in the definition of a fixed difference is provided by the t parameter. For example, t=0.05 means
#' that SNP allele frequencies of 95,5 and 5,95 percent will be regarded as fixed.
#'
#' @param gl -- name of the genlight object containing SNP genotypes or a genind object containing presence/absence data [required]
#' @param t -- threshold value for tollerance in when a difference is regarded as fixed [default 0]
#' @return Matrix of percent fixed differences (lower matrix), number of loci (upper matrix)
#' @import adegenet utils
#' @export
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' mat <- gl.fixed.diff(testset.gl, t=0.05)
#' @seealso \code{\link{is.fixed}}

gl.fixed.diff <- function(gl, t=0) {
x <- gl

  # Calculate percent allele frequencies
    gl.mat.sum <- gl.percent.freq(x)

  # GENERATE A MATRIX OF PAIRWISE FIXED DIFFERENCES
    
  # Report samples sizes for each population
    cat("Sample sizes")
    print(table(pop(gl)))
    cat("NOTE: fixed differences can arise through sampling error if sample sizes are small (N < 5)\n")

   cat("Calculating pairwise fixed differences\n")
   
  # Establish an array to hold the fixed differences and sample sizes
    npops<-nlevels(gl.mat.sum$popn)
    nloci<-nlevels(gl.mat.sum$locus)
    fixed <- array(-1, c(npops, npops))
    loc.count <- array(-1, c(npops, npops))

    # Set up the progress counter
    pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
    getTxtProgressBar(pb)

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
          # Compare and if not missing, increment
            cf <- is.fixed(gl.mat.sum$frequency[j],gl.mat.sum$frequency[k],t=0)
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
      setTxtProgressBar(pb, i/nloci)
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
    if (flag == 1) {
      cat("   Warning: Some comparisons involve sample sizes were less than 5.\n")
      cat("   Compounded Type I error rate may be high. Consider a priori amalgamation.\n")
    }

  # Convert missing values to NA
    fixed[fixed == -1] <- NA
    loc.count[loc.count == -1] <- NA
  # Convert to percentages if requested  
    fixed <- round(fixed*100/loc.count,4)
  # Tidy up adding row and column names
    rownames(fixed)<-levels(gl.mat.sum$popn)
    colnames(fixed)<-levels(gl.mat.sum$popn)
    fixed <- fixed[order(rownames(fixed)), order(colnames(fixed))]
    rownames(loc.count)<-levels(gl.mat.sum$popn)
    colnames(loc.count)<-levels(gl.mat.sum$popn)
    loc.count <- loc.count[order(rownames(loc.count)), order(colnames(loc.count))]
    
    # Put the sample sizes in the upper matrix, percent fixed differences in the lower matrix

    if (npops == 1) {
      cat("All populations amalgamated into one\n")
      fixed = 0
    } else {
      for (i in 1:npops-1) {
        for (j in (i+1):npops) {
          fixed[i,j] <- ind.count[i,j]
        }
      }
    }
    
  # Return the matrix
    return(fixed)
}
