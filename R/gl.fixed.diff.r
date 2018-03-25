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
#'
#' An absolute fixed difference is as defined above. However, one might wish to score fixed differences at some lower
#' level of allele frequency difference, say where percent allele fequencies are 95,5 and 5,95 rather than 100:0 and 0:100.
#' This adjustment can be done with the tloc parameter. For example, tloc=0.05 means that SNP allele frequencies of 
#' 95,5 and 5,95 percent will be regarded as fixed when comparing two populations at a locus.
#'
#' @param x -- name of the genlight object containing SNP genotypes [required]
#' @param tloc -- threshold defining a fixed difference (e.g. 0.05 implies 95:5 vs 5:95 is fixed) [default 0]
#' @param test -- if TRUE, calculate p values for the observed fixed differences [default FALSE]
#' @param reps -- number of replications to undertake in the simulation to estimate probability of false positives [default 1000]
#' @param delta -- the threshold value for the minor allele frequency to regard the difference between two populations to be fixed [default 0.02]
#' @param rm.global.monomorphs -- if TRUE, loci that are monomorphic across all individuals are removed before computation [default TRUE]
#' @param pb -- if TRUE, show a progress bar on time consuming loops [default FALSE]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A list containing the gl object x and the following square matricies
#'         [[1]] $gl -- the input genlight object;,
#'         [[2]] $fd -- raw fixed differences;,
#'         [[3]] $pcfd -- percent fixed differences;,
#'         [[4]] $nobs -- mean no. of individuals used in each comparison;,
#'         [[5]] $nloc -- total number of loci used in each comparison;,
#'         [[6]] $expobs -- if test=TRUE, the expected count of false positives for each comparison [by simulation],
#'         [[7]] $prob -- if test=TRUE, the significance of the count of fixed differences [by simulation])
#' @import adegenet utils
#' @export
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' fd <- gl.fixed.diff(testset.gl, tloc=0, test=TRUE, delta=0.02, reps=100, v=1 )
#' @seealso \code{\link{is.fixed}}

gl.fixed.diff <- function(x, tloc=0, test=FALSE, delta=0.02, reps=1000, rm.global.monomorphs=TRUE, pb=TRUE, v=2) {

  if (v > 0) {
    if (tloc > 0) {cat("Starting gl.fixed.diff: Comparing populations for fixed differences with tolerance",tloc,"\n")}
    if (tloc == 0) {cat("Starting gl.fixed.diff: Comparing populations for absolute fixed differences\n")}
  }
  
# Checking parameter values

  if(!(class(x) == "genlight")) {
    cat("Fatal Error: Specify a genlight object\n")
    stop("Execution terminated\n")
  }
  if (tloc > 0.5 || tloc < 0 ) {
    cat("Fatal Error: Parameter tloc should be positive in the range 0 to 0.5, reset to default of 0\n")
    stop("Execution terminated\n")
  }
  v <- as.integer(v)
  if (v < 0 || v > 5){
    cat("  Fatal Error: Parameter v must be between 0 and 5\n"); stop("Execution terminated\n")
  }
  
  # Checking count of populations
    if(nPop(x) < 2) {
      cat("Fatal Error: Distance calculation requires at least two populations, one or none present\n")
      stop("Execution terminated\n")
    }

  # Checking for and removing monomorphic loci
      x2 <- gl.filter.monomorphs(x,v=0)
      if (nLoc(x2) < nLoc(x)) {
        if(!rm.global.monomorphs) {
          if (v > 0) {cat("  Warning: Globally monomorphic loci retained, used in calculations\n")}
        } else {  
          if (v > 1) {cat("  Globally monomorphic loci removed\n")}
        x <- x2  
        }  
      }
      rm(x2)

  # Calculate percent allele frequencies
    ftable <- gl.percent.freq(x, v=v)

  # GENERATE A MATRIX OF PAIRWISE FIXED DIFFERENCES
    
  # Report samples sizes for each population
    if (v > 2){
      cat("Populations, aggregations and sample sizes")
      print(table(pop(x)))
    }
    if (min(table(pop(x))) < 10 && v > 1 ){
        cat("Warning: Fixed differences can arise through sampling error if sample sizes are small\n")
        cat("  Some sample sizes are small (N < 10, minimum in dataset =",min(table(pop(x))),")\n")
        if (!test) {cat("  Recommend manually amalgamating populations or setting test=TRUE to allow evaluation of statistical significance\n")}
    }

  # Establish an array to hold the fixed differences and sample sizes
    npops<-nlevels(ftable$popn)
    nloci<-nlevels(ftable$locus)
    fixed.matrix <- array(-1, c(npops, npops))
    exp.matrix <- array(-1, c(npops, npops))
    pcfixed.matrix <- array(-1, c(npops, npops))
    ind.count.matrix <- array(-1, c(npops, npops))
    loc.count.matrix <- array(-1, c(npops, npops))
    p.false.pos.matrix <- array(-1, c(npops, npops))

    # Set up the progress counter
    if (v > 1 & pb){
      progress <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
      getTxtProgressBar(progress)
      cat("\n")
    }
    
    # Cycle through the data to sum the fixed differences into a square matrix
    if (v > 1) {
      cat("Comparing populations pairwise -- this may take time. Please be patient\n")
    }
    for (popi in 1:(npops-1)){      # For each population
      for (popj in (popi+1):npops) { # For each other population
        # Pull out the data for each population
          p1 <- ftable[ftable$popn==levels(pop(x))[popi],]
          p2 <- ftable[ftable$popn==levels(pop(x))[popj],]
        # Calculate fixed differences
          fixed <- c(rep(0,nloci))
          for (i in 1:nloci) {       # For each locus
            fixed[i] <- is.fixed(p1$frequency[i],p2$frequency[i],tloc=tloc)
          }
        # Calculate stats across loci
          fixed.matrix[popi,popj] <- sum(fixed, na.rm=TRUE)
          loc.count.matrix[popi,popj] <- sum(!is.na(fixed))
          pcfixed.matrix[popi,popj] <- round(fixed.matrix[popi,popj]*100/sum(!is.na(fixed)),0)
          ind.count.matrix[popi,popj] <- round(mean(p1$nobs[p1$nobs>0])+mean(p2$nobs[p2$nobs>0]),1)
          # Make full matrix and add row and column names    
            fixed.matrix[popj,popi] <- fixed.matrix[popi,popj]
            fixed.matrix[popi,popi] <- 0; fixed.matrix[npops,npops] <- 0
            rownames(fixed.matrix)<-levels(pop(x))
            colnames(fixed.matrix)<-levels(pop(x))
            
            loc.count.matrix[popj,popi] <- loc.count.matrix[popi,popj]
            loc.count.matrix[popi,popi] <- NA; loc.count.matrix[npops,npops] <- NA
            rownames(loc.count.matrix)<-levels(pop(x))
            colnames(loc.count.matrix)<-levels(pop(x))
            
            pcfixed.matrix[popj,popi] <- pcfixed.matrix[popi,popj]
            pcfixed.matrix[popi,popi] <- 0; pcfixed.matrix[npops,npops] <- 0
            rownames(pcfixed.matrix)<-levels(pop(x))
            colnames(pcfixed.matrix)<-levels(pop(x))
            
            ind.count.matrix[popj,popi] <- ind.count.matrix[popi,popj]
            ind.count.matrix[popi,popi] <- NA; ind.count.matrix[npops,npops] <- NA
            rownames(ind.count.matrix)<-levels(pop(x))
            colnames(ind.count.matrix)<-levels(pop(x))
          
        # Calculate the probability that the observed differences are false positives
            if (test) {
              if (v > 1) {
                cat("    Testing for significance of fixed differences -- this may take a month of Sundays. Please be very patient\n")
              }
            outlist <- gl.utils.fdsim(x,c(levels(pop(x))[popi],levels(pop(x))[popj]),obs=fixed.matrix[popi,popj],delta=delta,reps=reps,v=v)
            p.false.pos.matrix[popi,popj] <- round(outlist$prob,4)
            exp.matrix[popi,popj] <- round(outlist$mnexpected,1)
          # Make full matrix and add row and column names    
            exp.matrix[popj,popi] <- exp.matrix[popi,popj]
            exp.matrix[popi,popi] <- 0; exp.matrix[npops,npops] <- 0
            rownames(exp.matrix)<-levels(pop(x))
            colnames(exp.matrix)<-levels(pop(x))
            
            p.false.pos.matrix[popj,popi] <- p.false.pos.matrix[popi,popj]
            p.false.pos.matrix[popi,popi] <- 0; p.false.pos.matrix[npops,npops] <- 0
            rownames(p.false.pos.matrix)<-levels(pop(x))
            colnames(p.false.pos.matrix)<-levels(pop(x))
            }
      }    
      if (v > 1 & pb){setTxtProgressBar(progress, popi/(npops-1))}
    }

  # Return the matricies
    if (v > 1) {
      cat("Returning a list containing the following square matricies:\n",
      "         [[1]] $gl -- input genlight object;\n",
      "         [[2]] $fd -- raw fixed differences;\n",
      "         [[3]] $pcfd -- percent fixed differences;\n",
      "         [[4]] $nobs -- mean no. of individuals used in each comparison;\n",
      "         [[5]] $nloc -- total number of loci used in each comparison;\n",
      "         [[6]] $expobs -- if test=TRUE, the expected count of false positives for each comparison [by simulation]\n",
      "         [[7]] $prob -- if test=TRUE, the significance of the count of fixed differences [by simulation]\n")
    }
    
    l <- list(gl=x,fd=fixed.matrix,pcfd=pcfixed.matrix,nobs=ind.count.matrix,nloc=loc.count.matrix,expobs=exp.matrix,pval=p.false.pos.matrix)

    if (v > 0) {
      cat("Completed gl.fixed.diff\n\n")
    }
    return(l)
}

