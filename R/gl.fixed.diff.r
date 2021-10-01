#' @name gl.fixed.diff
#' @title Generate a matrix of fixed differences and associated statistics for populations taken pairwise
#'
#' @description
#' This script takes SNP data or sequence tag P/A data grouped into populations in a genlight object (DArTSeq)
#' and generates a matrix of fixed differences between populations taken pairwise
#'
#' @details
#' A fixed difference at a locus occurs when two populations share no alleles or where all members 
#' of one population has a sequence tag scored, and all members of the other population has 
#' the sequence tag absent. The challenge with this approach is that when sample sizes are finite, 
#' fixed differences will occur through sampling error, compounded when many loci are examined. 
#' Simulations suggest that sample sizes of n1=5 and n2=5 are adequate to reduce the probability 
#' of [experiment-wide] type 1 error to negligible levels [ploidy=2]. A warning is issued if 
#' comparison between two populations involves sample sizes less than 5, taking into account 
#' allele drop-out.
#' 
#' Optionally, if test=TRUE, the script will test the fixed differences between final OTUs for statistical significance,
#' using simulation, and then further amalgamate populations that for which there are no significant fixed differences at 
#' a specified level of significance (alpha). To avoid conflation of true fixed differences with false positives in the
#' simulations, it is necessary to decide a threshold value (delta) for extreme true allele frequencies that will be considered
#' fixed for practical purposes. That is, fixed differences in the sample set will be considered to be positives (not false positives)
#' if they arise from true allele frequencies of less than 1-delta in one or both populations.  The parameter
#' delta is typically set to be small (e.g. delta = 0.02).
#' 
#' NOTE: The above test will only be calculated if tloc=0, that is, for analyses of absolute fixed differences. The
#' test applies in comparisons of allopatric populations only. For sympatric populations, use gl.pval.sympatry().
#'
#' An absolute fixed difference is as defined above. However, one might wish to score fixed differences at some lower
#' level of allele frequency difference, say where percent allele fequencies are 95,5 and 5,95 rather than 100:0 and 0:100.
#' This adjustment can be done with the tloc parameter. For example, tloc=0.05 means that SNP allele frequencies of 
#' 95,5 and 5,95 percent will be regarded as fixed when comparing two populations at a locus.
#'
#' @param x Name of the genlight object containing SNP genotypes or tag P/A data (SilicoDArT) 
#' or an object of class 'fd' [required]
#' @param tloc Threshold defining a fixed difference (e.g. 0.05 implies 95:5 vs 5:95 is fixed) 
#' [default 0]
#' @param test If TRUE, calculate p values for the observed fixed differences [default FALSE]
#' @param reps Number of replications to undertake in the simulation to estimate probability 
#' of false positives [default 1000]
#' @param delta Threshold value for the true population minor allele frequency (MAF) from which 
#' resultant sample fixed differences are considered true positives [default 0.02]
#' @param alpha Level of significance used to display non-significant differences between
#' populations as they are compared pairwise [default 0.05]
#' @param mono.rm If TRUE, loci that are monomorphic across all individuals are removed before 
#' beginning computations [default TRUE]
#' @param pb If TRUE, show a progress bar on time consuming loops [default FALSE]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 
#' 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return A list of Class "fd" containing the gl object and square matricies, as follows
#'         [[1]] $gl -- the output genlight object;
#'         [[2]] $fd -- raw fixed differences;
#'         [[3]] $pcfd -- percent fixed differences;
#'         [[4]] $nobs -- mean no. of individuals used in each comparison;
#'         [[5]] $nloc -- total number of loci used in each comparison;
#'         [[6]] $expfpos -- if test=TRUE, the expected count of false positives for each comparison [by simulation];
#'         [[7]] $sdfpos -- if test=TRUE, the standard deviation of the count of false positives for each comparison [by simulation];
#'         [[8]] $prob -- if test=TRUE, the significance of the count of fixed differences [by simulation])
#' @import utils
#' @export
#' @author Custodian: Arthur Georges -- Post to \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \donttest{
#' fd <- gl.fixed.diff(testset.gl, tloc=0, verbose=3 )
#' fd <- gl.fixed.diff(testset.gl, tloc=0, test=TRUE, delta=0.02, reps=100, verbose=3 )
#' }
#' @seealso \code{\link{is.fixed}}

gl.fixed.diff <- function(x, 
                          tloc = 0, 
                          test = FALSE, 
                          delta = 0.02, 
                          alpha = 0.05,
                          reps = 1000, 
                          mono.rm = TRUE, 
                          pb = FALSE, 
                          verbose = NULL) {

# SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
# FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func=funname,build="Jody", verbosity =verbose)
  
# CHECK DATATYPE 
  datatype <- utils.check.datatype(x,accept=c("genlight","SNP","SilicoDArT","fd"),verbose=verbose)
  if (datatype=="fd"){
    x <- x$gl
  }
  
# FUNCTION SPECIFIC ERROR CHECKING
  
  if (tloc > 0.5 || tloc < 0 ) {
    stop(error("Fatal Error: Parameter tloc should be positive in the range 0 to 0.5\n"))
  }  
  if(!(tloc == 0) & test==TRUE){
    cat(warn("  Warning: false positives can only be simulated for tloc=0, setting tloc to zero\n"))
    tloc <- 0
  }
  
  if ( verbose >= 2){
    if (tloc > 0) {cat(report("  Comparing populations for fixed differences with tolerance",tloc,"\n"))}
    if (tloc == 0) {cat(report("  Comparing populations for absolute fixed differences\n"))}
  }
  
  # Checking count of populations
    if(nPop(x) < 2) {
      stop(error("Fatal Error: Distance calculation requires at least two populations, one or none present\n"))
    }

  # Checking for and removing monomorphic loci
    if(!x@other$loc.metrics.flags$monomorphs){
      if (verbose >= 2) {cat(warn("  Warning: Monomorphic loci retained, used in calculations\n"))}
    } else {  
      if (verbose >= 2) {cat(report("  Monomorphic loci removed\n"))}
      x <- gl.filter.monomorphs(x,verbose=0)
    }  
#define dist2list function
  
dist2list <- function (dist) 
{
  if (!class(dist) == "dist") {
    stop(error("the input data must be a dist object."))
  }
  dat <- as.data.frame(as.matrix(dist))
  if (is.null(names(dat))) {
    rownames(dat) <- paste(1:nrow(dat))
  }
  value <- stack(dat)$values
  rnames <- rownames(dat)
  namecol <- expand.grid(rnames, rnames)
  colnames(namecol) <- c("col", "row")
  res <- data.frame(namecol, value)
  return(res)
}

# DO THE JOB      

  # Calculate percent allele frequencies
    ftable <- gl.percent.freq(x, verbose=0)

  # GENERATE A MATRIX OF PAIRWISE FIXED DIFFERENCES
    
  # Report samples sizes for each population
    if (verbose >= 3){
      cat("  Populations, aggregations and sample sizes")
      print(table(pop(x)))
    }
    if (min(table(pop(x))) < 10 && verbose >= 3 ){
        cat(warn("  Warning: Fixed differences can arise through sampling error if sample sizes are small\n"))
        cat(warn("    Some sample sizes are small (N < 10, minimum in dataset =",min(table(pop(x))),")\n"))
        if (!test) {cat(warn("    Recommend manually amalgamating populations or setting test=TRUE to allow evaluation of statistical significance\n"))}
    }

  # Establish an array to hold the fixed differences and sample sizes
    npops <- nlevels(ftable$popn)
    nloci <- nlevels(as.factor(ftable$locus))
    fixed.matrix <- array(-1, c(npops, npops))
    exp.matrix <- array(NA, c(npops, npops))
    sd.matrix <- array(NA, c(npops, npops))
    pcfixed.matrix <- array(-1, c(npops, npops))
    ind.count.matrix <- array(-1, c(npops, npops))
    loc.count.matrix <- array(-1, c(npops, npops))
    p.false.pos.matrix <- array(NA, c(npops, npops))

    # Set up the progress counter
    if (verbose >= 2 & pb){
      progress <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
      getTxtProgressBar(progress)
      cat("\n\n")
    }
    
    # Cycle through the data to sum the fixed differences into a square matrix
    if (verbose >= 2) {
      cat(report("  Comparing populations pairwise -- this may take time. Please be patient\n"))
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
            if(tloc != 0){
              cat(report("  Setting tloc to zero\n"))
              tloc.hold <- tloc
              tloc <- 0
            }
            outlist <- gl.fdsim(x,c(levels(pop(x))[popi],levels(pop(x))[popj]),obs=fixed.matrix[popi,popj],delta=delta,reps=reps,verbose=0)
            p.false.pos.matrix[popi,popj] <- round(outlist$prob,4)
            exp.matrix[popi,popj] <- round(outlist$mnexpected,1)
            sd.matrix[popi,popj] <- round(outlist$sdexpected,4)
            if (verbose >= 3) {
              if (p.false.pos.matrix[popi,popj] > alpha){
              cat("    ",popNames(x)[popi],"vs",
                  popNames(x)[popj],
                  " [p =", p.false.pos.matrix[popi,popj],",ns]\n")
              }
            }
          # Make full matrix and add row and column names    
            exp.matrix[popj,popi] <- exp.matrix[popi,popj]
            exp.matrix[popi,popi] <- 0; exp.matrix[npops,npops] <- 0
            rownames(exp.matrix)<-levels(pop(x))
            colnames(exp.matrix)<-levels(pop(x))
            
            sd.matrix[popj,popi] <- sd.matrix[popi,popj]
            sd.matrix[popi,popi] <- 0; sd.matrix[npops,npops] <- 0
            rownames(sd.matrix)<-levels(pop(x))
            colnames(sd.matrix)<-levels(pop(x))
            
            p.false.pos.matrix[popj,popi] <- p.false.pos.matrix[popi,popj]
            p.false.pos.matrix[popi,popi] <- 0; p.false.pos.matrix[npops,npops] <- 0
            rownames(p.false.pos.matrix)<-levels(pop(x))
            colnames(p.false.pos.matrix)<-levels(pop(x))
            }
      }    
      if (verbose >= 2 & pb){setTxtProgressBar(progress, popi/(npops-1))}
    }
 
  # Return the extreme low distances -- candidates for lack of significance
    if (verbose >= 5) {
      cat(report("\nDisplaying a list of 5% of pairs with smallest non-zero differences\n"))
      
      df <- dist2list(as.dist(fixed.matrix))
      
      df <- df[!df$col==df$row,]
      #df <- df[!df$value<=tpop,]
      pctile <- quantile(df$value[df$value>0],probs=0.05)
      df <- df[df$value<=pctile,]
      df <- df[order(df$value),]
      print(df)
      cat("\n")
    }

  # Return the matricies
    if (verbose >= 4) {
      if(pb){cat("\n")}
      cat(report("Returning a list containing the gl object and square matricies, as follows:\n",
      "         [[1]] $gl -- input genlight object;\n",
      "         [[2]] $fd -- raw fixed differences;\n",
      "         [[3]] $pcfd -- percent fixed differences;\n",
      "         [[4]] $nobs -- mean no. of individuals used in each comparison;\n",
      "         [[5]] $nloc -- total number of loci used in each comparison;\n",
      "         [[6]] $expfpos -- if test=TRUE, the expected count of false positives for each comparison [by simulation]\n",
      "         [[7]] $sdfpos -- if test=TRUE, the expected count of false positives for each comparison [by simulation]\n",
      "         [[8]] $prob -- if test=TRUE, the significance of the count of fixed differences [by simulation]\n"))
    }

    fixed.matrix <- as.dist(fixed.matrix)
    pcfixed.matrix <- as.dist(pcfixed.matrix)
    l <- list(gl=x,fd=fixed.matrix,pcfd=pcfixed.matrix,nobs=ind.count.matrix,nloc=loc.count.matrix,expfpos=exp.matrix,sdfpos=sd.matrix,pval=p.false.pos.matrix)
    class(l) <- "fd"
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
      cat(report("Completed:",funname,"\n"))
    }
    
    return(l)
}

