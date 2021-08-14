#' @name gl.assign.pa
#' @title Eliminate populations as possible source populations for an individual of unknown provenance, using private alleles
#' @description
#' This script eliminates from consideration as putative source populations, those populations for which the individual
#' has too many private alleles. The populations that remain are putative source populations, subject to further 
#' consideration. 
#' 
#' The algorithm identifies those target populations for which the individual has no private alleles or for which the
#' number of private alleles does not exceed a user specified threshold.
#' 
#' An excessive count of private alleles is an indication that the unknown does not belong to a target population (provided that
#' the sample size is adequate, say >=10). 
#' 
#' @param x -- name of the input genlight object [required]
#' @param unknown -- SpecimenID label (indName) of the focal individual whose provenance is unknown [required]
#' @param nmin -- minimum sample size for a target population to be included in the analysis [default 10]
#' @param threshold -- populations to retain for consideration; those for which the focal individual has 
#' less than or equal to threshold loci with private alleles [default 0]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return A genlight object containing the focal individual (assigned to population "unknown") and #' populations for which the focal individual is not distinctive (number of loci with private alleles less than or equal to thresold t.
#' 
#' @importFrom stats dnorm qnorm
#' @export
#' 
#' @author Custodian: Arthur Georges -- Post to \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' # Test run with a focal individual from the Macleay River (EmmacMaclGeor)
#'   x <- gl.assign.pa(testset.gl, unknown="UC_00146", nmin=10, threshold=1)
#'   
#' @seealso \code{\link{gl.assign.pca}}, \code{\link{gl.assign.mahandist}}

gl.assign.pa <- function (x, unknown, nmin=10, threshold=0, verbose=3) {

# SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
# FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func=funname,build="Jackson",v=verbose)
  
# CHECK DATATYPE 
  datatype <- utils.check.datatype(x,verbose=verbose)

# FUNCTION SPECIFIC ERROR CHECKING

  test <- unknown%in%indNames(x)
  if (!all(test,na.rm=FALSE)) {
    stop(error("Fatal Error: nominated focal individual (of unknown provenance) is not present in the dataset!\n"))
  }
  if (nmin <=0 & verbose >= 1){
    cat(warn("  Warning: the minimum size of the target population must be greater than zero, set to 10\n"))
    nmin <- 10
  }
  if (threshold < 0 & verbose >= 1){
    cat(warn("  Warning: the threshold for private alleles must be non-negative, set to 0\n"))
    threshold <- 0
  }

# DO THE JOB
  
  # Set a hard recommended minimum population size
  hard.min <- 10
  if (verbose >= 1 & (nmin < hard.min)){
    cat(warn("  Warning: Minimum sample size selected is less than 10 individuals\n"))
    cat(warn("    Risk of alleles present in the unknown being missed during sampling of populations with sample sizes less than 10\n"))
  }
  
  # Assign the unknown individual to population 'unknown'
  vec <- as.vector(pop(x))
  vec[indNames(x)==unknown] <- "unknown"
  pop(x) <- as.factor(vec)
  
  # Split the genlight object into one containing the unknown and one containing the remainding populations  
  unknowns <- x[pop(x)=="unknown",]
  knowns <- x[pop(x)!="unknown",]
  
  # Remove all known populations with less than nmin individuals
  pop.keep <- levels(pop(knowns))[table(pop(knowns)) >= nmin]
  # if (verbose >=2) {
  #   cat(report("  Retaining",length(pop.keep),"populations with sample size >=",
  #              nmin,":\n",paste(pop.keep,collapse=", "),"\n"))
  # }
  pop.toss <- levels(pop(knowns))[table(pop(knowns)) < nmin]
  if (verbose >=2) {
    cat(report("  Discarding",length(pop.toss),"populations with sample size <",nmin,":\n"))
        if(verbose >=3){cat(paste(pop.toss,collapse=", "),"\n")}
  }
  
  knowns <- knowns[pop(knowns) %in% pop.keep]
  
  # Warn of populations retained with sizes less than the hard wired miniumum    
  pop.warn <- levels(pop(knowns))[table(pop(knowns)) < hard.min]
  if (length(pop.warn >= 1)) {
    if (verbose >=1) {cat(warn("  Warning: Some retained populations have sample sizes less than",hard.min,":",pop.warn,"\n"))}
    if (verbose >=1) {cat(warn("    Substantial risk of private alleles arising as sampling error\n\n"))}
  }  
  
  # Report number of focal individuals (1) and number of target populations  
  n <- length(pop(unknowns))
  N <- length(levels(pop(knowns)))
  if (n != 1) {
    if (verbose >=0) {
      stop(error("Fatal Error: Number of unknown focal individuals > 1; population label 'unknown' already in use\n"))
    }  
  }
  if (verbose >=2) {cat(report("  Assigning",n,"unknown individual",unknown,"to",N,"target populations using",nLoc(x),"loci\n"))}
  
  
  # CALCULATE NUMBER OF LOCI WITH PRIVATE ALLELES
  
  # For each unknown individual
  for (i in 1:n) {
    # Grab the genotype of the unknown individual i
    unknown.ind <- as.matrix(unknowns[i])
    # Compare with each population
    count <- rep(0, N)
    count.NA <- rep(0, N)
    
    for (j in 1:N) {
      # Grab the genotypes for the known population j    
      known.pop <- as.matrix(knowns[pop(knowns)==levels(pop(knowns))[j]])
      
      # For each locus, count the number of non-private alleles  
      for (k in 1:nLoc(x)) {
        # Check to see if the unknown genotype is missing at locus k
        if (is.na(unknown.ind[k])) {
          count.NA[j] <- count.NA[j] + 1
        } else {  
          # Count the number of private alleles for each of unknown genotype 0,1,2
          # Where the unknown focal individual is homozygous reference allele [aa]
          if (unknown.ind[k]==0) {
            # If all individuals in the target population are bb [not aa or ab] then focal individual has private allele [a]
            if ( all(known.pop[,k]==2, na.rm=TRUE) )  { count[j] <- count[j] + 1}
          }
          # Where the unknown focal individual is homozygous for the alternate allele [bb]
          if (unknown.ind[k]==2) {
            # If all individuals in the target population are aa [not bb or ab] then focal individual has private allele [b]
            if ( all(known.pop[,k]==0, na.rm=TRUE) ) { count[j] <- count[j] + 1}
          }
          # Where the unknown focal individual is heterozgous [ab]
          if (unknown.ind[k]==1) {
            # If all individuals in the target population are aa, then [b] is private or if bb, then [a] is private
            if ( (all(known.pop[,k]==0, na.rm=TRUE) ) || (all(known.pop[,k]==2, na.rm=TRUE) ) ) { count[j] <- count[j] + 1}
          }
        }
      }
    }
    
    # Print out results 
    
    if (verbose >=3) {
      # cat(report("  Unknown individual:",unknown,"\n"))
      # cat(report("  Total number of SNP loci: ",nLoc(x),"\n"))
      cat("  Table showing populations against number of loci with private alleles\n")
      for (m in levels(as.factor(count))) {
        #cat("m=",m,"\n")
        if (as.numeric(as.character(m)) > threshold){
          cat(paste0("  >",threshold,"---"))
        }
        cat("  ",m,levels(pop(knowns))[count==m],"\n")
      } 
      #cat("\n")
    }
  }    
  # Save the data in a new gl object
  
  index <- ((pop(x) %in% levels(pop(knowns))[count<=threshold]) | (as.character(pop(x)) == "unknown"))
  gl <- x[index,]
  gl <- gl.filter.monomorphs(gl, verbose=0)

# Check that there is more than one population to assign (excluding 'unknown')
  if (verbose >= 2) {
  if (nPop(gl)==1) {
    if(verbose >= 2)cat(report("  There are no populations retained for assignment. Conclude that the unknown does not belong to one of the target populations or rerun with a higher threshold\n"))
  # } else if (nPop(gl)==2) {
  #   return(cat(report("  There are no further populations to compare for assignment.",levels(pop(gl))[1],"is the best assignment\n")))
  } else {
    if(verbose >= 2){cat(report("  Identified and retained",nPop(gl)-1,"putative source populations for",unknown,"\n"))}
  }
  }
# FLAG SCRIPT END

  if (verbose > 0) {
    cat(report("Completed:",funname,"\n"))
  }

  return(df)
}
