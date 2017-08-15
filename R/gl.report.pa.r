#' Report number of private alleles possessed by an individual of unknown provenance
#'
#' This script calculates the number of private alleles possessed by a focal individual of unknown
#' provenance when compared to a series of target populations.
#'
#' @param gl -- name of the input genlight object [required]
#' @param id -- identity label of the focal individual whose provenance is unknown [required]
#' @param nmin -- minimum sample size for a target population to be included in the analysis [default 10]
#' @param t -- populations to retain in the output genlight object; those for which the focal individual has less than or equal to t loci with private alleles [default 0]
#' @return A genlight object containing the focal individual (assigned to population "unknown") and 
#' populations for which the focal individual is not distinctive (number of loci with private alleles less than or equal to thresold t.
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' x <- testset.gl
#' # Test run with a focal individual from the Macleay River (EmmacMaclGeor)
#' x <- gl.report.pa(testset.gl, id="UC_00146", nmin=10, t=1)
#' 

gl.report.pa <- function (gl, id, nmin=10, t=0) {
x <- gl

  cat("IDENTIFYING LOCI WITH PRIVATE ALLELES\n\n")
  
# Set a recommended minimum population size
  hard.min <- 10

# Assign the unknown individual to population 'unknown'
  v <- as.vector(pop(x))
  v[indNames(x)==id] <- "unknown"
  pop(x) <- as.factor(v)
# Split the genlight object into one containing the unknown and one containing the remainding populations  
  unknowns <- x[pop(x)=="unknown",]
  knowns <- x[pop(x)!="unknown",]
# Remove all known populations with less than nmin individuals
  pop.keep <- levels(pop(knowns))[table(pop(knowns)) >= nmin]
  pop.toss <- levels(pop(knowns))[table(pop(knowns)) < nmin]

  cat("Discarding",length(pop.toss),"populations with sample size less than",nmin,":",pop.toss,"\n\n")
  cat("Retaining",length(pop.keep),"populations with sample size greater than or equal to",nmin,":",pop.keep,"\n\n")
  knowns <- knowns[pop(knowns) %in% pop.keep]
  pop.warn <- levels(pop(knowns))[table(pop(knowns)) < hard.min]
  if (length(pop.warn >= 1)) {
    cat("Warning: Some retained populations have sample sizes less than",hard.min,":",pop.warn,"\n")
    cat("  Risk of private alleles arising as sampling error\n\n")
  }  
# Calculate the number of individuals in the unknown population and the number of target populations 
  n <- length(pop(unknowns))
  N <- length(levels(pop(knowns)))
  if (n != 1) {

    cat("Fatal Error: Number of unknown focal individuals > 1; population label 'unknown' already in use\n")
    stop("Terminating execution")
  }
  cat("Assigning",n,"unknown individual(s) to",N,"target populations\n")
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
                # Count the number of non-private alleles for each of unknown genotype 0,1,2
                  if (unknown.ind[k]==0) {
                    if ( ( any(known.pop[,k]==0, na.rm=TRUE) || any(known.pop[,k]==1, na.rm=TRUE) ) ) { count[j] <- count[j] + 1}
                  }
                  if (unknown.ind[k]==2) {
                    if ( ( any(known.pop[,k]==2, na.rm=TRUE) || any(known.pop[,k]==1, na.rm=TRUE) ) ) { count[j] <- count[j] + 1}
                  }
                  if (unknown.ind[k]==1) {
                    count[j] <- count[j] + 1
                  }
              }
          }
      }
      # Calculate the number of loci with private alleles in the focal individual
        count <- nLoc(x)-count-count.NA
      # Print out results 

        cat("Unknown individual:",id,"\n")
        cat("Total number of SNP loci: ",nLoc(x),"\n\n")
        cat("Table showing number of loci with private alleles\n")
        for (m in levels(as.factor(count))) {
          cat("  ",m,levels(pop(knowns))[count==m],"\n")
        }  
      # Save the data in a new gl object
        cat("\n")
        index <- ((pop(x) %in% levels(pop(knowns))[count<=t]) | (as.character(pop(x)) == "unknown"))
        gl <- x[index,]
        gl <- gl.filter.monomorphs(gl, v=FALSE)
        cat("Data retained for the unknown individual and remaining candidate source populations (",t,"or less loci with private alleles)\n")
  }  
  return(gl)
}


