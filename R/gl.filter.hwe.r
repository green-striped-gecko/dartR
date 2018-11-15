#' Filters loci that show significant departure from Hardy-Weinberg Equilibrium
#' 
#' Calculates the probabilities of agreement with H-W equilibrium based on observed
#' frequencies of reference homozygotes, heterozygotes and alternate homozygotes. 
#' Uses the exact calculations contained in function prob.hwe() as developed by
#' Wigginton, JE, Cutler, DJ, and Abecasis, GR.
#' 
#' Input is a genlight {adegenet} object containing SNP genotypes (0 homozygous for reference SNP, 
#' 1 heterozygous, 2 homozygous for alternate SNP).
#' 
#' Loci are filtered if they show HWE departure in any one population.
#' Note that power to detect departures from HWE is affected by sample size and that
#' effective filtering may require substantial sample sizes (n > 20).
#' 
#' @param gl -- a genlight object containing the SNP genotypes [Required]
#' @param p -- level of significance (per locus) [Default 0.05]
#' @param basis -- basis for filtering out loci (any, HWE departure in any one population) [default basis="any"]
#' @param bon -- apply bonferroni correction to significance levels for filtering [default TRUE]  
#' @return a genlight object with the loci departing significantly from HWE removed
#' @author Arthur Georges (bugs? Post to https://groups.google.com/d/forum/dartr)
#' @export
#' @examples
#' list <- gl.filter.hwe(testset.gl,0.05, bon=TRUE)

gl.filter.hwe <- function(gl, p=0.05, basis="any", bon=TRUE) {
  
  x <- gl
  
  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required!\n"); stop()
  }
  
  if (basis == "any") {
    # Split the gl object into a list of populations
    poplist <- seppop(x)
  } else {
    cat("Fatal Error: basis parameter must be \"any\" or \"To be added\"\n")
    stop()
  }

# If HWE is to be calculated on one population
#  if(basis=="all"){
#    result <- utils.hwe(gl2, prob=p)
# else if HWE is to be calculated on each population    
#  } else {
    count <- 0
    for (i in poplist) {
      count <- count + 1
      if (count==1) {
        result <- utils.hwe(i, prob=p)
        Population <- rep(names(poplist)[count],nrow(result))
        result <- cbind(Population,result)
      } else {
        r <- utils.hwe(i, prob=p)
        Population <- rep(names(poplist)[count],nrow(r))
        r <- cbind(Population,r)
        result <- rbind(result, r)
      }
    }
#  }  
  rprob <- as.numeric(as.character(result$Prob))
  if (bon==TRUE) {
    result <- result[result$Bonsig=="*",]
  } else {
    result <- result[(rprob>0 & rprob<=p),]
  }
  failed.loci <- as.character(unique(result$Locus))
  cat("Loci examined:", nLoc(x),"\n")
  cat("  Deleted",length(failed.loci),"loci for significant departure from HWE\n")
  x <- x[,!locNames(x) %in% failed.loci]
  cat("  Loci retained:",nLoc(x))

  # Return the result
  return(x) 
}
