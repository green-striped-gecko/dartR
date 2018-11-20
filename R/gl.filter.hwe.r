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
#' @param x -- a genlight object containing the SNP genotypes [Required]
#' @param alpha -- level of significance (per locus) [Default 0.05]
#' @param basis -- basis for filtering out loci (any, HWE departure in any one population) [default basis="any"]
#' @param bon -- apply bonferroni correction to significance levels for filtering [default TRUE] 
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return a genlight object with the loci departing significantly from HWE removed
#' @author Arthur Georges (Post to https://groups.google.com/d/forum/dartr)
#' @export
#' @examples
#' list <- gl.filter.hwe(testset.gl, 0.05, bon=TRUE)

gl.filter.hwe <- function(x, alpha=0.05, basis="any", bon=TRUE, v=2) {
  
# ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required!\n"); stop()
  }
  
  if (v > 0) {
    cat("Starting gl.filter.hwe: Filtering on departure from HWE\n")
  }
  
  if (alpha < 0 | alpha > 1){
    cat("    Warning: level of significance per locus alpha must be an integer between 0 and 1, set to 0.05\n")
    alpha <- 0.05
  }
  
  if (basis != "any"){
    cat("    Warning: basis of assessment must be by population, other options not yet implemented, set to \'any\'\n")
    basis <- "any"
  }
  
  if (basis == "any") {
    # Split the gl object into a list of populations
    poplist <- seppop(x)
  } 

  count <- 0
    for (i in poplist) {
      count <- count + 1
      if (count==1) {
        result <- dartR:::utils.hwe(i, prob=alpha)
        Population <- rep(names(poplist)[count],nrow(result))
        result <- cbind(Population,result)
      } else {
        r <- dartR:::utils.hwe(i, prob=alpha)
        Population <- rep(names(poplist)[count],nrow(r))
        r <- cbind(Population,r)
        result <- rbind(result, r)
      }
    }
 
  rprob <- as.numeric(as.character(result$Prob))
  if (bon==TRUE) {
    result <- result[result$Bonsig=="*",]
  } else {
    result <- result[(rprob>0 & rprob<=alpha),]
  }
  failed.loci <- as.character(unique(result$Locus))

  if (v >= 2){
    cat("Loci examined:", nLoc(x),"\n")
    if (bon) {
      cat("  Deleted",length(failed.loci),"loci with significant departure from HWE, bonferroni corrected, at experiment-wide alpha =",alpha,"\n")
    } else {
      cat("  Deleted",length(failed.loci),"loci with significant departure from HWE at alpha =",alpha,"applied locus by locus\n")
    }  
  } 
  
  x <- x[,!locNames(x) %in% failed.loci]
  
  if (v >= 2){
    cat("  Loci retained:",nLoc(x),"\n")
  }
  
  if ( v > 0) {cat("gl.filter.hwe completed\n")}

  # Return the result
  return(x) 
}
