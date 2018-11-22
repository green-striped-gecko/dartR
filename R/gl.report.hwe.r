#' Reports departure from Hardy-Weinberg Equilibrium
#' 
#' Calculates the probabilities of agreement with H-W equilibrium based on observed
#' frequencies of reference homozygotes, heterozygotes and alternate homozygotes. 
#' Uses the exact calculations contained in function prob.hwe() as developed by
#' Wigginton, JE, Cutler, DJ, and Abecasis, GR.
#' 
#' #Input is a genlight {adegenet} object containing SNP genotypes (0 homozygous for reference SNP, 
#' #1 heterozygous, 2 homozygous for alternate SNP).
#' 
#' Tests are applied to all populations pooled (subset="all"), to each population treated separately
#' (subset="each") or to selected populations (subset=c("pop1","pop2")). Tests for Hwe are
#' only valid if there is no population substructure, and the tests have sufficient power
#' only when there is sufficient sample size (n>20).
#' 
#' @param gl -- a genlight object containing the SNP genotypes [Required]
#' @param p -- level of significance (per locus) [Default 0.05]
#' @param subset -- list populations to combine in the analysis | each | all [Default "all"] 
#' @return a dataframe containing loci, counts of reference SNP homozygotes, heterozygotes
#' and alternate SNP homozygotes; probability of departure from H-W equilibrium,
#' and per locus significance with and without Bonferroni Correction.
#' @author Arthur Georges (bugs? Post to \url{https://groups.google.com/d/forum/dartr})
#' @export
#' @examples
#' list <- gl.report.hwe(testset.gl,0.05, subset=c("EmmacMaclGeor", "EmmacCoopCully"))
#' list <- gl.report.hwe(testset.gl,0.05, subset="all")
#' list <- gl.report.hwe(testset.gl,0.05, subset="each")

gl.report.hwe <- function(gl, p=0.05, subset="each") {
  
  x <- gl
  flag <- 0
  
  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required!\n"); stop()
  }
  
  if (subset[1] == "all") {
    gl2 <- x
    flag <- 1
  } else if (subset[1] == "each") {
    # Split the gl object into a list of populations
    poplist <- seppop(x)
  } else if(nPop(x[pop(x) %in% subset])){
      flag <- 1
      gl2 <- x[pop(x) %in% subset]
      #gl2 <- gl.filter.monomorphs(gl2)
  } else {
    cat("Fatal Error: subset parameter must be \"each\", \"all\", or a list of populations\n")
    stop()
  }

# If HWE is to be calculated on one population
  if(flag==1){
    result <- utils.hwe(gl2, prob=p)
# else if HWE is to be calculated on each population    
  } else {
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
  }  
  rprob <- as.numeric(as.character(result$Prob))
  result <- result[(rprob>0 & rprob<=p),]
  result <- result[order(result$Locus),]
  cat("Reporting significant departures from Hardy-Weinberg Equilibrium\n")
  if (nrow(result)==0){
    cat("No significant departures\n")
  } else {
    cat("NB: Departures significant at the alpha level of",p,"are listed\n")
    if (p > 0.05) {
      cat("ns --",p,"< p < 0.05; * -- 0.05 < p < 0.01; ** -- 0.01 < p < 0.001; etc\n")
    } else {
      cat("ns -- p > 0.05; * -- 0.05 < p < 0.01; ** -- 0.01 < p < 0.001; etc\n")
    }
      cat("Critical values for significance of Bonferroni Corrected significance vary with sample size\n\n")
    print(result, row.names=FALSE)
  }  

  # Return the result
  return(result) 
}
