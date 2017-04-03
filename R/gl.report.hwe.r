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
#' @param gl -- a genlight object containing the SNP genotypes [Required]
#' @param p -- level of significance (per locus) [Default 0.05]
#' @param pop -- populations to combine in the analysis [Default "all"] 
#' @return a dataframe containing loci, counts of reference SNP homozygotes, heterozygotes
#' and alternate SNP homozygotes; probability of departure from H-W equilibrium,
#' and per locus significance.
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @export
#' @examples
#' gl.report.hwe(testset.gl,0.05, pop=c("EmmacMaclGeor", "EmmacCoopCully"))
#
# Amended 23-Oct-16

gl.report.hwe <- function(gl, p=0.05, pop="all") {
x <- gl
  
  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required!\n"); stop()
  }
    
  if (pop[1] != "all") {
    gl <- x[pop(x) %in% pop]
    gl <- gl.filter.monomorphs(gl)
  } else {
    gl <- x
  }  
  m <- as.matrix(gl)
    
# Intialize arrays
  hom.ref <- array(data=NA, dim=ncol(m))
  hom.snp <- array(data=NA, dim=ncol(m))
  het <- array(data=NA, dim=ncol(m))
  columns <- colnames(m)
  p.values <- array(data=NA, dim=ncol(gl))
  sig <- array(data=NA, dim=ncol(gl))
  
# For each locus, calculate counts of homozygotes and heterozygotes
# then probability of departure from H-W equilibrium.
  for (i in 1:ncol(gl)) {
    hom.ref[i] <- length(which(m[,i]==0))
    hom.snp[i] <- length(which(m[,i]==2))
    het[i] <- length(which(m[,i]==1))
    p.values[i] <- prob.hwe(het[i], hom.ref[i], hom.snp[i])
  }
# Locus-wise significance level 
  a <- p
# Population-wise significance level using Bonferoni Correction
  a.xpw <- p/ncol(m)
# Determine significance
  sig <- p.values < a  # per locus
  sig2 <- p.values < a.xpw  # experiment-wide
# Assemble results into a dataframe
  result <- data.frame(cbind(hom.ref, het, hom.snp, p.values))
  result <- cbind(columns, result, sig)
  names(result) <- c("Locus", "Hom_1", "Het", "Hom_2", "Prob", "Sig")

  # Report the results
  pc <- sum(sig,na.rm=TRUE)*100/length(sig)
  pc <- round(pc, digits=1)
  pc2 <- sum(sig2,na.rm=TRUE)*100/length(sig2)
  pc2 <- round(pc2, digits=1)
  cat("\nANALYSIS SUMMARY\n\n")
  cat("  Populations:",pop,"\n")
  cat(paste("  Number of loci examined:",ncol(as.matrix(x)),"\n"))
  cat(paste("  Number of polymorphic loci examined:",ncol(gl),"\n"))
  cat(paste("  Polymorphic loci that depart significantly from HWe:",pc,"%\n"))
  cat(paste("  Polymoprhic loci that depart significantly from HWe (Bonferroni Corrected):",pc2,"%\n"))

  # Return the result
  return(result) 
}
