#' Calculate diversity indices for SNPs
#'
#' This script takes a genlight object and calculates alpha and beta diversity for q=0:2 
#'
#' @param gl genlight object containing the SNP genotypes  [required]
#' 
#' @return a list of diversity indices
#' @export
#' @author Bernd Gruber (glbugs@@aerg.canberra.edu.au)

ttt <- testset.gl[c(1:7),c(2,17,111,33,144)]

gl.diversity <- function(gl) {

#split in pops
pops <- seppop(gl)  

  ##0Halpha (average number of alleles, ignoring missing values)
zero_H_alpha <- lapply(pops, function(x) mean(((colMeans(as.matrix(x), na.rm=T) %% 2)>0)+1-1, na.rm = T))


one_H_alpha <- lapply(pops, function(x) {
  p<- colMeans(as.matrix(x), na.rm = T)/2
  p <- p[!is.na(p)] #ignore loci with just missing data
  logp <- ifelse(!is.finite(log(p)), 0, log(p))
  log1_p <- ifelse(!is.finite(log(1-p)), 0, log(1-p))
  return(-mean(p*logp + (1-p)*log1_p))
})


return(list(zero_H_alpha= zero_H_alpha, one_H_alpha=one_H_alpha))
}



