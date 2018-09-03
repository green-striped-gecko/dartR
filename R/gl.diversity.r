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

dat <- list(toto=c(1,0,1,NA), titi=c(0,0,1,NA), tata=c(0,2,2, NA), tete=c(0,2,1,NA))
x2 <- new("genlight", dat, ploidy=2)
x2
pop(x2)<- factor(c("A","A","B","B"))


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

two_H_alpha <- lapply(pops, function(x) {
  p<- colMeans(as.matrix(x), na.rm = T)/2
  p <- p[!is.na(p)] #ignore loci with just missing data
  return(mean(1-(p*p+(1-p)*(1-p)))  )
})

zero_H_beta <- NA
one_H_beta <- NA
two_H_beta <- NA

npops <- length(pops)
pairs <- t(combn(npops,2))
#zero_H_beta
zero_H_beta <- apply(pairs,1, function(x)  {
  
  pop1 <- pops[[x[1]]]
  pop2 <- pops[[x[2]]]
  pp1 <- colMeans(as.matrix(pop1), na.rm = T)/2
  
  
  pp1 <- ifelse(pp1>0 & pp1<1, 0.5 , pp1)
  pp2 <- colMeans(as.matrix(pop2), na.rm = T)/2
  pp2 <- ifelse(pp2>0 & pp2<1, 0.5 , pp2)
  
  index <- !is.na(pp1) & !is.na(pp2)
  pp1 <- pp1[index]
  pp2 <- pp2[index]
  return(mean(abs(pp1-pp2)))
  
} )

mat_zero_H_beta <- matrix(NA, nrow = npops, ncol = npops)
mat_zero_H_beta[lower.tri(mat_zero_H_beta)] <- zero_H_beta

  
return(list(zero_H_alpha= zero_H_alpha, one_H_alpha=one_H_alpha,two_H_alpha=two_H_alpha ,zero_H_beta=mat_zero_H_beta))



}



