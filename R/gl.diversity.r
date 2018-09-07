#' Calculate diversity indices for SNPs
#'
#' !!Just an intro placeholder!! This script takes a genlight object and calculates alpha and beta diversity for q=0:2. Formulas are taken from Sherwin et al. 2017. The paper describes nicely the relationship between the different q levels and how they relate to population genetic processes such as dispersal and selection. For all indices the entropies (H) and corrosponding effective numbers Hill numbers (D), which reflect the amount of entities that are needed to get the observed valuea are calculated. In a nutshell the alpha indices between the different q-values should be similar if there are no deviation from expected allele frequencies and occurrences (e.g. all loci in HWE & equilibrium). If there is a deviation of an index this links to a process causing it such as dispersal, selection or strong drift. For a detailed explanation of all the indices, we recommend to resort to the literature provided below.
#'
#' @param gl genlight object containing the SNP genotypes  [required]
#' @param spectrumplot switch to provide a plot [TRUE]
#' 
#' @return a list of entropy indices for each level of q and equivalent numbers for alpha and beta diversity.
#' @export
#' @author Bernd Gruber (glbugs@@aerg.canberra.edu.au)
#' @references 
#' Sherwin, W.B., Chao, A., Johst, L., Smouse, P.E. (2017). Information Theory Broadens the Spectrum of Molecular Ecology and Evolution. TREE 32(12) 948-963. doi:10.1016/j.tree.2017.09.12
#' 
#' Chao et al. 2014
 

### To be done:

# adjust calculation of betas for population sizes (switch)


gl.diversity <- function(gl, spectrumplot=TRUE) {

if  (is.null(pop(gl)))  pop(gl)<- factor(rep("all", nInd(gl)))

#split in pops
pops <- seppop(gl)
     
  ##0Halpha (average number of alleles, ignoring missing values)

nlocpop <- lapply(pops, function(x) sum(!is.na(colMeans(as.matrix(x)))))
nlocpop$totalloc <- nLoc(gl)

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

mat_zero_H_beta <- NA
mat_one_H_beta <- NA
mat_two_H_beta <- NA

npops <- length(pops)

if (npops>1)
{
pairs <- t(combn(npops,2))

#missing loci 
#zero_H_beta
nlocpairpop <- apply(pairs,1, function(x)  {
  pop1 <- pops[[x[1]]]
  pop2 <- pops[[x[2]]]
  pp1 <- colMeans(as.matrix(pop1), na.rm = T)/2
  pp2 <- colMeans(as.matrix(pop2), na.rm = T)/2
  index <- !is.na(pp1) & !is.na(pp2)
  return(sum(index))
  
} )



mat_nloc_pops <- matrix(NA, nrow = npops, ncol = npops)
mat_nloc_pops[lower.tri(mat_nloc_pops)] <- nlocpairpop
colnames(mat_nloc_pops) <- rownames(mat_nloc_pops) <- names(pops)





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
colnames(mat_zero_H_beta) <- rownames(mat_zero_H_beta) <- names(pops)


#one_H_beta
# calculate one_H_alpha_all for combined pops
p<- colMeans(as.matrix(gl), na.rm = T)/2
p <- p[!is.na(p)] #ignore loci with just missing data
logp <- ifelse(!is.finite(log(p)), 0, log(p))
log1_p <- ifelse(!is.finite(log(1-p)), 0, log(1-p))
one_H_beta_all <- -mean(p*logp + (1-p)*log1_p)


one_H_beta <- apply(pairs,1, function(x) {
  return(one_H_beta_all-mean(c(one_H_alpha[[x[1]]], one_H_alpha[[x[2]]])))
})
mat_one_H_beta <- matrix(NA, nrow = npops, ncol = npops)
mat_one_H_beta[lower.tri(mat_one_H_beta)] <- one_H_beta
colnames(mat_one_H_beta) <- rownames(mat_one_H_beta) <- names(pops)

#two_H_beta
# calculate two_H_alpha_all for combined pops
p<- colMeans(as.matrix(gl), na.rm = T)/2
p <- p[!is.na(p)] #ignore loci with just missing data
two_H_beta_all <- mean(1-(p*p+(1-p)*(1-p)))  

two_H_beta <- apply(pairs,1, function(x) {
  m2Ha <- mean(c(two_H_alpha[[x[1]]], two_H_alpha[[x[2]]]))
#Johst-D
    return(  ((two_H_beta_all-m2Ha)/(1-m2Ha))*(npops/(npops-1)))
  #GST
 # return((two_H_beta_all-m2Ha)/two_H_beta_all  )
  
  
})
  mat_two_H_beta <- matrix(NA, nrow = npops, ncol = npops)
  mat_two_H_beta[lower.tri(mat_two_H_beta)] <- two_H_beta
  colnames(mat_two_H_beta) <- rownames(mat_two_H_beta) <- names(pops)

  mat_zero_D_beta <- mat_zero_H_beta+1
  mat_one_D_beta <-  exp(mat_one_H_beta)
  mat_two_D_beta <-  exp(mat_two_H_beta)

} # npops>1  

#Calculate hill numbers
zero_D_alpha <- unlist(zero_H_alpha)+1
one_D_alpha <- exp(unlist(one_H_alpha))
two_D_alpha <- 1/(1-unlist(two_H_alpha))




if (spectrumplot)
{

fs <- cbind(zero_D_alpha, one_D_alpha, two_D_alpha) 
cx <- max(1-(max(-12+nrow(fs),0)*0.025),0.5)
bb<-  barplot(fs, beside = T, names.arg = rep(rownames(fs),3), ylim=c(1,2.2), main="q-profile", col=rainbow(npops), las=2, xpd=FALSE, cex.names=cx)
 text(colMeans(bb), rep(2.1,3), labels = c("q=0", "q=1", "q=2")) 
}
return(list(nlocpop=unlist(nlocpop), nlocpairpop = mat_nloc_pops,
  zero_H_alpha= unlist(zero_H_alpha), one_H_alpha=unlist(one_H_alpha),two_H_alpha=unlist(two_H_alpha) ,
  zero_D_alpha= zero_D_alpha, one_D_alpha=one_D_alpha,two_D_alpha=two_D_alpha ,
  zero_H_beta=mat_zero_H_beta, one_H_beta=mat_one_H_beta, two_H_beta=mat_two_H_beta, zero_D_beta = mat_zero_D_beta, one_D_beta = mat_one_D_beta, two_D_beta = mat_two_D_beta))
}



