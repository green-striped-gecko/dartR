#' Simulates individuals based on the allele frequencies provided via a genlight object.
#'
#' This function simulates individuals based on the allele frequencies of a genlight object. The output is a genlight object with the same number of loci as the input genlight object. 
#'
#' @param gl -- name of the genlight object containing the SNP data 
#' @param n -- number of individuals that should be simulated
#' @param popname -- a population name for the simulated individuals [default Null]
#' @return a genlight object with n individuals. 
#' @details 
#' The function can be used to simulate populations for sampling designs or for power analysis. Check the example below where the effect of drift is explored, by simply simulating several generation a genlight object and putting in the allele frequencies of the previous generation. The beauty of the function is, that it is lightning fast.
#' 
#' @export
#' @author Bernd Gruber (bernd.gruber@@canberra.edu.au)
#' @examples
#' glsim <- gl.sim.ind(testset.gl, n=10, popname="sims")
#' glsim
#' ###Simulate drift over 10 generation 
#' assuming a bottleneck of only 10 individuals
#' [ignoring mating and mutation]
#' #Simulate 20 individuals with no structure
#' founder <- glSim(n.ind = 20, n.snp.nonstruc = 50, ploidy=2)
#' #number of fixed loci in the first generation
#' 
#' res <- sum(colMeans(as.matrix(founder), na.rm=TRUE) %%2 ==0)
#' simgl <- founder
#' #49 generations of only 10 individuals
#' for (i in 2:50) 
#' {
#'    simgl <- gl.sim.ind(simgl, n=10, popname="sims")
#'    res[i]<- sum(colMeans(as.matrix(simgl), na.rm=TRUE) %%2 ==0)
#' }
#' plot(1:50, res, type="b", xlab="generation", ylab="# fixed loci")
#' 
#' 

gl.sim.ind <- function(gl, n=50, popname=NULL)
{
  #allelefequency of the population
  p <- as.matrix(gl)
  alf <- colMeans(p, na.rm = T)/2
  alfinds <- matrix(rep(alf,n), nrow=n, byrow=T)
  simind<- apply(alfinds,c(1,2), function(x) sample(0:2,size = 1, prob =  c((1-x)^2, 2*x*(1-x),x^2)))
  #now create genlight objects.....
  
  glsim <- new("genlight", gen=simind, ploidy=2, ind.names=1:n, loc.names=locNames(gl), loc.all=gl@loc.all, position=position(gl), pop=rep(popname,n))
  return(glsim)
}

