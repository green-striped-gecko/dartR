#' Simulates a specified number of offsprings based on alleles provided by potential father(s) and mother(s)
#'
#' This takes a population (or a single individual) of fathers (provided as a genlight object) and mother(s) and simulates offsprings based on "random" mating. It can be used to simulate population dynamics and check the effect of those dynamics and allele frequencies, number of alleles. Another application is to simulate relatedness of siblings and compare it to actual relatedness found in the population to determine kinship.
#'
#' @param fathers -- genlight object of potential fathers
#' @param mothers -- genlight object of potential mothers simulated
#' @param noffpermother -- number of offsprings per mother
#' @param sexratio -- the sex ratio of simulated offsprings [females / females +males, 1 equals 100 percent females]
#' @return a genlight object with n individuals.
#' @importFrom stats runif 
#' @export
#' @author Bernd Gruber (bernd.gruber@@canberra.edu.au)
#' @examples
#' #Simulate 10 potential fathers
#' gl.fathers <- glSim(10, 20, ploidy=2)
#' #Simulate 10 potential mothers
#' gl.mothers <- glSim(10, 20, ploidy=2)
#' gl.sim.offspring(gl.fathers, gl.mothers, 2, sexratio=0.5)

gl.sim.offspring <- function(fathers, mothers, noffpermother, sexratio=0.5)
{
  
  noff <- nInd(mothers)*noffpermother
  mother <- sample(1:nInd(mothers),noff, replace = T)
  father <- sample(1:nInd(fathers), noff, replace = T)
  
  mmat <- as.matrix(mothers)[mother,]
  mhet <- sum(mmat==1)
  mother.half <- ifelse(mmat==1, sample(c(0,2), mhet, replace = T), mmat) 
  
  fmat <- as.matrix(fathers)[father,]
  fhet <- sum(fmat==1)
  father.half <- ifelse(fmat==1, sample(c(0,2), fhet, replace = T), fmat) 
  
  
  offmat <- (mother.half+father.half)/2
  gl2 <- new("genlight",gen=offmat,  ind.names = paste0("Po_",1:noff), loc.names=locNames(mothers), ploidy=rep(2, nrow(offmat)))
  
  #set sex ratio
  sr <- factor(ifelse(runif(nInd(gl2))<sexratio, "female","male"))
  gl2@other$sex <- sr
  return(gl2)
} 











