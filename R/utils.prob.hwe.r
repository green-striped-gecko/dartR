#' Exact SNP test of Hardy-Weinberg Equilibrium
#'
#' This code calculates an exact probability of departure from Hardy-Weinberg Equilibrium 
#' as described in Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of 
#' Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76:887-893.
#' 
#' Source code available from http://csg.sph.umich.edu/abecasis/Exact/r_instruct.html
#' 
#' Note: return code of -1.0 signals an error condition; return code of NA signals that 
#' all alleles are NA for a locus
#' 
#' @param obs_hets -- count of heterozygotes by locus
#' @param obs_hom1 -- count of homozygotes, reference state
#' @param obs_hom2 -- count of homozygotes, alternate state
#' @return Exact probability of agreement with HWE
#' @author Custodian: Luis Mijangos (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' hets <- 20
#' hom_1 <- 5
#' hom_2 <- 30
#' #p_value <- prob.hwe(hets, hom_1, hom_2)

utils.prob.hwe <- function(obs_hets, obs_hom1, obs_hom2) {
  
  if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0 ) {
    return(-1.0)
  }
  
  # Added by Arthur Georges Aug 8 2015 to avoid crash on prob.hwe(0,0,0)
  if (obs_hets==0 && obs_hom1==0 && obs_hom2 == 0) {
    return(-2.0)
  }
  
  # total number of genotypes
  N <- obs_hom1 + obs_hom2 + obs_hets
  
  # rarer homozygotes, more common homozygotes
  obs_homr <- min(obs_hom1, obs_hom2)
  obs_homc <- max(obs_hom1, obs_hom2)
  
  # number of rarer allele copies
  rare  <- obs_homr*2 + obs_hets
  
  # Initialize probability array
  probs <- rep(0, 1+rare)
  
  # Find midpoint of the distribution
  mid <- floor(rare * ( 2 * N - rare) / (2 * N))
  if ( (mid %% 2) != (rare %% 2) ) mid <- mid + 1
  
  probs[mid + 1] <- 1.0
  mysum <- 1.0
  
  # Calculate probabilities from midpoint down 
  curr_hets <- mid
  curr_homr <- (rare - mid) / 2
  curr_homc <- N - curr_hets - curr_homr
  
  while ( curr_hets >=  2) {
    probs[curr_hets - 1]  <- probs[curr_hets + 1] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0)  * (curr_homc + 1.0))
    mysum <- mysum + probs[curr_hets - 1]
    
    # 2 fewer heterozygotes -> add 1 rare homozygote, 1 common homozygote
    curr_hets <- curr_hets - 2
    curr_homr <- curr_homr + 1
    curr_homc <- curr_homc + 1
  }    
  
  # Calculate probabilities from midpoint up
  curr_hets <- mid
  curr_homr <- (rare - mid) / 2
  curr_homc <- N - curr_hets - curr_homr
  
  while ( curr_hets <= rare - 2) {
    probs[curr_hets + 3] <- probs[curr_hets + 1] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
    mysum <- mysum + probs[curr_hets + 3]
    
    # add 2 heterozygotes -> subtract 1 rare homozygote, 1 common homozygote
    curr_hets <- curr_hets + 2
    curr_homr <- curr_homr - 1
    curr_homc <- curr_homc - 1
  }    
  
  # P-value calculation
  target <- probs[obs_hets + 1]

  p <- min(1.0, sum(probs[probs <= target])/ mysum)
  
  return(p)
}
