get_tajima_D <- function(x){
  # Find allele frequencies (p1 and p2) for every locus in every population
  allele_freqs <- gl.percent.freq(x)
  names(allele_freqs)[names(allele_freqs) == "frequency"] <- "p1"
  allele_freqs$p1 <- allele_freqs$p1 / 100
  allele_freqs$p2 <- 1 - allele_freqs$p1
  
  # Get the names of all the populations
  pops <- unique(allele_freqs$popn)
  
  #split each population
  allele_freqs_by_pop <- split(allele_freqs, allele_freqs$popn)

  # Internal function to calculate pi
  calc_pi <- function(allele_freqs) {
    n = allele_freqs$nobs * 2  # vector of n values
    pi_sqr <- allele_freqs$p1 ^ 2 + allele_freqs$p2 ^ 2
    h = (n / (n - 1)) * (1 - pi_sqr) # vector of values of h
    sum(h,na.rm = T) # return pi, which is the sum of h across loci
  }
  
   get_tajima_D_for_one_pop <- function(allele_freqs_by_pop) {
    pi <- calc_pi(allele_freqs_by_pop)
    
    #Calculate number of segregating sites, ignoring missing data (missing data will not appear in the allele freq calculations)
    #S <- sum(!(allele_freqs_by_pop$p1 == 0 | allele_freqs_by_pop$p1 == 1))
    S <- sum(allele_freqs_by_pop$p1 >0 & allele_freqs_by_pop$p1 <1,na.rm = T)
    if(S == 0) {
      warning("No segregating sites")
      data.frame(pi = NaN, 
                 S = NaN, 
                 D = NaN, 
                 Pval.normal = NaN, 
                 Pval.beta = NaN)
    }
    
    n <- mean(allele_freqs_by_pop$nobs * 2 )
    
    tmp <- 1:(n - 1)
    a1 <- sum(1/tmp)
    a2 <- sum(1/tmp^2)
    b1 <- (n + 1)/(3 * (n - 1))
    b2 <- 2 * (n^2 + n + 3)/(9 * n * (n - 1))
    c1 <- b1 - 1/a1
    c2 <- b2 - (n + 2)/(a1 * n) + a2/a1^2
    e1 <- c1/a1
    e2 <- c2/(a1^2 + a2)
    
    # calculate D and do beta testing
    D <- (pi - S/a1) / sqrt(e1 * S + e2 * S * (S - 1))
    Dmin <- (2/n - 1/a1)/sqrt(e2)
    Dmax <- ((n/(2*(n - 1))) - 1/a1)/sqrt(e2)
    tmp1 <- 1 + Dmin * Dmax
    tmp2 <- Dmax - Dmin
    a <- -tmp1 * Dmax/tmp2
    b <- tmp1 * Dmin/tmp2
    p <- pbeta((D - Dmin)/tmp2, b, a)
    p <- ifelse(p < 0.5, 2 * p, 2 * (1 - p))
    
    data.frame(pi = pi, 
               S = S, 
               D = D, 
               Pval.normal = 2 * pnorm(-abs(D)), 
               Pval.beta = p)
  }
  
  output <- do.call("rbind", lapply(allele_freqs_by_pop, 
                                    get_tajima_D_for_one_pop))
  data.frame(population = rownames(output), output, row.names = NULL)
}
