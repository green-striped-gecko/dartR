matrix2gen <- function(snp_matrix, parallel = FALSE) {
  if (parallel) {
    i@gen <-
      parallel::mclapply(1:nrow(snp_matrix), function(i)
        new("SNPbin", as.integer(snp_matrix[i, ])), mc.silent = TRUE, mc.cleanup =
          TRUE, mc.preschedule = FALSE)
  } else {
    lapply(1:nrow(snp_matrix), function(i)
      new("SNPbin", as.integer(snp_matrix[i, ])))
  }
}

# function to sample alleles using allele frequencies as probability
s_alleles <- function(q_freq) {
  if (is.na(q_freq)) {
    return(NA)
  }
  alleles_sampled <-
    paste0(sample(
      c("a", "A"),
      size = 2,
      prob = c(q_freq, 1 - q_freq),
      replace = T
    ), collapse = "")
  
  if (alleles_sampled == "AA") {
    alleles_sam <- 0
  }
  
  if (alleles_sampled == "aA") {
    alleles_sam <- 1
  }
  
  if (alleles_sampled == "Aa") {
    alleles_sam <- 1
  }
  
  if (alleles_sampled == "aa") {
    alleles_sam <- 2
  }
  
  return(as.numeric(alleles_sam))
}

# function to sample genotypes based on Hardy-Weinberg equation
sample_genotype <- function(genotype_list = c(0, 1, 2), q_freq) {
  if (is.na(q_freq)) {
    return(NA)
  }
  #genotype probabilities based on Hardy-Weinberg equation
  # p^2 + 2pq + q^2 = 1
  geno_probs <- c(((1 - q_freq) ^ 2), # homozygote for the reference allele
                  (2 * (1 - q_freq) * q_freq), # heterozygote
                  (q_freq ^ 2)) # homozygote for the alternative allele)
  genotype_sampled <-
    sample(genotype_list, size = 1, prob = geno_probs)
  return(genotype_sampled)
}