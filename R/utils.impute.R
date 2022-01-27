matrix2gen <- function(snp_matrix, parallel) {
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
sample_alleles <- function(alleles_list = c("a", "A"), q) {
  if (is.na(q)) {
    return(NA)
  }
  alleles_sampled <-
    paste0(sample(
      alleles_list,
      size = 2,
      prob = c(q, 1 - q),
      replace = T
    ), collapse = "")
  alleles_sampled
  if (alleles_sampled == "AA") {
    alleles_sampled <- 0
  }
  if (alleles_sampled == "aA") {
    alleles_sampled <- 1
  }
  if (alleles_sampled == "Aa") {
    alleles_sampled <- 1
  }
  if (alleles_sampled == "aa") {
    alleles_sampled <- 2
  }
  return(alleles_sampled)
}

# function to sample genotypes based on Hardy-Weinberg equation
sample_genotype <- function(genotype_list = c(0, 1, 2), q) {
  if (is.na(q)) {
    return(NA)
  }
  #genotype probabilities based on Hardy-Weinberg equation
  # p^2 + 2pq + q^2 = 1
  geno_probs <- c(((1 - q) ^ 2), # homozygote for the reference allele
                  (2 * (1 - q) * q), # heterozygote
                  (q ^ 2)) # homozygote for the alternative allele)
  genotype_sampled <-
    sample(genotype_list, size = 1, prob = geno_probs)
  return(genotype_sampled)
}