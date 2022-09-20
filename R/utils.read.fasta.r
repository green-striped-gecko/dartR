

utils.read.fasta <-  function(file,
                              parallel = parallel,
                              n_cores = NULL) {
  ## find length of a genome
  NLOC <-
    nchar(scan(
      file,
      what = "character",
      sep = "\n",
      quiet = TRUE,
      skip = 1,
      nmax = 1
    ))
  
  POOL <- as.list(rep("-", NLOC))
  
  txt <- scan(file,
              what = "character",
              sep = "\n",
              quiet = TRUE)
  
  nb.ind <- length(grep("^>", txt))
  # find individuals' labels
  IND.LAB <- sub(">", "", txt[grep("^>", txt)])
  # split per individuals
  txt <- split(txt, rep(1:nb.ind, each = 2))
  if (parallel) {
    # each genome -> one vector
    txt <-
      parallel::mclapply(txt, function(e)
        strsplit(paste(e[-1], collapse = ""), split = ""),
        mc.cores = n_cores, mc.silent = TRUE, mc.cleanup =
          TRUE, mc.preschedule = FALSE)
  } else {
    # each genome -> one vector
    txt <-
      lapply(txt, function(e)
        strsplit(paste(e[-1], collapse = ""), split = ""))
  }
  
  ## POOL contains all alleles of each position
  # alleles current genomes
  temp <-
    as.list(apply(matrix(
      unlist(txt), byrow = TRUE, nrow = length(txt)
    ), 2, unique))
  # update global pool
  POOL <-
    mapply(function(x, y)
      unique(c(x, y)), POOL, temp, SIMPLIFY = FALSE)
  
  
  ## analyse pool of alleles
  letterOK <-
    c("a", "c", "t", "g", "A", "G", "C", "T", "M", "R", "W", "S", "Y", "K")
  # keep only proper letters
  POOL <-
    lapply(POOL, function(e)
      e[e %in% letterOK])
  ## POOL <- lapply(POOL, setdiff, "-")
  nb.alleles <- lengths(POOL)
  snp.posi <- which(nb.alleles == 2 | nb.alleles == 3)
  if (length(snp.posi) == 0) {
    cat(warn("  No polymorphism in the alignment - returning empty object./n"))
    return(new("genlight"))
  }
  
  txt <- scan(file,
              what = "character",
              sep = "\n",
              quiet = TRUE)
  
  ## read SNPs
  nb.ind <- length(grep("^>", txt))
  # split per individuals
  txt <- split(txt, rep(1:nb.ind, each = 2))
  if (parallel) {
    # each genome -> one SNP vector
    txt <-
      parallel::mclapply(txt, function(e)
        strsplit(paste(e[-1], collapse = ""), split = "")[[1]][snp.posi],
        mc.cores = n_cores, mc.silent = TRUE, mc.cleanup =
          TRUE, mc.preschedule = FALSE)
  } else {
    # each genome -> one SNP vector
    txt <-
      lapply(txt, function(e)
        strsplit(paste(e[-1], collapse = ""), split = "")[[1]][snp.posi])
  }
  
  txt2 <- lapply(txt, paste0, collapse = "")
  txt2 <- stringr::str_replace_all(txt2, "A", "1")
  txt2 <-  stringr::str_replace_all(txt2, "T", "2")
  txt2 <-  stringr::str_replace_all(txt2, "G", "3")
  txt2 <-  stringr::str_replace_all(txt2, "C", "4")
  txt2 <-  stringr::str_replace_all(txt2, "M", "5")
  txt2 <-  stringr::str_replace_all(txt2, "R", "5")
  txt2 <-  stringr::str_replace_all(txt2, "W", "5")
  txt2 <-  stringr::str_replace_all(txt2, "S", "5")
  txt2 <-  stringr::str_replace_all(txt2, "Y", "5")
  txt2 <-  stringr::str_replace_all(txt2, "K", "5")
  txt3 <- cbind(txt2, txt2)
  
  # make genotypes
  #to hack package checking...
  make_geno <- function() {
    
  }
  
  Rcpp::cppFunction(
    plugins = "cpp11",
    
    'List make_geno(StringMatrix mat) {
    int ind = mat.nrow();
    int loc = strlen(mat(0,0));
    List out(ind);
for (int i = 0; i < ind; i++) {
 std::string chr1 (mat(i,0));
 std::string chr2 (mat(i,1));
 StringVector temp(loc);
for (int j = 0; j < loc; j++) {
 StringVector geno = StringVector::create(chr1[j],chr2[j]);
    temp[j] = collapse(geno);
  }
      out[i] = temp;
    }
    return out;
  }'
  )
  
  plink_ped <- make_geno(txt3)
  plink_ped_2 <- lapply(plink_ped, function(x) {
    x[x == "55"] <- 1
    x[x == "11"] <- 0
    x[x == "22"] <- 2
    x[x == "33"] <- 0
    x[x == "44"] <- 2
    
    return(x)
    
  })
  
  res <- list()
  
  txt <-
    lapply(plink_ped_2, function(e)
      suppressWarnings(as.integer(e)))
  
  res <-
    c(res, lapply(txt, function(e)
      new(
        "SNPbin", snp = e, ploidy = 2L
      )))
  
  res <- new("genlight", res, ploidy = 2L)
  
  indNames(res) <- IND.LAB
  alleles(res) <- rep("C/G", nLoc(res))
  locNames(res) <-
    paste0(sub("\\..*", "", basename(file)), "_", snp.posi)
  return(res)
  
}

merge_gl_fasta <- function(gl_list, parallel = FALSE) {
  matrix_temp <- lapply(gl_list, function(y) {
    return(as.data.frame(cbind(ind_names = indNames(y), as.matrix(y))))
  })
  
  gl_temp <-
    Reduce(function(x, y) {
      merge(x, y, by = "ind_names", all = TRUE)
    }, matrix_temp)
  
  res_temp <-  matrix2gen(gl_temp[, 2:ncol(gl_temp)], parallel = parallel)
  
  res <- new("genlight", res_temp, ploidy = 2L)
  
  res_final <- gl.compliance.check(res, verbose = 0)
  
  return(res_final)
  
}
