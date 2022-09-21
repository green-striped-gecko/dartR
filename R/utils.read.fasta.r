

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
    cat(warn("  No polymorphism in the alignment - returning empty object.\n"))
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
  
  txt <- lapply(txt,stringr::str_replace_all,"A", "AA")
  txt <- lapply(txt,stringr::str_replace_all,"T" ,"TT")
  txt <- lapply(txt,stringr::str_replace_all, "G","GG")
  txt <- lapply(txt,stringr::str_replace_all, "C","CC")
  txt <- lapply(txt,stringr::str_replace_all, "M","AC")
  txt <- lapply(txt,stringr::str_replace_all, "R","AG")
  txt <- lapply(txt,stringr::str_replace_all, "W","AT")
  txt <- lapply(txt,stringr::str_replace_all, "S","CG")
  txt <- lapply(txt,stringr::str_replace_all, "Y","CT")
  txt <- lapply(txt,stringr::str_replace_all, "K","GT")
  
  t1 <- lapply(1:length(txt[[1]]),function(z){
    t2 <- unlist(lapply(txt,"[[",z))
    t3 <- table(t2)
    return(t3)
  })
  
  t1 <- lapply(t1,function(x){
    
    drop_allele <- which(names(x) == "-" | 
            names(x)=="V" |
            names(x)=="H" |
            names(x)=="D" | 
            names(x)=="B" |
            names(x)=="N")
    
    if(length(drop_allele)>0){
      return(x[-drop_allele])
    }else{
      return(x)
    }
  })
  
  myRef <- strsplit(names(unlist(lapply(t1,which.max))),"")
  myAlt <- strsplit(names(unlist(lapply(t1,which.min))),"")
  
  ref_alt <- lapply(1:length(myRef),function(x){
    return(c(myRef[[x]],myAlt[[x]]))
  })
  
  ref_alt <- lapply(ref_alt,unique)
  
  loc.all <- unlist(lapply(ref_alt,function(x){
    paste0(x[1],"/",x[2])
          }))

  txt2 <- lapply(1:length(txt[[1]]), function(x) {
    hom_ref <- paste0(ref_alt[[x]][1],ref_alt[[x]][1])
    hom_alt <- paste0(ref_alt[[x]][2],ref_alt[[x]][2])

    res <- lapply(txt,function(y){
      if(y[x]==hom_ref){
        return(0)
      } else if(y[x]==hom_alt){
        return(2)
      } else {
        return(1)
      }
    })
  })
  
  txt3 <- as.data.frame(Reduce(cbind,txt2))
  
  txt3[] <- lapply(txt3, as.integer)
  
  res <- list()
  
  res <- c(res, apply(txt3,1, function(e)
      new(
        "SNPbin", snp = e, ploidy = 2L
      )))
  
  res <- new("genlight", res, ploidy = 2L)
  
  indNames(res) <- IND.LAB
  alleles(res) <- loc.all
  locNames(res) <-  paste0(sub("\\..*", "", basename(file)), "_", snp.posi)
  return(res)
  
}

merge_gl_fasta <- function(gl_list, parallel = FALSE) {
  
  matrix_temp <- lapply(gl_list, function(y) {
    return(as.data.frame(cbind(ind_names = indNames(y), as.matrix(y))))
  })
  
  loc_names <- lapply(matrix_temp,function(x){
  colnames(x)[-1]  
  })

  loc_names <- Reduce("c",loc_names)
  
  loc_all <- lapply(gl_list, function(y) {
    return(y$loc.all)
  })
  
  loc_all <- Reduce("c",loc_all)
  
  gl_temp <-
    Reduce(function(x, y) {
      merge(x, y, by = "ind_names", all = TRUE)
    }, matrix_temp)
  
  res_temp <-  matrix2gen(gl_temp[, 2:ncol(gl_temp)], parallel = parallel)
  
  res <- new("genlight", res_temp, ploidy = 2L)
  
  res$loc.names <- loc_names
  alleles(res) <- loc_all
  res$ind.names <- gl_temp$ind_names
  
  res_final <- gl.compliance.check(res, verbose = 0)
  
  return(res_final)
  
}
