utils.read.ped <- function (file, 
                            n, 
                            snps, 
                            which, 
                            split = "\t| +",
                            sep = ".", 
                            na.strings = "0", 
                            lex.order = FALSE,
                            show_warnings=T){
  r0 <- as.raw(0)
  r1 <- as.raw(1)
  r2 <- as.raw(2)
  r3 <- as.raw(3)
  con <- gzfile(file)
  open(con)
  if (missing(n)) {
    n <- 0
    repeat {
      line <- readLines(con, n = 1)
      if (length(line) == 0) 
        break
      n <- n + 1
    }
    if (n == 0) 
      stop("Nothing read")
    seek(con, 0)
  }
  gen <- missing(snps)
  map <- NULL
  if (!gen) {
    m <- length(snps)
    if (m == 1) {
      map <- read.table(snps, comment.char = "")
      m <- nrow(map)
      if (missing(which)) {
        which <- 1
        repeat {
          snps <- map[, which]
          if (!any(duplicated(snps))) 
            break
          if (which == ncol(map)) 
            stop("No unambiguous snp names found on file")
          which <- which + 1
        }
      } else {
        snps <- map[, which]
      }
    }
  } else {
    line <- readLines(con, n = 1)
    fields <- strsplit(line, split)[[1]]
    nf <- length(fields)
    if (nf%%2 != 0) 
      stop("Odd number of fields")
    m <- (nf - 6)/2
    seek(con, 0)
  }
  nf <- 6 + 2 * m
  result <- matrix(raw(n * m), nrow = n)
  ped <- character(n)
  mem <- character(n)
  pa <- character(n)
  ma <- character(n)
  sex <- numeric(n)
  aff <- numeric(n)
  rownms <- character(n)
  a1 <- a2 <- rep(NA, m)
  a1m <- a2m <- rep(TRUE, m)
  mallelic <- rep(FALSE, m)
  for (i in 1:n) {
    line <- readLines(con, n = 1)
    fields <- strsplit(line, "\t| +")[[1]]
    to.na <- fields %in% na.strings
    fields[to.na] <- NA
    ped[i] <- fields[1]
    mem[i] <- fields[2]
    pa[i] <- fields[3]
    ma[i] <- fields[4]
    sex[i] <- suppressWarnings(as.numeric(fields[5]))
    aff[i] <- as.numeric(fields[6])
    alleles <- matrix(fields[7:nf], byrow = TRUE, ncol = 2)
    one <- two <- rep(FALSE, m)
    for (k in 1:2) {
      ak <- alleles[, k]
      akm <- is.na(ak)
      br1 <- !akm & a1m
      a1[br1] <- ak[br1]
      a1m[br1] <- FALSE
      br2 <- !akm & (a1 == ak)
      one[br2] <- TRUE
      br3 <- !akm & !a1m & (a1 != ak)
      br4 <- br3 & a2m
      a2[br4] <- ak[br4]
      a2m[br4] <- FALSE
      br5 <- br3 & (a2 == ak)
      two[br5] <- TRUE
      mallelic <- mallelic | !(akm | one | two)
    }
    gt <- rep(r0, m)
    gt[one & !two] <- r1
    gt[one & two] <- r2
    gt[two & !one] <- r3
    result[i, ] <- gt
  }
  close(con)
  if (any(a1m & show_warnings==T)){ 
    warning("no data for ", sum(a1m), " loci")
  }
  mono <- (a2m & !a1m)
  if (any(mono & show_warnings==T)){ 
    warning(sum(mono), " loci were monomorphic")
  }
  if (any(mallelic & show_warnings==T)) {
    result[, mallelic] <- r0
    warning(sum(mallelic), " loci were multi-allelic --- set to NA")
  }
  if (gen){ 
    snps <- paste("locus", 1:m, sep = sep)
  }
  if (any(duplicated(ped))) {
    if (any(duplicated(mem))) {
      rnames <- paste(ped, mem, sep = sep)
      if (any(duplicated(rnames))){ 
        stop("could not create unique subject identifiers")
      }
    } else{
      rnames <- mem
    }
  }  else {
    rnames <- ped
  }
  dimnames(result) <- list(rnames, snps)
  result <- new("SnpMatrix", result)
  if (lex.order) {
    swa <- (!(is.na(a1) | is.na(a2)) & (a1 > a2))
    switch.alleles(result, swa)
    a1n <- a1
    a1n[swa] <- a2[swa]
    a2[swa] <- a1[swa]
    a1 <- a1n
  }
  fam <- data.frame(row.names = rnames, pedigree = ped, member = mem, 
                    father = pa, mother = ma, sex = sex, affected = aff)
  if (is.null(map)){ 
    map <- data.frame(row.names = snps, snp.name = snps, 
                      allele.1 = a1, allele.2 = a2)
  } else {
    map$allele.1 <- a1
    map$allele.2 <- a2
    names(map)[which] <- "snp.names"
  }
  list(genotypes = result, fam = fam, map = map)
}