library(dartR)
library(vcfR)
library(pegas)

source("gl2fasta_fix.r")
test_4 <- testset.gl
test_4 <- gl.filter.callrate(test_4,threshold = 1)
gl2fasta(test_4, method=4, outfile="test.fasta",verbose=3,outpath = getwd())
test_fasta <- fasta2genlight("test.fasta")
test_fasta$other$loc.metrics <- as.data.frame(matrix(nrow = nLoc(test_fasta),ncol = nInd(test_fasta)))
test_fasta <- utils.reset.flags(test_fasta)
pop(test_fasta) <- as.factor(rep(1,nInd(test_fasta)))
res_1 <- get_tajima_D(test_fasta)

test_pop <- testset.gl
pop(test_pop) <- as.factor(rep(1,nInd(test_pop)))
res_2 <- get_tajima_D(test_pop)


vcf_1 <- gl.read.vcf("mdp_genotype.vcf")
vcf_2 <- vcf_1[,which(vcf_1$chromosome=="1")]
res_1 <- get_tajima_D(test_fasta)

vcf_2 <- read.vcfR("mdp_genotype_chr1.vcf")
getCHROM(vcf_2)
DNAbin_1 <-  vcfR2DNAbin(vcf_2,consensus = TRUE , extract.haps = FALSE)
seg.sites(DNAbin_1)
res_2 <- tajima.test(DNAbin_1)


./vcftools --vcf in.vcf --out tajimasd --TajimaD 10000

./bgzip /Users/s441489/OneDrive\ -\ University\ of\ Canberra/dartR_dev/Data_Example_DiploidUnphased.vcf
./tabix -p vcf /Users/s441489/OneDrive\ -\ University\ of\ Canberra/dartR_dev/Data_Example_DiploidUnphased.vcf.gz


##### PopGenome
vcf_2 <- read.big.fasta("test.DnaSP.fas")
vcf_2 <- readVCF("2_5.vcf.gz",numcols = 135,tid = "un",frompos = 1,topos = 682719)
GENOME.class <- neutrality.stats(vcf_2)
GENOME.class@Tajima.D

##### PopGenome
# 1.254947
##### DnaSP
# 1.254947

### strataG

test_stratag <- read.fasta("test.fasta")
res_stratag <- tajimasD(test_stratag, CI = 0.95)


# 
# 25.	Compress file using htslib
# 25.1.	cd htslib	
# 25.2.	./bgzip ~/filtered_variants_output.vcf
# 25.3.	Output: filtered_variants.vcf.gz
# 26.	Index file using htslib
# 26.1.	cd htslib
# 26.2.	./tabix -p vcf ~/filtered_variants.vcf.gz	
# 26.3.	Ouput: filtered_variants.vcf.gz.tbi

x <- seppop(testset.gl)
x[[1]]

gl2fasta(vcf_1,outpath=getwd())

fasta1 <- fasta2DNAbin("output.fasta")
fasta1


res <- get_tajima_D(testset.gl)

source("gl2genind.R")
library(hierfstat)
library(rlist)

data(dolph.seqs)

tajimasD(dolph.seqs)

test <- get_tajima_D(testset.gl)

getLocusNames_2

x<- testset.gl
gtype <- genlight2gtypes_2(x)
res_tajimas <- tajimasD(gtype, CI = 0.95)

x <- testset.gl
x2 <- genlight2gtypes_2(x)
# x3 <- genind2gtypes(x2)
res_tajimas <- tajimasD_2(x2, CI = 0.95)

genlight2gtypes_2 <- function(x){
  if (!inherits(x, "genlight"))
    stop("'x' must be a genlight object")
  genotypes <- list(c("A", "A"), c("A", "G"), c("G", "G"),c(NA,NA))
  x2 <- as.data.frame(x)
  x2[is.na(x2)] <- 3
  
  gen.mat <- list.cbind(lapply(x2, function(num.alt) {
    do.call(rbind, genotypes[num.alt + 1])
  })[])
  loci <- x@loc.names
  # if (is.null(loci)) 
  #   loci <- paste0("L", 1:adegenet::nLoc(x))
  colnames(gen.mat) <- paste0(rep(loci, each = 2), ".", 1:2)
  # has.pop <- !is.null(x@pop)
  gen.mat <- cbind(id = as.character(1:nrow(gen.mat)), strata = as.character(x@pop), as.data.frame(gen.mat))
  df2gtypes(x = gen.mat, ploidy = 2, other = list(genind = adegenet::other(x)))
}

tajimasD_2 <- function (x, CI = 0.95) 
{
  x.to.D <- function(x, Dmin, Dmax) x * (Dmax - Dmin) + Dmin
  beta.D <- function(tajima.D, alpha, beta, Dmin, Dmax) {
    gamma(alpha + beta) * (Dmax - tajima.D)^(alpha - 1) * 
      (tajima.D - Dmin)^(beta - 1)/(gamma(alpha) * gamma(beta) * 
                                      (Dmax - Dmin)^(alpha + beta - 1))
  }
  # x <- if (inherits(x, "gtypes")) {
  #   getSequences(x, as.haplotypes = FALSE, as.multidna = TRUE)
  # }
  # else {
  #   as.multidna(x)
  # }
  x <- getSequences_2(x, simplify = FALSE)
  purrr::map(x, function(dna) {
    dna <- as.matrix(dna)
    num.sites <- ncol(dna)
    pws.diff <- ape::dist.dna(dna, model = "N", pairwise.deletion = TRUE, 
                              as.matrix = TRUE)
    pi <- mean(pws.diff[lower.tri(pws.diff)])
    S <- ncol(variableSites(dna)$site.freqs)
    n <- nrow(dna)
    n.vec <- 1:(n - 1)
    a1 <- sum(1/n.vec)
    a2 <- sum(1/n.vec^2)
    b1 <- (n + 1)/(3 * (n - 1))
    b2 <- 2 * (n^2 + n + 3)/(9 * n * (n - 1))
    c1 <- b1 - 1/a1
    c2 <- b2 - (n + 2)/(a1 * n) + a2/a1^2
    e1 <- c1/a1
    e2 <- c2/(a1^2 + a2)
    D_obs <- (pi - S/a1)/sqrt(e1 * S + e2 * S * (S - 1))
    if (is.nan(D_obs)) {
      warning("D cannot be computed (division by zero in final equation)")
      return(c(D = NA, p.value = NA))
    }
    DMin <- (2/n - 1/a1)/e2^(1/2)
    DMax <- if (n%%2 == 0) {
      (n/(2 * (n - 1)) - 1/a1)/e2^(1/2)
    }
    else {
      ((n + 1)/(2 * n) - 1/a1)/e2^(1/2)
    }
    Alpha <- -(1 + DMin * DMax) * DMax/(DMax - DMin)
    Beta <- (1 + DMin * DMax) * DMin/(DMax - DMin)
    lci.p <- (1 - CI)/2
    LCI <- x.to.D(stats::qbeta(lci.p, Beta, Alpha), DMin, 
                  DMax)
    UCI <- x.to.D(stats::qbeta(1 - lci.p, Beta, Alpha), DMin, 
                  DMax)
    tryCatch({
      prob <- stats::integrate(beta.D, lower = ifelse(D_obs < 
                                                        0, DMin, D_obs), upper = ifelse(D_obs < 0, D_obs, 
                                                                                        DMax), alpha = Alpha, beta = Beta, Dmin = DMin, 
                               Dmax = DMax)
      tibble::tibble(D = D_obs, p.value = prob$value, LCI = LCI, 
                     UCI = UCI)
    }, error = function(e) {
      warning("error in Tajima's D integration, NA returned")
      tibble::tibble(D = NA, p.value = NA, LCI = NA, UCI = NA)
    })
  }) %>% dplyr::bind_rows() %>% dplyr::mutate(locus = names(x)) %>% 
    as.data.frame()
}


getSequences_2 <- function(x, loci = NULL, ids = NULL, simplify = TRUE,
                           exclude.gap.only = TRUE, ...) {
  # check that object isn't empty
  # if(is.null(x@dna)) {
  #   warning("'x' is empty. NULL returned.", call. = FALSE)
  #   return(NULL)
  # }
  
  if(is.character(loci) | is.null(loci)) {
    loci <- 1: dplyr::n_distinct(x@data$locus)
    #.checkLocusNames(x, loci)
  } else if(is.numeric(loci) | is.logical(loci)) {
    loci <- getLocusNames_2(x)[loci]
  } else stop("'loci' must be a character, numeric, or logical")
  
  # loop through loci
  new.dna <- sapply(loci, function(this.locus) {
    # extract this DNAbin object
    dna <- as.list(x@data[this.locus])
    if(exclude.gap.only) dna <- dna[!.isGapOnly(dna)]
    # return sequences for IDs which are present
    locus.ids <- .checkIDs(dna, ids)
    if(is.null(locus.ids)) NULL else dna[locus.ids]
  }, simplify = T)
  new.dna <- new.dna[!sapply(new.dna, is.null)]
  
  if(length(new.dna) == 1 & simplify) new.dna[[1]] else new.dna
}


#
## Internal functions: not documented, not exported
##
## Thibaut Jombart, April 2015
##

## compute the number of missing sequences for each gene in a list of DNAbin matrices or phyDat objects 'x'
.nMissingSequences <- function(x){
  ## only keep non-empty matrices
  x <- lapply(x, as.character)
  x <- x[sapply(x, nrow)>0]
  out <- sapply(x, function(e) sum(apply(e=="-",1,all)))
  return(sum(out,na.rm=TRUE))
}


# Compare requested locus names with actual names in multidna object
.checkLocusNames <- function(x, loci = NULL) {
  # check that locus names can be found
  if(!is.null(loci)) {
    missing <- setdiff(loci, getLocusNames_2(x))
    if(length(missing) > 0) {
      missing <- paste(missing, collapse = ", ")
      warning(paste("The following loci could not be found:", missing), call. = FALSE)
    }
    loci <- intersect(loci, getLocusNames_2(x))
  } else if(is.logical(loci) | is.numeric(loci)) {
    loci <- getLocusNames_2(x)[loci]
  }
  # set to all locus names if none specified
  if(is.null(loci)) loci <- getLocusNames_2(x)
  # return NULL if no locus names match
  if(length(loci) == 0) return(NULL)
  loci
}

# Compare requested ids with actual ids in DNAbin object
.checkIDs <- function(dna, ids = NULL) {
  # check that id names can be found
  if(!is.null(ids)) {
    missing <- setdiff(ids, labels(dna))
    if(length(missing) > 0) {
      missing <- paste(missing, collapse = ", ")
      warning(paste("The following ids could not be found:", missing), call. = FALSE)
    }
    ids <- intersect(ids, labels(dna))
  } else if(is.logical(ids) | is.numeric(ids)) {
    ids <- labels(dna)[ids]
  }
  # set to all locus names if none specified
  if(is.null(ids)) ids <- labels(dna)
  # return NULL if no locus names match
  if(length(ids) == 0) return(NULL)
  ids
}


# Return a vector of logicals denoting if each DNAbin sequence contains only gaps
.isGapOnly <- function(dna) {
  sapply(as.character(as.list(dna)), function(this.seq) all(this.seq == "-"))
}
# the same for phyDat objects
.isGapOnlyPhyDat <- function(dna, gap="-") {
  allLevels <- attr(dna, "allLevels")
  ind <- match(gap, allLevels)
  sapply(dna, function(this.seq) all(this.seq == ind))
}

getLocusNames_2 <- function(x, ...){getLociNames(x)}

