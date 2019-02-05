#' Convert a genlight object to nexus format PAUP SVDquartets
#'
#' The output nexus file contains the snp data in one of two forms, depending upon what you regard
#' as most appropriate. One form, that used by Chifman and Kubatko, has two lines per individual,
#' one providing the reference SNP the second providing the alternate SNP (method=1). 
#' A second form, recommended by Dave Swofford, has
#' a single line per individual, resolving heterozygous SNPs by replacing them with standard
#' ambiguity codes (method=2).
#'
#' Reference: Chifman, J. and L. Kubatko. 2014. Quartet inference from SNP data under the coalescent, Bioinformatics, 30: 3317-3324
#' 
#' 
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param outfile -- file name of the output file (including extension).
#' @param outpath -- path where to save the output file (set to tempdir by default)
#' @param method -- method = 1, nexus file with two lines per individual; method = 2, nexus
#' file with one line per individual, ambiguity codes [default 2]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return NULL
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl2svdquartets(testset.gl[1:100,])

gl2svdquartets <- function(x, outfile="svd.nex", outpath=tempdir(), method=2, v=2) {

  outfile <- file.path(outpath, outfile)
  
  if (v > 0) {cat(paste("Starting gl2svdquartets: Create nexus file\n\n"))}
  if (v > 1) {cat(paste("    Extacting SNP bases and creating records for each individual\n"))}

# Extract the reference base and the alternate base for each locus
  #snp <- as.character(x@other$loc.metrics$SNP)
  #snp <- strsplit(snp,":")
  #snp <- unlist(lapply(snp, `[`, 2))
  #snp <- strsplit(snp,">")
  ref <- unname(sapply(x@loc.all, function(x) strsplit(x, split = "/")[1][[1]][1]))
  #ref <- unlist(lapply(snp, `[`, 1))
  alt <-  unname(sapply(x@loc.all, function(x) strsplit(x, split = "/")[1][[1]][2]))
  #alt <- unlist(lapply(snp, `[`, 2))
  
# Sort the data on population
  df <- data.frame(as.matrix(x))
  df <- cbind(indNames(x),pop(x),df)
  df <- df[order(df$pop),]
  indlabels <- df[,1]
  poplabels <- df[,2]
  m <- df[,3:(nLoc(x)+2)]
  
  # Create a lookup table for the ambiguity codes
  #     A  T  G  C
  #  A  A  W  R  M)
  #  T  W  T  K  Y
  #  G  R  K  G  S
  #  C  M  Y  S  C  
  
  conversion <- matrix(c("A","W", "R","M", "W", "T", "K", "Y", "R", "K", "G", "S", "M", "Y", "S", "C"),nrow=4, ncol=4)
  colnames(conversion) <- c("A","T","G","C")
  rownames(conversion) <- colnames(conversion)

  # Add individual labels to sequences
  refseq <- array(data=NA,dim=nInd(x))
  altseq <- array(data=NA,dim=nInd(x))
  ambseq <- array(data=NA,dim=nInd(x))
  if (method == 1) {
    for (ind in 1:nInd(x)) {
      refseq[ind] <- paste0(indlabels[ind],"_1   ")
      altseq[ind] <- paste0(indlabels[ind],"_2   ")
    }
  } else {
    for (ind in 1:nInd(x)) {
      ambseq[ind] <- paste0(indlabels[ind],"   ")
    }
  } 

# progressively add the bases  
  for (ind in 1:nInd(x)) {
    for (loc in 2:nLoc(x)){
    if (is.na(m[ind,loc])) {
      refseq[ind] <- paste0(refseq[ind],"?")
      altseq[ind] <- paste0(altseq[ind],"?")
      ambseq[ind] <- paste0(ambseq[ind],"?")
    } else if (m[ind,loc] == 0) {
      refseq[ind] <- paste0(refseq[ind],ref[loc])
      altseq[ind] <- paste0(altseq[ind],ref[loc])
      ambseq[ind] <- paste0(ambseq[ind],ref[loc])
    } else if (m[ind,loc] == 1) {
      refseq[ind] <- paste0(refseq[ind],ref[loc])
      altseq[ind] <- paste0(altseq[ind],alt[loc])
      ambseq[ind] <- paste0(ambseq[ind],conversion[ref[loc], alt[loc]])
    } else if (m[ind,loc] == 2) {
      refseq[ind] <- paste0(refseq[ind],alt[loc])
      altseq[ind] <- paste0(altseq[ind],alt[loc])
      ambseq[ind] <- paste0(ambseq[ind],alt[loc])
    } else {
      cat ("Fatal Error: SNP must be scored as 0, 1, 2 or NA\n"); stop()
    }
    }
  }
  
# Create the taxpartition (popname : 25-60)
  a <- array(data=NA,dim=length(poplabels))
  b <- array(data=NA,dim=length(poplabels))
  a[1] <- 1
  b <- table(poplabels)
  for (i in 2:length(b)) {
    b[i] = b[i] + b[i-1]
    a[i] = b[i-1] +1
  }
  plabels <- unique(poplabels)
  
# Create the svd file
  if (v > 1) {cat(paste("    Writing results to nexus file",outfile,"\n"))}
  sink(outfile)
  cat("#nexus\n")
  cat("BEGIN DATA;\n")
  cat(paste0("     dimensions ntax = ",2*nInd(x)," nchar = ",nLoc(x)," ;\n"))
  cat("     format datatype = dna gap = - ;\n\n")
  cat("matrix\n")
  if (method == 1) {
    for (i in 1:nInd(x)){
      cat(paste0(refseq[i],"\n"))
      cat(paste0(altseq[i],"\n"))
    }
  } else {
    for (i in 1:nInd(x)){
      cat(paste0(ambseq[i],"\n"))
    }
  }  
  cat(";\n")
  cat("end;\n\n")
  cat("begin sets;\n")
  cat("    taxpartition pops =\n")
  for (i in 1:(length(plabels)-1)) {
    cat(paste0("        ",plabels[i]," : ",a[i],"-",b[i],",\n"))
  }
  cat("       ",paste0(plabels[length(plabels)]," : ",a[length(plabels)],"-",b[length(plabels)],";\n"))
  cat("end;\n\n")
  cat("begin paup;\n")
  cat("log file=svd.txt;\n")
  cat("lset nthreads=3;\n")
  cat("SVDQuartets\n")
  cat("    evalQuartets=random\n")
  cat("    speciesTree=yes\n")
  cat("    partition=pops\n")
  cat("    bootstrap=standard\n")
  cat("    nreps=10000\n")
  cat("    ambigs=distribute\n")
  cat("    treeFile=svd.tre;\n")
  cat("savetrees file=svd_boot.tre from=1 to=1 maxdecimals=2;\n")
  cat("log stop;\n")
  cat("quit;\n")
  cat("end;\n")
  
  sink()
  
  if (v > 2) {cat(paste("    Records written to",outfile,":",nInd(x),"\n"))}
  if (v > 0) {cat("gl2svdquartets Completed\n")}
  
  return(NULL)

}

