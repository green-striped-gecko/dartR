#' Convert a genlight object to nexus format PAUP SVDquartets
#'
#' Package SNPRelate relies on a bit-level representation of a SNP dataset that competes with \{adegenet\} genlight
#' objects and associated files. This function saves a genlight object to a gds format file.
#'
#' Reference: Chifman, J. and L. Kubatko. 2014. Quartet inference from SNP data under the coalescent, Bioinformatics, 30: 3317-3324
#' 
#' 
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param outfile -- file name of the output file (including extension).
#' @return NULL
#' @export
#' @importFrom
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' gl2svdquartets(testset.gl)

gl2svdquartets <- function(x, outfile="svd.nex") {
  
  if (v > 0) {cat(paste("Starting gl2svdquartets: Create nexus file suitable for svdquartets\n\n"))}

# Extract the reference base and the alternate base for each locus
  v <- as.character(x@other$loc.metrics$SNP)
  v <- strsplit(v,":")
  v <- unlist(lapply(v, `[`, 2))
  v <- strsplit(v,">")
  ref <- unlist(lapply(v, `[`, 1))
  alt <- unlist(lapply(v, `[`, 2))
  
# Sort the data on population
  df <- data.frame(as.matrix(x))
  df <- cbind(indNames(x),pop(x),df)
  df <- df[order(df$pop),]
  indlabels <- df[,1]
  poplabels <- df[,2]
  m <- df[,3:(nLoc(x)+2)]

  #Initialize The first base of the sequence for each individual  
  refseq <- array(data=NA,dim=nInd(x))
  altseq <- array(data=NA,dim=nInd(x))
  for (ind in 1:nInd(x)) {
    if (is.na(m[ind,1])) {
      refseq[ind] <- "?"; altseq[ind] <- "?"
    } else if (m[ind,1] == 0) {
      refseq[ind] <- ref[1]; altseq[ind] <- ref[1]
    } else if (m[ind,1] == 1) {
      refseq[ind] <- ref[1]; altseq[ind] <- alt[1]
    } else if (m[ind,1] == 2) {
      refseq[ind] <- alt[1]; altseq[ind] <- alt[1]
    } else {
      cat ("Fatal Error: SNP must be scored as 0, 1, 2 or NA\n"); stop()
    }
  }
  
# progressively add the bases  
  for (ind in 1:nInd(x)) {
    for (loc in 2:nLoc(x)){
    if (is.na(m[ind,loc])) {
      refseq[ind] <- paste0(refseq[ind],"?"); altseq[ind] <- paste0(altseq[ind],"?")
    } else if (m[ind,loc] == 0) {
      refseq[ind] <- paste0(refseq[ind],ref[loc]); altseq[ind] <- paste0(altseq[ind],ref[loc])
    } else if (m[ind,loc] == 1) {
      refseq[ind] <- paste0(refseq[ind],ref[loc]); altseq[ind] <- paste0(altseq[ind],alt[loc])
    } else if (m[ind,loc] == 2) {
      refseq[ind] <- paste0(refseq[ind],alt[loc]); altseq[ind] <- paste0(altseq[ind],alt[loc])
    } else {
      cat ("Fatal Error: SNP must be scored as 0, 1, 2 or NA\n"); stop()
    }
    }
  }
  
# Add individual names to sequences (should have done this first)
  for (ind in 1:nInd(x)) {
    refseq[ind] <- paste0(indlabels[ind],"_1   ",refseq[ind])
    altseq[ind] <- paste0(indlabels[ind],"_2   ",altseq[ind])
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
  sink(outfile)
  cat("#nexus\n")
  cat("BEGIN DATA;\n")
  cat(paste0("     dimensions ntax = ",2*nInd(x)," nchar = ",nLoc(x)," ;\n"))
  cat("     format datatype = dna gap = - ;\n\n")
  cat("matrix\n")
  for (i in 1:nInd(x)){
    cat(paste0(refseq[i],"\n"))
    cat(paste0(altseq[i],"\n"))
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
  
  if (v > 0) {cat("gl2svdquartets Completed\n")}
  
  return(NULL)

}

