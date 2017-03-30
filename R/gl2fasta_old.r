#' Export DArT genlight object \{adegenet\} to fastA format
#'
#' Concatenates DArT trimmed sequences after replacing allelic states at heterozygous loci with standard ambiguity codes.
#'
#' This script generates haplotypes from SNP genotypes in a DArT genlight object in preparation for analyses
#' that expect only a single character state for each individual at each locus. For example, many 
#' phylogenetic analyses have difficulty accommodating heterozygosity. Under this approach, heterozytotes are represented
#' by ambiguity codes then the sequence tags are concatenated across loci to generate
#' a single haplotype to be used in the analysis.
#'
#' The script writes out the composite haplotypes for each individual as a fastA file. Requires
#' trimmed sequences to be among the locus metrics.
#'
#' @param gl -- genlight object
#' @param outfile -- filename of the output fastA file [default genlight.fasta]
#' @return NULL
#' @author Bernd Gruber (glbugs@aerg.canberra.edu.au)
#' @import adegenet
#' @importFrom seqinr write.fasta
#' @importFrom utils combn edit flush.console getTxtProgressBar read.csv setTxtProgressBar txtProgressBar write.csv write.table
#' @importFrom graphics axis barplot box image lines text
#' @importFrom methods new
#' @importFrom stats dist nobs optimize pchisq variable.names

gl2fasta_old <- function(gl, outfile = "genlight.fasta") {

# Error checks
  if(class(gl) == "genlight") {
    cat("Analysing a genlight object\n")
  } else {
    cat("Fatal Error: Specify a genlight object\n")
    stop()
  }
  if(length(gl@other$loc.metrics$TrimmedSequence) == 0) {
    cat("Fatal Error: Data must include Trimmed Sequences\n"); stop()
  }
    
#  ptm <- proc.time()[3]

  allnames <- locNames(gl)
  snpdesc = gl@other$loc.metrics$SNP
  trimmed <- gl@other$loc.metrics$TrimmedSequence
  snpmatrix <- as.matrix(gl)

  conversion <- matrix(c("A","W", "R","M", "W", "T", "K", "Y", "R", "K", "G", "S", "M", "Y", "S", "C"),nrow=4, ncol=4)
  colnames(conversion)<- c("A","T","G","C")
  rownames(conversion)<- colnames(conversion)
  
#     A  T  G  C
#  A  A  W  R  M
#  T  W  T  K  Y
#  G  R  K  G  S
#  C  M  Y  S  C  

  allelepos = as.numeric(gsub("(\\d{1,3}):(.)>(.)", "\\1", snpdesc, perl=T))
  allel1 =gsub("(\\d{1,3}):(.)>(.)", "\\2", snpdesc, perl=T)
  allel2 = gsub("(\\d{1,3}):(.)>(.)", "\\3", snpdesc, perl=T)

  sequences<-NA

  for (i in 1:nInd(gl)) {
    indseq<- NA
    for (j in 1:nLoc(gl)) {
      if (is.na(snpmatrix[i,j])) {
        code <- "N"
      } else {
        if (snpmatrix[i,j]==0)  {a1=allel1[j]; a2=allel1[j] }
        if (snpmatrix[i,j]==1)  {a1=allel1[j]; a2=allel2[j] }
        if (snpmatrix[i,j]==2)  {a1=allel2[j]; a2=allel2[j] }
        code <- conversion[a1, a2]
      }
      snppos <- allelepos[j]
      if (code !="N") {
        indseq[j] <- paste0( substr(gl@other$loc.metrics$TrimmedSequence[j],1,snppos),code,substr(gl@other$loc.metrics$TrimmedSequence[j],snppos+2, 500))
      } else {
        indseq[j] <- paste(rep("N", nchar(as.character(gl@other$loc.metrics$TrimmedSequence[j]))), collapse = "")
      }
    }
    indseqall = paste(indseq, collapse="")
    write.fasta(indseqall, gl@other$ind.metrics$phylo.label[i], file.out=outfile, open="a")
#    cat(paste("Individual:", i,"Took: ", round(proc.time()[3]-ptm),"seconds\n") )
  }
  return(NULL)
}
