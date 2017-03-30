#' Export DArT genlight object \{adegenet\} to fastA format
#'
#' @param gl -- genlight object
#' @param outfile -- filename of the output fastA file [default genlight.fasta]
#' @return NULL
#' @export
#' @author Bernd Gruber (glbugs@aerg.canberra.edu.au)
#' @import adegenet
#' @importFrom utils combn edit flush.console getTxtProgressBar read.csv setTxtProgressBar txtProgressBar write.csv write.table
#' @importFrom graphics axis barplot box image lines text
#' @importFrom methods new
#' @importFrom stats dist nobs optimize pchisq variable.names

dart2fasta <- function(gl, outfile = "genlight.fasta"){

  if(class(gl) == "genlight") {
    cat("Filtering a genlight object\n")
  } else {
    cat("Fatal Error: Specify a genlight object\n")
    cat("  Syntax: gl2fasta(gl, outfile=\"test.fasta\"\n")
    stop()
  }
  if (length(gl@other.metrics$TrimmedSequence) >= 0) {
    cat("Trimmed sequences located\n")
  } else {
    cat("Fatal Error: genlight object does not contain trimmed sequences.\n")
    stop()
  }
    
  ptm <- proc.time()[3]

  allnames <- locNames(gl)
  snpdesc = gl@other$metrics$SNP
  allseq <- gl@other$metrics$TrimmedSequence
  snpmatrix <- as.matrix(gl)

  conversion <- matrix(c("A","W", "R","M", "W", "T", "K", "Y", "R", "K", "G", "S", "M", "Y", "S", "C"),nrow=4, ncol=4)
  colnames(conversion)<- c("A","T","G","C")
  rownames(conversion)<- colnames(conversion)

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
        indseq[j] <- paste0( substr(gl@other$metrics$TrimmedSequence[j],1,snppos),code,substr(gl@other$metrics$TrimmedSequence[j],snppos+2, 500))
      } else {
        indseq[j] <- paste(rep("N", nchar(as.character(gl@other$metrics$TrimmedSequence[j]))), collapse = "")
      }
    }
    indseqall = paste(indseq, collapse="")
    write.fasta(indseqall, gl@other$covariates$phylo.label[i], file.out=paste0(outfile, ".fasta"), open="a")
    cat(paste("Individual:", i,"Took: ", round(proc.time()[3]-ptm),"seconds\n") )
  }
  return(NULL)
}
