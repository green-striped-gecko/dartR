#' Summary of base pair frequencies
#'
#' This script calculates the frequencies of the four bases, and the frequency of transitions and
#' transversions in a DArT genlight object.
#'
#' The script checks if trimmed sequences are included in the locus metadata, and if so, tallies up
#' the numbers of A,T,G and C bases. Only the reference state at the SNP locus is counted. Counts of transitions
#' and transversions assume that there is no directionality, that is C>T is the same as T>C, because
#' the reference state is arbitrary.
#' 
#' @param gl -- name of the DArT genlight object [required]
#' @return Matrix containing the percent frequencies of each base (A,C,T,G) and the transition and 
#' transversion frequencies.
#' @export
#' @import adegenet stringr
#' @author Arthur Georges (bugs? Post to https://groups.google.com/d/forum/dartr)
#' @examples
#' lst <- gl.report.bases(testset.gl)
#' lst

gl.report.bases <- function(gl) {
x <- gl
  
  if(any(names(x@other$loc.metrics) == "TrimmedSequence")) {
    cat("Key variable TrimmedSequence found.\n\n")
  } else {
    cat("Fatal Error: Dataset does not include variable TrimmedSequence!\n"); stop()
  }
  
# Count up the number of bases, and the number of each of ATGC, and other
  A <- sum(str_count(x@other$loc.metrics$TrimmedSequence, "A"))
  G <- sum(str_count(x@other$loc.metrics$TrimmedSequence, "G"))
  C <- sum(str_count(x@other$loc.metrics$TrimmedSequence, "C"))
  T <- sum(str_count(x@other$loc.metrics$TrimmedSequence, "T"))
  total <- sum(str_length(x@other$loc.metrics$TrimmedSequence))
  total.ATGC <- sum(A,G,C,T)
  if (total != total.ATGC) {
    cat("Warning: Codes other than A, T, G and C present\n")
  }
  other <- total - total.ATGC
  other <- other*100/total
  A <- A*100/total
  G <- G*100/total
  T <- T*100/total
  C <- C*100/total
  
# Calculate the fragment lengths  
  mn <- mean(str_length(x@other$loc.metrics$TrimmedSequence))
  mx <- max(str_length(x@other$loc.metrics$TrimmedSequence))
  mi <- min(str_length(x@other$loc.metrics$TrimmedSequence))

# Extract the SNPs  
  matrix <- str_split_fixed(x@other$loc.metrics$SNP,":",2)
  state.change <- matrix[,2]
  
# Sum the transitions and transversions  
  tv <- sum(str_count(state.change, "A>C")) +
    sum(str_count(state.change, "C>A")) +
    sum(str_count(state.change, "G>T")) +
    sum(str_count(state.change, "T>G")) +
    sum(str_count(state.change, "A>T")) +
    sum(str_count(state.change, "T>A")) +
    sum(str_count(state.change, "G>C")) +
    sum(str_count(state.change, "C>G"))
  
  ts <- sum(str_count(state.change, "A>G")) +
    sum(str_count(state.change, "G>A")) +
    sum(str_count(state.change, "C>T")) +
    sum(str_count(state.change, "T>C"))
  
  if (ts+tv != length(x@other$loc.metrics$TrimmedSequence)) {
    cat("Warning: Sum of transitions plus transversions does not equal number of loci.\n")
  }
  ts <- ts*100/length(x@other$loc.metrics$TrimmedSequence)
  tv <- tv*100/length(x@other$loc.metrics$TrimmedSequence)
  ratio <- ts/tv

# Report
  cat(paste("Average trimmed sequence length:",round(mn,digits=1),"(",mi,"to",mx,")"),"\n")
  cat(paste("Total number of trimmed sequences:",length(x@other$loc.metrics$TrimmedSequence)),"\n")
  cat("Base frequencies (%)\n")
  cat(paste("  A:",round(A,2)),"\n")
  cat(paste("  G:",round(G,2)),"\n")
  cat(paste("  T:",round(T,2)),"\n")
  cat(paste("  C:",round(C,2)),"\n\n")
  cat(paste("Transitions  :",round(ts,2),"\n"))
  cat(paste("Transversions:",round(tv,2),"\n"))
  cat(paste("tv/ts ratio:", round(ratio,4),"\n\n"))

# Create return matrix
  col1 <- c("A","G","T","C","tv","ts")
  col2 <- c(round(A,2),round(G,2),round(T,2),round(C,2),round(tv,2),round(ts,2))
  matrix <- cbind(col1, col2)
  
  return(matrix)
  
}
