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
#' For presence/absence data (SilicoDArT), it is not possible to count transitions and tranversions or tv/ts ratio
#' because the SNP data is not available, only a single sequence tag.
#' 
#' @param x -- name of the genlight object containing the SNP or presence/absence data [required]
#' @param plot -- if TRUE, histograms of base composition are produced [default TRUE]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @param silent -- if FALSE, function returns an object, otherwise NULL [default TRUE]
#' @return If silent==TRUE, returns NULL; otherwise a matrix containing the percent frequencies of each base (A,C,T,G) and the transition and transversion frequencies.
#' @export
#' @import stringr
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' lst <- gl.report.bases(testset.gl,silent=FALSE)
#' lst

gl.report.bases <- function(x, plot=TRUE, silent=TRUE, verbose = 2) {

# TIDY UP FILE SPECS

  build <- "Jacob"
  funname <- match.call()[[1]]
  # Note does not draw upon or modify the loc.metrics.flags
  
  # FLAG SCRIPT START
  
  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 3\n")
    verbose <- 3
  }
  
  cat("Starting",funname,"[ Build =",build,"]\n")

# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }
  
    if (all(x@ploidy == 1)){
      cat("  Processing Presence/Absence (SilicoDArT) data\n")
    } else if (all(x@ploidy == 2)){
      cat("  Processing a SNP dataset\n")
    } else {
      cat ("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)"); stop("Terminating Execution!")
    }

# FUNCTION SPECIFIC ERROR CHECKING

    if(!any(names(x@other$loc.metrics) == "TrimmedSequence")) {
      cat("  Fatal Error: Dataset does not include variable TrimmedSequence!\n"); stop()
    }

# DO THE JOB
  
# Count up the number of bases, and the number of each of ATGC, and other
  if (verbose >= 2){
    cat("  Counting the bases\n")
  }
  A <- sum(str_count(x@other$loc.metrics$TrimmedSequence, "A"))
  G <- sum(str_count(x@other$loc.metrics$TrimmedSequence, "G"))
  C <- sum(str_count(x@other$loc.metrics$TrimmedSequence, "C"))
  T <- sum(str_count(x@other$loc.metrics$TrimmedSequence, "T"))
  total <- sum(str_length(x@other$loc.metrics$TrimmedSequence))
  total.ATGC <- sum(A,G,C,T)
  if (verbose >= 2){
  if (total != total.ATGC) {
    cat("  Warning: Codes other than A, T, G and C present\n")
  }
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

  if (verbose >= 2){
    cat("  Counting Transitions and Transversions\n")
  }  
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
  
  if (verbose >= 2){
  if (ts+tv != length(x@other$loc.metrics$TrimmedSequence)) {
    cat("  Warning: Sum of transitions plus transversions does not equal number of loci.\n")
  }
  }  
  ts <- ts*100/length(x@other$loc.metrics$TrimmedSequence)
  tv <- tv*100/length(x@other$loc.metrics$TrimmedSequence)
  ratio <- ts/tv

# Report
  cat(paste("  Average trimmed sequence length:",round(mn,digits=1),"(",mi,"to",mx,")"),"\n")
  cat(paste("  Total number of trimmed sequences:",length(x@other$loc.metrics$TrimmedSequence)),"\n")
  cat("  Base frequencies (%)\n")
  cat(paste("    A:",round(A,2)),"\n")
  cat(paste("    G:",round(G,2)),"\n")
  cat(paste("    T:",round(T,2)),"\n")
  cat(paste("    C:",round(C,2)),"\n\n")
  if (all(ploidy(x)==1)) {
    cat("  Presence Absence data (SilicoDArT), transition/transversions cannot be calculated\n")
    tv <- NA
    ts <- NA
  } else {
    cat(paste("  Transitions  :",round(ts,2),"\n"))
    cat(paste("  Transversions:",round(tv,2),"\n"))
    cat(paste("  tv/ts ratio:", round(ratio,4),"\n\n"))
  }
  
  if (plot) {
    par(mfrow = c(2, 1),pty="s")
    df <- cbind(A,C,T,G)
    if (all(x@ploidy==2)){
      title <- paste0("SNP data (DArTSeq)\nBase Frequencies")
    } else {
      title <- paste0("Fragment P/A data (SilicoDArT)\nBase Frequencies")
    }  
    barplot(df, main=title, col=rainbow(1), width=c(.1,.1,.1,.1))
    if (all(x@ploidy!=1)){
      df <- cbind(ts,tv)
      title <- paste("Transitions and Transversion Rates\n (ts/tv ratio =",round(ratio,2),")")
      barplot(df, main=title, col=rainbow(1))
    }
  }  
  
# Create return matrix
  if (verbose >= 2){
    cat("  Creating an output matrix to return\n")
  }
  col1 <- c("A","G","T","C","tv","ts")
  col2 <- c(round(A,2),round(G,2),round(T,2),round(C,2),round(tv,2),round(ts,2))
  matrix <- cbind(col1, col2)
  #colnames(matrix) <- c("Base","Pcent")
  
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

  if(silent==TRUE){
    return(NULL)
  } else{
      return(matrix)
  }  
  
}
