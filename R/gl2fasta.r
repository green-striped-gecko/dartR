#' Concatenates DArT trimmed sequences and outputs a fastA file.
#'
#' Concatenated sequence tags are useful for phylogenetic methods where information
#' on base frequencies and transition and transversion ratios are required (for
#' example, Maximum Liklihood methods). Where relevant, heterozygous loci are resolved 
#' before concatenation by either assigning ambiguity codes or by random allele 
#' assignment. 
#'
#' Four methods are employed
#' 
#' Method 1 -- heterozygous positions are replaced by
#' the standard ambiguity codes. The resultant sequence fragments are concatenated across loci to generate
#' a single combined sequence to be used in subsequent ML phylogenetic analyses.
#'
#' Method=2 -- the heterozyous state is resolved by randomly assigning one or the other SNP variant
#' to the individual. The resultant sequence fragments are concatenated across loci to generate
#' a single composite haplotype to be used in subsequent ML phylogenetic analyses.
#'
#' Method 3 -- heterozygous positions are replaced by
#' the standard ambiguity codes. The resultant SNP bases are concatenated across loci to generate
#' a single combined sequence to be used in subsequent MP phylogenetic analyses.
#'
#' Method=4 -- the heterozyous state is resolved by randomly assigning one or the other SNP variant
#' to the individual. The resultant SNP bases are concatenated across loci to generate
#' a single composite haplotype to be used in subsequent MP phylogenetic analyses.
#' 
#' Trimmed sequences for which the SNP has been trimmed out, rarely, by adaptor mis-identity are deleted.
#'
#' The script writes out the composite haplotypes for each individual as a fastA file. Requires
#' 'TrimmedSequence' to be among the locus metrics (\code{@other$loc.metrics}) and information of the type of alleles (slot loc.all e.g. "G/A") and the position of the SNP in slot position of the ```genlight``` object (see testset.gl@position and testset.gl@loc.all for how to format these slots.)
#'  
#' @param x -- name of the DArT genlight object [required]
#' @param method -- 1 | 2 | 3 | 4. Type method=0 for a list of options  [method=1]
#' @param outfile -- name of the output file (fasta format) [output.fasta]
#' @param outpath -- path where to save the output file (set to tempdir by default)
#' @param probar -- if TRUE, a progress bar will be displayed for long loops [default = TRUE]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A new gl object with all loci rendered homozygous
#' @export
#' @importFrom seqinr write.fasta
#' @importFrom utils combn edit flush.console getTxtProgressBar read.csv setTxtProgressBar txtProgressBar write.csv write.table
#' @importFrom graphics axis barplot box image lines text
#' @importFrom methods new
#' @importFrom stats dist nobs optimize pchisq variable.names
#' @import stringr
#' @author Bernd Gruber and Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl <- gl.filter.repavg(testset.gl,t=1)
#' gl <- gl.filter.overshoot(gl,verbose=3)
#' gl <- gl.filter.callrate(testset.gl,t=.98)
#' gl <- gl.filter.monomorphs(gl)
#' gl2fasta(gl, method=1, outfile="test.fasta",verbose=3)


gl2fasta <- function(x, method=1, outfile="output.fasta", outpath=tempdir(), probar=FALSE, verbose=NULL) {

# TIDY UP FILE SPECS
  
  build <- "Jacob"
  funname <- match.call()[[1]]
  hold <- x
  outfile <- file.path(outpath, outfile)
  # Note draws upon and modifies the loc.metrics.flags for Call Rate, and modifies the flags for all other metrics if method='ind'.
  
# FLAG SCRIPT START
  # set verbosity
  if (is.null(verbose) & !is.null(x@other$verbose)) verbose=x@other$verbose
  if (is.null(verbose)) verbose=2
 
  
  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }
  
  if (verbose >= 1){
    cat("Starting",funname,"\n")
  }
  
# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    stop("  Fatal Error: genlight object required!\n")
  }
  
  if (verbose >= 2){
    if (all(x@ploidy == 1)){
      stop("Fatal Error: Detected Presence/Absence (SilicoDArT) data. Please provide a SNP dataset with trimmed sequences\n")
    } else if (all(x@ploidy == 2)){
      cat("  Processing a SNP dataset\n")
    } else {
      stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
    }
  }
  
  # Check monomorphs have been removed up to date
  if (x@other$loc.metrics.flags$monomorphs == FALSE){
    if (verbose >= 2){
      cat("  Warning: Dataset contains monomorphic loci which will be included in the output fastA file\n")
    }  
  }
  
# FUNCTION SPECIFIC ERROR CHECKING
  
  if(length(x@other$loc.metrics$TrimmedSequence) != nLoc(x)) {
    stop("Fatal Error: Data must include Trimmed Sequences for each loci in a column called 'TrimmedSequence' in the @other$loc.metrics slot.\n")
  }
  if(length(x@position) != nLoc(x)) {
    stop("Fatal Error: Data must include position information for each loci in the @position slot.\n")
  }
  if(length(x@loc.all) != nLoc(x)) {
    stop("Fatal Error: Data must include type of alleles in the @loc.all slot.\n")
  }
  
  if(method==1 && verbose >= 2){
    cat("Assigning ambiguity codes to heterozygote SNPs, concatenating trimmed sequence\n")
  }else if (method==2 && verbose >= 2) {
    cat("Randomly allocating heterozygotes (1) to homozygotic state (0 or 2), concatenating trimmed sequence\n")
  }else if (method==3 && verbose >= 2) {
    cat("Assigning ambiguity codes to heterozygote SNPs, concatenating SNPs\n")
  }else if (method==4 && verbose >= 2) {
    cat("Randomly allocating heterozygotes (1) to homozygotic state (0 or 2), concatenating SNPs\n")
  }else{
    if(verbose >= 2){
      cat("Method not properly specified\n")
      cat("Replace score for heterozygotic loci with"):
      cat("  method=1 -- ambiguity codes, concatenate fragments) [default]\n")
      cat("  method=2 -- random assignment to one of the two homogeneous states, concatenate fragments\n")
      cat("  method=3 -- ambiguity codes, concatenate SNPs only\n")
      cat("  method=4 -- random assignment to one of the two homogeneous states, concatenate SNPs only\n")
    }
    stop("Fatal Error: Parameter method out of range.\n")
    
  }
  
# DO THE JOB

# METHOD = AMBIGUITY CODES

if (method==1 || method==3) {

  allnames <- locNames(x)
  snp = as.character(x@loc.all)
  trimmed <- as.character(x@other$loc.metrics$TrimmedSequence)
  snpmatrix <- as.matrix(x)
  
# Create a lookup table for the ambiguity codes
#     A  T  G  C
#  A  A  W  R  M)
#  T  W  T  K  Y
#  G  R  K  G  S
#  C  M  Y  S  C  

  conversion <- matrix(c("A","W", "R","M", "W", "T", "K", "Y", "R", "K", "G", "S", "M", "Y", "S", "C"),nrow=4, ncol=4)
  colnames(conversion) <- c("A","T","G","C")
  rownames(conversion) <- colnames(conversion)

# Extract alleles 1 and 2
  allelepos = x@position
  allele1 =gsub("(.)/(.)", "\\1", snp, perl=T)
  allele2 = gsub("(.)/(.)", "\\2", snp, perl=T)

  
  if(verbose >= 2){
      cat(paste("  Removing loci for which snp position is outside the length of the trimmed sequences\n" ) )
  }
  x <- gl.filter.overshoot(x,verbose=0)
  
  sequences <- NA

# Prepare the output fastA file
  if(verbose >= 2){
    cat("Generating haplotypes ... This may take some time\n")
  }

  sink(outfile)

  for (i in 1:nInd(x)) {
    seq <- NA
    for (j in 1:nLoc(x)) {
      if (is.na(snpmatrix[i,j])) {
        code <- "N"
      } else {
        if (snpmatrix[i,j]==0)  {a1=allele1[j]; a2=allele1[j] }
        if (snpmatrix[i,j]==1)  {a1=allele1[j]; a2=allele2[j] }
        if (snpmatrix[i,j]==2)  {a1=allele2[j]; a2=allele2[j] }
        code <- conversion[a1, a2]
      }
      snppos <- allelepos[j]
      if(method==1){
        if (code !="N") {
          seq[j] <- paste0( substr(as.character(x@other$loc.metrics$TrimmedSequence[j]),1,snppos),code,substr(x@other$loc.metrics$TrimmedSequence[j],snppos+2, 500))
        } else {
          seq[j] <- paste(rep("N", nchar(as.character(x@other$loc.metrics$TrimmedSequence[j]))), collapse = "")
        }
      } else if(method==3){
          seq[j] <- code
      }
    }
    # Join all the trimmed sequence together into one long "composite" haplotype
      result <- paste(seq,sep="",collapse="")
    # Write the results to file in fastA format
      cat(paste0(">", indNames(x)[i],"_",pop(x)[i], "\n"))
      cat(result, " \n")

    # cat(paste("Individual:", i,"Took: ", round(proc.time()[3]-ptm),"seconds\n") )
  }

# Close the output fastA file
  sink()
}

# METHOD = RANDOM ASSIGNMENT

if (method==2 || method==4) {

# Randomly allocate heterozygotes (1) to homozygote state (0 or 2)
  matrix <- as.matrix(x)
  #cat("Randomly allocating heterozygotes (1) to homozygote state (0 or 2)\n")
  #pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
  #getTxtProgressBar(pb)
  r <- nrow(matrix)
  c <- ncol(matrix)
  for (i in 1:r) {
    for (j in 1:c) {
      if (matrix[i,j] == 1 && !is.na(matrix[i,j])) {
        # Score it 0 or 2
        matrix[i,j] <- (sample(1:2, 1)-1)*2
      }
    }
  #  setTxtProgressBar(pb, i/r)
  }
  
  if(verbose >= 2){
    cat(paste("  Removing loci for which snp position is outside the length of the trimmed sequences\n" ) )
  }
  x <- gl.filter.overshoot(x,verbose=0)
  
# Prepare the output fastA file
  if(verbose >= 2){cat("Generating haplotypes ... This may take some time\n")}

  sink(outfile)

# For each individual, and for each locus, generate the relevant haplotype 
  seq <- rep(" ", c)
  #pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
  #getTxtProgressBar(pb)
  for (i in 1:r) {
    for (j in 1:c) {
      # Reassign some variables
        trimmed <- as.character(x@other$loc.metrics$TrimmedSequence[j])
        snp <- x@loc.all[j]
        snpos <- x@position[j]
      # Shift the index for snppos to start from 1 not zero
        snpos <- snpos +1
      
      # If the score is homozygous for the reference allele
        if(method==2){
          if (matrix[i,j] == 0 && !is.na(matrix[i,j])) {
            seq[j] <- trimmed
          }
        } else if(method==4){
            seq[j] <- str_sub(trimmed, start=(snpos), end=(snpos))
        }
      # If the score is homozygous for the alternate allele
        if (matrix[i,j] == 2 && !is.na(matrix[i,j])) {
          # Split the trimmed into a beginning sequence, the SNP and an end sequences
            start <- str_sub(trimmed, end=snpos-1)
            snpbase <- str_sub(trimmed, start=(snpos), end=(snpos))
            end <- str_sub(trimmed, start=snpos+1)
        # Extract the SNP transition bases (e.g. A and T)
          state1 =gsub("(.)/(.)", "\\1", snp, perl=T)
          state2 = gsub("(.)/(.)", "\\2", snp, perl=T)
        # Change the SNP state to the alternate
          if (snpbase == state1) {
            snpbase <- state2
          } else {
            snpbase <- state1
          }

          if(method==2){
          # Paste back to form the alternate fragment
            target <- paste0(start,snpbase,end)
          # Remove adaptors and save the trimmed alternate sequence
            seq[j] <- str_sub(target, start=1, end=nchar(trimmed))
          }else if(method==4){
            seq[j] <- snpbase
          }
        }

      # If the SNP state is missing, assign NNNNs  
        if (is.na(matrix[i,j])) {
          seq[j] <- "N"
          if(method==2){
            seq[j] <- str_pad(seq[j], nchar(trimmed), side = c("right"), pad = "N")
          }
        }
    }
    #setTxtProgressBar(pb, i/r)

    # Join all the trimmed sequence together into one long "composite" haplotype
    result <- paste(seq,sep="",collapse="")
    # Write the results to file in fastA format
    cat(paste0(">", indNames(x)[i],"_",pop(x)[i], "\n"))
    cat(result, " \n")
  
  } # Select the next individual and repeat
  
  # Close the output fastA file
    sink()

  }
  
# FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat("Completed:",funname,"\n")
  }

return(TRUE)
  
}


