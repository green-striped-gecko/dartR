#' Concatenates DArT trimmed sequences and outputs a fastA file.
#'
#' Heterozygous loci are resolved by either assigning ambiguity codes or by random allele assignment. 
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
#' The script writes out the composite haplotypes for each individual as a fastA file. Requires
#' trimmed sequences to be among the locus metrics.
#' 
#' @param gl -- name of the DArT genlight object [required]
#' @param method -- 1 | 2 | 3 | 4. Type method=0 for a list of options  [method=1]
#' @param outfile -- name of the output file (fasta format) [output.fasta]
#' @return A new gl object with all loci rendered homozygous
#' @export
#' @import adegenet
#' @importFrom seqinr write.fasta
#' @importFrom utils combn edit flush.console getTxtProgressBar read.csv setTxtProgressBar txtProgressBar write.csv write.table
#' @importFrom graphics axis barplot box image lines text
#' @importFrom methods new
#' @importFrom stats dist nobs optimize pchisq variable.names
#' @import stringr
#' @author Bernd Gruber and Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' \dontrun{
#' gl2fasta(gl, method=1, outfile="test.fasta")
#' }

gl2fasta <- function(gl, method=1, outfile="output.fasta") {
  
  if(class(gl) == "genlight") {
    cat("Analysing a genlight object\n")
  } else {
    cat("Fatal Error: Specify a genlight object\n")
    stop()
  }
  if(length(gl@other$loc.metrics$TrimmedSequence) == 0) {
    cat("Fatal Error: Data must include Trimmed Sequences\n"); stop()
  }
  if(method==1){
    cat("Assigning ambiguity codes to heterozygote SNPs, concatenating trimmed sequence\n")
  }else if (method==2) {
    cat("Randomly allocating heterozygotes (1) to homozygotic state (0 or 2), concatenating trimmed sequence\n")
  }else if (method==3) {
    cat("Assigning ambiguity codes to heterozygote SNPs, concatenating SNPs\n")
  }else if (method==4) {
    cat("Randomly allocating heterozygotes (1) to homozygotic state (0 or 2), concatenating SNPs\n")
  }else{
    cat("Fatal Error: Parameter method out of range.\n")
    cat("Replace score for heterozygotic loci with"):
    cat("  method=1 -- ambiguity codes, concatenate fragments) [default]\n")
    cat("  method=2 -- random assignment to one of the two homogeneous states, concatenate fragments\n")
    cat("  method=3 -- ambiguity codes, concatenate SNPs only\n")
    cat("  method=4 -- random assignment to one of the two homogeneous states, concatenate SNPs only\n")
    stop()
  }

# METHOD = AMBIGUITY CODES

if (method==1 || method==3) {

  allnames <- locNames(gl)
  snp = as.character(gl@other$loc.metrics$SNP)
  trimmed <- as.character(gl@other$loc.metrics$TrimmedSequence)
  snpmatrix <- as.matrix(gl)
  
# Create a lookup table for the ambiguity codes
#     A  T  G  C
#  A  A  W  R  M
#  T  W  T  K  Y
#  G  R  K  G  S
#  C  M  Y  S  C  

  conversion <- matrix(c("A","W", "R","M", "W", "T", "K", "Y", "R", "K", "G", "S", "M", "Y", "S", "C"),nrow=4, ncol=4)
  colnames(conversion) <- c("A","T","G","C")
  rownames(conversion) <- colnames(conversion)

# Extract alleles 1 and 2
  allelepos = as.numeric(gsub("(\\d{1,3}):(.)>(.)", "\\1", snp, perl=T))
  allele1 =gsub("(\\d{1,3}):(.)>(.)", "\\2", snp, perl=T)
  allele2 = gsub("(\\d{1,3}):(.)>(.)", "\\3", snp, perl=T)

  sequences <- NA

# Prepare the output fastA file
  cat("Generating haplotypes ... This may take some time\n")

  sink(outfile)

  for (i in 1:nInd(gl)) {
    seq <- NA
    for (j in 1:nLoc(gl)) {
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
          seq[j] <- paste0( substr(gl@other$loc.metrics$TrimmedSequence[j],1,snppos),code,substr(gl@other$loc.metrics$TrimmedSequence[j],snppos+2, 500))
        } else {
          seq[j] <- paste(rep("N", nchar(as.character(gl@other$loc.metrics$TrimmedSequence[j]))), collapse = "")
        }
      } else if(method==3){
          seq[j] <- code
      }
    }
#    seqall = paste(seq, collapse="")
#    write.fasta(seqall, gl@other$ind.metrics$phylo.label[i], file.out=outfile, open="a")

    # Join all the trimmed sequence together into one long "composite" haplotype
      result <- paste(seq,sep="",collapse="")
    # Write the results to file in fastA format
      cat(paste0(">", indNames(gl)[i],"|",pop(gl)[i], "\n"))
      cat(result, " \n")

#    cat(paste("Individual:", i,"Took: ", round(proc.time()[3]-ptm),"seconds\n") )
  }

# Close the output fastA file
  sink()
}

# METHOD = RANDOM ASSIGNMENT

if (method==2 || method==4) {

# Randomly allocate heterozygotes (1) to homozygote state (0 or 2)
  matrix <- as.matrix(gl)
  cat("Randomly allocating heterozygotes (1) to homozygote state (0 or 2)\n")
  pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
  getTxtProgressBar(pb)
  r <- nrow(matrix)
  c <- ncol(matrix)
  for (i in 1:r) {
    for (j in 1:c) {
      if (matrix[i,j] == 1 && !is.na(matrix[i,j])) {
        # Score it 0 or 2
        matrix[i,j] <- (sample(1:2, 1)-1)*2
      }
    }
    setTxtProgressBar(pb, i/r)
  }
  
# Check that the sequences are all the same length
  
  mx <- max(str_length(gl@other$loc.metrics$AlleleSequence))
  mn <- min(str_length(gl@other$loc.metrics$AlleleSequence))
  if (mn == mx) {
    cat(paste("\nAllelic sequences (incl. adaptors) are all the same length:", mx),"\n")
  } else {
    cat("Allele sequences not all the same length, padding with NNNNs\n")
    gl@other$loc.metrics$AlleleSequence <- str_pad(gl@other$loc.metrics$AlleleSequence, mx, side = c("right"), pad = "N")
  }
  
# Prepare the output fastA file
  cat("Generating haplotypes ... This may take some time\n")

  sink(outfile)

# For each individual, and for each locus, generate the relevant haplotype 
  seq <- rep(" ", c)
  for (i in 1:r) {
    for (j in 1:c) {
      # Reassign some variables
        fragment <- as.character(gl@other$loc.metrics$AlleleSequence[j])
        trimmed <- as.character(gl@other$loc.metrics$TrimmedSequence[j])
        snp <- as.character(gl@other$loc.metrics$SNP[j])
        snpos <- gl@other$loc.metrics$SnpPosition[j]
      # Shift the index for snppos to start from 1 not zero
        snpos <- snpos +1
      
      # If the score is homozygous for the reference allele
        if(method==2){
          if (matrix[i,j] == 0 && !is.na(matrix[i,j])) {
            seq[j] <- trimmed
          }
        } else if(method==4){
            seq[j] <- str_sub(fragment, start=(snpos), end=(snpos))
        }

      # If the score is homozygous for the alternate allele
        if (matrix[i,j] == 2 && !is.na(matrix[i,j])) {
          # Split the fragment into a beginning sequence, the SNP and an end sequences
            start <- str_sub(fragment, end=snpos-1)
            snpbase <- str_sub(fragment, start=(snpos), end=(snpos))
            end <- str_sub(fragment, start=snpos+1)
        # Extract the SNP transition bases (e.g. A and T)
          state.change <- str_split_fixed(snp,":",2)
          state1 <- str_sub(state.change[2], start=1, end=1)
          state2 <- str_sub(state.change[2], start=3, end=3)
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
            #seq[j] <- gl.utils.1(trimmed,target)
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

    # Join all the trimmed sequence together into one long "composite" haplotype
    result <- paste(seq,sep="",collapse="")
    # Write the results to file in fastA format
    cat(paste0(">", indNames(gl)[i],"|",pop(gl)[i], "\n"))
    cat(result, " \n")
  
  } # Select the next individual and repeat
  
  # Close the output fastA file
    sink()

  }

return(TRUE)
  
}


