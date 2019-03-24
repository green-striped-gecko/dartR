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
#' @param gl -- name of the DArT genlight object [required]
#' @param method -- 1 | 2 | 3 | 4. Type method=0 for a list of options  [method=1]
#' @param outfile -- name of the output file (fasta format) [output.fasta]
#' @param outpath -- path where to save the output file (set to tempdir by default)
#' @param probar -- if TRUE, a progress bar will be displayed for long loops [default = TRUE]
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
#' gl <- gl.filter.callrate(testset.gl,t=.98)
#' gl2fasta(gl, method=1, outfile="test.fasta")


gl2fasta <- function(gl, method=1, outfile="output.fasta", outpath=tempdir(), probar=TRUE) {
  outfile <- file.path(outpath, outfile)
  
  if(class(gl) != "genlight") {
    stop("Fatal Error: Specify a genlight object\n")
  }
  if(length(gl@other$loc.metrics$TrimmedSequence) != nLoc(gl)) {
    stop("Fatal Error: Data must include Trimmed Sequences for each loci in a column called 'TrimmedSequence' in the @other$loc.metrics slot.\n")
  }
  if(length(gl@position) != nLoc(gl)) {
     stop("Fatal Error: Data must include position information for each loci in the @position slot.\n")
  }
  if(length(gl@loc.all) != nLoc(gl)) {
    stop("Fatal Error: Data must include type of alleles in the @loc.all slot.\n")
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
  snp = as.character(gl@loc.all)
  trimmed <- as.character(gl@other$loc.metrics$TrimmedSequence)
  snpmatrix <- as.matrix(gl)
  
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
  allelepos = gl@position
  allele1 =gsub("(.)/(.)", "\\1", snp, perl=T)
  allele2 = gsub("(.)/(.)", "\\2", snp, perl=T)

  
  lenTrim <- nchar(as.character(gl@other$loc.metrics$TrimmedSequence))
  
  index <- lenTrim>allelepos
  if (sum(index)!=nLoc(gl)) 
  {
    cat(paste("Not all snp position are within the length of the trimmed sequences. Those loci will be deleted (",sum(!index),")."  ) )
  }
  
  sequences <- NA

# Prepare the output fastA file
  cat("Generating haplotypes ... This may take some time\n")

  sink(outfile)

  for (i in 1:nInd(gl)) {
    seq <- NA
    for (j in 1:nLoc(gl)) {
      if (index[j]) {
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
          seq[j] <- paste0( substr(as.character(gl@other$loc.metrics$TrimmedSequence[j]),1,snppos),code,substr(gl@other$loc.metrics$TrimmedSequence[j],snppos+2, 500))
        } else {
          seq[j] <- paste(rep("N", nchar(as.character(gl@other$loc.metrics$TrimmedSequence[j]))), collapse = "")
        }
      } else if(method==3){
          seq[j] <- code
      }
      }#run only if index is true
    }
    # seqall = paste(seq, collapse="")
    # write.fasta(seqall, gl@other$ind.metrics$phylo.label[i], file.out=outfile, open="a")
    # Join all the trimmed sequence together into one long "composite" haplotype
      result <- paste(seq,sep="",collapse="")
    # Write the results to file in fastA format
      cat(paste0(">", indNames(gl)[i],"_",pop(gl)[i], "\n"))
      cat(result, " \n")

    # cat(paste("Individual:", i,"Took: ", round(proc.time()[3]-ptm),"seconds\n") )
  }

# Close the output fastA file
  sink()
}

# METHOD = RANDOM ASSIGNMENT

if (method==2 || method==4) {

# Randomly allocate heterozygotes (1) to homozygote state (0 or 2)
  matrix <- as.matrix(gl)
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
  
  lenTrim <- nchar(as.character(gl@other$loc.metrics$TrimmedSequence))
  index <- lenTrim > gl@position
  if (sum(index)!=nLoc(gl)) 
  {
    cat(paste("Not all SNP positions are within the length of the trimmed sequences. Those loci will be deleted (",sum(!index),")."  ) )
  }
  
# Prepare the output fastA file
  cat("Generating haplotypes ... This may take some time\n")

  sink(outfile)

# For each individual, and for each locus, generate the relevant haplotype 
  seq <- rep(" ", c)
  #pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
  #getTxtProgressBar(pb)
  for (i in 1:r) {
    for (j in 1:c) {
      if (index[j]) {
      # Reassign some variables
        trimmed <- as.character(gl@other$loc.metrics$TrimmedSequence[j])
        snp <- gl@loc.all[j]
        snpos <- gl@position[j]
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
    }
    #setTxtProgressBar(pb, i/r)

    # Join all the trimmed sequence together into one long "composite" haplotype
    result <- paste(seq,sep="",collapse="")
    # Write the results to file in fastA format
    cat(paste0(">", indNames(gl)[i],"_",pop(gl)[i], "\n"))
    cat(result, " \n")
  
  } # Select the next individual and repeat
  
  # Close the output fastA file
    sink()

  }

return(TRUE)
  
}


