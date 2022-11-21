#' Converts a genlight object into a format suitable for input to the BPP 
#' program
#' 
#' This function generates the sequence alignment file and the Imap file. The 
#' control file should produced by the user. 
#' 
#' If method = 1, heterozygous positions are replaced by standard ambiguity 
#' codes.
#' 
#' If method = 2, the heterozygous state is resolved by randomly assigning one 
#' or the other SNP variant to the individual.
#'
#' Trimmed sequences for which the SNP has been trimmed out, rarely, by adapter
#'  mis-identity are deleted.
#'
#' This function requires 'TrimmedSequence' to be among the locus metrics
#' (\code{@other$loc.metrics}) and information of the type of alleles (slot
#' loc.all e.g. 'G/A') and the position of the SNP in slot position of the
#'  ```genlight``` object (see testset.gl@position and testset.gl@loc.all for
#'  how to format these slots.)
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param method One of 1 | 2, see details [default = 1].
#' @param outfile Name of the sequence alignment file ["output_bpp.txt"].
#' @param imap Name of the Imap file ["Imap.txt"].
#' @param outpath Path where to save the output file (set to tempdir by default)
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log; 3, progress and results summary; 5, full report 
#'  [default 2 or as specified using gl.set.verbosity].
#' @details
#' It's important to keep in mind that analyses based on coalescent theory, 
#' like those done by the programme BPP, are meant to be used with sequence
#'  data. In this type of data, large chunks of DNA are sequenced, so when we
#'   find polymorphic sites along the sequence, we know they are all on the same
#'    chromosome. This kind of data, in which we know which chromosome each 
#'    allele comes from, is called "phased data." Most data from reduced 
#'    representation genome-sequencing methods, like DArTseq, is unphased, 
#'    which means that we don't know which chromosome each allele comes from. 
#'    So, if we apply coalescence theory to data that is not phased, we will get
#'     biased results. As in Ellegren et al., one way to deal with this is to 
#'     "haplodize" each genotype by randomly choosing one allele from 
#'     heterozygous genotypes (2012) by using method = 2.
#'     
#'     Be mindful that there is little information in the literature on the
#'      validity of this method. 
#' @return NULL
#' @references
#' \itemize{
#' \item Ellegren, Hans, et al. "The genomic landscape of species divergence in 
#'Ficedula flycatchers." Nature 491.7426 (2012): 756-760.
#' \item Flouri T., Jiao X., Rannala B., Yang Z. (2018) Species Tree Inference 
#'with BPP using Genomic Sequences and the Multispecies Coalescent. Molecular
#' Biology and Evolution, 35(10):2585-2593. doi:10.1093/molbev/msy147
#'}
#' @export
#' @author Custodian: Luis Mijangos (Post to
#'  \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' require(dartR.data)
#' test <- platypus.gl
#' test <- gl.filter.callrate(test,threshold = 1)
#' test <- gl.filter.monomorphs(test)
#' test <- gl.subsample.loci(test,n=25)
#' gl2bpp(x = test)

gl2bpp <- function(x,
                   method = 1,
                   outfile = "output_bpp.txt", 
                   imap = "Imap.txt",
                   outpath = tempdir(),
                   verbose = NULL) {
  
  outfilespec <- file.path(outpath, outfile)
  outfilespec_imap <- file.path(outpath, imap)

  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jackson",
                   verbosity = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
  # CHECK IF PACKAGES ARE INSTALLED
  pkg <- "seqinr"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }
  
  # Check monomorphs have been removed up to date
  if (x@other$loc.metrics.flags$monomorphs == FALSE) {
    if (verbose >= 2) {
      cat(
        warn(
          "  Warning: Dataset contains monomorphic loci which will be included in the output fasta file\n"
        )
      )
    }
  }
  
  if (length(x@other$loc.metrics$TrimmedSequence) != nLoc(x)) {
    stop(
      error(
        "Fatal Error: Data must include Trimmed Sequences for each loci in a column called 'TrimmedSequence' in the @other$loc.metrics slot.\n"
      )
    )
  }
  if (length(x@position) != nLoc(x)) {
    stop(
      error(
        "Fatal Error: Data must include position information for each loci in the @position slot.\n"
      )
    )
  }
  if (length(x@loc.all) != nLoc(x)) {
    stop(error(
      "Fatal Error: Data must include type of alleles in the @loc.all slot.\n"))
  }
  
  if (verbose >= 2) {
    cat(
      report(
        "Assigning ambiguity codes to heterozygote SNPs\n"
      )
    )
  } 
  
  # DO THE JOB
  if (verbose >= 2) {
    cat(report(
      paste(
        "  Removing loci for which SNP position is outside the length of the trimmed sequences\n"
      )
    ))
  }
  x <- gl.filter.overshoot(x, verbose = 0)
  
  # METHOD = AMBIGUITY CODES
  
  if (method == 1) {
  
    allnames <- locNames(x)
    snp <- as.character(x@loc.all)
    trimmed <- as.character(x@other$loc.metrics$TrimmedSequence)
    snpmatrix <- as.matrix(x)
    
    # Create a lookup table for the ambiguity codes A T G C A A W R M) T W T K Y G R K G S C M Y S C
    
    conversion <-
      matrix(
        c(
          "A",
          "W",
          "R",
          "M",
          "W",
          "T",
          "K",
          "Y",
          "R",
          "K",
          "G",
          "S",
          "M",
          "Y",
          "S",
          "C"
        ),
        nrow = 4,
        ncol = 4
      )
    colnames(conversion) <- c("A", "T", "G", "C")
    rownames(conversion) <- colnames(conversion)
    
    # Extract alleles 1 and 2
    allelepos <- x@position
    allele1 <- gsub("(.)/(.)", "\\1", snp, perl = T)
    allele2 <- gsub("(.)/(.)", "\\2", snp, perl = T)
    
    # Prepare the output fastA file
    if (verbose >= 2) {
      cat(report("Generating haplotypes ... This may take some time\n"))
    }
    
    sink(outfilespec)
    
    for (j in 1:nLoc(x)) {
      seq <- NA
      for (i in 1:nInd(x)) {
        if (is.na(snpmatrix[i, j])) {
          code <- "N"
        } else {
          if (snpmatrix[i, j] == 0) {
            a1 <- allele1[j]
            a2 <- allele1[j]
          }
          if (snpmatrix[i, j] == 1) {
            a1 <- allele1[j]
            a2 <- allele2[j]
          }
          if (snpmatrix[i, j] == 2) {
            a1 <- allele2[j]
            a2 <- allele2[j]
          }
          code <- conversion[a1, a2]
        }
        snppos <- allelepos[j]
          if (code != "N") {
            seq[i] <-
              paste0(
                substr(
                  as.character(
                    x@other$loc.metrics$TrimmedSequence[j]
                  ),
                  1,
                  snppos
                ),
                code,
                substr(
                  x@other$loc.metrics$TrimmedSequence[j],
                  snppos + 2,
                  500
                )
              )
          } else {
            seq[i] <-
              paste(rep("N", nchar(
                as.character(
                  x@other$loc.metrics$TrimmedSequence[j]
                )
              )), collapse = "")
          }
      }

      cat(paste0(nInd(x)," ",nchar(seq[1])),"\n")
      cat(paste0(locNames(x)[j],"^",indNames(x)," ",seq,"\n"))

    }
    
    # Close the output file
    sink()
  }
    
    # METHOD = RANDOM ASSIGNMENT
    
    if (method == 2 ) {
      # Randomly allocate heterozygotes (1) to homozygote state (0 or 2)
      matrix <- as.matrix(x)
      # cat('Randomly allocating heterozygotes (1) to homozygote state (0 or 2)\n') pb <- txtProgressBar(min=0, max=1, style=3,
      # initial=0, label='Working ....') getTxtProgressBar(pb)
      r <- nrow(matrix)
      c <- ncol(matrix)
      for (i in 1:r) {
        for (j in 1:c) {
          if (matrix[i, j] == 1 && !is.na(matrix[i, j])) {
            # Score it 0 or 2
            matrix[i, j] <- (sample(1:2, 1) - 1) * 2
          }
        }
      }
      
      if (verbose >= 2) {
        cat(report("Generating haplotypes ... This may take some time\n"))
      }
      
      sink(outfilespec)
      
      # For each individual, and for each locus, generate the relevant haplotype

      for (j in 1:c) {
        seq <- NA
        # Reassign some variables
        trimmed <- as.character(x@other$loc.metrics$TrimmedSequence[j])
        snp <- x@loc.all[j]
        snpos <- x@position[j]
        # Shift the index for snppos to start from 1 not zero
        snpos <- snpos + 1
        for (i in 1:r) {
          # If the score is homozygous for the reference allele
            if (matrix[i, j] == 0 && !is.na(matrix[i, j])) {
              seq[i] <- trimmed
            }
   
          # If the score is homozygous for the alternate allele
          if (matrix[i, j] == 2 && !is.na(matrix[i, j])) {
            # Split the trimmed into a beginning sequence, the SNP and an end sequences
            start <- stringr::str_sub(trimmed, end = snpos - 1)
            snpbase <- stringr::str_sub(trimmed,
                               start = (snpos),
                               end = (snpos))
            end <- stringr::str_sub(trimmed, start = snpos + 1)
            # Extract the SNP transition bases (e.g. A and T)
            state1 <- gsub("(.)/(.)", "\\1", snp, perl = TRUE)
            state2 <- gsub("(.)/(.)", "\\2", snp, perl = TRUE)
            # Change the SNP state to the alternate
            if (snpbase == state1) {
              snpbase <- state2
            } else {
              snpbase <- state1
            }
            
              # Paste back to form the alternate fragment
              target <- paste0(start, snpbase, end)
              # Remove adaptors and save the trimmed alternate sequence
              seq[i] <-  stringr::str_sub(target,
                                 start = 1,
                                 end = nchar(trimmed))
          }
          
          # If the SNP state is missing, assign NNNNs
          if (is.na(matrix[i, j])) {
            seq[i] <- "N"
              seq[i] <-
                stringr::str_pad(
                  seq[i],
                  nchar(trimmed),
                  side = c("right"),
                  pad = "N"
                )
          }
        }

        cat(paste0(nInd(x)," ",nchar(seq[1])),"\n")
        cat(paste0(locNames(x)[j],"^",indNames(x)," ",seq,"\n"))
        
      }  
      
      # Close the output fastA file
      sink()
      
    }
    
    # Imap file
    sink(outfilespec_imap)
    cat(paste(indNames(x), pop(x),"\n"))
    sink()
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  return(NULL)
  
}
