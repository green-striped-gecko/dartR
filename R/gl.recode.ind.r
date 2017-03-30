#' Recode individual (=specimen) labels in a genelight or genind object \{adegenet\}
#'
#' This script recodes individual labels and/or deletes individuals from a DaRT genlight SNP file or a SilicoDArT genind file
#' based on information provided in a csv file.
#'
#' Renaming individuals may be required when there have been errors in labelling arising
#' in the process from sample to DArT files. There may be occasions where renaming
#' individuals is required for preparation of figures. Caution needs to be exercised
#' because of the potential for breaking the "chain of evidence" between the samples themselves
#' and the analyses. Recoding individuals can be done with a recode table (csv).
#' 
#' The script, having deleted individuals, identifies resultant monomorphic loci or loci
#' with all values missing and deletes them (using gl.filter.monomorphs.r)
#'
#' @param gl -- name of the genlight object containing SNP genotypes or a genind object containing presence/absence data [required]
#' @param ind.recode -- name of the csv file containing the individual relabelling [required]
#' @return A genlight or genind object with the recoded and reduced data
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' \dontrun{
#'    gl <- gl.recode.ind(gl, ind.recode="ind_recode_table_0.csv")
#'    glind <- gl.recode.ind(glind, ind.recode="ind_recode_table_0.csv")
#' }
#' @seealso \code{\link{gl.filter.monomorphs}}
#' 
#'

gl.recode.ind <- function(gl, ind.recode){
x <- gl

  if(class(x)!="genind" & class(x)!="genlight") {
    cat("Fatal Error: genind or genlight object required for gl.recode.ind.r!\n"); stop()
  }

# RELABEL INDIVIDUALS
  cat("Processing",class(x),"object\n")
  cat("  Relabelling individuals (=specimens) as per ", ind.recode, "\n")
  recode.table <- read.csv(ind.recode, stringsAsFactors=FALSE, header=FALSE);
# Error check
  if(length(unique(indNames(x))) != length(unique(recode.table[,1]))) {
    cat("Fatal Error: Individual names in data file are not the same as in the recode table\n"); stop()
  }
# Apply the recode to the individuals
  ind.list <- as.character(indNames(x))
  ntr <- length(recode.table[,1])
  for (i in 1:nInd(x)) {
    for (j in 1:ntr) {
      if (ind.list[i]==recode.table[j,1]) {ind.list[i] <- recode.table[j,2]}
    }
  }
  indNames(x) <- ind.list

# Remove rows flagged for deletion
  cat("  Removing entities flagged for deletion in ", ind.recode, "\n")
  x2 <- x[!x$ind.names=="delete" & !x$ind.names=="Delete"]
  
  x2 <- gl.filter.monomorphs(x2)

# REPORT A SUMMARY
  cat("Summary of recoded dataset\n")
  cat(paste("  No. of loci:",nLoc(x2),"\n"))
  cat(paste("  No. of individuals:", nInd(x2),"\n"))
  cat(paste("  No. of populations: ", length(levels(factor(pop(x2)))),"\n"))

    return(x2)
}

