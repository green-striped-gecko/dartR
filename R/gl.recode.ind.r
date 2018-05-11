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
#' The script, having deleted individuals, optionally identifies resultant monomorphic loci or loci
#' with all values missing and deletes them (using gl.filter.monomorphs.r). The script also optionally
#' recalculates statistics made redundant by the deletion of individuals from the dataset.
#' 
#' The script returns a genlight object with the new individual labels and the recalculated locus metadata.
#'
#' @param gl -- name of the genlight object containing SNP genotypes or a genind object containing presence/absence data [required]
#' @param ind.recode -- name of the csv file containing the individual relabelling [required]
#' @param recalc -- Recalculate the locus metadata statistics if any individuals are deleted in the filtering [default TRUE]
#' @param mono.rm -- Remove monomorphic loci [default TRUE]
#' @param v -- verbosity: 0, silent; 1, brief; 2, verbose [default 1]
#' @return A genlight or genind object with the recoded and reduced data
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' \donttest{
#'    gl <- gl.recode.ind(testset.gl, ind.recode="testset_pop_recode.csv")
#' }
#' @seealso \code{\link{gl.filter.monomorphs}}
#' 
#'

gl.recode.ind <- function(gl, ind.recode, recalc=TRUE, mono.rm=TRUE, v=1){
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

  # If there are individuals to be deleted, then recalculate relevant locus metadata and remove monomorphic loci
  
  if ("delete" %in% x$ind.names | "Delete" %in% x$ind.names) {
    # Remove rows flagged for deletion
      cat("Deleting individuals or samples flagged for deletion\n")
      x2 <- x[!x$ind.names=="delete" & !x$ind.names=="Delete"]
    #  Remove monomorphic loci
      if (mono.rm) {x2 <- gl.filter.monomorphs(x2,v=v)}
      if (recalc) {
    # Recalculate statistics
        x2 <- utils.recalc.avgpic(x2,v=v)
        x2 <- utils.recalc.callrate(x2,v=v)
        x2 <- utils.recalc.freqhets(x2,v=v)
        x2 <- utils.recalc.freqhomref(x2,v=v)
        x2 <- utils.recalc.freqhomsnp(x2,v=v)
      }
  }

# REPORT A SUMMARY
  if (v==2) {
    cat("Summary of recoded dataset\n")
    cat(paste("  No. of loci:",nLoc(x2),"\n"))
    cat(paste("  No. of individuals:", nInd(x2),"\n"))
    cat(paste("  No. of populations: ", length(levels(factor(pop(x2)))),"\n"))
    if (!recalc) {cat("Note: Locus metrics not recalculated\n")}
    if (!mono.rm) {cat("note: Resultant monomorphic loci not deleted\n")}
  }
  if (v==1) {  
    if (!recalc) {cat("Note: Locus metrics not recalculated\n")}
    if (!mono.rm) {cat("note: Resultant monomorphic loci not deleted\n")}
  }

    return(x2)
}

