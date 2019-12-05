#' Convert a genlight object to gds format
#'
#' Package SNPRelate relies on a bit-level representation of a SNP dataset that competes with \{adegenet\} genlight
#' objects and associated files. This function saves a genlight object to a gds format file.
#' 
#' 
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param outfile -- file name of the output file (including extension) [default gl2gds.gds]
#' @param outpath -- path where to save the output file [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath="." when calling this function to direct output files to your working directory.
#' @param verbose -- specify the level of verbosity: 0, silent, fatal errors only; 1, flag function begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return NULL
#' @export
#' @importFrom SNPRelate snpgdsCreateGeno snpgdsOpen snpgdsSummary snpgdsClose
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' \donttest{
#' gl2gds(testset.gl)
#' }

# Last amended 3-Feb-19

gl2gds <- function(x, outfile="gl2gds.gds", outpath=tempdir(), verbose=2) {

# TIDY UP FILE SPECS

  outfilespec <- file.path(outpath, outfile)
  funname <- match.call()[[1]]

# FLAG SCRIPT START

  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }

  if (verbose > 0) {
    cat("Starting",funname,"\n")
  }

# STANDARD ERROR CHECKING
  
  if(!is(x, "genlight")) {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }

  # Work around a bug in adegenet if genlight object is created by subsetting
      if (nLoc(x)!=nrow(x@other$loc.metrics)) { stop("The number of rows in the loc.metrics table does not match the number of loci in your genlight object!")  }

  # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nInd(x))
      pop(x) <- as.factor(pop(x))
    }

  # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x, verbose=0)
    if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {cat("  Warning: genlight object contains monomorphic loci\n")}

# DO THE JOB

# Shift snp position by 1 (DArT starts at position 0; SNPRelate starts at position 1)
  snp.pos <- x@position + 1
  
# Convert any NA values to 0 (genlight objects have NA for missing; SNPRelate has 0 in this instance)  
  for (i in 1:nLoc(x)) {
    if(is.na(snp.pos[i])) {snp.pos[i] <- 0}
  }

# Create the gds file
  if (verbose >= 2) {cat("  Converting SNP data to gds format\n")}
  SNPRelate::snpgdsCreateGeno(gds.fn=outfilespec,
                         genmat = as.matrix(x),
                         sample.id = indNames(x),
                         snp.id = locNames(x),
                         snp.rs.id = x@other$loc.metrics$AlleleID,
                         snp.chromosome = x@chromosome,
                         snp.position = snp.pos,
                         snp.allele = x@loc.all,
                         other.vars = x@other,
                         snpfirstdim=FALSE)
  
# Open the GDS file, which will print out a summary of contents
  if (verbose >= 2 ) {cat("  Writing data to file",outfilespec,"\n")}
  genofile <- SNPRelate::snpgdsOpen(outfilespec)
  cat("Structure of gds file\n\n")
  SNPRelate::snpgdsSummary(genofile)
  print(genofile)

# Close the GDS file
  if (verbose >= 2 ) {cat("  Closing file",outfilespec,"\n")}
  SNPRelate::snpgdsClose(genofile)

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

  return(NULL)

}

