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

gl2gds <- function(x, outfile="gl2gds.gds", outpath=tempdir(), verbose=NULL) {

# TIDY UP FILE SPECS

  outfilespec <- file.path(outpath, outfile)
  funname <- match.call()[[1]]

# FLAG SCRIPT START
  # set verbosity
  if (is.null(verbose) & !is.null(x@other$verbose)) verbose=x@other$verbose
  if (is.null(verbose)) verbose=2
 

  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }

  if (verbose > 0) {
    cat("Starting",funname,"\n")
  }

# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }

  if (verbose >= 2){
    if (all(x@ploidy == 1)){
      stop("Fatal Error: Detected Presence/Absence (SilicoDArT) data. Please provide a SNP dataset\n")
    } else if (all(x@ploidy == 2)){
      cat("  Processing a SNP dataset\n")
    } else {
      stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
    }
  }

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

