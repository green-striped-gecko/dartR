#' Convert a genlight object to nexus format PAUP SVDquartets
#'
#' Package SNPRelate relies on a bit-level representation of a SNP dataset that competes with \{adegenet\} genlight
#' objects and associated files. This function saves a genlight object to a gds format file.
#' 
#' 
#' @param gl -- name of the genlight object containing the SNP data [required]
#' @param outfile -- file name of the output file (including extension).
#' @return NULL
#' @export
#' @importFrom SNPRelate snpgdsCreateGeno snpgdsOpen snpgdsSummary snpgdsClose
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' \dontrun{
#' gl2gds(testset.gl)
#' }

gl2gds <- function(gl, outfile="gl2gds.gds") {
  
  cat(paste("Converting gl object to gds formatted file", outfile, "\n\n"))

# Load the R packages: gdsfmt and SNPRelate
#
## try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("SNPRelate")
# library(gdsfmt)
# library(SNPRelate)
# library(adegenet)
  
# Shift snp position by 1 (DArT starts at position 0; SNPRelate starts at position 1)
  snp.pos <- gl@position + 1
  
# Convert any NA values to 0 (genlight objects have NA for missing; SNPRelate has 0 in this instance)  
  for (i in 1:nLoc(gl)) {
    if(is.na(snp.pos[i])) {snp.pos[i] <- 0}
  }

# Create the gds file
  snpgdsCreateGeno(gds.fn=outfile,
                         genmat = as.matrix(gl),
                         sample.id = indNames(gl), 
                         snp.id = locNames(gl),
                         snp.rs.id = gl@other$loc.metrics$CloneID,
                         snp.chromosome = gl@chromosome,
                         snp.position = snp.pos,
                         snp.allele = gl@loc.all,
                         other.vars = gl@other,
                         snpfirstdim=FALSE)
  
# Open the GDS file, which will print out a summary of contents
  genofile <- snpgdsOpen(outfile)
  cat("Structure of gds file\n\n")
  snpgdsSummary(genofile)
  print(genofile)

# Close the GDS file
  snpgdsClose(genofile)

  return(NULL)

}

