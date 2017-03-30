#' Write out data from a gl object \link{adegenet} to csv file
#'
#' This script writes to file the SNP genotypes with specimens as entities (columns) and
#' loci as attributes (rows). Each row has associated locus metadata. Each column, with header
#' of specimen id, has population in the first row. 
#' 
#' The data coding differs from the DArT 1 row format in that 0 = reference homozygous, 2 =
#' alternate homozygous, 1 = heterozyous, and NA = missing SNP assignment. 
#'
#' @param gl -- name of the genlight object containing SNP genotypes [required]
#' @param outfile -- name of the csv file to write the data to [required]
#' @return saves a glenlight object to csv
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' \dontrun{
#' gl.write.csv(gl, outfile="SNP_1row.csv")
#' }

gl.write.csv <- function(gl, outfile=NULL) {
x <- gl

# Add individual names and population names to rows
  x1 <- cbind.data.frame(pop(x),as.matrix(x))
# Transpose to have id as columns, loci as rows  
  x1 <- t(x1)
# Create two filler rows to bring number of data rows and number of locus.metadata rows together
  filler1 <- rep("*",length(x@other$loc.metrics[1,]))
# Bind the filler rows to the locus metadata
  x2 <- rbind(filler1, as.matrix(x@other$loc.metrics))
# Bind the locus metadata to the data  
  x3 <- cbind.data.frame(x2,x1)
# Output   
  write.table(x3, file=outfile,sep=",",row.names=FALSE)
  
  a <- NULL
  return(a)
}
