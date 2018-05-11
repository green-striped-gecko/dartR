#' Create a proforma recode_pop_table file for reassigning population names
#'
#' Renaming populations may be required when there have been errors in assignment arising
#' in the process from sample to DArT files or when one wishes to amalgamate populations, or delete populations.
#' Recoding populations can also be done with a recode table (csv).
#'
#' This script facilitates 
#' the construction of a recode table by producing a proforma file with
#' current population names in two identical columns. Edit the second
#' column to reassign populations. Use keyword Delete to delete a population.
#' 
#' Apply the recoding using gl.recode.pop(). 
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param outfile -- name of the new proforma file [default recode_pop_table.csv]
#' @param outpath -- path where to save the output file (set to tempdir by default)
#' @return A vector containing the new population names
#' @export
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' \donttest{
#' result <- gl.make.recode.pop(testset.gl, outfile="Emmac_recode_pop.csv")
#' }

 gl.make.recode.pop <- function(x, outfile="recode_pop_table.csv", outpath=tempdir()) {

 mat <- cbind(levels(pop(x)),levels(pop(x)))
 write.table(mat, file=file.path(outpath, outfile), sep=",", row.names=FALSE, col.names=FALSE)
 
 cat(paste("Proforma recode table written to:",outfile,"\n"))
 
 return(levels(pop(x)))
 
 }
