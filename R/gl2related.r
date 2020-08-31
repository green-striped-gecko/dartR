#' Convert a genlight object to format suitable to be run with Coancestry
#'
#' The output txt file contains the snp data and an additional column with the names if the individual. The file then can be used and loaded into coancestry or - if installed - run with the related package. Be aware the related package was crashing in previous versions, but in general is using the same code as coancestry and therefore should have identical results. Also running coancestry with thousands of SNPs via the GUI seems to be not reliable and therefore for comparisons between coancestry and related we suggest to use the command line version of coancestry.
#' 
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param outfile -- file name of the output file (including extension) [default 'related.txt']
#' @param outpath -- path where to save the output file [default tempdir()]
#' @param save -- a switch if you want to save the file or not. This might be useful for someone who wants to use the coancestry function to calculate relatedness and not export to coancestry. See the example below. [default TRUE]
#' @return a data.frame that can be used to run with the related package
#' @export
#' @author Bernd Gruber (bugs? Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gtd <- gl2related(bandicoot.gl[1:10,1:20], save=FALSE)
#' \dontrun{
#' ##running with the related package
#' #install.packages("related", repos="http://R-Forge.R-project.org")
#' library(related)
#' coan <- coancestry(gtd, wang=1) 
#' head(coan$relatedness)
#' ##check ?coancestry for information how to use the function.
#' }
#' @references Jack Pew, Jinliang Wang, Paul Muir and Tim Frasier (2014). 
#' related: related: an R package for analyzing pairwise relatedness
#'  data based on codominant molecular markers. 
#'  R package version 0.8/r2.
#'   \url{https://R-Forge.R-project.org/projects/related/}




gl2related <- function(x, outfile="related.txt", outpath=tempdir(), save=TRUE)
{
gd <- as.matrix(x)

mm <- matrix(NA, nrow=nrow(gd), ncol =ncol(gd)*2 )
mm[, seq(1,nLoc(x)*2,2)] <- gd
mm <- ifelse(mm==1,0,mm) 

mm[, seq(2,nLoc(x)*2,2)] <- gd
mm <- ifelse(mm==1,2,mm) 

mm <- mm+1
mm <- ifelse(is.na(mm),0,mm)

gtd <- data.frame(V1=as.character(indNames(x)), mm)
gtd$V1 <- as.character(gtd$V1)

fn <- file.path(outpath, outfile)
#export as file
if (save) write.table(gtd, file = fn, sep = "\t",col.names = F, row.names = F)

return(gtd)

}
