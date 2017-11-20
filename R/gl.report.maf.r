#' Report minor allele frequency for loci
#'
#' Minor allele frequencies are reported for every locus. Also the number of individuals that have been used to calculate 
#' allele frequencies are reported to enable one to discard frequencies of loci where only a few individuals are available.
#'
#' @param gl -- name of the genlight object containing the SNP data [required]
#' @param minind -- minimum number of individuals required to calculate the minor allele frequency. Please note allelefrequncies for such loci are still reported, but a warning is issued if the number is below minind.
#' @param plot -- switch if a histogram of the distribution of allele frequencies should be provided
#' @return the functions returns a list with two slots (vectors). maf contains the minor allele frequency and n the number of individuals the minor allele frequency is based on.
#' @export
#' @author Bernd Gruber (glbugs@aerg.canberra.edu.au)
#' @examples
#' frequ <- gl.report.maf(testset.gl)
#' summary(frequ$maf)
#' summary(frequ$n)


gl.report.maf <- function(gl, minind=5, plot=TRUE)
{
  af <- (colMeans(as.matrix(gl), na.rm = T)/2 )
  af <- ifelse(af>0.5, 1-af, af)
  minn <- colSums(!is.na(as.matrix(gl)))
  
  if (sum(is.na(af))>0)
  {
    cm <- which(is.na(af))
    
    warning(paste(length(cm), "loci have just missing data."))
    cat(names(cm))
  }
  
  if (sum(minn<minind)>0) warning(paste(sum(minn<minind),"Allele frequencies are based on less than",minind,"individuals"))
  
  if (plot) hist(af, breaks=seq(0,0.5,0.05), col=rainbow(11), main="", xlab="Minor Allele Frequency")
  return(list(maf=af, n=minn))
}
