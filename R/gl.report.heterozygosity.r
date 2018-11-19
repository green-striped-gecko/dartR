#' Reports hetrozygosity
#' 
#' Calculates the observed heterozygisities by population from a genlight object
#' and plots as a barchart ordered on heterozygosity.
#' 
#' #Input is a genlight {adegenet} object containing SNP genotypes (0 homozygous for reference SNP, 
#' #1 heterozygous, 2 homozygous for alternate SNP).
#' 
#' @param gl -- a genlight object containing the SNP genotypes [Required]
#' @return a dataframe containing population labels, observed heterozygosities and sample sizes
#' @author Bernd Gruber (bugs? Post to https://groups.google.com/d/forum/dartr)
#' @export
#' @importFrom grDevices rainbow
#' @importFrom graphics par
#' @examples
#' gl.report.heterozygosity(testset.gl)
#
# Amended 9-Mar-17

gl.report.heterozygosity <- function(gl) {
x <- gl
  
  # Split the genlight object into a list of populations
  sgl <- seppop(x)
  # Calculate heterozygosity for each population in the list
  hs <- data.frame(lapply(sgl, function(x) mean(colMeans(as.matrix(x, na.rm=TRUE)==1), na.rm=TRUE) ))
  # Calculate sample sizes
  sums <- data.frame(lapply(sgl, function(x) mean(colSums(as.matrix(x, na.rm=TRUE)==1), na.rm=TRUE) ))
  n <- t(sums/hs)
  n <- cbind(row.names(n),n)
  # Join the sample sizes with the heteozygosities
  d1 <- data.frame(names=names(hs), hs=as.numeric(hs[1,]))
  d2 <- data.frame(n)
  names(d2)<- c("names","freq")
  dd <- join(d1,d2, by="names")
  # Order by heterozygosity
  dd <- dd[order(dd$hs),]
  # Plot the results
  par(mai=c(2.5,1,0.5,0.2))
  barplot(dd$hs, names.arg=paste(dd$names, dd$freq, sep=" | "), las=2, cex.names=1, space=0, border=F, col=rainbow(nrow(dd)), main="Observed Heterozygosity")

  # Return the result
  return(dd) 
}

