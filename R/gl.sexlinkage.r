#' Identify loci that are sex linked in specimens in a genlight \{adegenet\} object
#'
#' Alleles unique to the Y or W chromosome and monomorphic on the X chromosomes will appear in the SNP dataset as 
#' genotypes that are heterozygotic in all individuals of the heterogametic sex and homozygous in all individuals 
#' of the homogametic sex.
#' 
#' This script will identify loci with alleles that behave in this way, as putative sex specific SNP markers.
#'
#' @param gl -- name of the genlight object containing the SNP data [required]
#' @param tm -- tolerance, that is tm=0.95 means that only 95% of males need be heterozygous in an XY system (females
#' in a ZW system) [default 0]
#' @param tf -- tolerance, that is tf=0.95 means that only 95% of females need be homozygous in an XY system (males
#' in a ZW system) [default 0]
#' @param v -- verbosity: 0, silent; 1, brief; 2, verbose [default 1]
#' @return The list of sex specific loci
#' @export
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' result <- gl.sexlinkage(testset.gl)

gl.sexlinkage <- function(gl, tm=0, tf=0, v=1) {
  x <- gl

# Extract the data for the females
  matf <- as.matrix(x[x@other$ind.metrics$Sex=="F"])
  # For each individual
    f <- array(data=NA, dim=c(ncol(matf),3))
    for (i in 1:ncol(matf)) {
      for (j in 1:3) {
        f[i,j] <- length(which(matf[,i]==(j-1)))
      }  
    }
    dff <- data.frame(f)
    row.names(dff) <- locNames(x)
    colnames(dff) <- c("F0","F1","F2")

# Extract the data for the males
  matm <- as.matrix(x[x@other$ind.metrics$Sex=="M"])
# For each individual
  m <- array(data=NA, dim=c(ncol(matm),3))
  for (i in 1:ncol(matm)) {
    for (j in 1:3) {
      m[i,j] <- length(which(matm[,i]==(j-1)))
    }  
  }
  dfm <- data.frame(m)
  row.names(dfm) <- locNames(x)
  colnames(dfm) <- c("M0","M1","M2")

# Combine the two files
df <- cbind(dff,dfm)

# Check for hets in all males, homs in all females (XY); ditto for ZW
zw <- df[df$F1/(df$F0+df$F1+df$F2)==(1-tf) && df$M1/(df$M0+df$M1+df$M2)==(0+tm),]
xy <- df[df$F1/(df$F0+df$F1+df$F2)==(0+tf) && df$M1/(df$M0+df$M1+df$M2)==(1-tm),]

if (nrow(zw) == 0){
  cat("No sex linked markers consistent with female heterogamety (ZZ/ZW)\n")
} else {
  cat("Sex linked loci consistent with female heterogamety (ZZ/ZW)\n")
  cat("  0 = homozygous reference; 1 = heterozygous; 2 = homozygous alternate\n")
  print(zw)
}
if (nrow(zw) == 0){
  cat("No sex linked markers consistent with male heterogamety (XX/XY)\n")
} else {
  cat("Sex linked loci consistent with male heterogamety (XX/XY)\n")
  cat("  0 = homozygous reference; 1 = heterozygous; 2 = homozygous alternate\n")
  print(xy)
}

list <- c(zw,xy)

return(list)

}