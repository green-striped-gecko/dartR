#' Generate percentage allele frequencies by locus and population
#'
#' This is a support script, to take SNP data or SilocoDArT presence/absence data grouped into populations in a genelight or genind object \{adegenet\}
#' and generate a table of allele frequencies for each population and locus
#'
#' @param gl -- name of the genlight containing the SNP genotypes or genind object containing the presence/absence data [required]
#' @return A matrix with allele frequencies (genelight) or presence/absence frequencies (genind) broken down by population and locus
#' @export
#' @import adegenet
#' @importFrom data.table melt
#' @importFrom plyr ddply
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' m <-  gl.percent.freq(testset.gl)
#' m

gl.percent.freq<- function(gl) {
x <- gl

# Determine data type
  if(class(x)=="genlight"){
    cat("Using SNP data from a genlight object\n")
    m <- as.matrix(x)
  } else if(class(x)=="genind"){
    cat("Using RFP data (SilicoDArT) from a genind object\n")
    m <- as.matrix(x)
  } else {
    cat("Fatal Error: Input data must be a genlight or genind object\n")
  }
  
# Assign population names to the individuals, discarding specimen numbers
  rownames(m) <- pop(x)
    
# Convert the SNP data to long format
  m.long <- reshape2::melt(m, na.rm=FALSE)
  colnames(m.long) <- c("popn", "locus", "snp")
    
# Calculate sums and counts broken down by population and locus
  cat("  Tallying allele frequencies, this may take some time\n")
  m.summed<-plyr::ddply(
            m.long,
            .(popn,locus),
            summarize,
            sums=sum(as.numeric(as.character(snp)),na.rm=TRUE),
            count=sum(!is.na(snp)),
            missing=sum(is.na(snp))
            )
  names(m.summed)<-c("popn","locus","sum","nobs","nmissing")
  
# Calculate some new variables, and in particular, percentage frequencies
  if(class(x)=="genlight") {m.summed$frequency<-m.summed$sum*100/(2*m.summed$nobs)}
  if(class(x)=="genind") {m.summed$frequency<-m.summed$sum*100/m.summed$nobs}
  m.summed$frequency[is.nan(m.summed$frequency)]<-NA
  m.summed$n<-m.summed$nobs+m.summed$nmissing
  
# Sort the data on locus and population in preparation for analysis
  m.summed <- m.summed[order(m.summed$locus, m.summed$popn),]

  cat("  Calculation of allele frequencies complete\n")
  return(m.summed)
  
}
