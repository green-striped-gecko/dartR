#' Convert a genlight object to nexus format for PhyloNet network analysis of biallelic markers
#'
#' The output nexus file contains the snp data in the format expected of PhyloNet in order to undertake a network analysis.
#' Note that you will need to consider and adjust the parameters of the output file in consultation with the PhyloNet documentation. 
#' 
#' Reference: Zhu, J, When, D., Yu, Y., Meudt, H.M., Nakhley, L., and Matsen, F.A. 2018. Bayesian inference of phylogenetic networks from bi-allelic genetic markers, PLOS Computational Biology, 14: e1005932
#' 
#' 
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param outfile -- file name of the output file (including extension).
#' @param cl -- PhyloNet parameter to include in the outfile: length of the MCMC chain [default 500,000]
#' @param bl -- PhyloNet parameter to include in the outfile: number of iterations in burn-in period [default 200,000]
#' @param pl -- PhyloNet parameter to include in the outfile: number of cores available for the analysis [default 4]
#' @param mr -- PhyloNet parameter to include in the outfile:The maximum number of reticulation nodes in the sampled phylogenetic networks [default 4]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return NULL
#' @export
#' @importFrom
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' gl2svdquartets(testset.gl)

gl2phylonet <- function(x, outfile="phylonet.nex", cl=500000, bl=200000, pl=4, mr=4, v=2) {
  
  if (v > 0) {cat(paste("Starting gl2phlyonet: Create nexus file\n\n"))}
  if (v > 1) {cat(paste("    Extacting SNP states (0,1,2,?) and creating records for each individual\n"))}


# Extract the SNP data, and order on population
  m <- as.matrix(x)
  m <- m[order(pop(x)),]
  pop(x) <- pop(x)[order(pop(x))]
  m <- replace(m, is.na(m), "?")
  a <- paste0(indNames(x),"  ")
  m <- cbind(a,m)
  
# Construct the taxon list (list of individuals = taxa)
  taxalist <- NULL 
  for (i in 1:length(indNames(x))){ 
    taxalist <- paste0(taxalist,indNames(x)[i],",")
  }  
  taxalist <- gsub(".$", "", taxalist) 
  
# Construct the population list (=tm) 
  for (i in 1:length(indNames(x))){
    if (i == 1){
      tm <- paste0("<",as.character(pop(x)[i]),":",indNames(x)[i])
    } else {
      if (pop(x)[i] == pop(x)[i-1]) {
        tm <- paste0(tm,",",indNames(x)[i])
      } else {
        tm <- paste0(tm,"; ",as.character(pop(x)[i]),":",indNames(x)[i])
      }
    }
  }  
  tm <- paste0(tm,">")

# Create the nexus file
  if (v > 1) {cat(paste("    Writing results to nexus file",outfile,"\n"))}
  sink(outfile)
  cat("#nexus\n")
  cat("BEGIN DATA;\n")
  cat(paste0("Dimensions ntax=",nInd(x)," nchar=",nLoc(x),";\n"))
  cat("Format datatype=dna symbols=\"012\" missing=? gap=-;\n")
  cat("Matrix\n\n")
  
  for (i in 1:nInd(x)){
    cat(paste0(m[i,],collapse=""))
    cat("\n")
  }
  cat(";\n")
  cat("End;\n\n")
  
  cat("BEGIN PHYLONET;\n") 
  cat("MCMC_BiMarkers\n")
  cat("-cl",format(cl, scientific=FALSE),"\n")
  cat("-bl",format(bl, scientific=FALSE),"\n")
  cat("-sf 500\n") 
  cat("-prebl 5000\n") 
  cat("-premc3 (2.0,4.0)\n") 
  cat("-premr 1\n") 
  cat("-diploid\n")
  cat("-op\n") 
  cat("-varytheta\n") 
  cat("-pp 1.5\n") 
  cat("-ee 2.0\n") 
  cat("-mr 2\n") 
  cat("-pl",pl,"\n") 
  cat("-ptheta 0.1\n") 
  cat("-sd 12345678\n") 
  cat(paste0("-taxa (",taxalist,")\n")) 
  cat(paste0("-tm ",tm,";\n"))
  cat("END;\n")
  
  sink()
  
  if (v > 2) {cat(paste("    Records written to",outfile,":",nInd(x),"\n"))}
  if (v > 2) {cat(paste("    Log written to",logfile,"\n"))}
  if (v > 2) {cat(paste("    Trees written to svd_boot.tre\n"))}
  if (v > 0) {
    cat("gl2phylonet Completed\n")
    cat("You should edit",outfile,"to amend options as necessary.\n")
    cat("Note: the PhyloNet nexus file is very pernickity about case and indenting (Matrix not matrix)\n")
  }
  
  return(NULL)

}


