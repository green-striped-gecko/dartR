#' Convert a genlight object to nexus format suitable for phylogenetic analysis by SNAPP (via BEAUti)
#'
#' The output nexus file contains the snp data and relevant PAUP command lines suitable for BEAUti.
#' 
#' Reference: Bryant, D., Bouckaert, R., Felsenstein, J., Rosenberg, N.A. and RoyChoudhury, A. (2012). Inferring species trees directly from biallelic genetic markers: bypassing gene trees in a full coalescent analysis. Molecular Biology and Evolution 29:1917-1932.
#' 
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param outfile -- file name of the output file (including extension).
#' @param outpath -- path where to save the output file [default tempdir()]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return NULL
#' @export
#' @importFrom
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' gl2snapp(testset.gl)

gl2snapp <- function(x, outfile="snapp.nex", outpath=tempdir(), v=2) {
  
  if (v > 0) {cat(paste("Starting gl2snapp: Create nexus file\n\n"))}
  if (v > 1) {cat(paste("    Extacting SNP data and creating records for each individual\n"))}

# Extract the reference base and the alternate base for each locus (excuse the contortion)
  
  m <- as.matrix(x)
  m[is.na(m)] <- "?"
  colnames(m) <- NULL
  df <- data.frame(m)
  df <- cbind(indNames(x),pop(x),df)
  df <- df[order(df$pop),]
  indlabels <- df[,1]
  poplabels <- df[,2]
  
# Create the snapp file
  if (v > 1) {cat(paste("    Writing results to nexus file",outfile,"\n"))}

  sink(outfile)
  
  cat("#nexus\n")
  cat("BEGIN DATA;\n")
  cat(paste0("     dimensions ntax = ",nInd(x)," nchar = ",nLoc(x)," ;\n"))
  cat("     format datatype=integerdata missing=? symbols=\"012\";\n")
  cat("matrix\n")
    for (i in 1:nInd(x)){
      cat(paste0(poplabels[i],"_",indlabels[i]))
      cat("  ")
      cat(m[i,], sep="")
      cat("\n")
    }
  cat(";\n")
  cat("end;\n\n")

  sink()
  
  if (v > 2) {cat(paste("    Records written to",outfile,":",nInd(x),"\n"))}
  if (v > 0) {cat("\ngl2snapp Completed\n")}
  
  return(NULL)

}


