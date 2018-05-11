#' Converts genlight objects to STRUCTURE formated files
#'
#' This function exports genlight objects to STRUCTURE formatted files (be aware there is a gl2faststruture version as well). It is based on the code provided by Lindsay Clark (see \url{https://github.com/lvclark/R_genetics_conv})  and this function is basically a wrapper around her numeric2structure function. See also: Lindsay Clark. (2017, August 22). lvclark/R_genetics_conv: R_genetics_conv 1.1 (Version v1.1). Zenodo: \url{http://doi.org/10.5281/zenodo.846816}.
#' 
#' @param gl -- genlight containing lat longs [required]
#' @param indNames -- switch if individuals names should be added (defult to indNames in gl)
#' @param addcolumns -- additional columns to be added  before genotypes
#' @param ploidy -- defaults to 2
#' @param exportMarkerNames -- switch if loci names should be included (locNames(gl))
#' @param outfile -- name (path) of the output shape file
#' @param outpath -- path of the output file. Default is to tempdir(). If to be saved in the current working directory change to "."
#' @param v -- verbosity: if v=0 no output, v=1 reports name and path of output file. default 1
#' @export
#' @author Bernd Gruber (wrapper) and Lindsay V. Clark [lvclark@illinois.edu] 
#' @examples
#' \donttest{
#' gl2structure(testset.gl)
#'}


gl2structure <- function(gl, indNames =NULL, addcolumns = NULL, ploidy = 2,exportMarkerNames = TRUE, outfile="gl.str", outpath=tempdir(), v=1)
 {
   if(!"genlight" %in% class(gl)){
     stop("Function was designed for genlight objects.")
   }
nInd <- nInd(gl)
if (is.null(indNames)) indNames=indNames(gl)
if(length(indNames) != nInd){
stop("Number of individuals does not match between indNames and genmat.")
}

if (!is.null(addcolumns) && is.null(dim(addcolumns))) addcolumns <- data.frame(pop=addcolumns)

if(!is.null(addcolumns) && nrow(addcolumns) != nInd){
stop("Number of individuals does not match between addColumns and gl.")
}
genmat <- as.matrix(gl)
if(!all(genmat %in% c(0:ploidy,NA))){
stop("genmat must only contain 0, 1, 2... ploidy and NA")
}
if(length(outfile) != 1 || !is.character(outfile)){
stop("file must be a single character string.")
}
if(length(ploidy) != 1 || !is.numeric(ploidy)){
stop("ploidy must be a single number")
}
if(!exportMarkerNames %in% c(TRUE, FALSE)){
stop("exportMarkerNames must be TRUE or FALSE")
}
  
  # make sets of possible genotypes
  G <- list()
  for(i in 0:ploidy){
  G[[i + 1]] <- c(rep(1, ploidy - i), rep(2, i))
  }
  G[[ploidy + 2]] <- rep(-9, ploidy) # for missing data
  
  # set up data frame for Structure
  StructTab <- data.frame(ind = rep(indNames, each = ploidy))
  # add any additional columns
  if(!is.null(addcolumns)){
  
  for(i in 1:dim(addcolumns)[2]){
 StructTab <- data.frame(StructTab, rep(addcolumns[,i], each = ploidy))
 if(!is.null(dimnames(addcolumns)[[2]])){
 names(StructTab)[i + 1] <- dimnames(addcolumns)[[2]][i]
 } else {
 names(StructTab)[i + 1] <- paste("X", i, sep = "")
 }
  }
  }
  
  # add genetic data
  for(i in 1:dim(genmat)[2]){
  thesegen <- genmat[,i] + 1
  thesegen[is.na(thesegen)] <- ploidy + 2
  StructTab[[dimnames(genmat)[[2]][i]]] <- unlist(G[thesegen])
  }
  
  # add marker name header
  if(exportMarkerNames){
  cat(paste(locNames(gl), collapse = "\t"), sep = "\n", file = outfile)
  }
  
  # export all data
  write.table(StructTab, row.names = FALSE, col.names = FALSE, append = TRUE,
 sep = "\t", file = outfile, quote = FALSE)
  if (v==1)  cat(paste("Structure file saved as:", outfile,"\nin folder:",outpath))
}

