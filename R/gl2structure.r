#' Converts genlight objects to STRUCTURE formated files
#'
#' This function exports genlight objects to STRUCTURE formatted files (be aware there is a gl2faststructure version as well). It is based on the code provided by Lindsay Clark (see \url{https://github.com/lvclark/R_genetics_conv})  and this function is basically a wrapper around her numeric2structure function. See also: Lindsay Clark. (2017, August 22). lvclark/R_genetics_conv: R_genetics_conv 1.1 (Version v1.1). Zenodo: doi.org/10.5281/zenodo.846816.
#' 
#' @param x -- name of the genlight object containing the SNP data and location data, lat longs [required]
#' @param indNames -- specify individuals names to be added [if NULL, defaults to indNames(x)]
#' @param addcolumns -- additional columns to be added before genotypes [default NULL]
#' @param ploidy -- set the ploidy [defaults 2]
#' @param exportMarkerNames -- if TRUE, locus names locNames(x) will be included [default TRUE]
#' @param outfile -- file name of the output file (including extension) [default gl.str]
#' @param outpath -- path where to save the output file [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath="." when calling this function to direct output files to your working directory.
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @export
#' @author Bernd Gruber (wrapper) and Lindsay V. Clark [lvclark@illinois.edu] 
#' @examples
#' \donttest{
#' gl2structure(testset.gl)
#'}


gl2structure <- function(x, indNames=NULL, addcolumns=NULL, ploidy=2, exportMarkerNames=TRUE, outfile="gl.str", outpath=tempdir(), verbose=NULL){

# TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  outfilespec <- file.path(outpath, outfile)
  
# SET VERBOSITY
  
  if (is.null(verbose)){ 
    if(!is.null(x@other$verbose)){ 
      verbose <- x@other$verbose
    } else { 
      verbose <- 2
    }
  } 
  
  if (verbose < 0 | verbose > 5){
    cat(paste("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n"))
    verbose <- 2
  }
  
# FLAG SCRIPT START
  
  if (verbose >= 1){
    if(verbose==5){
      cat("Starting",funname,"[ Build =",build,"]\n")
    } else {
      cat("Starting",funname,"\n")
    }
  }

# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    stop("  Fatal Error: genlight object required!\n")
  }
  
  if (verbose >= 2){
    if (all(x@ploidy == 1)){
      stop("Fatal Error: Detected Presence/Absence (SilicoDArT) data. Please provide a SNP dataset\n")
    } else if (all(x@ploidy == 2)){
      cat("  Processing a SNP dataset\n")
    } else {
      stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
    }
  }

# FUNCTION SPECIFIC ERROR CHECKING

  nInd <- nInd(x)
  if (is.null(indNames)) {indNames=indNames(x)}
  if(length(indNames) != nInd){
    stop("Fatal Error: No. of individuals listed in user-specified indNames and no. of individuals in supplied genlight object x do not match\n")
  }

  if (!is.null(addcolumns) && is.null(dim(addcolumns))) addcolumns <- data.frame(pop=addcolumns)
  if(!is.null(addcolumns) && nrow(addcolumns) != nInd){
    stop("Fatal Error: No. of individuals in user-specified addColumns and no. of individuals in supplied genlight object x does not match\n")
  }

  genmat <- as.matrix(x)
  if(!all(genmat %in% c(0:ploidy,NA))){
    stop("Fatal Error: genmat must only contain 0, 1, 2... ploidy and NA\n")
  }

  if(length(outfile) != 1 || !is.character(outfile)){
    stop("Fatal Error: output file must be a single character string\n")
  }

  if(length(ploidy) != 1 || !is.numeric(ploidy)){
    stop("Fatal Error: ploidy must be a single number\n")
  }

  if(!exportMarkerNames %in% c(TRUE, FALSE)){
    stop("Fatal Error: exportMarkerNames must be TRUE or FALSE\n")
  }

# DO THE JOB

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
    cat(paste(locNames(x), collapse = "\t"), sep = "\n", file = outfile)
  }
  
# export all data
  write.table(StructTab, row.names = FALSE, col.names = FALSE, append = TRUE,
  sep = "\t", file = outfilespec, quote = FALSE)
  if (verbose >=2 )  cat(paste("Structure file saved as:", outfile,"\nin folder:",outpath))

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

  return(NULL)
}

