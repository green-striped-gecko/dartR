#' Export DArT genlight object \{adegenet\} to faststructure format (to run faststructure elsewhere)
#'
#' Recodes in the quite specific faststructure format (e.g first six columns need to be there, but are ignored...check faststructure documentation (if you find any :-( )))
#'
#' The script writes out the a file in faststructure format. 
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param outfile -- file name of the output file (including extension) [default gl.str]
#' @param outpath -- path where to save the output file [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath="." when calling this function to direct output files to your working directory.
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @param probar switch to show/hide progress bar
#' @return NULL
#' @export
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @importFrom utils getTxtProgressBar setTxtProgressBar txtProgressBar 

gl2faststructure <- function(x, outfile="gl.str", outpath=tempdir(), probar=FALSE, verbose=NULL){

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
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
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

# DO THE JOB

  x <- as.matrix(x)
  #add six dummy colums
  nc <- ncol(x)+6
  if (probar) pb <- txtProgressBar(min=0, 1, style=3, initial=NA)
  zz <- file(outfilespec, "w")
  for (i in 1:nrow(x))
  {
    dummy <- rbind(x[i,], x[i,])
    index <- colSums(dummy, na.rm=T)==2
    dummy[, index] <- c(0,2)
    dummy <- ifelse(is.na(dummy), -9, dummy)
    dummy <- ifelse(dummy==0,1, dummy)
    dummy <- cbind(i,i,i,i,i,i,dummy)
    write(t(dummy), file=outfilespec, sep="\t", ncolumns = nc, append=TRUE)  
    if (probar)   setTxtProgressBar(pb, i/nrow(x))  
  }
  close(zz)
  if (probar) close(pb)
  #gi <-read.structure("project_data.str", ask=F, n.ind=50, n.loc=100, col.pop=0, col.lab = 1, onerowperind = TRUE, sep="\t")
  if (verbose >= 2){cat(paste0("Saved faststructure file: ", outfilespec, "\n") )}
  if (verbose >= 3){cat(paste("Consists of", nrow(x), "individuals and ", ncol(x), "loci."))}

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

  return(NULL)
}
