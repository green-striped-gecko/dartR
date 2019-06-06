#' Export DArT genlight object \{adegenet\} to faststructure format (to run faststructure elsewhere)
#'
#' Recodes in the quite specific faststructure format (e.g first six columns need to be there, but are ignored...check faststructure documentation (if you find any :-( )))
#'
#' The script writes out the a file in faststructure format. 
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param outfile -- file name of the output file (including extension) [default gl.str]
#' @param outpath -- path where to save the output file [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath="." when calling this function to direct output files to your working directory.
#' @param verbose -- specify the level of verbosity: 0, silent, fatal errors only; 1, flag function begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @param probar switch to show/hide progress bar
#' @return NULL
#' @export
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @importFrom utils getTxtProgressBar setTxtProgressBar txtProgressBar 

gl2faststructure <- function(x, outfile="gl.str", outpath=tempdir(), probar=TRUE, verbose=2){

# TIDY UP FILE SPECS

  outfilespec <- file.path(outpath, outfile)
  funname <- match.call()[[1]]

# FLAG SCRIPT START

  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }

  if (verbose > 0) {
    cat("Starting",funname,"\n")
  }

# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }


  # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nInd(x))
      pop(x) <- as.factor(pop(x))
    }

  # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x, verbose=0)
    if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {cat("  Warning: genlight object contains monomorphic loci\n")}

# FUNCTION SPECIFIC ERROR CHECKING

# DO THE JOB

  x <- as.matrix(x)
  #add six dummy colums
  nc <- ncol(x)+6
  if (probar) pb <- txtProgressBar(min=0, 1, style=3, initial=NA)
  zz <- file(outfile, "w")
  for (i in 1:nrow(x))
  {
    dummy <- rbind(x[i,], x[i,])
    index <- colSums(dummy, na.rm=T)==2
    dummy[, index] <- c(0,2)
    dummy <- ifelse(is.na(dummy), -9, dummy)
    dummy <- ifelse(dummy==0,1, dummy)
    dummy <- cbind(i,i,i,i,i,i,dummy)
    write(t(dummy), file=outfile, sep="\t", ncolumns = nc, append=TRUE)  
    if (probar)   setTxtProgressBar(pb, i/nrow(x))  
  }
  close(zz)
  if (probar) close(pb)
  #gi <-read.structure("project_data.str", ask=F, n.ind=50, n.loc=100, col.pop=0, col.lab = 1, onerowperind = TRUE, sep="\t")
  if (verbose >= 2){cat(paste0("Saved faststructure file: ",getwd(),"/", outfile, "\n") )}
  if (verbose >= 3){cat(paste("Consists of", nrow(x), "individuals and ", ncol(x), "loci."))}

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

  return(NULL)
}
