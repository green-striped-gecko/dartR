#' Converts a genind object to genlight object
#' 
#' @param gi -- a genind object
#' @param parallel -- switch to deactivate parallel version. Default set to FALSE. Might not be worth to run it parallel most of the times.
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A genlight object, with all slots filled.
#' @export
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @details Be aware due to ambiguity which one is the reference allele a combination of gi2gl(gl2gi(gl)) does not return an identical object (but in terms of analysis this conversions are equivalent)

gi2gl <- function(gi, parallel=FALSE, verbose=NULL){
  
# TRAP COMMAND, SET VERSION

  funname <- match.call()[[1]]
  build <- "Jacob"
  
# SET VERBOSITY
   # set verbosity
  if (is.null(verbose) & !is.null(gi@other$verbose)) verbose=gi@other$verbose
  if (is.null(verbose)) verbose=2
  if (verbose < 0 | verbose > 5){
      cat(paste("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to default [2] \n"))
      verbose <- 2
  }

  if (verbose >= 1){
    cat("Starting",funname,"\n")
  }
  
# STANDARD ERROR CHECKING
  
  if(class(gi)!="genind") {
    stop("  Fatal Error: genind object required!\n")
  }
  
# DO THE JOB
  
  locna <- gi@loc.n.all
  ccc<-1
  for (i in 2:length(locna)) 
  {
    if (locna[i-1]==1)  ccc[i] <- ccc[i-1]+1 else ccc[i]<- ccc[i-1]+2
  }
  gl <-new("genlight", gi@tab[,ccc], pop = pop(gi), other=gi@other, ploidy=2, loc.names=locNames(gi), ind.names=indNames(gi), parallel=parallel)
  if (is.null(gl@other$loc.metrics.flags$monomorphs)) gl@other$loc.metrics.flags$monomorphs <- FALSE
# FLAG SCRIPT END
  
  # ADD TO HISTORY
  
  nh <- length(gl@other$history)
  gl@other$history[[nh + 1]] <- match.call()   
  
  if (verbose >= 1) {
    cat("Completed:",funname,"\n")
  }
  
  return(gl)
}

