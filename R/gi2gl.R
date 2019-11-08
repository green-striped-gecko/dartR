#' Converts a genind object to genlight object
#' 
#' @param gi -- a genind object
#' @param parallel -- switch to deactivate parallel version. Default set to TRUE. Only for testing purpose.
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A genlight object, with all slots filled.
#' @export
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @details Be aware due to ambiguity which one is the reference allele a combination of gi2gl(gl2gi(gl)) does not return an identical object (but in terms of analysis this conversions are equivalent)

gi2gl <- function(gi, parallel=TRUE, verbose=2){
  
# TRAP COMMAND, SET VERSION

  funname <- match.call()[[1]]
  build <- "Jacob"
  
# SET VERBOSITY
  
  if (is.null(verbose)){
    verbose <- x@other$verbose
  }
  if (verbose < 0 | verbose > 5){
      cat(paste("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to default",x@other$verbose,"\n"))
      verbose <- x@other$verbose
  }

# FLAG SCRIPT START

  if (verbose >= 1){
    cat("Starting",funname,"[ Build =",build,"]\n")
  }
  
# STANDARD ERROR CHECKING
  
  if(class(x)!="genind") {
    stop("  Fatal Error: genind object required!\n")
  }
  
# DO THE JOB
  
  locna <- x@loc.n.all
  ccc<-1
  for (i in 2:length(locna)) 
  {
    if (locna[i-1]==1)  ccc[i] <- ccc[i-1]+1 else ccc[i]<- ccc[i-1]+2
  }
  gl <-new("genlight", x@tab[,ccc], pop = pop(x), other=x@other, ploidy=2, loc.names=locNames(x), ind.names=indNames(x), parallel=parallel)
  
# FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat("Completed:",funname,"\n")
  }
  
  return(gl)
}

