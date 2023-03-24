#'@name gl.sample
#'
#'@title Samples individuals from populations
#'
#'@description This is a convenience function to prepare a bootstrap approach in dartR. For a bootstrap approach it is often desirable to sample a defined number of individuals for each of the populations in a genlight object and then calculate a certain quantity for that subset (redo a 1000 times)
#'
#'@param x genlight object containing SNP/silicodart genotypes
#'@param nsample the number of individuals that should be sampled
#'@param replace a switch to sample by replacement (default).
#'@param verbose set verbosity
#'@details This is convenience function to facilitate a bootstrap approach
#'@return returns a genlight object with nsample samples from each populations.
#'
#'@author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#'@examples 
#'\dontrun{
#' #bootstrap for 2 possums populations to check effect of sample size on fixed alleles
#' gl.set.verbosity(0)
#' pp <- possums.gl[1:60,]
#' nrep <- 1:10
#' nss <- seq(1,10,2)
#' res <- expand.grid(nrep=nrep, nss=nss)
#' for (i in 1:nrow(res)) {
#' dummy <- gl.sample(pp, nsample=res$nss[i], replace=TRUE)
#' pas <- gl.report.pa(dummy, plot.out = F)
#' res$fixed[i] <- pas$fixed[1]
#' }
#' boxplot(fixed ~ nss, data=res)
#'}
#'@family base dartR
#'@export 
#'
gl.sample <- function(x,
                  nsample = min(table(pop(x))),
                  replace = TRUE,
                  verbose = NULL) {
   #remove metadata to speed up  
  #if (!is.null(x@other$loc.metrics))  x@other$loc.metrics<- NULL
  #if (!is.null(x@other$ind.metrics))  x@other$ind.metrics<- NULL
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  # FLAG SCRIPT START
 funname <- match.call()[[1]]
  #utils.flag.start(func=funname,build="Jody",verbose=verbose)
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose=verbose)
  # FUNCTION SPECIFIC ERROR CHECKING
  
  # DO THE JOB
  #find samples
  ss <- sapply(1:nPop(x), function(z) which(pop(x)==levels(pop(x))[z]), simplify = F) 
  samps <- unlist(lapply(ss, function(x) sample(x, nsample, replace=replace))) 
  #subset x by samples
    ns <- ceiling(length(samps)/nInd(x))
    ff <- rep(1:ns,nInd(x))[1:length(samps)] 
    sp <- split(samps, ff)
    px <- lapply(sp, function(z) x[z,] )
    xx <- do.call(rbind, px)
  n10 <- nchar(as.character(nInd(xx)))
  lzs <- paste0("%0",as.character(n10),"d")
  indNames(xx)<- paste0(sprintf(lzs,1:nInd(xx)),"_",indNames(xx))
  return(xx)
}
    
