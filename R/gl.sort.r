#'@name gl.sort
#'
#'@title re-sorts genlight objects
#'
#'@description Often it is desirable to have the genlight object sorted individuals by population names, indiviual name, for example to have a more informative gl.smearplot (showing banding patterns for populations). Also sorting by loci can be informative in some instances. This function provides the ability to sort individuals of a genlight object by providing the order of individuals or populations and also by loci metric providing the order of locis. See examples below for specifics.
#'
#'@param x genlight object containing SNP/silicodart genotypes
#'@param sort.by either "ind", "pop". Default is pop
#'@param order.by that is used to order individuals or loci. Depening on the order.by parameter, this needs to be a vector of length of nPop(genlight) for populations or  nInd(genlight) for individuals. If not specified alphabetical order of populations or individuals is used. For sort.by="ind" order.by can be also a vector specifying the order for each individual (for example another ind.metrics)
#'@param verbose set verbosity
#'@details This is convenience function to facilitate sorting of individuals within the genlight object. For example if you want to visualise the "band" of population in a gl.smearplot then the order of individuals is important. Also
#'@return Returns a reordered genlight object. Sorts also the ind/loc.metrics and coordinates accordingly
#'
#'@author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#'@examples 
#'#sort by populations
#'bc <- gl.sort(bandicoot.gl)
#'#sort from West to East
#'bc2 <- gl.sort(bandicoot.gl, sort.by="pop" ,
#'order.by=c("WA", "SA", "VIC", "NSW", "QLD"))
#'#sort by missing values
#'miss <- rowSums(is.na(as.matrix(bandicoot.gl)))
#'bc3 <- gl.sort(bandicoot.gl, sort.by="ind", order.by=miss)
#'gl.smearplot(bc3)
#'@family base dartR
#'@export 
#'
gl.sort <- function(x,
                  sort.by = "pop",
                  order.by = NULL,
                  verbose = NULL) {
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func=funname,build="Jody",v=verbose)
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose=verbose)
  # FUNCTION SPECIFIC ERROR CHECKING
  if (is.na(pmatch(sort.by,c("pop","ind"))))
    stop(error("sort.by is not either pop, ind or loci. Please specify one of those options."))
  if (!is.null(order.by)) {
    if (sort.by=="pop") len <- nPop(x)
    if (sort.by=="ind") len <- nInd(x)
    #if (sort.by=="other ") len <- nLoc(x)
    if (len!=length(order.by)) stop(error("Length of order.by does not match length of",sort.by, ". Check your input parameters and make sure length of your sort.by selection matches the length of the vector provided by order.by"))
  }
  # DO THE JOB
  
  if (sort.by=="pop") {
    if (is.null(order.by)) index <- order(pop(x)) else {
      if (!(sum(pop(x) %in% order.by))==nInd(x)) stop(error("order.by does not contain all levels of sort.by."))
      index  <- order(factor(pop(x), levels=order.by))
    } 
    xx <- x[index,]
    xx@other$latlon <- x@other$latlon[index,]
    xx@other$ind.metrics <- x@other$ind.metrics[index,]
    if (!is.null(order.by)) pop(xx)<- factor(pop(xx), levels=order.by)
    }
  if (sort.by=="ind") {
    if (is.null(order.by)) index <- order(indNames(x)) else {
      if (!(length(order.by)==nInd(x))) stop(error("order.by does not contain all levels of sort.by."))
      index  <- order(order.by)
    } 
      xx <- x[index,]
    xx@other$latlon <- x@other$latlon[index,]
    xx@other$ind.metrics <- x@other$ind.metrics[index,]
  }
  # ADD TO HISTORY
  if (is(xx,"genlight")) {
    nh <- length(xx@other$history)

    xx@other$history[[nh + 1]] <- c(match.call())
  } 
  return(xx)
}
    