#' Isolation by distance
#' 
#' This functions performs an isolation by distance analysis based on a mantel test and also produces an isolation by distance plot. If a genlight object with coordinates is provided) then a Euclidean and genetic distance matrix are calculated (currently. Currently only pairwise Fst between population is implemented. Coordinates are expected as lat long and converted to Google Earth Mercator projection. If coordinates are already projected, set projected=TRUE. If such an object is provided an isolation by distance analysis and plot is performed on log(Euclidean distance) against population based pairwise Fst/1-Fst (see  Rousseau's distance measure. Genetics April 1, 1997 vol. 145 no. 4 1219-1228)
#' You can provide also your own genetic and Euclidean distance matrix. The function is based on the code provided by the adegenet tutorial (\url{http://adegenet.r-forge.r-project.org/files/tutorial-basics.pdf}), using the functions  \link[vegan]{mantel} (package vegan), \link[StAMPP]{stamppFst} (package StAMPP) and Mercator in package dismo.
#' 
#' @importFrom vegan mantel
#' @importFrom MASS kde2d
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics abline title points
#' @importFrom stats as.dist lm
#' @importFrom StAMPP stamppFst
#' @param x genlight object. If provided a standard analysis on Fst/1-Fst and log(distance) is performed
#' @param distance type of distance that is calculated and used for the analysis. Can be either population based "Fst" [\link[StAMPP]{stamppFst}], "D" [\link[StAMPP]{stamppNeisD}] or individual based "propShared", [gl.propShared], "euclidean" [gl.dist.ind, method="Euclidean"].
#' @param coordinates Can be either "latlon", "xy" or a two column data.frame with column names "lat","lon", "x", "y")  Coordinates are provided via \code{gl@other$latlon} ['latlon'] or via \code{gl@other$xy} ['xy']. If latlon data will be projected to meters using Mercator system [google maps] or if xy then distance is directly calculated on the coordinates.
#' @param logdist TRUE/FALSE switch if log of distance should be used [default is set to TRUE].
#' @param logoffset 0. If you have individuals/populations with zero distances between each other log transfromation results in zero distances, hence you need to an offset to the distances.
#' @param Dgen genetic distance matrix if no genlight object is provided
#' @param Dgeo Euclidean distance matrix if no genlight object is provided
#' @param permutations number of permutations in the mantel test
#' @param plot should an isolation by distance plot be returned. Default is plot=TRUE
#' @param paircols should pairwise dots colored by population/individual pairs
#' @param plot_theme Theme for the plot. See details for options [default theme_dartR()].
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]

#' @return returns a list of the following components: Dgen (the genetic distance matrix), Dgeo (the Euclidean distance matrix), mantel (the statistics of the mantel test)
#' @export
#' @author Bernd Gruber (bugs? Post to \url{https://groups.google.com/d/forum/dartr})
#' @seealso \link[vegan]{mantel}, \link[StAMPP]{stamppFst}
#' @references 
#' Rousset (1997) Genetic Differentiation and Estimation of Gene Flow from F-Statistics Under Isolation by Distancenetics 145(4), 1219-1228.
#' @examples 
#' \donttest{
#' ibd <- gl.ibd(bandicoot.gl)
#' ibd <- gl.ibd(bandicoot.gl,plot = FALSE)
#' }


gl.ibd <-  function(x=NULL, 
           distance="Fst",
           coordinates="latlon",
           logdist=TRUE,
           logoffset = 0,
           Dgen=NULL,
           Dgeo=NULL,
           permutations=999, 
           plot=TRUE,
           paircols=TRUE,
           plot_theme = theme_dartR(),
           verbose=options()$dartR_verbose) {
    
    # CHECK IF PACKAGES ARE INSTALLED
    if (!(requireNamespace("dismo", quietly = TRUE))) {
      stop("Package dismo needed for this function to work. Please install it.") } else { 
    
    # TRAP COMMAND
    funname <- match.call()[[1]]
    
    # GENERAL ERROR CHECKING
    verbose <- gl.check.verbosity(verbose)
    x <- utils.check.gl(x, verbose=0)
    
    
    #### SETTING DATA TYPE ####
    if (all(x@ploidy == 1)){
      if (verbose>0) cat(report("  Processing Presence/Absence (SilicoDArT) data\n"))
      datatype <- "SilicoDArT"
    } else if (all(x@ploidy == 2)){
      if (verbose>0) cat(report("  Processing a SNP dataset\n"))
      datatype <- "SNP"
    } else {
      stop (error("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)"))
    }

 #specific error checks       

    
    

 if (!is.null(Dgen) & !is.null(Dgeo)) 
  {
   if (verbose>0)  cat(report("Analysis performed using provided genetic and Euclidean distance matrices. If a genlight object is provided, it is ignored."))
   ta="dgendgeo"
}

if (is(x,"genlight")) 
{
  if (verbose>0) cat(report("SAnalysis performed on the genlight object."))
  
 ta="genlight" 
  
}
  
#check coordinates
coords <- NULL
if (is(coordinates,"character")) {
if (coordinates=="latlon")  {
  if (is.null(x@other$latlon))  stop(error("Cannot find coordinates in x@other$latlon")) 
  coords <- dismo::Mercator(x@other$latlon[,c("lon","lat")]) 
  if (verbose>0) cat(report("Coordinates transformed to Mercator (google) projection to calculate distances in meters.\n"))
}

if (coordinates=="xy") {
  if (is.null(x@other$xy))  stop(error("Cannot find coordinates in x@other$xy"))
  coords <- x@other$xy
}
    
}
    
if (is(coordinates,"data.frame")) {
      if (length(setdiff(colnames(coordinates),c("lat","lon")))==0)
        coords <- dismo::Mercator(coordinates[,c("lon","lat")])   
      
      if (length(setdiff(colnames(coordinates),c("x","y")))==0) 
        coords <- coordinates[,c("x","y")]
      if (is.null(coords)) stop(error("No valid coordinates provided. check the provided data.frame and its format."))
    }

if (is.null(coords)) stop(error("No valid coordinates provided!"))
    
#make sure coordinates have the correct length
if (nrow(coords)!=nInd(x) &  ta=="genlight") stop(error("Cannot find coordinates for each individual in slot @other$latlon"))

typedis=NULL
if (distance=="Fst" | distance=="D") typedis="pop"
if (distance=="propShared" | distance=="euclidean") typedis="ind"

#population based distances

if (is.null(Dgeo) & typedis=="pop")
{
  if (nPop(x)>1) {
  pop.xy <- apply(coords, 2, function(a) tapply(a, pop(x), mean, na.rm=T) )
  Dgeo <- dist(pop.xy)
  } else {stop(error("Less than 2 populations provided, therefore no pairwise distances can be calculated"))}
}

if (is.null(Dgeo) & typedis=="ind")
{
  if (nInd(x)>1) {
    
    Dgeo <- dist(coords)
  } else {stop(error("Less than 2 individuals provided, therefore no pairwise distances can be calculated"))}
}

#apply logarithm to distance

if (logdist) {
   if (sum(Dgeo==0,na.rm=T)>0 & logoffset==0) stop(error("Cannot log transform distances due to zero euclidean distances between pairs. Set logoffset to different from zero!")) else  Dgeo <- log(Dgeo+logoffset) 
}

if (is.null(Dgen) & distance=="Fst")
  Dgen <- as.dist(StAMPP::stamppFst(x, nboots=1))
if (is.null(Dgen) & distance=="D")
  Dgen <- as.dist(StAMPP::stamppNeisD(x,pop=TRUE))
if (is.null(Dgen) & distance=="propShared")
  Dgen <- as.dist(1-gl.propShared(x))
if (is.null(Dgen) & distance=="euclidean")
  Dgen <- as.dist(gl.dist.ind(x, method = "Euclidean",verbose=0, plot=FALSE))


### order both matrices to be alphabetically as levels in genlight (ind or pop)
if (class(x)=="genlight") {
if (typedis=="pop") ordering <- levels(pop(x)) else ordering <- order((indNames(x)))
Dgen <- as.dist(as.matrix(Dgen)[ordering, ordering])
Dgeo <- as.dist(as.matrix(Dgeo)[ordering, ordering])
}

#make sure both matrices are distance objects
Dgen <- as.dist(Dgen)
Dgeo <- as.dist(Dgeo)


if (is.null(Dgeo)) stop(error("Cannot calculate distance matrix or no distance matrix provided!"))
if (is.null(Dgen)) stop(error("Cannot calculate genetic distance matrix or no genetic distance matrix provided!"))


manteltest <- vegan::mantel(Dgen, Dgeo, na.rm=TRUE, permutations = permutations)
print(manteltest)

if (plot) 
  {
    plot(Dgeo, Dgen)
    abline(lm(as.numeric(Dgen)~as.numeric(Dgeo)))
    title("Isolation by distance")
  }
  
####### Printing outputs, using package patchwork
res <- data.frame(Dgen=as.numeric(Dgen), Dgeo=as.numeric(Dgeo))
if (!paircols) {
  p1 <- ggplot(res, aes(x=Dgeo, y=Dgen))+geom_point()+plot_theme}    else {
  
  cols <- which(lower.tri(as.matrix(Dgen)),arr.ind = T)
  c1 <- cols[,2]
  c2 <- cols[,1]
  res <- data.frame(Dgen=as.numeric(Dgen), Dgeo=as.numeric(Dgeo), col1=colnames(as.matrix(Dgen))[c1], col2=colnames(as.matrix(Dgen))[c2])
  p1 <- ggplot(res)+geom_point(aes(Dgeo, Dgen, col=col1), size=4)+geom_point(aes(Dgeo, Dgen, col=col2), size=1.5, shape=15)+plot_theme
  
  }

if (plot) print(p1)
if (verbose>0) print(manteltest)

# creating temp file names
temp_plot <- tempfile(pattern =paste0("dartR_plot",paste0(names(match.call()),"_",as.character(match.call()),collapse = "_")))
temp_table <- tempfile(pattern = paste0("dartR_table",paste0(names(match.call()),"_",as.character(match.call()),collapse = "_"),"_"))

# saving to tempdir
saveRDS(p1, file = temp_plot)
if(verbose>=2){cat(report("  Saving the plot in ggplot format to the tempfile as",temp_plot,"using saveRDS\n"))}
saveRDS(manteltest, file = temp_table)
if(verbose>=2){cat(report("  Saving the result of the manteltest to the tempfile as",temp_table,"using saveRDS\n"))}
if(verbose>=2){cat(report("  NOTE: Retrieve output files from tempdir using gl.list.reports() and gl.print.reports()\n"))}

# FLAG SCRIPT END

if (verbose >= 1) {
  cat(report("\nCompleted:", funname, "\n\n"))
}

out <- list(Dgen=Dgen, Dgeo=Dgeo, mantel=manteltest)
return(out)
}
}

