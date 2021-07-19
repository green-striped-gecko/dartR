#' @name gl.plot.structure
#'
#' @title Plot a STRUCTURE analysis using a genlight object
#'
#' @description 
#' This function takes a structure run object (output from
#'  \code{\link{gl.runstructure}}) and plots the typical structure bar
#'   plots on a spatial map, providing a barplot for each subpopulation. Therefore it requires coordinates from a genlight object. This kind of plots should support the interpretation of the spatial structure of a population, but in principle is not different from \code{\link{gl.plot.structure}}
#' @param sr structure run object from \code{\link{gl.run.structure}} [required].
#' @param k the number for k the q matrix should be based on. Needs to
#'  be within you simulated range of k's in your sr structure run object.
#'  @param x 
#' @return an sr object (structure.result list output). Each list entry is a single structure
#'run output (there are k.range * num.k.rep number of runs). For example the
#' summary output of the first run can be accessed via \code{sr[[1]]$summary} 
#' or the q-matrix of the third run via \code{sr[[3]]$q.mat}. To conveniently
#' summarise the outputs across runs (clumpp) you need to run
#' gl.plotstructure on the returned sr object. For evanno plots run gl.evanno on your sr object.
#'
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#' #@examples
#' @export
#' @seealso \code{\link{gl.run.structure}},  \link[strataG]{clumpp},
#' @references 
#' Pritchard, J.K., Stephens, M., Donnelly, P. (2000) Inference of population structure using multilocus genotype data. Genetics 155, 945-959.
#' 
#' Archer, F. I., Adams, P. E. and Schneiders, B. B. (2016) strataG: An R package for manipulating, summarizing and analysing population genetic data. Mol Ecol Resour. doi:10.1111/1755-0998.12559
#' 
#' Evanno, G., Regnaut, S., and J. Goudet. 2005. Detecting the number of clusters of individuals using the software STRUCTURE: a simulation study. Molecular Ecology 14:2611-2620.


gl.map.structure <- function(sr, k ) {}

# cex =1.5
# 
# sx <- abs(diff(range(centers[,"lon"])))/(10*10)*cex
# sy <- sx/max(ff)* cex
# 
# ii=1
# m1 <- m
# for (ii in 1:5) {
#   for ( i in 1:5) {
#     oo <- (i-1)*sx
#     
#     m1 <- m1   %>% addRectangles(cx[ii]+oo, cy[ii], cx[ii]+oo+sx, cy[ii]+ff[ii,i]*sy, opacity = 0, color =  rainbow(5)[i], fillOpacity = 0.5)
#     
#   }
# }
# 
# m1 %>% leaflet::addProviderTiles(provider) %>% addLegend(labels=colnames(ff), colors=rainbow(5),position ="topright" )
