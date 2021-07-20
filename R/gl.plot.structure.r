#' @name gl.plot.structure
#'
#' @title Plot a STRUCTURE analysis using a genlight object
#'
#' @description 
#' This function takes a structure run object (output from
#'  \code{\link{gl.run.structure}}) and plots the typical strcture bar
#'   plot that visiualise the q matrix of a structure run. 
#'
#' @param sr structure run object from \code{\link{gl.run.structure}} [required].
#' @param k the number for k the q matrix should be based on. Needs to
#'  be within you simulated range of k's in your sr structure run object.
#' @param sort how q matrix is sorted (by population, group proportion etc.)
#'  [default is by the population definition and group proportion from group
#'   1 to k-1, meaning first sorted by group, then by proportion of group 1,
#'    then group 2,.... then group k-1].
#' @param CLUMPP path to the clumpp executable file. Windows: CLUMPP.exe,
#'  macos and linux: CLUMPP (no exe).
#' @param ... additional parameter passed to the clumpp function within
#'  package strataG (\link[strataG]{clumpp}).
#' @param plot_theme Theme for the plot. See details for options [default
#'  theme_dartR()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log ; 3, progress and results summary; 5, full report [default
#'   NULL, unless specified using gl.set.verbosity]
#'
#' @details The function outputs a barplot which is the typical output of
#'  structure. This function needs the use of clumpp, which is an external
#'   program that needs to be installed For a evanno plot use gl.evanno. 
#' 
#' @return a barplot in ggplot format
#'
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' 
#' @examples
#' \dontrun{
#' #CLUMPP needs to be installed to be able to run the example
#' #only the first 100 loci
#' #bc <- bandicoot.gl[,1:100]
#' #sr <- gl.run.structure(bc, k.range = 2:5, num.k.rep = 3, exec = "./structure.exe")
#' #qmat <- gl.plot.structure(sr, k=3, CLUMPP="d:/structure/")
#' #gl.map.structure(qmat, bc, scalex=1, scaley=0.5)
#' }
#' @export
#' @seealso \code{\link{gl.run.structure}},  \link[strataG]{clumpp}, \code{\link{gl.plot.structure}}
#' @references 
#' Pritchard, J.K., Stephens, M., Donnelly, P. (2000) Inference of population structure using multilocus genotype data. Genetics 155, 945-959.
#' Archer, F. I., Adams, P. E. and Schneiders, B. B. (2016) strataG: An R package for manipulating, summarizing and analysing population genetic data. Mol Ecol Resour. doi:10.1111/1755-0998.12559
#' 
#' Mattias Jakobsson and Noah A. Rosenberg. 2007. CLUMPP: a cluster matching and permutation program for dealing with label switching and multimodality in analysis of population structure. Bioinformatics 23(14):1801-1806. Available at \href{http://web.stanford.edu/group/rosenberglab/clumppDownload.html}{clumpp}
#' 

###@importFrom strataG genind2gtypes structureRun

gl.plot.structure <- function(sr, k, sort=NULL, CLUMPP="./" ,... , plot_theme,verbose ){
  
#IS strataG INSTALLED?
  
  if (!requireNamespace("strataG", quietly = TRUE)) #not already installed?
  {
    ap <- available.packages()  #check CRAN
    
    oncran <- match("strataG", ap)
    if (is.na(oncran)) {
      warning("package strataG needs to be installed. It is currently not on CRAN, hence I try to install it manually via Github using devtools:\n  devtools::install_github('EricArcher/strataG'")
      devtools::install_github('EricArcher/strataG')
      
    }
  } else {
# DO THE JOB
#run clump
  
if (!is(sr,"structure.result")) stop("sr is not a structure result object returned from gl.run.structure.")

#change range of simulated ks in structure object
ks <- range((lapply(sr, function(x) x$summary[1])))
  
if (is.null(k) | k<ks[1] | k > ks[2]) stop("No k provided of k is not in range of the simulated k's in your structure run object")
if (Sys.info()['sysname']=="Windows") clumppfile <- "CLUMPP.exe" else clumppfile <- "CLUMPP"
if (!file.exists(file.path(CLUMPP, clumppfile))) stop("Cannot find clumpp executable. Please provide full path.") 
owd <- getwd()
setwd(CLUMPP)
q.mat <- strataG::clumpp(sr, k = k)
setwd(owd)

qq <- q.mat[,4:(k+3)]

if (is.null(sort)) {
ll <- data.frame(cbind(as.numeric(factor(q.mat$orig.pop)), qq[,]))
zz <- do.call(order,  unname(as.list(ll)))

bb <- t(qq[zz,])
colnames(bb)<- q.mat$id
narg <- paste(q.mat$id,q.mat$orig.pop[zz], sep="_")

} else {
  bb <- t(qq[sort,])
  colnames(bb)<- q.mat$id[sort]
  narg <- q.mat$id[sort]
  }

#bgg <- reshape2::melt(bb)
#bgg <- transform(bgg, Var1=factor(Var1, rownames(bb)), Var2=factor(Var2, colnames(bb)))

bbpp <- barplot(bb, col = 1:k,las = 2, main = paste0("K=",k), border=1:k, space = 0, names.arg=narg)
#ggplot(bgg, aes(x=Var2, y=value, fill=Var1), )+geom_bar(stat="identity", width = 1)
  }
return(q.mat)
}
