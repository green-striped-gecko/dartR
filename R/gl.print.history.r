#' Prints history of a genlight object
#' 
#' @param x -- a genlight object (with history [optional])
#' @param history -- [optional] either a link to a history slot (gl\@other$history), or a vector indicating which part of the history of x is used [c(1,3,4) uses the first, third and forth entry from x\@other$history]. If no history is provided the complete history of x is used (recreating the identical object x).
#' @return prints a table with all history records. Currently the style cannot be changed.
#' @export
#' @importFrom gridExtra grid.table ttheme_default
#' @author Bernd Gruber (bugs? Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples 
#' dartfile <- system.file("extdata","testset_SNPs_2Row.csv", package="dartR")
#' metadata <- system.file("extdata","testset_metadata.csv", package="dartR")
#' gl <- gl.read.dart(dartfile, ind.metafile = metadata, probar=FALSE) 
#' gl2 <- gl.filter.callrate(gl, method="loci", threshold=0.9)
#' gl3 <- gl.filter.callrate(gl2, method="ind", threshold=0.95)
#' #Now "replay" part of the history "onto" another genlight object
#' bc.fil <- gl.play.history(gl.compliance.check(bandicoot.gl), 
#' history=gl3@other$history[c(2,3)], verbose=1)
#' gl.print.history(bc.fil)




gl.print.history<- function(x = NULL, history = NULL)
{
  if (class(x) == "genlight")
    if (is.null(history))
      hist2 <- x@other$history
    else
      hist2 <- x@other$history[history]
    
    if (is.null(x) & is.list(history))
      hist2 <- history
    
    
    nh <- length(hist2)
    if (nh==0) warning("You did not specify a history correctly. Check your genlight object.") else {
      
      #for (i in 1:length(hist2)) {
      #  hist2[[i]]$x <- "gl"
      #}
      
      dd <- data.frame(nr=1:nh, history=as.character(hist2))
      
      #max width
      dd$history = sapply(lapply(dd$history, strwrap, width=80), paste, collapse="\n")
      
      dd[nh+1,] <- c("->",as.character(substitute(x)) )
      #set table theme
      tt <- ttheme_default()
      tt$rowhead$fg_params$x=0
      tt$core$fg_params$fontsize=11
      tt$core$fg_params$hjust=0
      tt$core$fg_params$x=c(rep(0.5, nh),0.2, rep(0.01, nh+1))
      tt$core$fg_params$fontfamily="mono"
      tt$core$fg_params$fontface="bold"
      plot(0, type="n", xlab="", ylab="", axes=F)
      grid.table(dd, theme=tt)
    } 
    
}


