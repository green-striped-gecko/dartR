
#' @importFrom gridExtra grid.table ttheme_default


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


