#' Prints dartR reports saved in tempdir
#' 
#'@return Prints a table with all reports saved in tempdir. Currently the style cannot be changed.
#'@export
#'@importFrom gridExtra grid.table ttheme_default
#'@author Bernd Gruber & Luis Mijangos (bugs? Post to \url{https://groups.google.com/d/forum/dartr})
#'@examples 
#' reports <- gl.list.reports()

gl.list.reports <- function(){
  
  files_tempdir <- list.files(tempdir())
  files_tempdir <- files_tempdir[which(str_match(files_tempdir,"dartR")=="dartR")]
  
  nh <- length(files_tempdir)
  if (nh==0){
    cat(warn("There are no dartR reports saved in the current tempdir.") )
  } else {
  
      dd <- data.frame(nr=1:nh, reports=files_tempdir)
    
    #max width
    dd$reports = sapply(lapply(dd$reports, strwrap, width=80), paste, collapse="\n")
    
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


