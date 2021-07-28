#' Prints dartR reports saved in tempdir
#' 
#'@return Prints a table with all reports saved in tempdir. Currently the style cannot be changed.
#'@export
#'@importFrom knitr kable

#'@author Bernd Gruber & Luis Mijangos (bugs? Post to \url{https://groups.google.com/d/forum/dartr})
#'@examples 
#'\dontrun{
#' gl.report.callrate(testset.gl)
#' gl.list.reports()
#'}
gl.list.reports <- function(){
  
  files_tempdir <- list.files(tempdir())
  files_tempdir <- files_tempdir[which(str_match(files_tempdir,"dartR")=="dartR")]
  
  nh <- length(files_tempdir)
  if (nh==0){
    cat(warn("There are no dartR reports saved in the current tempdir.") )
  } else {
  
      dd <- data.frame(nr=1:nh, reports=files_tempdir)
      dd$function_call <- NA
        for (i in 1:length(files_tempdir)) {
          dd$function_call[i] <- readRDS(paste0(tempdir(),"/", files_tempdir[i]))[[1]]
        }
    dd$function_call  <-substr(dd$function_call,2,999)
    
    dd$nr <- crayon::green(dd$nr)
    dd$reports <- crayon::blue(dd$reports)
    dd$function_call <- crayon::yellow(dd$function_call)
    
      
    #colnames(dd) <- c("Report number","Report","Function call") 
    colnames(dd)<- crayon::bold(colnames(dd))
    print(knitr::kable(dd, align = c("c","l","l"), format="simple",width=80))

    
    #  if (max(nchar(dd$function_call))>30) 
    #   dd$function_call  <-gsub('(.{1,30})', '\\1\n',dd$function_call)
     
    # #set table theme
    # tt <- ttheme_default()
    # tt$rowhead$fg_params$x=0
    # tt$core$fg_params$fontsize=11
    # tt$core$fg_params$hjust=0.0
    # tt$core$fg_params$x=c(rep(0.9, nh),rep(0.01,2*nh))
    # tt$core$fg_params$fontfamily="mono"
    # tt$core$fg_params$fontface="bold"
    # plot(0, type="n", xlab="", ylab="", axes=F)
    # grid.table(dd, theme=tt)
    
  }
}


