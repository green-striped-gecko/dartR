#' Prints dartR reports saved in tempdir
#' 
#'@return Prints a table with all reports saved in tempdir. Currently the style cannot be changed.
#'@export
#'@author Bernd Gruber & Luis Mijangos (bugs? Post to \url{https://groups.google.com/d/forum/dartr})
#'@examples 
#'\dontrun{
#' gl.report.callrate(testset.gl,save2tmp=TRUE)
#' gl.list.reports()
#'}
gl.list.reports <- function(){
  
  files_tempdir_temp <- list.files(tempdir())
  files_tempdir_plot <- files_tempdir_temp[which(str_match(files_tempdir_temp,"Plot")=="Plot")]
  files_tempdir_table <- files_tempdir_temp[which(str_match(files_tempdir_temp,"Table")=="Table")]
  files_tempdir_blast <- files_tempdir_temp[which(str_match(files_tempdir_temp,"Blast")=="Blast")]
  files_tempdir_files <- files_tempdir_temp[which(str_match(files_tempdir_temp,"File")=="File")]
  
  files_tempdir <- c(files_tempdir_plot,files_tempdir_table,files_tempdir_blast)
  
  nh <- length(files_tempdir)
  if (nh==0){
    cat(warn("There are no dartR reports saved in the current tempdir.") )
  } else {
  
      dd <- data.frame(nr=1:nh, reports=files_tempdir,time=file.info(paste0(tempdir(),"/",files_tempdir))$atime)
      dd$function_call <- NA
        for (i in 1:length(files_tempdir)) {
          dd$function_call[i] <- readRDS(paste0(tempdir(),"/", files_tempdir[i]))[[1]]
        }
    dd$function_call  <- substr(dd$function_call,2,999)
    
   
    dd$function_call <- crayon::yellow(dd$function_call)
    dd$time <- crayon::magenta(dd$time)
    tables <- which(str_match(dd$reports,"Table")=="Table")
    plots <- which(str_match(dd$reports,"Plot")=="Plot")
    dd[tables,"reports"] <- crayon::blue(dd[tables,"reports"])
    dd[plots,"reports"] <- crayon::cyan(dd[plots,"reports"])
    dd <- dd[order(dd$time),]
    dd$nr <- 1:nrow(dd)
    dd$nr <- crayon::green(dd$nr)
      
    #colnames(dd) <- c("Report number","Report","Function call") 
    colnames(dd)<- crayon::bold(colnames(dd))
    print(knitr::kable(dd, align = c("c","l","l"), format="simple",width=80,row.names=FALSE))

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


