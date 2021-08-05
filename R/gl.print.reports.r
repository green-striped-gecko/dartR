#' Prints dartR reports saved in tempdir
#' @param print_report number of report from \code{\link{gl.list.reports}} that is to be printed
#'@return Prints reports that were saved in tempdir. 
#'@export
#'@author Bernd Gruber & Luis Mijangos (bugs? Post to \url{https://groups.google.com/d/forum/dartr})
#'@examples 
#' \dontrun{
#' reports <- gl.print.reports(1)
#' }

gl.print.reports <- function(print_report){
  
  files_tempdir_temp <- list.files(tempdir())
  files_tempdir_plot <- files_tempdir_temp[which(str_match(files_tempdir_temp,"Plot")=="Plot")]
  files_tempdir_table <- files_tempdir_temp[which(str_match(files_tempdir_temp,"Table")=="Table")]
  files_tempdir_blast <- files_tempdir_temp[which(str_match(files_tempdir_temp,"Blast")=="Blast")]
  files_tempdir_files <- files_tempdir_temp[which(str_match(files_tempdir_temp,"File")=="File")]
  
  files_tempdir <- c(files_tempdir_plot,files_tempdir_table,files_tempdir_blast)
  
  nh <- length(files_tempdir)
  if (print_report > nh){
    cat(warn("The report requested is not in the current tempdir."))
  } else {
    
    dd <- data.frame(nr=1:nh, reports=files_tempdir,time=file.info(paste0(tempdir(),"/",files_tempdir))$atime)
    dd <- dd[order(dd$time),]
    dd$nr <- 1:nrow(dd)

    report_to_print <- readRDS(paste0(tempdir(),"/", dd$reports[print_report]))
    suppressWarnings(print(report_to_print[[2]]))
  }
}
