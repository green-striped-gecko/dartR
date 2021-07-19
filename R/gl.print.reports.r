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
  
  files_tempdir <- list.files(tempdir())
  files_tempdir <- files_tempdir[which(str_match(files_tempdir,"dartR")=="dartR")]
  
  nh <- length(files_tempdir)
  if (print_report > nh){
    cat(warn("The report requested is not in the current tempdir.") )
  } else {
    report_to_print <- readRDS(paste0(tempdir(),"/", files_tempdir[print_report]))
    suppressWarnings(print(report_to_print[[2]]))
  }
}
