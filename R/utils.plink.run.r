#' @name utils.plink.run
#' @title Runs PLINK from within R
#'
#' @description Runs PLINK from within R. 
#' @param dir.in The path where the data files are
#' @param   pink.cmd The 'name' to call plink. This will depend on the file name
#'   (without the extension '.exe' if on windows) or the name of the PATH variable
#' @param   plink.path The path where the executable is. If plink is listed in
#'   the PATH then there is no need for this. This is what the option "path"
#'   means
#' @param fn The root of the data file name
#' @param out The root of the output file name
#' @param  syntax the flags to pass to plink call
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report [default NULL].
#' @return A character vector with the command used for PLINK.
#' @details
#' PLINK needs to be installed on the 
#'   machine and syntax used need to be appropriate for the version installed.
#' @references
#' Purcell, Shaun, et al. 'PLINK: a tool set for whole-genome association and
#' population-based linkage analyses.' The American journal of human genetics
#' 81.3 (2007): 559-575.
#' @export
#' @author Custodian: Carlo Pacioni and Luis Mijangos (Post to
#'  \url{https://groups.google.com/d/forum/dartr})

utils.plink.run <- function(dir.in, 
                      plink.cmd="plink", 
                      plink.path="path", 
                      out="hapmap1", 
                      syntax, 
                      verbose=NULL) {
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jody",
                   verbosity = verbose)
  
  # DO THE JOB
  old.wd <- getwd()
  on.exit(setwd(old.wd))
  # I set wd to where the user wants the output file and copy there the executable
  # because providing the path to PLINK doesn't work for me on Windows and I couldn't 
  # find a better solution
  setwd(dir.in) 
  if(!plink.path == "path") {
    is.win <- Sys.info()['sysname'] == "Windows" # Possibly is not needed to 
    # specifically search for the .exe on windows, it is just a safety feature
    exe <- list.files(plink.path, 
                      pattern = if(is.win) ".exe$" else paste0("^", plink.cmd), 
                      full.names = TRUE, include.dirs = TRUE)
    file.copy(from = exe, to = file.path(dir.in, basename(exe)))
  }
  cmd <- paste(plink.cmd, syntax, "--out", out, collapse=" ")
  system(cmd)
  file.remove(from = exe, to = file.path(dir.in, basename(exe)))
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }
  
  return(cmd)
}
