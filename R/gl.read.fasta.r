#' @name gl.read.fasta
#' @title Reads FASTA files and converts them to genlight object
#' @description
#' The following IUPAC Ambiguity Codes are taken as heterozygotes:
#'  W, S, M, K, R, Y. 
#'  
#'  Take into account that the function will take a little bit longer the first
#'   time you run it because C++ functions must be compiled.
#' @param fasta_files Fasta files to read [required].
#' @param parallel A logical indicating whether multiple cores -if available-
#'  should be used for the computations (TRUE), or not (FALSE); requires the
#'   package parallel to be installed [default FALSE].
#' @param n_cores If parallel is TRUE, the number of cores to be used in the
#'  computations; if NULL, then the maximum number of cores available on the
#'   computer is used [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @details
#' Ambiguity characters are often used to code heterozygotes. However, using
#'  heterozygotes as ambiguity characters may bias many estimates. See more
#'   information in the link below:
#' https://evodify.com/heterozygotes-ambiguity-characters/
#' @return A genlight object.
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \dontrun{
#' folder_samples <- getwd()
#' file_names <- list.files(path = folder_samples, pattern = "*.fas",
#' full.names = TRUE)
#' gl_list <- gl.read.fasta(file_names)
#' }
#' @family reading functions
#' @export


gl.read.fasta <- function(fasta_files,
                          parallel = FALSE,
                          n_cores = NULL,
                          verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jody",
                   verbosity = verbose)
  
  # DO THE JOB
  
  gl_list <- lapply(fasta_files,
                    utils.read.fasta,
                    parallel = parallel,
                    n_cores = n_cores)
  
  fin_res <- merge_gl_fasta(gl_list, parallel = parallel)
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # RETURN
  return(invisible(fin_res))
  
}