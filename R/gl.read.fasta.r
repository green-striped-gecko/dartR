#' @name gl.read.fasta
#' @title Reads FASTA files and converts them to genlight object
#' @description
#' The following IUPAC Ambiguity Codes are taken as heterozygotes:
#' \itemize{
#'  \item M is heterozygote for	AC and CA	
#'  \item R is heterozygote	for AG and GA	
#'  \item W is heterozygote	for AT and TA	
#'  \item S is heterozygote	for CG and GC
#'  \item Y is heterozygote	for CT and TC
#'  \item K is heterozygote	for GT and TG
#'  }
#'  
#' The following IUPAC Ambiguity Codes are taken as missing data:
#' \itemize{
#'  \item V
#'  \item H 
#'  \item D 
#'  \item B
#'  \item N 
#'  }
#'  
#'  The function can deal with missing data in individuals, e.g. when FASTA 
#'  files have different number of individuals due to missing data.
#'  
#'  The allele with the highest frequency is taken as the reference allele.
#'  
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
#' \url{https://evodify.com/heterozygotes-ambiguity-characters/}
#' @return A genlight object.
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#'  # Folder where the fasta files are located. 
#'  folder_samples <- system.file('extdata', package ='dartR')
#'  # listing the FASTA files, including their path. Files have an extension
#'  # that contains "fas".
#'  file_names <- list.files(path = folder_samples, pattern = "*.fas", 
#'                           full.names = TRUE)
#'  # reading fasta files
#'   obj <- gl.read.fasta(file_names)
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
                    n_cores = n_cores,
                    verbose = verbose)
  
  x <- merge_gl_fasta(gl_list, parallel = parallel,verbose = verbose)
  
  x <- gl.compliance.check(x, verbose = verbose)
 
   # add history
  x@other$history <- list(match.call())
  x <- gl.recalc.metrics(x, verbose = 0)
  
  if (verbose >= 2) {
    cat(
      important(
        "  Genlight object does not have individual metrics. You need to add them 'manually' to the @other$ind.metrics slot.\n"
      )
    )
  }
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # RETURN
  return(invisible(x))
  
}