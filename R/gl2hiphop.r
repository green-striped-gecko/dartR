#' Converts a genlight objects into hiphop format
#'
#' This function exports genlight objects to the format used by the parentage
#' assignment R package hiphop. Hiphop can be used for paternity and maternity
#' assignment and outperforms conventional methods where closely related
#' individuals occur in the pool of possible parents. The method compares the
#' genotypes of offspring with any combination of potentials parents and scores
#' the number of mismatches of these individuals at bi-allelic genetic markers
#'  (e.g. Single Nucleotide Polymorphisms).
#'
#' @references
#' Cockburn, A., Penalba, J.V.,Jaccoud, D.,Kilian, A., Brouwer, L.,Double, M.C.,
#'  Margraf, N., Osmond, H.L., van de Pol, M. and Kruuk, L.E.B.(in revision).
#'  HIPHOP: improved paternity assignment among close relatives using a simple
#'  exclusion method for bi-allelic markers. Molecular Ecology Resources, DOI to
#'  be added upon acceptance
#' @param x Name of the genlight object containing the SNP data [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @return Dataframe containing all the genotyped individuals (offspring and
#'  potential parents) and their genotypes scored using bi-allelic markers.
#' @export
#' @importFrom dplyr %>% mutate_all
#' @author Luis Mijangos (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' \donttest{
#' result <- gl2hiphop(testset.gl)
#' }

gl2hiphop <- function(x,
                      verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # DO THE JOB
    
    x <- as.matrix(x[,])
    x[x == 1] <- "het"
    x[x == 2] <- 1
    x[x == "het"] <- 2
    x <- as.data.frame(x)
    x <- x %>%
        dplyr::mutate_all(as.numeric)
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(x)
}
