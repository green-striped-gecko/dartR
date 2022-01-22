#' @name gl.impute
#' @title Imputation of missing data
#' @description
#' This function imputes genotypes on a population-by-population basis, where 
#' populations can be considered panmictic.
#' @param x Name of the genlight object containing the SNP data [required].
#' @param method Imputation method, either "frequency" or "HW" [default "HW"].
#' @param parallel A logical indicating whether multiple cores -if available- 
#' should be used for the computations (TRUE), or not (FALSE); requires the 
#' package parallel to be installed [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, 
#' progress log ; 3, progress and results summary; 5, full report 
#' [default 2 or as specified using gl.set.verbosity].
#' @details 
#' We advice that imputation should be performed on sampling locations, before 
#' any aggregation. The imputation is achieved by replacing missing values using
#' either of two methods:
#' \itemize{
#' \item "frequency" where allele frequencies for the population to which 
#' the respective individual belongs are used as probabilities to sample 2 
#' alleles for each missing data in each individual. 
#' \item "HW" where the Hardy-Weinberg equation (p^2 + 2pq + q^2 = 1) within 
#' each population is used to sample genotypes for each missing data in each 
#' individual.
#' }
#' 
#' Note that loci that have all the values missing are not imputed. Consider 
#' using the function \code{\link{gl.filter.allna}}.
#' @return A genlight object with the missing data imputed.
#' @export
#' @author Custodian: Luis Mijangos (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl <- gl.impute(platypus.gl)
#'

# this is a a function to convert a matrix to a format suitable for the slot 
# @gen in a genlight object

gl.impute <-  function(x,
                       method = "HW",
                       parallel = FALSE,
                       verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jody",
                   verbosity = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # DO THE JOB
  
  #separating populations
  pop_list_temp <- seppop(x)
  pop_list <- list()
  
  for (y in pop_list_temp) {
    loci_all_nas <- sum(glNA(y) > nInd(y))
    nas_number <- sum(glNA(y)) / 2
    number_imputations <- nas_number - (loci_all_nas * nInd(y))
    if (method == "frequency") {
      if (verbose >= 2) {
        cat(
          important(
            "  Population",
            popNames(y),
            "has",
            loci_all_nas,
            "loci with all missing values.\n"
          )
        )
        cat(
          important(
            "  ",
            number_imputations,
            "genotypes were imputed by sampling alleles using their frequency in the population as probalility.\n"
          )
        )
      }
      q_allele <- glMean(y)
      pop_matrix <- as.matrix(y)
      loc_na <- which(is.na(pop_matrix), arr.ind = T)
      pop_matrix[loc_na] <-
        unname(unlist(lapply(q_allele[loc_na[, 2]], function(x) {
          sample_alleles(q = x)
        })))
      y@gen <- matrix2gen(pop_matrix, parallel = parallel)
      pop_list <- c(pop_list, y)
    }
    if (method == "HW") {
      if (verbose >= 2) {
        cat(
          important(
            "  Population",
            popNames(y),
            "has",
            loci_all_nas,
            "loci with all missing values.\n"
          )
        )
        cat(
          important(
            "  ",
            number_imputations,
            "genotypes were imputed by sampling genotypes using the Hardy-Weinberg equation.\n"
          )
        )
      }
      q_allele <- glMean(y)
      pop_matrix <- as.matrix(y)
      loc_na <- which(is.na(pop_matrix), arr.ind = T)
      pop_matrix[loc_na] <-
        unname(unlist(lapply(q_allele[loc_na[, 2]], function(x) {
          sample_genotype(q = x)
        })))
      y@gen <- matrix2gen(pop_matrix, parallel = parallel)
      pop_list <- c(pop_list, y)
    }
  }
  
  # if more than 1 population
  if (length(pop_list) > 1) {
    x3 <- NULL
    # merge back populations
    for (pop in pop_list) {
      x3 <- rbind (x3, pop)
    }
  }
  # if 1 population
  if (length(pop_list) == 1) {
    x3 <- pop_list[[1]]
  }
  
  x3$chromosome <- x@chromosome
  x3$position <- x$position
  x3$ploidy <- x$ploidy
  x3$strata <- x$strata
  x3$hierarchy <- x$hierarchy
  x3$other <- x$other
  
  x3 <- gl.compliance.check(x3)
  
  # ADD TO HISTORY
  x3@other$history <- x@other$history
  nh <- length(x3@other$history)
  x3@other$history[[nh + 1]] <- match.call()
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # RETURN
  return(x3)
}
