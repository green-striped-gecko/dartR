#' @name gl.impute
#' @title Imputation of missing data
#' @description
#' This function imputes genotypes on a population-by-population basis, where
#' populations can be considered panmictic, or imputes the state for
#' presence-absence data.
#' @param x Name of the genlight object containing the SNP or presence-absence
#' data [required].
#' @param method Imputation method, either "frequency" or "HW" or "neighbour" 
#' or "random" [default "HW"].
#' @param parallel A logical indicating whether multiple cores -if available-
#' should be used for the computations (TRUE), or not (FALSE); requires the
#' package parallel to be installed [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log ; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @details
#' We recommend that imputation be performed on sampling locations, before
#' any aggregation. Imputation is achieved by replacing missing values using
#' either of two methods:
#' \itemize{
#' \item If "frequency", genotypes scored as missing at a locus in an individual
#'  are imputed using the average allele frequencies at that locus in the 
#'  population from which the individual was drawn.
#' \item If "HW", genotypes scored as missing at a locus in an individual are 
#' imputed by sampling at random assuming Hardy-Weinberg equilibrium. Applies 
#' only to genotype data.
#' \item If "neighbour", substitute the missing values for the focal individual
#'  with the values taken from the nearest neighbour. Repeat with next nearest
#'  and so on until all missing values are replaced.
#' \item if "random, missing data are substituted by random values (0, 1 or 2). 
#' }
#'
#'   The nearest neighbour is the one with the smallest Euclidean distance in 
#'   all the dataset.
#'
#'   The advantage of this approach is that it works regardless of how many
#'   individuals are in the population to which the focal individual belongs,
#'   and the displacement of the individual is haphazard as opposed to:
#'
#'   (a) Drawing the individual toward the population centroid (HW and Frequency).
#'
#'   (b) Drawing the individual toward the global centroid (glPCA).
#'
#' Note that loci that are missing for all individuals in a population are not 
#' imputed with method 'frequency' or 'HW'. Consider using the function 
#' \code{\link{gl.filter.allna}} with by.pop=TRUE to remove them first.
#'
#' @return A genlight object with the missing data imputed.
#' @export
#' @author Custodian: Luis Mijangos 
#' (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl <- gl.impute(gl.filter.allna(platypus.gl,verbose=0),method="neighbour")

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
  
  if (method == "frequency" | method == "HW") {
    pop_list_temp <- seppop(x)
    pop_list <- list()
    
    for (y in pop_list_temp) {
      loci_all_nas <- sum(glNA(y) > nInd(y))
      nas_number <- sum(glNA(y)) / 2
      number_imputations <- nas_number - (loci_all_nas * nInd(y))
      
      if (method == "frequency") {
        if (verbose >= 2 & loci_all_nas >= 1) {
          cat(
            warn(
              "  Warning: Population ",
              popNames(y),
              " has ",
              loci_all_nas,
              " loci with all missing values.\n"
            )
          )
          if (verbose >= 3) {
            cat(
              report(
                "  Method= 'frequency':",
                number_imputations,
                "values imputed.\n"
              )
            )
          }
        }
        
        q_allele <- glMean(y)
        pop_matrix <- as.matrix(y)
        loc_na <- which(is.na(pop_matrix), arr.ind = T)
pop_matrix[loc_na] <- unname(unlist(lapply(q_allele[loc_na[, 2]], function(x) {
            return(as.numeric(s_alleles(q_freq = x)))
          })))
        y@gen <- matrix2gen(pop_matrix, parallel = parallel)
        pop_list <- c(pop_list, y)
      }
      
      if (method == "HW") {
        if (verbose >= 2 & loci_all_nas >= 1) {
          cat(
            warn(
              "  Warning: Population ",
              popNames(y),
              " has ",
              loci_all_nas,
              " loci with all missing values.\n"
            )
          )
          if (verbose >= 3) {
            cat(report(
              "  Method= 'HW':",
              number_imputations,
              "values imputed.\n"
            ))
          }
        }
        
        q_allele <- glMean(y)
        pop_matrix <- as.matrix(y)
        loc_na <- which(is.na(pop_matrix), arr.ind = T)
        pop_matrix[loc_na] <-
          unname(unlist(lapply(q_allele[loc_na[, 2]], function(x) {
            return(sample_genotype(q_freq = x))
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
        x3 <- rbind(x3, pop)
      }
    }
    
    # if 1 population
    if (length(pop_list) == 1) {
      x3 <- pop_list[[1]]
    }
    
  }
  
  if (method == "neighbour") {
    pop_list_temp <- seppop(x)
    pop_list <- list()
    
    for (y in pop_list_temp) {
      loci_all_nas <- sum(glNA(y) > nInd(y))
      nas_number <- sum(glNA(y)) / 2
      number_imputations <- nas_number - (loci_all_nas * nInd(y))
    }
    
    if (verbose >= 2 & loci_all_nas >= 1) {
      cat(
        warn(
          "  Warning: Population ",
          popNames(y),
          " has ",
          loci_all_nas,
          " loci with all missing values.\n"
        )
      )
      if (verbose >= 3) {
        cat(report(
          "  Method= 'neighbour':",
          number_imputations,
          "values imputed.\n"
        ))
      }
    }
    
    x3 <- x
    
    eucl_dis <-
      gl.dist.ind(x,
                  method = "Euclidean",
                  verbose = 0,
                  plot.out = FALSE)
    pw_dis <- as.data.frame(as.table(as.matrix(eucl_dis)))
    
    x_matrix <- as.matrix(x)
    
    for (ind in 1:nInd(x)) {
      ind_imp <- x_matrix[ind, ]
      pw_dis_2 <- pw_dis[which(ind == pw_dis$Var1), ]
      pw_dis_3 <- pw_dis_2[order(pw_dis_2$Freq), ]
      pw_dis_4 <- pw_dis_3[-(pw_dis_3 == 0), ]
      
      while (sum(is.na(ind_imp)) > 0) {
        if (nrow(pw_dis_4) == 0) {
          cat(important(
            "  No more individuals left to impute individual",
            ind,
            "\n"
          ))
          break()
        }
        
        neig <- as.numeric(pw_dis_4[1, "Var2"])
        neig_matrix <- as.matrix(x[neig])
        
        loc_na <- unname(which(is.na(ind_imp)))
        
        ind_imp[loc_na] <- neig_matrix[loc_na]
        
        x_matrix[ind, ] <- ind_imp
        pw_dis_4 <- pw_dis_4[-1, ]
      }
      
    }
    
    x3@gen <- matrix2gen(x_matrix, parallel = parallel)
    
  }
  
  if (method == "random") {
    pop_list_temp <- seppop(x)
    pop_list <- list()
    
    for (y in pop_list_temp) {
      loci_all_nas <- sum(glNA(y) > nInd(y))
      nas_number <- sum(glNA(y)) / 2
      number_imputations <- nas_number - (loci_all_nas * nInd(y))
    }
    
    if (verbose >= 2 & loci_all_nas >= 1) {
      cat(
        warn(
          "  Warning: Population ",
          popNames(y),
          " has ",
          loci_all_nas,
          " loci with all missing values.\n"
        )
      )
      if (verbose >= 3) {
        cat(report(
          "  Method= 'random':",
          number_imputations,
          "values imputed.\n"
        ))
      }
    }

    x3 <- x
    
    x_matrix <- as.matrix(x)
    loc_na <- which(is.na(x_matrix), arr.ind = T)
    x_matrix[loc_na] <- sample(c(0:2),size=nrow(loc_na),replace = T)
    x3@gen <- matrix2gen(x_matrix, parallel = parallel)
    
  }
  
  x3$chromosome <- x@chromosome
  x3$position <- x$position
  x3$ploidy <- x$ploidy
  x3$strata <- x$strata
  x3$hierarchy <- x$hierarchy
  x3$other <- x$other
  
  x3 <- gl.compliance.check(x3, verbose = 0)
  
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
