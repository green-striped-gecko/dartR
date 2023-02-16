#' @name utils.assignment_2
#' @title Population assignment probabilities
#' @description
#' This function takes one individual and estimates
#' their probability of coming from individual populations
#' from multilocus genotype frequencies.
#
#' @param x Name of the genlight object containing the SNP data [required].
#' @param unknown Name of the individual to be assigned to a population [required].
# @param inbreeding_par The inbreeding parameter [default 0].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @details
#' This function is a re-implementation of the function multilocus_assignment
#'  from package gstudio.
#'  Description of the method used in this function can be found at:
#' https://dyerlab.github.io/applied_population_genetics/population-assignment.html
#' @return A \code{data.frame} consisting of assignment probabilities for each
#'  population.
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' require("dartR.data")
#' res <- utils.assignment_2(platypus.gl,unknown="T27")
#' @export

utils.assignment_2 <- function(x,
                             unknown,
                             # inbreeding_par = 0,
                             verbose = 2) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jody",
                   verbosity = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  if (unknown %in% indNames(x) == FALSE) {
    stop(error(
      paste("  Individual", unknown, "is not in the genlight object\n")
    ))
  }
  
  # DO THE JOB
  
  # filtering loci with all missing data by population
  # x <- gl.filter.allna(x, by.pop = TRUE, verbose = 0)
  unknown_pop <- gl.keep.ind(x, ind.list = unknown, verbose = 0)
  unknown_pop <- data.frame(gl2alleles(unknown_pop))
  
  x <- gl.drop.ind(x, ind.list = unknown, verbose = 0)
  
  pop_names <- popNames(x)
  
  pop_list <- seppop(x)
  gl_alleles <- do.call(rbind, strsplit(x$loc.all, "/"))
  
  frequencies <- lapply(pop_list, function(y) {
    freq_allele <- gl.alf(y)
    freqs_gl <-
      data.frame(
        Allele1 = gl_alleles[, 1],
        Allele2 = gl_alleles[, 2],
        count1 = freq_allele[, 1] * (nInd(y)*2),
        count2 = freq_allele[, 2] * (nInd(y)*2)
      )
    return(freqs_gl)
  })
  
  ret <- data.frame(Population = pop_names, Likelihood = 0)
  
  for (popx in 1:nPop(x)) {
    
    # alpha_ = k in Baudouin and Lebrun (2000)
    alpha_ <- 2
    # the total number of different allelic states at this locus over all 
    # reference populations
    k <- 2
    # the total number of genes to be assigned
    m <- nLoc(x) 
    # the total number of genes in the reference population sample
    n <- nLoc(x) 
    
    term1 <- lgamma(m+1)
    term2 <- lgamma(n+alpha_)
    

    if (verbose >= 2) {
      cat(
        report(
          "  Assigning individual",
          unknown,
          "against population",
          pop_names[popx],
          "\n"
        )
      )
    }
    
    popfreq <- frequencies[[popx]]
    
    loc <-
      as.data.frame(do.call(rbind, strsplit(unname(
        unlist(unknown_pop)
      ), ":")))
    colnames(loc) <- c("a1", "a2")
    
    df_assign <- cbind(loc, popfreq)
    
    df_assign$a1_count <- NA
    df_assign[which(df_assign$a1 == df_assign$a2 & df_assign$a1 == df_assign$Allele1),"a1_count"] <- 2
    df_assign[which(df_assign$a1 == df_assign$a2 & df_assign$a1 == df_assign$Allele2),"a1_count"] <- 0
    df_assign[which(df_assign$a1 != df_assign$a2),"a1_count"] <- 1
    
    df_assign$a2_count <- NA
    df_assign[which(df_assign$a1 == df_assign$a2 & df_assign$a1 == df_assign$Allele2),"a2_count"] <- 2
    df_assign[which(df_assign$a1 == df_assign$a2 & df_assign$a1 == df_assign$Allele1),"a2_count"] <- 0
    df_assign[which(df_assign$a1 != df_assign$a2),"a2_count"] <- 1
    
    term3_1 <- lgamma(( df_assign$a1_count +  df_assign$count1  + (alpha_/k) ))
    term3_2 <- lgamma(( df_assign$a2_count +  df_assign$count2  + (alpha_/k) ))
    
    term3 <- sum(term3_1,term3_2,na.rm = TRUE)
    
    term4_1 <- lgamma(df_assign$a1_count + 1)
    term4_2 <- lgamma(df_assign$a2_count + 1)
    
    term4 <- sum(term4_1,term4_2,na.rm = TRUE)
    
    term5_1 <- lgamma(df_assign$a1_count + (alpha_/k) )
    term5_2 <- lgamma(df_assign$a2_count + (alpha_/k) )
    
    term5 <- sum(term5_1,term5_2,na.rm = TRUE)
    
    term6 <- lgamma(m+n+alpha_)
    
    log_L <- term1 + term2 + term3 - term4 - term5 - term6
    
    # assign Likelihood
    ret[popx, "Likelihood"] <- -log_L
  }
  
  ret <- ret[order(ret$Likelihood,decreasing = TRUE),]
  ret$score <- ret$Likelihood / sum(ret$Likelihood)
  ret$score <- round(ret$score, 5)
  ret$Likelihood <- round(ret$Likelihood, 5)
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # RETURN
  
  return(invisible(ret))
  
}

gl2alleles <- function (gl) {
  x <- as.matrix(gl[, ])
  homs1 <-
    paste(substr(gl@loc.all, 1, 1), "/", substr(gl@loc.all, 1, 1), sep = "")
  hets <- gl@loc.all
  homs2 <-
    paste(substr(gl@loc.all, 3, 3), "/", substr(gl@loc.all, 3, 3), sep = "")
  xx <- matrix(NA, ncol = ncol(x), nrow = nrow(x))
  
  for (i in 1:nrow(x)) {
    for (ii in 1:ncol(x)) {
      inp <- x[i, ii]
      if (!is.na(inp)) {
        if (inp == 0)
          xx[i, ii] <- homs1[ii]
        else if (inp == 1)
          xx[i, ii] <- hets[ii]
        else if (inp == 2)
          xx[i, ii] <- homs2[ii]
      } else {
        xx[i, ii] <- NA
      }
    }
  }
  xx <- gsub("/", ":", xx)
  return(xx)
}
