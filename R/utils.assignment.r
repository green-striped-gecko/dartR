#' @name utils.assignment
#' @title Population assignment probabilities
#' @description
#' This function takes one individual and estimates
#' their probability of coming from individual populations
#' from multilocus genotype frequencies.
#
#' @param x Name of the genlight object containing the SNP data [required].
#' @param ind Name of the individual to be assigned to a population [required].
#' @param inbreeding_par The inbreeding parameter [default 0]. 
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
#' res <- utils.assignment(platypus.gl,ind="T27")
#' @export

utils.assignment <- function(x,
                            ind,
                            inbreeding_par = 0,
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
  
  # DO THE JOB
  
  # filtering loci with all missing data by population 
  x <- gl.filter.allna(x, by.pop = TRUE, verbose = 0)
  
  pop_list <- seppop(x)
  gl_alleles <- do.call(rbind, strsplit(x$loc.all, "/"))
  loc_names <- locNames(x)
  
  frequencies <- lapply(pop_list, function(y) {
    freq_allele <- gl.alf(y)
    pop_name <- popNames(y)
    freqs_gl_1 <-
      data.frame(
        Population = pop_name,
        Locus = loc_names,
        Allele = gl_alleles[, 1],
        Frequency = freq_allele[, 1]
      )
    freqs_gl_2 <-
      data.frame(
        Population = pop_name,
        Locus = loc_names,
        Allele = gl_alleles[, 2],
        Frequency = freq_allele[, 2]
      )
    freqs_gl <- rbind(freqs_gl_1, freqs_gl_2)
    return(freqs_gl)
  })
  
  frequencies <- rbindlist(frequencies)
  
  individual <- gl.keep.ind(x, ind.list = ind, verbose = 0)
  individual <- data.frame(gl2alleles(individual))
  colnames(individual) <- loc_names
  
  loci <- loc_names
  
  pops <- sort(unique(frequencies$Population))
  ret <- data.frame(Population = pops, Probability = 0)
  
  for (pop in pops) {
    prob <- 1
    
    if (verbose >= 2) {
      cat(report(
        "  Assigning individual",
        ind,
        "against population",
        pop,
        "\n"
      ))
    }
    
    for (locus in loci) {
      popfreq <- frequencies[frequencies$Population == pop &
                               frequencies$Locus == locus ,]
      
      loc <- individual[[locus]]
      
      if (!is.na(loc)) {
        all_alleles <- strsplit(loc, ":")[[1]]
        
        f <-
          prod(unlist(lapply(all_alleles, function(z)
            return(
              popfreq$Frequency[popfreq$Allele == z]
            ))))
        
        if (all_alleles[1] != all_alleles[2]) {
          f <- f * 2
        }
        
        if (inbreeding_par > 0) {
          f <- f * (1 - inbreeding_par)
          if (all_alleles[1] == all_alleles[2]) {
            f <-
              f + (popfreq$Frequency[popfreq$Allele == all_alleles[1]] * inbreeding_par)
          }
        }
        
        prob <- prob * f
      } else{
        prob <- prob
        
      }
      
    }
    
    # assign probability
    ret$Probability[ret$Population == pop] <- prob
  }
  
  ret <- ret[order(-ret$Probability),]
  
  ret$Posterior <- ret$Probability / sum(ret$Probability)
  
  ret$Posterior <- round(ret$Posterior, 5)
  ret$Probability <- round(ret$Probability, 5)
  
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
      }
      else
        xx[i, ii] = NA
    }
  }
  xx <- gsub("/", ":", xx)
  return(xx)
}
