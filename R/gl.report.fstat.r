#' @name gl.report.fstat
#' @title
#' Reports various statistics of genetic differentiation between
#' populations
#' @family matched reports
#' @description
#'
#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param nboots Number of bootstraps [default 0].
#' @param conf Confident intervals [default 0.95].
#' @param boot.type A vector of character strings representing the type of
#' intervals required. The value should be any subset of the values
#' c("norm","basic", "stud", "perc") [default "norm"].
#' @param parallel The type of parallel operation to be used (if any). If
#' missing, the default is taken from the option "boot.parallel" (and if that is
#'  not set, "no").
#' @param ncpus Number of processes to be used in parallel operation: typically
#' one would chose this to the number of available CPUs [default 1].
#' @param digits Number of Digits to use [default 4].
#' @param plot.stat Statistic to plot. One of "Ho","Hs","Ht","Dst","Htp","Dstp",
#' "Fst","Fstp","Fis","Dest","Gst_max","Gst_H" [default "Fst"].
#' @param plot.display If TRUE, histograms of base composition are displayed in
#' the plot window
#' [default TRUE].
#' @param palette.divergent A divergent palette for the distance values
#'  [default gl.colors("div")].
#' @param plot.dir Directory in which to save files [default = working directory]
#' @param plot.file Name for the RDS binary file to save (base name only,
#' exclude extension) [default NULL]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default NULL, unless specified using gl.set.verbosity]
#' @param ... Parameters passed to function \link[gplots]{heatmap.2},
#'  such as width and height, when the ggplot is to be saved.
#' @details
#'
#'  Even though Fst and its relatives can predict evolutionary processes
#'  (Holsinger & Weir, 2009), they are not true measures of genetic
#'   differentiation in the sense that they are dependent on the diversity
#'    within populations (Meirmans & Hedrick, 2011), the number of populations
#'    analysed (Alcala & Rosenberg, 2017) and are not monotonic
#'    (Sherwin et al., 2017). Recent approaches have been developed to
#'    accommodate these mathematical restrictions (G'ST and Jost's D; Hedrick,
#'    2005; Jost, 2008). More recently, novel approaches based on information
#'    theory (Mutual Information; Sherwin et al., 2017) and allele frequencies
#'    (Allele Frequency Difference; Berner, 2019) have distinct properties that
#'     make them valuable resources to interpret genetic differentiation.
#'
#'     Note that each measure of genetic differentiation has advantages and
#'     drawbacks, and the decision of using a particular measure is usually
#'     based on the research question.
#'      \itemize{
#'      \item
#'     "Ho" - Observed heterozygosity corrected for sample size
#'     (Nei, 1987, pp. 164–165) is calculated as:
#'
#'     \figure{Ho_equation.jpg}
#'
#'     where Pkii represents the proportion of homozygote i in sample k and np
#'     the number of samples.
#'
#'     \item
#'     "Hs" - Expected heterozygosity corrected for sample size
#'     (Nei, 1987, pp. 164–165) is calculated as:
#'
#'     \figure{He_equation.jpg}
#'
#'     \figure{He_equation_2.jpg}
#'
#'     \item
#'     "Ht" - Overall heterozygosity corrected for sample size
#'     (Nei, 1987, pp. 164–165) is calculated as:
#'
#'     \figure{Ht_equation.jpg}
#'
#'     \figure{Ht_equation_2.jpg}
#'
#'     \item
#'     "Dst" - Amount of heterozygosity among samples (Nei, 1987, pp. 164–165)
#'     is calculated as:
#'
#'     \figure{Dst_equation.jpg}
#'
#'     \item
#'     "Htp" - Overall heterozygosity corrected for sample size
#'     (Nei, 1987, pp. 164–165) is calculated as:
#'
#'     \figure{Htp_equation.jpg}
#'
#'     \item
#'     "Dstp" - Amount of heterozygosity among samples corrected for sample size
#'      (Nei, 1987, pp. 164–165) is calculated as:
#'
#'     \figure{Dstp_equation.jpg}
#'
#'     \item
#'     "Fst" - Nei's Gst (Nei, 1987, pp. 164–165) is calculated as:
#'
#'     \figure{Fst_equation.jpg}
#'
#'     \item
#'     "Fstp" - Fst corrected for sample size (Nei, 1987, pp. 164–165) is
#'     calculated as:
#'
#'     \figure{Fstp_equation.jpg}
#'
#'     \item
#'     "Fis" - Inbreeding coefficient is calculated as:
#'
#'     \figure{Fis_equation.jpg}
#'
#'     \item
#'     "Dest" - Jost’s D (Jost, 2008) is calculated as:
#'
#'     \figure{Dest_equation.jpg}
#'
#'     \item
#'     "Gst_max" - The maximum level that Gst can obtain for the observed amount
#'      of genetic variation (Hedrick 2005) is calculated as:
#'
#'     \figure{GstMax_equation.jpg}
#'
#'     where k is the number of subpopulations.
#'
#'     \item
#'     "Gst_H" - Gst standardized by the maximum level that it can obtain for the
#'     observed amount of genetic variation (Hedrick 2005) is calculated as:
#'
#'     \figure{Gst_H.jpg}
#'
#'     }
#'
#'  \strong{Bootstrapping}
#'
#'  Nonparametric Bootstrap Confidence Intervals (CI) are calculated using the
#'  function \link[boot]{boot.ci} (package rrBLUP). The type of intervals is
#'  controlled by the parameter "boot.type".
#'
#' Efron (1979) suggested that Bootstraps should be at least 200.
#'
#'  \strong{Plotting}
#'
#'  The plot can be customised by including any parameter from the function
#'  \link[gplots]{heatmap.2}
#'
#'  A color vector can be obtained with  \code{\link{gl.select.colors}} and then
#'   passed to the function with the plot.colors parameter.
#'
#' If a plot.file is given, the plot arising from this function is saved as an
#'  "RDS" binary file using saveRDS(); can be reloaded with readRDS(). A file
#'   name must be specified for the plot to be saved.
#'
#'  If a plot directory (plot.dir) is specified, the ggplot binary is saved to that
#'  directory; otherwise to the tempdir().
#'
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' res <- gl.report.fstat(platypus.gl)
#'
#' @references
#' \itemize{
#' \item
#' Alcala, N., & Rosenberg, N. A. (2017). Mathematical constraints on FST:
#' Biallelic markers in arbitrarily many populations. Genetics (206), 1581-1600.
#' \item
#' Berner, D. (2019). Allele frequency difference AFD–an intuitive alternative
#' to FST for quantifying genetic population differentiation. Genes, 10(4), 308.
#' \item
#' Efron, B. (1979). Bootstrap methods: Another look at the jackknife. Annals of
#' Statistics 7, 1–26.
#' \item
#' Hedrick, P. W. (2005). A standardized genetic differentiation measure.
#' Evolution, 59(8), 1633-1638.
#' \item
#' Holsinger, K. E., & Weir, B. S. (2009). Genetics in geographically structured
#'  populations: defining, estimating and interpreting FST. Nature Reviews
#'  Genetics, 10(9), 639- 650.
#'  \item
#'  Jost, L. (2008). GST and its relatives do not measure differentiation.
#'  Molecular Ecology, 17(18), 4015-4026.
#'  \item
#'  Meirmans, P. G., & Hedrick, P. W. (2011). Assessing population structure:
#'  FST and related measures. Molecular Ecology Resources, 11(1), 5-18.
#'  \item
#'  Nei, M. (1987). Molecular evolutionary genetics: Columbia University Press.
#'  \item
#'  Sherwin, W. B., Chao, A., Jost, L., & Smouse, P. E. (2017). Information
#'  theory broadens the spectrum of molecular ecology and evolution. Trends in
#'   Ecology & Evolution, 32(12), 948-963.
#' }
#' @export
#' @return A list containing matrices with genetic statistics taken pairwise by
#' population
#'
# ----------------------
# Function

gl.report.fstat <- function(x,
                            nboots = 0,
                            conf = 0.95,
                            boot.type = "norm",
                            parallel = "no",
                            ncpus = 1,
                            digits = 4,
                            plot.stat = "Fst",
                            plot.display = TRUE,
                            palette.divergent = gl.colors("div"),
                            plot.file = NULL,
                            plot.dir = NULL,
                            verbose = NULL,
                            ...) {
  # PRELIMINARIES -- checking ----------------
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.2",
                   verbosity = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # DO THE JOB
  
  # keeping populations with more than 1 individuals
  pop_names <- popNames(x)[which(table(pop(x)) > 1)]
  
  x <- gl.keep.pop(x, pop.list = pop_names)
  class(x) <- "dartR"
  
  # bootstrapping function
  pop.diff <- function(x, indices) {
    x2 <- x[, indices]
    pop_diff <- utils.basic.stats(x2, digits = digits)
    return(pop_diff$overall)
  }
  
  pops <- seppop(x)
  npops <- length(pops)
  pairs_pops <- t(combn(npops, 2))
  pairs_pops_names <- apply(pairs_pops, 1, function(y) {
    paste0(names(pops)[y[1]], "_vs_", names(pops)[y[2]])
  })
  
  ### pairwise
  if (npops > 2) {
    
    pairpop_res <- apply(pairs_pops, 1, function(y) {
      tpop <- rbind.dartR(pops[[y[1]]], pops[[y[2]]])
      res_tmp <- utils.basic.stats(tpop, digits = digits)$overall
      return(res_tmp)
    })
    
    if (nboots > 0) {
      # bootstrapping
      pairpop_boot <- apply(pairs_pops, 1, function(y) {
        tpop <- rbind.dartR(pops[[y[1]]], pops[[y[2]]])
        res_boots <- boot::boot(
          data = tpop,
          statistic = pop.diff,
          R = nboots,
          parallel = parallel,
          ncpus = ncpus
        )
        return(res_boots)
      })
      
      # confidence intervals
      res_CI <- replicate(length(pairpop_boot),
                          as.data.frame(matrix(nrow = 12, ncol = 2)),
                          simplify = FALSE)
      
      for (boot_n in 1:length(pairpop_boot)) {
        for (stat_n in 1:12) {
          t0_tmp <- pairpop_boot[[boot_n]]$t0[stat_n]
          t_tmp <- pairpop_boot[[boot_n]]$t[, stat_n]
          res_1 <- boot:::basic.ci(t0 = t0_tmp,
                                   t = t_tmp,
                                   conf = conf)[4:5]
          res_CI[[boot_n]][stat_n, ] <- (res_1)
        }
      }
    }
    
  } else{
    
    tpop <- rbind.dartR(pops[[1]], pops[[2]])
    
    pairpop_res <- utils.basic.stats(tpop, digits = digits)$overall
    
    if (nboots > 0) {
      # bootstrapping
      pairpop_boot <- boot::boot(
        data = tpop,
        statistic = pop.diff,
        R = nboots,
        parallel = parallel,
        ncpus = ncpus
      )
      
      # confidence intervals
      res_CI <- as.data.frame(matrix(nrow = 12, ncol = 2))
      
        for (stat_n in 1:12) {
          t0_tmp <- pairpop_boot$t0[stat_n]
          t_tmp <- pairpop_boot$t[, stat_n]
          res_1 <- boot:::basic.ci(t0 = t0_tmp,
                                   t = t_tmp,
                                   conf = conf)[4:5]
          res_CI[stat_n, ] <- (res_1)
        }
      
    }
    
  }
    
    if (npops > 2) {
      stat_pop <- asplit(pairpop_res, 2)
      stat_pop <- lapply(stat_pop,function(x){
        as.data.frame(x)
      })
      CI <- Map(cbind, stat_pop, res_CI)
      CI <- lapply(CI, function(y) {
        colnames(y) <-  c("Value", "HCI", "LCI")
        return(y)
      })
      names(CI) <- pairs_pops_names

    } else{
      
      stat_pop <- pairpop_res
      CI <- cbind(stat_pop, res_CI)
      colnames(CI) <-  c("Value", "HCI", "LCI")

    }

      mat_pops <-
        rep(list(matrix(
          NA, nrow = npops, ncol = npops
        )), 12)
      
      if (npops > 2) {
    for (i in 1:length(mat_pops)) {
      mat_pops[[i]][lower.tri(mat_pops[[i]])] <- pairpop_res[i, ]
      colnames(mat_pops[[i]]) <-
        rownames(mat_pops[[i]]) <- names(pops)
      mat_pops[[i]][upper.tri(mat_pops[[i]])] <-
        t(mat_pops[[i]])[rev(lower.tri(mat_pops[[i]]))]
    }
      }else{
        
        for (i in 1:length(mat_pops)) {
          mat_pops[[i]][lower.tri(mat_pops[[i]])] <- pairpop_res[i]
          colnames(mat_pops[[i]]) <-
            rownames(mat_pops[[i]]) <- names(pops)
          mat_pops[[i]][upper.tri(mat_pops[[i]])] <-
            t(mat_pops[[i]])[rev(lower.tri(mat_pops[[i]]))]
        }
        
      }
    
     if (npops > 2) {
      names(mat_pops) <- rownames(pairpop_res)
     } else{
       names(mat_pops) <- names(pairpop_res)
     }
    
  # }
  
  # Printing outputs -----------
      print(list(Stat_matrices=mat_pops, Confident_Intervals=CI))

  # solution to print gplots https://stackoverflow.com/a/19191951
  create_heatmap <- function(...) {
    plot_heatmap <- function()
      gl.plot.heatmap(...)
  }
  
  p3 <- create_heatmap(mat_pops[[plot.stat]],
                       palette.divergent = palette.divergent)
  
  # PLOT THE RESULTS -----------------
  if (plot.display & npops > 2) {
    p3()
  }
  
  if (plot.display & npops <= 2) {
    report(cat("No plot is displayed if only two populations are analysed\n"))
  }
  
  # Optionally save the plot ---------------------
  
  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(p3,
                           dir = plot.dir,
                           file = plot.file,
                           verbose = verbose)
  }
  
  # FINISH UP -------------------
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  # ----------------------
  
  # RETURN
  
  if (nboots > 0) {
    return(list(Stat_matrices=mat_pops, Confident_Intervals=CI))
  } else{
    return(mat_pops)
  }
  
}
