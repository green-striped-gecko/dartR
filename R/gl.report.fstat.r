#' @name gl.report.fstat
# PRELIMINARIES -- Set parameters --------------
#' @title Reports various statistics of genetic differentiation between populations 
#' @family matched reports

#' @description

#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param nboots Number of bootstraps [default 0].
#' @param conf Confident intervals [default 0.95].
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
#' @param plot.theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot.colors List of two color names for the borders and fill of the
#'  plots [default gl.select.colors(library="brewer",palette="Blues",select=c(7,5))].
#' @param plot.dir Directory in which to save files [default = working directory]
#' @param plot.file Name for the RDS binary file to save (base name only, exclude extension) [default NULL]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default NULL, unless specified using gl.set.verbosity]
#' @param ... Parameters passed to function \link[gplots]{heatmap.2},
#'  such as width and height, when the ggplot is to be saved.
#' @details
#' 
#'  Even though FST and its relatives can predict evolutionary processes 
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
#'  \figure{Ho_equation.jpg}   
#' 
#'  \strong{Bootstrapping}
#'  
#' Efron (1979) suggested that Bootstraps should be at least 200.
#' 
#'  A color vector can be obtained with gl.select.colors() and then passed to the function
#'  with the plot.colors parameter.

#' Themes can be obtained from in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }

#' If a plot.file is given, the ggplot arising from this function is saved as an "RDS"
#' binary file using saveRDS(); can be reloaded with readRDS(). A file name must be
#' specified for the plot to be saved.

#'  If a plot directory (plot.dir) is specified, the ggplot binary is saved to that
#'  directory; otherwise to the tempdir().

#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

#' @examples
#' 
#' @references
#' \itemize{
#' \item
#' Holsinger, K. E., & Weir, B. S. (2009). Genetics in geographically structured
#'  populations: defining, estimating and interpreting F ST. Nature Reviews 
#'  Genetics, 10(9), 639- 650.
#'  
#' Efron, B. (1979). Bootstrap methods: Another look at the jackknife.
#' Annals of Statistics 7, 1–26.
#' }
#' @export
#' @return 
# ----------------------
# Function

gl.report.fstat <- function(x,
                        nboots = 0,
                        conf = 0.95,
                        parallel = "no",
                        ncpus = 1,
                        digits = 4,
                        plot.stat = "Fst", 
                        plot.display = TRUE,
                        plot.theme = theme_dartR(),
                        plot.colors = gl.select.colors(
                          library = "brewer",
                          palette = "Blues",
                          select = c(7, 5),
                          verbose = 0
                        ),
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
  
  x <- gl.keep.pop(x,pop.list = pop_names)
  
  # bootstrapping function
  pop.diff <- function(x, indices) {
    x2 <- x[, indices]
    pop_diff <- utils.basic.stats(x2, digits = digits)
    
    return(pop_diff$overall)

  }
  
  pops <- seppop(x)
  npops <- length(pops)
  pairs_pops <- t(combn(npops, 2))
  pairs_pops_names <- apply(pairs_pops,1,function(y){
    paste0(names(pops)[y[1]],"_vs_",names(pops)[y[2]])
  })
  
  ### pairwise
  if (npops > 2) {
    pairpop <- apply(pairs_pops, 1, function(y) {
      tpop <- rbind(pops[[y[1]]], pops[[y[2]]])
      
      if (nboots > 0) {
        res <- boot::boot(data = tpop,
                    statistic = pop.diff,
                    R = nboots)
      } else{
        res <- utils.basic.stats(tpop, digits = digits)$overall
      }
      
      return(res)
      
    })
    
  } else{
    tpop <- rbind(pops[[1]], pops[[2]])
    
    if (nboots > 0) {
      pairpop <- boot::boot(data = tpop,
                      statistic = pop.diff,
                      R = nboots)
      
    } else{
      pairpop <- utils.basic.stats(tpop, digits = digits)$overall
    }
  }
  
  if (nboots > 0) {
    mat_pops <- rep(list(matrix(NA, nrow = npops, ncol = npops)), 12)
    
    # if more than 2 pops
    if (npops > 2) {
      pairpop_2 <- Reduce(cbind, lapply(pairpop, "[[", 1))
    } else{
      pairpop_2 <- as.data.frame(pairpop[[1]])
    }
    
    for (i in 1:length(mat_pops)) {
      mat_pops[[i]][lower.tri(mat_pops[[i]])] <- pairpop_2[i, ]
      colnames(mat_pops[[i]]) <- rownames(mat_pops[[i]]) <- names(pops)
    }
    
    names(mat_pops) <- rownames(pairpop_2)
    
    # calculation of confident intervals are based on the function Bootstrap.CI
    # from Anne Chao's package spadeR
    if (npops > 2) {
      pairpop_boot <- lapply(pairpop, "[[", "t")
      pairpop.mean <- lapply(pairpop_boot, apply, 2, mean)
      LCI <- as.list(1:length(pairpop.mean))
      UCI <- as.list(1:length(pairpop.mean))
    } else{
      pairpop_boot <- pairpop$t
      pairpop.mean <- colMeans(pairpop_boot)
      LCI <- NULL
      UCI <- NULL
    }
    
    for (i in 1:length(pairpop.mean)) {
      if (npops > 2) {
        LCI_tmp <-
          -apply(pairpop_boot[[i]], 2, function(x)
            quantile(x, probs = (1 - conf) / 2)) + pairpop.mean[[i]]
        UCI_tmp <-
          apply(pairpop_boot[[i]], 2, function(x)
            quantile(x, probs = 1 - (1 - conf) / 2)) - pairpop.mean[[i]]
        LCI[[i]] <- pairpop_2[, i] - LCI_tmp
        UCI[[i]] <- pairpop_2[, i] + UCI_tmp
      } else{
        # LCI_tmp <- -apply(pairpop_boot,2,function(x)quantile(x,probs = (1-conf)/2)) + pairpop.mean[i]
        # UCI_tmp <-  apply(pairpop_boot,2,function(x)quantile(x,probs = 1-(1-conf)/2)) - pairpop.mean[i]
        # LCI <- pairpop_2 - LCI_tmp
        # UCI <- pairpop_2 + UCI_tmp
        LCI_tmp <-
          apply(pairpop_boot, 2, function(x)
            quantile(x, probs = (1 - conf) / 2))
        UCI_tmp <-
          apply(pairpop_boot, 2, function(x)
            quantile(x, probs = 1 - (1 - conf) / 2))
        LCI <- LCI_tmp
        UCI <- UCI_tmp
      }
    }
    
    if (npops > 2) {
      stat_pop <- lapply(pairpop, "[[", "t0")
      CI <- Map(cbind, stat_pop, LCI, UCI)
      CI <- lapply(CI,function(y){
        colnames(y) <-  c("Mean","HCI","LCI")
        return(y)
      })
      names(CI) <- pairs_pops_names
      
    } else{
      stat_pop <- pairpop$t0
      CI <- Reduce(cbind, list(stat_pop, LCI, UCI))
      colnames(CI) <-  c("Mean","HCI","LCI")
      names(CI) <- pairs_pops_names
    }
    
  } else{
    if (npops > 2) {
      mat_pops <-
        rep(list(matrix(
          NA, nrow = npops, ncol = npops
        )), nrow(pairpop))
    } else{
      mat_pops <-
        rep(list(matrix(
          NA, nrow = npops, ncol = npops
        )), length(pairpop))
    }
    
    for (i in 1:length(mat_pops)) {
      mat_pops[[i]][lower.tri(mat_pops[[i]])] <- pairpop[i,]
      colnames(mat_pops[[i]]) <- rownames(mat_pops[[i]]) <- names(pops)
      mat_pops[[i]][upper.tri(mat_pops[[i]])] <-
        t(mat_pops[[i]])[rev(lower.tri(mat_pops[[i]]))]
      
    }
    
    if (npops > 2) {
      names(mat_pops) <- rownames(pairpop)
    } else{
      names(mat_pops) <- names(pairpop)
    }
    
  }

  # Printing outputs -----------

  # solution to print gplots https://stackoverflow.com/a/19191951
  create_heatmap <- function(...) {
    plot_heatmap <- function() gl.plot.heatmap(...)
  }
  
  p3 <- create_heatmap(mat_pops[[plot.stat]],
                       cexRow = 0.75,
                       cexCol = 0.75)

  # PLOT THE RESULTS -----------------
  
  if (plot.display) {
    p3()
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
    return(list(mat_pops, CI))
  } else{
    return(mat_pops)
  }

  }
