#' @name gl.spatial.autoCorr
#' @title Spatial autocorrelation following Smouse and Peakall 1999
#'
#' @description  Global spatial autocorrelation is a multivariate approach
#' combining all loci into a single analysis. The autocorrelation coefficient
#'  "r" is calculated for each pair of individuals in each specified distance
#'  class. For more information see Smouse and Peakall 1999, Peakall et al. 2003
#'   and Smouse et al. 2008.
#'
#' @details This function executes a modified version
#'  of \code{spautocorr} from the package \code{PopGenReport}. Differently 
#' from \code{PopGenReport}, this function also computes the 95\% confidence 
#' intervals around the r via bootstraps, the 95% confidence interval around the 
#' null hypothesis of no spatial structure and the one-tail test via permutation, 
#' and the correction factor described by Peakall et al 2003.
#' 
#' The input can be i) a genlight object (which has to have the latlon slot 
#' populated), ii) a pair of \code{Dgeo} and \code{Dgen}, which have to be either
#'  \code{matrix} or \code{dist} objects, or iii) a \code{list} of the 
#'  \code{matrix} or \code{dist} objects if the 
#'  analysis needs to be carried out for multiple populations (in this case, 
#'  all the elements of the \code{list} have to be of the same class (i.e. 
#'  \code{matrix} or \code{dist}) and the population order in the two lists has 
#'  to be the same. 
#' 
#' If the input is a genlight object, the function calculates the linear distance
#' for \code{Dgeo} and the relevant \code{Dgen} matrix (see \code{Dgen_method}) 
#' for each population. 
#' When the method selected is a genetic similarity matrix (e.g. "simple" 
#' distance), the matrix is internally transformed with \code{1 - Dgen} so that 
#' positive values of autocorrelation coefficients indicates more related 
#' individuals similarly as implemented in GenAlEx. If the user provide the 
#' distance matrices, care must be taken in interpreting the results because
#' similarity matrix will generate negative values for closely related 
#' individuals.
#' 
#' If \code{max(Dgeo)>1000} (e.g. the geographic distances are in thousands of 
#' metres), values are divided by 1000 (in the example before these would then 
#' become km) to facilitate readability of the plots.
#'
#' If \code{bins} is of length = 1 it is interpreted as the number of (even)
#' bins to use. In this case the starting point is always the minimum value in 
#' the distance matrix, and the last is the maximum. If it is a numeric vector 
#' of length>1, it is interpreted as the breaking points. In this case, the 
#' first has to be the lowest value, and the last has to be the highest. There 
#' are no internal checks for this and it is user responsibility to ensure that
#' distance classes are properly set up. If that is not the case, data that fall
#' outside the range provided will be dropped. The number of bins will be 
#' \code{length(bins) - 1}.
#'
#' The permutation constructs the 95\% confidence intervals around the null
#' hypothesis of no spatial structure (this is a two-tail test). The same data
#' are also used to calculate the probability of the one-tail test (See 
#' references below for details).
#'
#' Bootstrap calculations are skipped and \code{NA} is returned when the number 
#' of possible combinations given the sample size of any given distance class is
#' < \code{reps}.
#' 
#' Methods available to calculate genetic distances for SNP data:
#' \itemize{
#' \item "propShared" using the function \code{\link{gl.propShared}}.
#' \item "grm" using the function \code{\link{gl.grm}}.
#' \item "Euclidean" using the function \code{\link{gl.dist.ind}}.
#' \item "Simple" using the function \code{\link{gl.dist.ind}}.
#' \item "Absolute" using the function \code{\link{gl.dist.ind}}.
#' \item "Manhattan" using the function \code{\link{gl.dist.ind}}.
#' }
#' 
#' Methods available to calculate genetic distances for SilicoDArT data:
#' \itemize{
#' \item "Euclidean" using the function \code{\link{gl.dist.ind}}.
#' \item "Simple" using the function \code{\link{gl.dist.ind}}.
#' \item "Jaccard" using the function \code{\link{gl.dist.ind}}.
#' \item "Bray-Curtis" using the function \code{\link{gl.dist.ind}}.
#' }
#'   
#'  Plots and table are saved to the temporal directory (tempdir) and can be
#'  accessed with the function \code{\link{gl.print.reports}} and listed with
#'  the function \code{\link{gl.list.reports}}. Note that they can be accessed
#'  only in the current R session because tempdir is cleared each time that the
#'   R session is closed.
#'
#' Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'  
#' @param x Genlight object [default NULL].
#' @param Dgen Genetic distance matrix if no genlight object is provided
#' [default NULL].
#' @param Dgeo Geographic distance matrix if no genlight object is provided. This
#'  is typically an Euclidean distance but it can be any meaningful 
#'  (geographical) distance metrics [default NULL].
#' @param coordinates Can be either 'latlon', 'xy' or a two column data.frame
#' with column names 'lat','lon', 'x', 'y')  Coordinates are provided via
#' \code{gl@other$latlon} ['latlon'] or via \code{gl@other$xy} ['xy']. If latlon
#' data will be projected to meters using Mercator system [google maps] or if
#' xy then distance is directly calculated on the coordinates [default .
#' @param Dgen_method Method to calculate genetic distances. See details
#'  [default "Euclidean"].
#' @param Dgeo_trans Transformation to be used on the geographic distances. See
#' Dgen_trans [default "Dgeo"].
#' @param Dgen_trans You can provide a formula to transform the genetic
#' distance. The transformation can be applied as a formula using Dgen as the
#'  variable to be transformed. For example: \code{Dgen_trans = 'Dgen/(1-Dgen)'.
#'   Any valid R expression can be used here 
#'   [default 'Dgen', which is the identity function.]}
#' @param bins The number of bins for the distance classes
#' (i.e. \code{length(bins) == 1)} or a vectors with the break points. See 
#' details [default 5].
#' @param reps The number to be used for permutation and bootstrap analyses
#' [default 100].
#' @param plot.pops.together Plot all the populations in one plot. Confidence 
#' intervals from permutations are not shown [default FALSE].
#' @param permutation Whether permutation calculations for the null hypothesis 
#' of no spatial structure should be carried out [default TRUE].
#' @param bootstrap Whether bootstrap calculations to compute the 95\% confidence
#' intervals around r should be carried out [default TRUE].
#' @param plot_theme Theme for the plot. See details [default NULL].
#' @param plot_colors_pop A color palette for populations or a list with
#' as many colors as there are populations in the dataset [default NULL].
#' @param CI_color Color for the shade of the 95\% confidence intervals around 
#' the r estimates [default "red"].
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session
#' temporary directory (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log ; 3, progress and results summary; 5, full report [default
#'   NULL, unless specified using gl.set.verbosity].
#' @return Returns a data frame with the following columns:
#' \enumerate{
#' \item Bin  The distance classes
#' \item N The number of pairwise comparisons within each distance class
#' \item r.uc The uncorrected autocorrelation coefficient
#' \item Correction the correction
#' \item r The corrected autocorrelation coefficient
#' \item L.r The corrected autocorrelation coefficient lower limit
#' (if \code{bootstap = TRUE})
#' \item U.r The corrected autocorrelation coefficient upper limit
#' (if \code{bootstap = TRUE})
#' \item L.r.null.uc The uncorrected lower limit for the null hypothesis of no 
#' spatial autocorrelation (if \code{permutation = TRUE})
#' \item U.r.null.uc  The uncorrected upper limit for the null hypothesis of no 
#' spatial autocorrelation (if \code{permutation = TRUE})
#' \item L.r.null The corrected lower limit for the null hypothesis of no 
#' spatial autocorrelation (if \code{permutation = TRUE})
#' \item U.r.null The corrected upper limit for the null hypothesis of no 
#' spatial autocorrelation (if \code{permutation = TRUE})
#' \item p.one.tail The p value of the one tail statistical test
#' }
#'
#' @author Carlo Pacioni, Bernd Gruber & Luis Mijangos 
#' (Post to \url{https://groups.google.com/d/forum/dartr})
#' @references
#' \itemize{
#' \item Smouse PE, Peakall R. 1999. Spatial autocorrelation analysis of
#' individual multiallele and multilocus genetic structure. Heredity 82:
#' 561-573.
#' \item Double, MC, et al. 2005. Dispersal, philopatry and infidelity: 
#' dissecting local genetic structure in superb fairy-wrens (Malurus cyaneus). 
#' Evolution 59, 625-635.
#' \item Peakall, R, et al. 2003. Spatial autocorrelation analysis offers new
#' insights into gene flow in the Australian bush rat, Rattus fuscipes.
#' Evolution 57, 1182-1195.
#' \item Smouse, PE, et al. 2008. A heterogeneity test for fine-scale genetic
#' structure. Molecular Ecology 17, 3389-3400.
#' \item Gonzales, E, et al. 2010. The impact of landscape disturbance on 
#' spatial genetic structure in the Guanacaste tree, Enterolobium
#' cyclocarpum (Fabaceae). Journal of Heredity 101, 133-143.
#' \item Beck, N, et al. 2008. Social constraint and an absence of sex-biased
#' dispersal drive fine-scale genetic structure in white-winged choughs.
#' Molecular Ecology 17, 4346-4358.
#' }
#' @examples
#' require("dartR.data")
#' res <- gl.spatial.autoCorr(platypus.gl, bins=seq(0,10000,2000))
#' # using one population, showing sample size
#' test <- gl.keep.pop(platypus.gl,pop.list = "TENTERFIELD")
#' res <- gl.spatial.autoCorr(test, bins=seq(0,10000,2000),CI_color = "green")
#' @importFrom tidyr pivot_wider
#' @export

gl.spatial.autoCorr <- function(x = NULL,
                                Dgeo = NULL,
                                Dgen = NULL,
                                coordinates = "latlon", 
                                Dgen_method = "Euclidean",
                                Dgeo_trans = "Dgeo",
                                Dgen_trans = "Dgen",
                                bins = 5,
                                reps = 100,
                                plot.pops.together = FALSE,
                                permutation = TRUE,
                                bootstrap = TRUE,
                                plot_theme = NULL,
                                plot_colors_pop = NULL,
                                CI_color = "red",
                                plot.out = TRUE,
                                save2tmp = FALSE,
                                verbose = NULL) {
  
  # CHECK IF PACKAGES ARE INSTALLED
  pkg <- "dismo"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }

  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jackson",
                   verbosity = verbose)
  
  # CHECK DATATYPE
  if (!is.null(x)) {
    dt <- utils.check.datatype(x, verbose = 0)
  }
  
  # specific error checks
  if (!is.numeric(bins)) {
    stop(error("  The argument 'bins' should be a numeric vector\n"))
  }
  
  if (!is.null(Dgen) & !is.null(Dgeo)) {
    if (verbose > 0)
      cat(
        report(
          "  Analysis performed using provided genetic and Euclidean distance matrices. If a genlight object is provided, it is ignored.\n"
        )
      )
    ta <-"dgendgeo"
  }
  
  if (is(x, "genlight")) {
    if (verbose > 0)
      cat(report("  Analysis performed on the genlight object.\n"))
    ta <-"genlight"
  }
  
  # avoid global binding error
  Bin <-
    r <-
    L.r <-
    U.r <- L.r.null <- U.r.null <- Freq <- Var1 <-  NULL 
  
  # DO THE JOB #
  
  #### if a genlight object is provided ####
  if (!is.null(x) & is(x, "genlight")) {
 
  pop_list <- seppop(x)
  
  Dgen_list <- list()
  Dgeo_list <- list()
  
  for (i in seq_along(pop_list)) {
    
    if (verbose > 2) {
      cat(
        report(paste("  Analysing population",names(pop_list[i]),"\n")
        )
      )
    }
    
    x_temp <- pop_list[[i]]
    
  # check coordinates (if no Dgen and Dgeo is provided)
    coords <- NULL
    if (is(coordinates, "character")) {
      if (coordinates == "latlon") {
        if (is.null(x_temp@other$latlon))
          stop(error(
            "Cannot find coordinates in x@other$latlon"
          ))
        coords <- dismo::Mercator(x_temp@other$latlon[, c("lon", "lat")])
        coordstring <-"x@other$latlon (Mercator transformed)"
      }
      
      if (coordinates == "xy") {
        if (is.null(x_temp@other$xy))
          stop(error("Cannot find coordinates in x@other$xy"))
        coords <- x_temp@other$xy
        coordstring <-"x@other$xy"
      }
    }
    
    if (is(coordinates, "data.frame")) {
      if (length(setdiff(colnames(coordinates), c("lat", "lon"))) == 0) {
        coords <- dismo::Mercator(coordinates[, c("lon", "lat")])
        coordstring <-"data.frame lat/lon (Mercator transformed)"
      }
      
      if (length(setdiff(colnames(coordinates), c("x", "y"))) == 0) {
        coords <- coordinates[, c("x", "y")]
        coordstring <-"data.frame x/y"
      }
      
      if (is.null(coords)) {
        stop(
          error(
            "No valid coordinates provided. Check the provided data.frame and its format.\n"
          )
        )
      }
    }
    
    if (is.null(coords)) {
      stop(error("No valid coordinates provided!\n"))
    }
    
    # make sure coordinates have the correct length
    if (nrow(coords) != nInd(x_temp) & ta == "genlight") {
      stop(error(
        "Cannot find coordinates for each individual in slot @other$latlon.\n"
      ))
    }
    
      if (nInd(x_temp) > 1) {
        Dgeo <- dist(coords)
      } else {
        stop(
          error(
            "Less than 2 individuals provided, therefore no pairwise distances can be calculated.\n"
          )
        )
      }
    
    # calculate genetic distances
    if (Dgen_method == "propShared") {
      Dgen <- as.dist(gl.propShared(x_temp))
    } else {
      if (Dgen_method == "grm") {
        Dgen <- as.dist(gl.grm(x_temp, plotheatmap=FALSE, verbose = 0))
          } else {
            Dgen <- gl.dist.ind(x_temp, method = Dgen_method, plot.out = FALSE,
                          verbose = 0)
          }
    }
    if ((dt == "SNP" &
        Dgen_method == "propShared" |
        Dgen_method == "Euclidean" |
        Dgen_method == "Simple" |
        Dgen_method == "Absolute") |
      (dt == "SilicoDArT" & Dgen_method == "Euclidean")) {
      
      # Reverse genetic distance matrix so that correlated values
      # indicated more similar individuals as we are used to see plots in GenAleEx
      Dgen <- 1 - Dgen
    }
    
    distance <- Dgen_method

  # convert matrices to distance objects
    Dgen_list[[i]] <- as.dist(Dgen)
    Dgeo_list[[i]] <- as.dist(Dgeo)
      } # Close for(i in 1:length(pop_list))
  pop.names <- popNames(x)
  } # close if a genlight object is provided


  #### if distances are provided ####
  if (is.null(x)) {
    chk.D <- function(D, name.D) {
      if(!(is(D, "dist") | is.matrix(D))) 
        stop(error(paste0("  ", name.D, 
                          " is neither a list, a matrix nor a distance\n")))
      if(is.matrix(D)) D <- as.dist(D)
      return(list(D))
    }
    
    if(!is(Dgeo, "list")) {
      Dgeo_list <- chk.D(D=Dgeo, name.D = "Dgeo")
    } else {
      Dgeo_list <- Dgeo
    }
    if(!is(Dgen, "list")) {
      Dgen_list <- chk.D(D=Dgen, name.D = "Dgen")
    }  else {
      Dgen_list <- Dgen
    }

    # now distances are a list
    if(length(Dgeo_list) != length(Dgen_list))
      stop(error( "  The arguments Dgen and Dgeo should be of same length\n"))

    if(is.null(names(Dgeo_list))) pop.names <- paste0("Pop", seq_along(Dgen_list))
    
    chk.D.list <- function(D.list, name.D) {
      if(length(unique(sapply(D.list, class))) != 1) {
        stop(error(paste0("  ", name.D, 
              " is a list, but its elements are of different classes. These should be either all matrices or distances\n")))
      }
      
      if(!(is(D.list[[1]], "dist") | is.matrix(D.list[[1]]))) {
        stop(error(paste0("  ", name.D, 
              " is a list, but its element are neither all matrices nor distances\n")))
      }
      
      if(is.matrix(D.list[[1]])) D.list <- lapply(D.list, as.dist)
      return(D.list)
    }
    
    Dgeo_list <- chk.D.list(Dgeo_list, name.D = "Dgeo")
    Dgen_list <- chk.D.list(Dgen_list, name.D = "Dgen")
      
    len.elements.Dgeo <- sapply(Dgeo_list, length)
    len.elements.Dgen <- sapply(Dgen_list, length)
    
    if (is.character(all.equal(len.elements.Dgeo, len.elements.Dgen))) {
      stop(error("  The arguments Dgen and Dgeo should have identical dimensions\n"))
    }
    coordstring = "Dgeo provided."
    distance = "Dgen provided"
    typedis = "ind"
  } # Close if matrices are provided
  #----------------------------------------------------------------------------#

  #### Apply transformations ####
  apply.transformation <- function(D, transFUN, name.D) {
    assign(name.D, value = D)
    new.obj <- eval(parse(text = transFUN))
    return(new.obj)
  }

  Dgen_list <- lapply(Dgen_list, apply.transformation, transFUN=Dgen_trans, name.D="Dgen")
  Dgeo_list <- lapply(Dgeo_list, apply.transformation, transFUN=Dgeo_trans, name.D="Dgeo")
      
  lapply(Dgeo_list, function(D) {
    if(sum(is.infinite(D)) > 0) {
      stop(
        error(
          "Most likely some pairwise individual distances were zero and the transformation created missing values [e.g. log(Dgeo)]. Consider adding a suitable tranformation e.g. an offset to your Dgeo transformation if using a log transformation [e.g. Dgeo_trans='log(Dgeo+1)'] or adding some 'noise' to the coordinates.\n"
        )
      )
    }

  } 
  )
  
  convert2matrix <- function(D) {
    D <- as.matrix(D)
    diag(D) <- 0
    return(D)
  }
  
  Dgen_list <- lapply(Dgen_list, convert2matrix) 
  Dgeo_list <- lapply(Dgeo_list, convert2matrix)
  
  #### Execute utils.spautocorr on a list ####
      res <- list()
      
      for(z in seq_along(Dgeo_list)) {
        
        Dgeo <- Dgeo_list[[z]]
        Dgen <- Dgen_list[[z]]
        
        sample.size <- nrow(Dgeo)
        crt <- 1 / (sample.size - 1) # correction
        nbins <- if (length(bins) == 1) {
          bins
        } else {
          length(bins) - 1
        }
        
        splist <-
          utils.spautocor(Dgen,
                          Dgeo,
                          permutation = FALSE,
                          bins = bins,
                          reps = reps)
        
        if (permutation) {
          bssplist <- replicate(reps,
                                utils.spautocor(
                                  Dgen,
                                  Dgeo,
                                  permutation = TRUE,
                                  bins = bins,
                                  reps = reps
                                ))
          
          #convert the output into a matrix
          bs <- matrix(
            unlist(bssplist),
            nrow = reps,
            ncol = nbins,
            byrow = TRUE
          )
          bs.l <- apply(bs, 2, quantile, probs = 0.025, na.rm = TRUE)
          bs.u <- apply(bs, 2, quantile, probs = 0.975, na.rm = TRUE)
          
          p.one.tail <-
            sapply(seq_along(splist$r.uc), function(i, r.rc, r, crt = crt) {
              if (is.na(r[i])) {
                NA
              } else{
                if (r[i] >= 0) {
                  sum(r.rc[, i] >= r[i]) / length(r.rc[, i])
                } else{
                  sum(r.rc[, i] <= r[i]) / length(r.rc[, i])
                }
              }
            }, r = splist$r.uc + crt,  r.rc = bs + crt)
          
        }
        
        if (bootstrap) {
          errors <-
            replicate(reps,
                      utils.spautocor(
                        Dgen,
                        Dgeo,
                        bootstrap = TRUE,
                        bins = bins,
                        reps = reps
                      ))
          errors <-
            matrix(unlist(errors),
                   nrow = reps,
                   ncol = nbins,
                   byrow = TRUE)
          err.l <- apply(errors, 2, quantile, probs = 0.025, na.rm = TRUE)
          err.u <- apply(errors, 2, quantile, probs = 0.975, na.rm = TRUE)
        }
        
        res_temp <- cbind(splist, Correction = crt, r = splist$r.uc + crt)
        if (bootstrap) {
          res_temp <- cbind(res_temp, L.r = err.l + crt, U.r = err.u + crt)
        }
        
        if (permutation) {
          res_temp <- cbind(
            res_temp,
            L.r.null.uc = bs.l,
            U.r.null.uc = bs.u,
            L.r.null = bs.l + crt,
            U.r.null = bs.u + crt,
            p.one.tail = p.one.tail
          )
        }
        res[[z]] <- res_temp
      }
      
      names(res) <- pop.names
  #-------- Close Execute utils.spautoCorr ------------------------------------#

#### PRINTING OUTPUTS ####

  if (plot.out) {
    
    if (is.null(plot_theme)) {
      plot_theme <- theme_dartR()
    }
    
      if (is.null(plot_colors_pop)) {
        plot_colors_pop <- discrete_palette
      }
      
      if (is(plot_colors_pop, "function")) {
        if(is.null(x)) n.pop <- length(Dgeo_list) else
          n.pop <- nPop(x)
        plot_colors_pop <- plot_colors_pop(n.pop)
      }
      
      if (!is(plot_colors_pop, "function")) {
        plot_colors_pop <- plot_colors_pop
      }
      
      spa_multi <-data.table::rbindlist(res, use.names = TRUE, 
                                        fill = TRUE, idcol = "Population")
      if(spa_multi [, max(Bin)] > 1000) {
        lbls <- round(spa_multi$Bin/1000, 0) 
      } else {
        lbls <- spa_multi$Bin
      }
      
      p3 <- ggplot(spa_multi, aes_string("Bin", "r", col="Population")) +
        geom_line(size=1) +
        geom_point(size=2) +
        geom_hline(yintercept = 0, col = "black", size=1) +
        scale_color_manual(values = plot_colors_pop) +
        scale_x_continuous(breaks = spa_multi$Bin,
                      labels = lbls) +
        ylab("Autocorrelation (r)") + 
        xlab("Distance class") + 
        plot_theme
      
      if (bootstrap) {
        p3 <- p3 +   
          geom_errorbar(aes(ymin=L.r, ymax=U.r), 
                        width=spa_multi[, mean(tail(Bin, -1) - head(Bin, -1))]/10) 
      }

      if (permutation & plot.pops.together == FALSE) {
        p3 <- p3 +  
          geom_ribbon(aes(ymin=L.r.null,ymax=U.r.null), fill = CI_color, 
                      alpha=0.25,show.legend = FALSE) + 
          geom_line(aes(y = L.r.null), col = "black", linetype = "dashed") +
          geom_point(aes(y = L.r.null), col = "black") +
          geom_line(aes(y = U.r.null), col = "black", linetype = "dashed") +
          geom_point(aes(y = U.r.null), col = "black") +
          facet_wrap(~Population, nrow = length(Dgen_list), scales = "free_y") +
           theme(
          #   strip.text.x = element_text(size = 12),
          #   axis.text.x = element_text(
          #     size = 12
          #   ), 
             legend.position = "none") 
      }
      
      if(length(Dgen_list) == 1) {
        p3 <- p3 + 
        scale_x_continuous(breaks = spa_multi$Bin,
                            labels = lbls,
                            sec.axis = sec_axis(
                                  trans = ~ .,
                                  breaks = spa_multi$Bin,
                                  labels = spa_multi$N)) +
        theme(strip.text = element_blank(), legend.position = "none")
      }
    suppressWarnings(
    suppressMessages(print(p3))
    )
    }
  
  if (verbose > 0) {
    cat(report("  Coordinates used from:", coordstring, "\n"))
    cat(report("  Transformation of Dgeo:", Dgeo_trans, "\n"))
    cat(report("  Genetic distance:", distance, "\n"))
    cat(report("  Tranformation of Dgen: ", Dgen_trans, "\n"))
    print(res)
  }
  
  # SAVE INTERMEDIATES TO TEMPDIR
  
  # creating temp file names
  if (save2tmp) {
    if (plot.out) {
      temp_plot <- tempfile(pattern = "Plot_")
      match_call <-
        paste0(names(match.call()),
               "_",
               as.character(match.call()),
               collapse = "_")
      # saving to tempdir
      saveRDS(list(match_call, p3), file = temp_plot)
      
      if (verbose >= 2) {
        cat(report("  Saving the ggplot to session tempfile\n"))
      }
    }
    
    temp_table <- tempfile(pattern = "Table_")
    saveRDS(list(match_call, res), file = temp_table)
    if (verbose >= 2) {
      cat(report("  Saving tabulation to session tempfile\n"))
      cat(
        report(
          "  NOTE: Retrieve output files from tempdir using gl.list.reports() and gl.print.reports()\n"
        )
      )
    }
  }
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # RETURN
  return(invisible(res))
  
}
