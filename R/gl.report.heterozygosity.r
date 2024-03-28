#' @name gl.report.heterozygosity
#' @title Reports observed, expected and unbiased heterozygosities and FIS
#' (inbreeding coefficient) by population or by individual from SNP data
#' @family unmatched report
#'
#' @description Calculates the observed, expected and unbiased expected (i.e.
#' corrected for sample size) heterozygosities and FIS (inbreeding coefficient)
#' for each population or the observed heterozygosity for each individual in a
#' genlight object.

#' @param x Name of the genlight object containing the SNP [required].
#' @param method Calculate heterozygosity by population (method='pop') or by
#' individual (method='ind') [default 'pop'].
#' @param n.invariant An estimate of the number of invariant sequence tags used
#' to adjust the heterozygosity rate [default 0].
#' @param nboots Number of bootstrap replicates to obtain confident intervals
#' [default 0].
#' @param conf The confidence level of the required interval  [default 0.95].
#' @param CI.type Method to estimate confident intervals. One of
#' "norm", "basic", "perc" or "bca" [default "bca"].
#' @param ncpus Number of processes to be used in parallel operation. If ncpus
#' > 1 parallel operation is activated,see "Details" section [default 1].
#' @param plot.display Specify if plot is to be produced [default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot.colors.pop A color palette for population plots or a list with
#' as many colors as there are populations in the dataset
#' [default gl.colors("dis")].
#' @param plot.colors.ind List of two color names for the borders and fill of
#' the plot by individual [default gl.colors(2)].
#' @param error.bar statistic to be plotted as error bar either "SD" (standard 
#' deviation) or "SE" (standard error) or "CI" (confident intervals)
#'  [default "SD"].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session
#' temporary directory (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, unless specified using gl.set.verbosity].
#' @details
#' Observed heterozygosity for a population takes the proportion of
#' heterozygous loci for each individual then averages over the individuals in
#' that population. The calculations take into account missing values.

#' Expected heterozygosity for a population takes the expected proportion of
#' heterozygotes, that is, expected under Hardy-Weinberg equilibrium, for each
#' locus, then averages this across the loci for an average estimate for the
#' population.
#'
#' Expected heterozygosity is calculated using the correction for sample size
#' following equation 2 from Nei 1978.
#'
#' Observed heterozygosity for individuals is calculated as the proportion of
#' loci that are heterozygous for that individual.
#'
#' Finally, the loci that are invariant across all individuals in the dataset
#' (that is, across populations), is typically unknown. This can render
#' estimates of heterozygosity analysis specific, and so it is not valid to
#' compare such estimates across species or even across different analyses. This
#' is a similar problem faced by microsatellites. If you have an estimate of the
#' number of invariant sequence tags (loci) in your data, such as provided by
#' \code{\link{gl.report.secondaries}}, you can specify it with the n.invariant
#' parameter to standardize your estimates of heterozygosity.
#'
#' \strong{NOTE}: It is important to realise that estimation of adjusted
#' heterozygosity requires that secondaries not to be removed.
#'
#' Heterozygosities and FIS (inbreeding coefficient) are calculated by locus
#' within each population using the following equations:
#' \itemize{
#' \item Observed heterozygosity (Ho) = number of homozygotes / n_Ind,
#' where n_Ind is the number of individuals without missing data.
#' \item Observed heterozygosity adjusted (Ho.adj) <- Ho * n_Loc /
#'  (n_Loc + n.invariant),
#' where n_Loc is the number of loci that do not have all missing data  and
#' n.invariant is an estimate of the number of invariant loci to adjust
#' heterozygosity.
#' \item Expected heterozygosity (He) = 1 - (p^2 + q^2),
#' where p is the frequency of the reference allele and q is the frequency of
#' the alternative allele.
#' \item Expected heterozygosity adjusted (He.adj) = He * n_Loc /
#' (n_Loc + n.invariant)
#' \item Unbiased expected heterozygosity (uHe) = He * (2 * n_Ind /
#' (2 * n_Ind - 1))
#' \item Inbreeding coefficient (FIS) = 1 - (mean(Ho) / mean(uHe))
#' }

#'\strong{ Function's output }

#' Output for method='pop' is an ordered barchart of observed heterozygosity,
#' unbiased expected heterozygosity and FIS (Inbreeding coefficient) across populations
#' together with a table of mean observed and expected heterozygosities and FIS
#' by population and their respective standard deviations (SD).

#' In the output, it is also reported by population: the number of loci used to
#'  estimate heterozygosity(n.Loc), the number of polymorphic loci (polyLoc),
#'  the number of monomorphic loci (monoLoc) and loci with all missing data
#'   (all_NALoc).

#' Output for method='ind' is a histogram and a boxplot of heterozygosity across
#' individuals.

#'  Plots and table are saved to the session temporary directory (tempdir)

#'  Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'  
#'   \strong{Error bars}
#'  
#'  The best method for presenting or assessing genetic statistics depends on 
#'  the type of data you have and the specific questions you're trying to 
#'  answer. Here's a brief overview of when you might use each method:
#'  
#'   \strong{1. Confidence Intervals ("CI"):}
#'   
#'  - Usage: Often used to convey the precision of an estimate.
#'  
#'  - Advantage: Confidence intervals give a range in which the true parameter 
#'  (like a population mean) is likely to fall, given the data and a specified 
#'  probability (like 95%).
#'  
#'  - In Context: For genetic statistics, if you're estimating a parameter,
#'   a 95% CI gives you a range in which you're 95% confident the true parameter
#'    lies.
#'  
#'   \strong{2. Standard Deviation ("SD"):}
#'   
#'  - Usage: Describes the amount of variation from the average in a set of data.
#'  
#'  - Advantage: Allows for an understanding of the spread of individual data
#'   points around the mean.
#'   
#'  - In Context: If you're looking at the distribution of a quantitative trait 
#'  (like height) in a population with a particular genotype, the SD can 
#'  describe how much individual heights vary around the average height.
#'  
#'   \strong{3. Standard Error ("SE"):}
#'   
#'  - Usage: Describes the precision of the sample mean as an estimate of the 
#'  population mean.
#'  
#'  - Advantage: Smaller than the SD in large samples; it takes into account 
#'  both the SD and the sample size. 
#'  
#'  - In Context: If you want to know how accurately your sample mean represents
#'   the population mean, you'd look at the SE.
#'   
#'    \strong{Recommendation:}
#'    
#'   - If you're trying to convey the precision of an estimate, confidence 
#'   intervals are very useful.
#'   
#'   - For understanding variability within a sample, standard deviation is key.
#'   
#'   - To see how well a sample mean might estimate a population mean, consider 
#'   the standard error.
#'   
#'   In practice, geneticists often use a combination of these methods to 
#'   analyze and present their data, depending on their research questions and 
#'   the nature of the data.
#'   
#'  \strong{Confident Intervals}
#'
#' The uncertainty of a parameter, in this case the mean of the statistic, can
#' be summarised by a confidence interval (CI) which includes the true parameter
#' value with a specified probability (i.e. confidence level; the parameter
#' "conf" in this function).
#'
#' In this function, CI are obtained using Bootstrap which is an inference
#' method that samples with replacement the data (i.e. loci) and calculates the
#'  statistics every time.
#'
#'  This function uses the function \link[boot]{boot} (package boot) to perform
#'  the bootstrap replicates and the function \link[boot]{boot.ci}
#'  (package boot) to perform the calculations for the CI.
#'
#'  Four different types of nonparametric CI can be calculated
#'   (parameter "CI.type" in this function):
#'   \itemize{
#'    \item First order normal approximation interval ("norm").
#'    \item Basic bootstrap interval ("basic").
#'    \item Bootstrap percentile interval ("perc").
#'    \item Adjusted bootstrap percentile interval ("bca").
#'    }
#'
#' The studentized bootstrap interval ("stud") was not included in the CI types
#'  because it is computationally intensive, it may produce estimates outside
#'  the range of plausible values and it has been found to be erratic in
#'  practice, see for example the "Studentized (t) Intervals" section in:
#'
#'    \url{https://www.r-bloggers.com/2019/09/understanding-bootstrap-confidence-interval-output-from-the-r-boot-package}
#'
#'     Nice tutorials about the different types of CI can be found in:
#'
#'     \url{https://www.datacamp.com/tutorial/bootstrap-r}
#'
#'     and
#'
#'    \url{https://www.r-bloggers.com/2019/09/understanding-bootstrap-confidence-interval-output-from-the-r-boot-package}
#'
#'      Efron and Tibshirani (1993, p. 162) and Davison and Hinkley
#'      (1997, p. 194) suggest that the number of bootstrap replicates should
#'      be between 1000 and 2000.
#'
#'  \strong{It is important} to note that unreliable confident intervals will be
#'   obtained if too few number of bootstrap replicates are used.
#'   Therefore, the function \link[boot]{boot.ci} will throw warnings and errors
#'    if bootstrap replicates are too few. Consider increasing then number of
#'    bootstrap replicates to at least 200.
#'
#'    The "bca" interval is often cited as the best for theoretical reasons,
#'    however it may produce unstable results if the bootstrap distribution
#'     is skewed or has extreme values. For example, you might get the warning
#'     "extreme order statistics used as endpoints" or the error "estimated
#'     adjustment 'a' is NA". In this case, you may want to use more bootstrap
#'     replicates or a different method or check your data for outliers.
#'
#'    The error "estimated adjustment 'w' is infinite" means that the estimated
#'    adjustment ‘w’ for the "bca" interval is infinite, which can happen when
#'    the empirical influence values are zero or very close to zero. This can
#'    be caused by various reasons, such as:
#'
#'    The number of bootstrap replicates is too small, the statistic of interest
#'     is constant or nearly constant across the bootstrap samples, the data
#'     contains outliers or extreme values.
#'
#'     You can try some possible solutions, such as:
#'
#' Increasing the number of bootstrap replicates, using a different type of
#' bootstrap confidence interval or removing or transforming the outliers or
#'  extreme values.
#'  
#'  \strong{Parallelisation}
#'
#'  If the parameter ncpus > 1, parallelisation is enabled. In Windows, parallel
#'   computing employs a "socket" approach that starts new copies of R on each
#'    core. POSIX systems, on the other hand (Mac, Linux, Unix, and BSD),
#'    utilise a "forking" approach that replicates the whole current version of
#'     R and transfers it to a new core.
#'
#'     Opening and terminating R sessions in each core involves a significant
#'     amount of processing time, therefore parallelisation in Windows machines
#'    is only quicker than not usung parallelisation when nboots > 1000-2000.
#' @author Custodian: Luis Mijangos (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#'
#' @references
#' Nei, M. (1978). Estimation of average heterozygosity and genetic distance
#' from a small number of individuals. Genetics, 89(3), 583-590.
#'
#' @examples
#'  \donttest{
#' require("dartR.data")
#' df <- gl.report.heterozygosity(platypus.gl)
#' df <- gl.report.heterozygosity(platypus.gl,method='ind')
#' n.inv <- gl.report.secondaries(platypus.gl)
#' gl.report.heterozygosity(platypus.gl, n.invariant = n.inv[7, 2])
#' }
#' df <- gl.report.heterozygosity(platypus.gl)

#' @seealso \code{\link{gl.filter.heterozygosity}}

#' @export
#' @return A dataframe containing population labels, heterozygosities, FIS,
#' their standard deviations and sample sizes

gl.report.heterozygosity <- function(x,
                                     method = "pop",
                                     n.invariant = 0,
                                     nboots = 0,
                                     conf = 0.95,
                                     CI.type = "bca",
                                     ncpus = 1,
                                     plot.display = TRUE,
                                     plot.theme = theme_dartR(),
                                     plot.colors.pop = gl.colors("dis"),
                                     plot.colors.ind = gl.colors(2),
                                     error.bar = "SD",
                                     save2tmp = FALSE,
                                     verbose = NULL) {
  # standard error function
  std.error <- function(x) {
    res <- sd(x, na.rm = TRUE) / sqrt(length(x))
    return(res)
  }
  
  # bootstrapping function
  pop.het <- function(x,
                      indices,
                      n.invariant) {
    pop.het_fun <- function(df,
                            n.invariant) {
      Ho.loc <- colMeans(df == 1, na.rm = TRUE)
      n_loc <- apply(df, 1, function(y) {
        sum(!is.na(y))
      })
      Ho.adj.loc <- Ho.loc * n_loc / (n_loc + n.invariant)
      q_freq <- colMeans(df, na.rm = TRUE) / 2
      p_freq <- 1 - q_freq
      He.loc <- 2 * p_freq * q_freq
      n_ind <- apply(df, 2, function(y) {
        sum(!is.na(y))
      })
      ### CP ### Unbiased He (i.e. corrected for sample size) hard
      # coded for diploid
      uHe.loc <-
        (2 * as.numeric(n_ind) / (2 * as.numeric(n_ind) - 1)) * He.loc
      Hexp.adj.loc <- He.loc * n_loc / (n_loc + n.invariant)
      
      FIS.loc <- 1 - (Ho.loc / He.loc)
      
      all.res <- c(
        Ho.loc = mean(Ho.loc, na.rm = TRUE),
        Ho.adj.loc = mean(Ho.adj.loc, na.rm = TRUE),
        He.loc = mean(He.loc, na.rm = TRUE),
        uHe.loc = mean(uHe.loc, na.rm = TRUE),
        Hexp.adj.loc = mean(Hexp.adj.loc, na.rm = TRUE),
        FIS.loc = mean(FIS.loc, na.rm = TRUE)
      )
      
      return(all.res)
    }
    
    df <- x[, indices]
    
    res <- pop.het_fun(df,
                       n.invariant = n.invariant)
    
    return(res)
    
  }
  
  # Counting individuals function
  ind.count <- function(x) {
    # the loci that are completely missing
    loci.na <-
      which(colSums(is.na(as.matrix(x))) == nrow(as.matrix(x)))
    # the number of samples in the matrix the number of non-genotyped
    # samples remove the loci that are completely missing
    if (length(loci.na) > 0) {
      nind <-
        mean(nrow(as.matrix(x)) -
               colSums(is.na(as.matrix(x)))[-loci.na])
      # the number of samples in the matrix the number of
      # non-genotyped samples
    } else {
      nind <- mean(nrow(as.matrix(x)) - colSums(is.na(as.matrix(x))))
    }
    
    return(nind)
  }
  
  # setting parallel
  # if (ncpus > 1) {
  if (grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
    parallel <- "multicore"
  }
  ## if windows
  if (!grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
    parallel <- "snow"
  }
  # }
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jackson",
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype <-
    utils.check.datatype(x, accept = "SNP", verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
  if (!(method == "pop" | method == "ind")) {
    cat(
      warn(
        "Warning: Method must either be by population or by individual,
                set to method='pop'\n"
      )
    )
    method <- "pop"
  }
  
  if (n.invariant < 0) {
    cat(warn(
      "Warning: Number of invariant loci must be non-negative, set to
            zero\n"
    ))
    n.invariant <- 0
    if (verbose == 5) {
      cat(report(
        "  No. of invariant loci can be esimated using
                    gl.report.secondaries\n"
      ))
    }
  }
  
  if (any(grepl(x@other$history, pattern = "gl.filter.secondaries") == TRUE) &
      n.invariant > 0) {
    cat(
      warn(
        "  Warning: Estimation of adjusted heterozygosity requires that
                secondaries not to be removed. A gl.filter.secondaries call was
                found in the history. This may cause the results to be
                incorrect\n"
      )
    )
  }
  
  if( nboots == 0 & error.bar == "CI" ){
    cat(error(
      "  Number of boostraps ('nboots' parameter) must be > 0 to calculate confident 
   intervals \n"))
    stop()
  }
  
  # DO THE JOB
  
  ########### FOR METHOD BASED ON POPULATIONS
  
  if (method == "pop") {
    # Set a population if none is specified (such as if the genlight object
    # has been generated manually)
    if (is.null(pop(x)) |
        is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2) {
        cat(
          warn(
            "  No population assignments detected,
                             individuals assigned to a single population
                        labelled 'pop1'\n"
          )
        )
      }
      pop(x) <- array("pop1", dim = nInd(x))
      pop(x) <- as.factor(pop(x))
    }
    
    # Split the genlight object into a list of populations
    sgl <- seppop(x)
    
    # OBSERVED HETEROZYGOSITY
    if (verbose >= 2) {
      cat(
        report(
          "  Calculating Observed Heterozygosities, averaged across
                    loci, for each population\n"
        )
      )
    }
    
    # Calculate heterozygosity for each population in the list
    # CP = Carlo Pacioni CP ###
    Ho.loc <-
      lapply(sgl, function(x)
        colMeans(as.matrix(x) == 1, na.rm = TRUE))
    
    Ho <-
      unlist(lapply(sgl, function(x)
        mean(
          colMeans(as.matrix(x) == 1, na.rm = TRUE), na.rm = TRUE
        )))
    
    ### CP ### observed heterozygosity standard deviation
    HoSD <-
      unlist(lapply(sgl, function(x)
        sd(
          colMeans(as.matrix(x) == 1, na.rm = TRUE), na.rm = TRUE
        )))
    
    HoSE <- unlist(lapply(sgl, function(x)
      std.error(colMeans(
        as.matrix(x) == 1, na.rm = TRUE
      ))))
    
    ##########
    
    # Calculate the number of loci that are not all NAs CP ###
    n_loc <-
      unlist(lapply(sgl, function(x)
        sum(!(
          colSums(is.na(as.matrix(x))) == nrow(as.matrix(x))
        ))))
    ##########
    
    # calculate the number of polymorphic and monomorphic loci by population
    poly_loc <- NULL
    mono_loc <- NULL
    all_na_loc <- NULL
    
    for (y in 1:length(sgl)) {
      y_temp <- sgl[[y]]
      hold <- y_temp
      mono_tmp <- gl.allele.freq(y_temp, simple = TRUE, verbose = 0)
      loc.list <- rownames(mono_tmp[which(mono_tmp$alf1 == 1 |
                                            mono_tmp$alf1 == 0), ])
      loc.list_NA <- which(colSums(is.na(as.matrix(y_temp)))==nInd(y_temp))
      # rownames(mono_tmp[which(is.na(mono_tmp$alf1)), ])
      
      # Remove NAs from list of monomorphic loci and loci with all NAs
      # loc.list <- loc.list[!is.na(loc.list)]
      
      # remove monomorphic loci and loci with all NAs
      if (length(loc.list) > 0) {
        y_temp <- gl.drop.loc(y_temp, loc.list = loc.list, verbose = 0)
      }
      
      poly_loc <-  c(poly_loc, nLoc(y_temp))
      mono_loc <- c(mono_loc, (nLoc(hold) - nLoc(y_temp) - length(loc.list_NA)))
      all_na_loc <- c(all_na_loc, length(loc.list_NA))
      
    }
    
    # Apply correction CP ###
    Ho.adj <- Ho * n_loc / (n_loc + n.invariant)
    # Manually compute SD for Ho.adj sum of the square of differences from
    # the mean for polymorphic sites plus sum of the square of differences
    #(which is the Ho.adj because Ho=0) from the mean for invariant sites
    Ho.adjSD <-
      sqrt((
        mapply(function(x, Mean)
          sum((x - Mean) ^ 2, na.rm = TRUE), Ho.loc, Mean = Ho.adj) +
          n.invariant * Ho.adj ^
          2
      ) / (n_loc +
             n.invariant - 1))
    
    Ho.adjSE <-  Ho.adjSD / sqrt(poly_loc+mono_loc)
    
    n_ind <- sapply(sgl, ind.count)
    
    ##########
    
    # EXPECTED HETEROZYGOSITY
    if (verbose >= 2) {
      cat(report("  Calculating Expected Heterozygosities\n\n"))
    }
    
    Hexp <- array(NA, length(sgl))
    HexpSD <- array(NA, length(sgl))
    HexpSE <- array(NA, length(sgl))
    
    uHexp <- array(NA, length(sgl))
    uHexpSD <- array(NA, length(sgl))
    uHexpSE <- array(NA, length(sgl))
    
    Hexp.adj <- array(NA, length(sgl))
    Hexp.adjSD <- array(NA, length(sgl))
    Hexp.adjSE <- array(NA, length(sgl))
    
    FIS <- array(NA, length(sgl))
    FISSD <- array(NA, length(sgl))
    FISSE <- array(NA, length(sgl))
    
    # dataframes to store data for boxplots
    uHe_df <- as.list(rep(NA, length(sgl)))
    Fis_df <- as.list(rep(NA, length(sgl)))
    
    # For each population
    for (i in 1:length(sgl)) {
      gl <- sgl[[i]]
      gl <- utils.recalc.freqhomref(gl, verbose = 0)
      gl <- utils.recalc.freqhomsnp(gl, verbose = 0)
      gl <- utils.recalc.freqhets(gl, verbose = 0)
      p <- gl@other$loc.metrics$FreqHomRef
      q <- gl@other$loc.metrics$FreqHomSnp
      hets <- gl@other$loc.metrics$FreqHets
      p <- (2 * p + hets) / 2
      q <- (2 * q + hets) / 2
      H <- 1 - (p ^ 2 + q ^ 2)
      
      Hexp[i] <- mean(H, na.rm = TRUE)
      HexpSD[i] <- sd(H, na.rm = TRUE)
      HexpSE[i] <- std.error(H)
      
      ### CP ### Unbiased He (i.e. corrected for sample size) hard
      # coded for diploid
      uH <- (2 * as.numeric(n_ind[i]) / (2 * as.numeric(n_ind[i]) - 1)) * H
      uHe_df[[i]] <-  uH
      ### CP ###
      uHexp[i] <- mean(uH, na.rm = TRUE)
      uHexpSD[i] <- sd(uH, na.rm = TRUE)
      uHexpSE[i] <- std.error(uH)
      
      Hexp.adj[i] <- Hexp[i] * n_loc[i] / (n_loc[i] + n.invariant)
      Hexp.adjSD[i] <-
        sqrt((sum((H - Hexp.adj[i]) ^ 2, na.rm = TRUE) + n.invariant * Hexp.adj[i] ^ 2) /
               (n_loc[i] + n.invariant - 1))
      Hexp.adjSE[i] <- Hexp.adjSD[i] / sqrt(poly_loc[i]+mono_loc[i])
      
      FIS_temp <- (uH - Ho.loc[[i]]) /  uH 
      Fis_df[[i]] <- FIS_temp
      FIS[i] <- mean(FIS_temp, na.rm = TRUE)
      FISSD[i] <- sd(FIS_temp, na.rm = TRUE)
      FISSE[i] <- std.error(FIS_temp)
      
    }
    
    # bootstrapping
    # creating matrices to store CI
    npops <- nPop(x)
    res_CI <- replicate(npops,
                        as.data.frame(matrix(nrow = 6, ncol = 2)),
                        simplify = FALSE)
    
    if (npops > 1 & nboots > 0) {
      pop_boot <- lapply(sgl, function(y) {
        df <- as.data.frame(as.matrix(y))
        
        res_boots <- boot::boot(
          data = df,
          statistic = pop.het,
          n.invariant = n.invariant,
          R = nboots,
          parallel = parallel,
          ncpus = ncpus
        )
        return(res_boots)
      })
      
      # confidence intervals
      pop_res <- rbind(Ho, Ho.adj, Hexp, uHexp,Hexp.adj, FIS)
      
      for (pop_n in 1:length(sgl)) {
        for (stat_n in 1:6) {
          res_CI_tmp <- boot::boot.ci(
            boot.out = pop_boot[[pop_n]],
            conf = conf,
            type = CI.type,
            index = stat_n,
            t0 =  pop_res[stat_n, pop_n],
            t = pop_boot[[pop_n]]$t[, stat_n]
          )
          
          res_CI[[pop_n]][stat_n,] <-
            tail(as.vector(res_CI_tmp[[4]]), 2)
          
        }
      }
      
      
    }
    
    if (npops == 1 & nboots > 0) {
      # observed values
      pop_res <- rbind(Ho, Ho.adj, Hexp, uHexp,Hexp.adj, FIS)
      res_CI <- as.data.frame(matrix(nrow = 6, ncol = 2))
      df <- as.data.frame(as.matrix(sgl[[1]]))
      
      # bootstrapping
      pop_boot <- boot::boot(
        data = df,
        statistic = pop.het,
        R = nboots,
        n.invariant = n.invariant,
        parallel = parallel,
        ncpus = ncpus
      )
      
      # confidence intervals
      for (stat_n in 1:6) {
        res_CI_tmp <-
          boot::boot.ci(
            boot.out = pop_boot,
            conf = conf,
            type = CI.type,
            index = stat_n,
            t0 =  pop_res[stat_n],
            t = pop_boot$t[, stat_n]
          )
        
        res_CI[stat_n,] <- tail(as.vector(res_CI_tmp[[4]]), 2)
        
      }
    }
    
    if (npops > 2 & nboots > 0) {
      names(res_CI) <- popNames(x)
      
      CI <- lapply(res_CI, function(y) {
        y <-  c(y[,1],y[,2])
        return(y)
      })
      
      stat_list <- as.list(1:12)
      for(stat_i in 1:12){
        stat_list[[stat_i]] <- unlist(lapply(CI,"[",stat_i))
      }
      
    }
    
    if (npops == 1 & nboots > 0) {
      CI <- c(res_CI[,1],res_CI[,2])
      
      stat_list <- as.list(1:12)
      for(stat_i in 1:12){
        stat_list[[stat_i]] <- CI[stat_i]
      }
      
    }
    
    ### CP ###
    if(nboots > 0){
      
      df <-
        data.frame(
          pop = popNames(x),
          n.Ind = round(n_ind, 6),
          n.Loc = n_loc,
          n.Loc.adj = n_loc / (n_loc + n.invariant),
          polyLoc = poly_loc ,
          monoLoc = mono_loc ,
          all_NALoc = all_na_loc,
          
          Ho = round(as.numeric(Ho),6),
          HoSD = round(HoSD,6),
          HoSE = round(HoSE, 6),
          HoLCI = round(stat_list[[1]], 6),
          HoHCI = round(stat_list[[7]], 6),
          
          Ho.adj = round(as.numeric(Ho.adj),6),
          Ho.adjSD = round(Ho.adjSD,6),
          Ho.adjSE = round(Ho.adjSE, 6),
          Ho.adjLCI = round(stat_list[[2]], 6),
          Ho.adjHCI = round(stat_list[[8]], 6),
          
          He = round(Hexp, 6),
          HeSD = round(HexpSD, 6),
          HeSE = round(HexpSE,6),
          HeLCI = round(stat_list[[3]], 6),
          HeHCI = round(stat_list[[9]], 6),
          
          uHe = round(uHexp, 6),
          uHeSD = round(uHexpSD, 6),
          uHeSE = round(uHexpSE ,6),
          uHeLCI = round(stat_list[[4]], 6),
          uHeHCI = round(stat_list[[10]], 6),
          
          He.adj = round(Hexp.adj, 6),
          He.adjSD = round(Hexp.adjSD, 6),
          He.adjSE = round(Hexp.adjSE, 6),
          He.adjLCI = round(stat_list[[5]], 6),
          He.adjHCI = round(stat_list[[11]], 6),
          
          FIS = round(FIS,6),
          FISSD = round(FISSD,6),
          FISSE = round(FISSE , 6),
          FISLCI = round(stat_list[[6]], 6),
          FISHCI = round(stat_list[[12]], 6)
        )
    }else{
      
      df <-
        data.frame(
          pop = popNames(x),
          n.Ind = round(n_ind, 6),
          n.Loc = n_loc,
          n.Loc.adj = n_loc / (n_loc + n.invariant),
          polyLoc = poly_loc ,
          monoLoc = mono_loc ,
          all_NALoc = all_na_loc,
          
          Ho = round(as.numeric(Ho),6),
          HoSD = round(HoSD,6),
          HoSE = round(HoSE, 6),
          HoLCI = NA,
          HoHCI = NA,
          
          Ho.adj = round(as.numeric(Ho.adj),6),
          Ho.adjSD = round(Ho.adjSD,6),
          Ho.adjSE = round(Ho.adjSE, 6),
          Ho.adjLCI = NA,
          Ho.adjHCI = NA,
          
          He = round(Hexp, 6),
          HeSD = round(HexpSD, 6),
          HeSE = round(HexpSE,6),
          HeLCI = NA,
          HeHCI = NA,
          
          uHe = round(uHexp, 6),
          uHeSD = round(uHexpSD, 6),
          uHeSE = round(uHexpSE ,6),
          uHeLCI = NA,
          uHeHCI = NA,
          
          He.adj = round(Hexp.adj, 6),
          He.adjSD = round(Hexp.adjSD, 6),
          He.adjSE = round(Hexp.adjSE, 6),
          He.adjLCI = NA,
          He.adjHCI = NA,
          
          FIS = round(FIS,6),
          FISSD = round(FISSD,6),
          FISSE = round(FISSE , 6),
          FISLCI = NA,
          FISHCI = NA
        )
    }
    ##########
    
    if (plot.display) {
      
      error_L <- error_H <- value <- color <- variable <- He.adj <- NULL
      
      # printing plots and reports assigning colors to populations
      if (is(plot.colors.pop, "function")) {
        colors_pops <- plot.colors.pop(length(levels(pop(x))))
      }
      
      if (!is(plot.colors.pop, "function")) {
        colors_pops <- plot.colors.pop
      }
      
      if (n.invariant == 0) {
        
        pop_list_plot <- df
        pop_list_plot$pop <- as.factor(pop_list_plot$pop)
        pop_list_plot$color <- colors_pops
        
        pop_list_plot_stat <- pop_list_plot[,c("Ho", "uHe", "FIS", "n.Ind",  "pop",  "color")]
        pop_list_plot_stat <- reshape2::melt(pop_list_plot_stat, id = c("pop", "color", "n.Ind"))
        
        if(error.bar=="SD"){
          pop_list_plot_error <- pop_list_plot[,c("HoSD", "uHeSD", "FISSD","pop")]
          pop_list_plot_error <- reshape2::melt(pop_list_plot_error,id = c("pop"))
          colnames(pop_list_plot_error) <- c("pop","variable","error")
          pop_list_plot_error <- pop_list_plot_error[,c("pop","error")]
          pop_list_plot_stat <- cbind(pop_list_plot_stat,error=pop_list_plot_error$error)
        }
        
        if(error.bar=="SE"){
          pop_list_plot_error <- pop_list_plot[,c("HoSE","uHeSE","FISSE","pop")]
          pop_list_plot_error <- reshape2::melt(pop_list_plot_error,id = c("pop"))
          colnames(pop_list_plot_error) <- c("pop","variable","error")
          pop_list_plot_error <- pop_list_plot_error[,c("pop","error")]
          pop_list_plot_stat <- cbind(pop_list_plot_stat,error=pop_list_plot_error$error)
        }
        
        if(error.bar=="CI"){
          pop_list_plot_error_L <- pop_list_plot[,c("HoLCI","uHeLCI","FISLCI","pop")]
          pop_list_plot_error_H <- pop_list_plot[,c("HoHCI","uHeHCI","FISHCI","pop")]
          pop_list_plot_error_L <- reshape2::melt(pop_list_plot_error_L,id = c("pop"))
          pop_list_plot_error_H <- reshape2::melt(pop_list_plot_error_H,id = c("pop"))
          colnames(pop_list_plot_error_L) <- c("pop","variable_L","error_L")
          colnames(pop_list_plot_error_H) <- c("pop","variable_H","error_H")
          pop_list_plot_error <- cbind(pop_list_plot_error_L,pop_list_plot_error_H)
          pop_list_plot_stat <- cbind(pop_list_plot_stat,
                                      error_L = pop_list_plot_error$error_L,
                                      error_H = pop_list_plot_error$error_H)
          
        }
        
        p3 <-
          ggplot(data = pop_list_plot_stat, aes(x = pop, 
                                                y = value,
                                                fill = pop)) +
          geom_bar(stat = "identity", 
                   color = "black", 
                   position = position_dodge())+ 
          facet_wrap(~variable, nrow=1,scales = "free") +
          scale_fill_manual(values = pop_list_plot_stat$color) +
          scale_x_discrete(labels = paste(pop_list_plot_stat$pop,
                                          round(pop_list_plot_stat$n.Ind,
                                                0),
                                          sep = " | ")) +
          plot.theme +
          theme(
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(
              angle = 90,
              hjust = 1,
              face = "bold",
              size = 12
            ),
            axis.title.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "none"
          ) 
        
        if(error.bar=="SD"){
          p3 <- p3 + 
            geom_errorbar(aes(ymin = value, 
                              ymax = value + error), 
                          width=0.5)+
            ggtitle(label = "Heterozygosities and FIS by Population",
                    subtitle = "Error bars show Standard Deviation")
        }
        
        if(error.bar=="SE"){
          p3 <- p3 + 
            geom_errorbar(aes(ymin = value - error, 
                              ymax = value + error), 
                          width=0.5) +
            ggtitle(label = "Heterozygosities and FIS by Population",
                    subtitle = "Error bars show Standard Error")
        }
        
        if(error.bar=="CI"){
          p3 <- p3 + 
            geom_errorbar(aes(ymin = error_L, 
                              ymax = error_H), 
                          width=0.5)+
            ggtitle(label = "Heterozygosities and FIS by Population",
                    subtitle = "Error bars show Confident Intervals")
        }
        
      } else {
        
        df.ordered <- df
        df.ordered$color <- colors_pops
        df.ordered <- df.ordered[order(df.ordered$Ho.adj), ]
        df.ordered$pop <- factor(df.ordered$pop, levels = df.ordered$pop)
        p1 <-
          ggplot(df.ordered, aes(
            x = pop,
            y = Ho.adj,
            fill = pop
          )) + geom_bar(position = "dodge",
                        stat = "identity",
                        color = "black") +
          scale_fill_manual(values = df.ordered$color) +
          scale_x_discrete(labels = paste(df.ordered$pop,
                                          round(df.ordered$n.Ind, 0),
                                          sep = " | ")) + plot.theme + theme(
                                            axis.ticks.x = element_blank(),
                                            axis.text.x = element_blank(),
                                            axis.title.x = element_blank(),
                                            axis.ticks.y = element_blank(),
                                            axis.title.y = element_blank(),
                                            legend.position = "none"
                                          ) +
          labs(fill = "Population") +
          ggtitle("Adjusted Observed Heterozygosity by Population")
        
        p2 <-
          ggplot(df.ordered, aes(
            x = pop,
            y = He.adj,
            fill = pop
          )) + geom_bar(position = "dodge",
                        stat = "identity",
                        color = "black") +
          scale_fill_manual(values = df.ordered$color) +
          scale_x_discrete(labels = paste(df.ordered$pop,
                                          round(df.ordered$n.Ind, 0),
                                          sep = " | ")) + plot.theme + theme(
                                            axis.ticks.x = element_blank(),
                                            axis.text.x = element_text(
                                              angle = 90,
                                              hjust = 1,
                                              face = "bold",
                                              size = 12
                                            ),
                                            axis.title.x = element_blank(),
                                            axis.ticks.y = element_blank(),
                                            axis.title.y = element_blank(),
                                            legend.position = "none"
                                          ) +
          labs(fill = "Population") +
          ggtitle("Adjusted Expected Heterozygosity by Population")
        
        p3 <- (p1 / p2)
      }
    }
    
    # OUTPUT REPORT
    if (verbose >= 3) {
      cat("  Reporting Heterozygosity by Population\n")
      cat("\n  No. of loci =", nLoc(x), "\n")
      cat("  No. of individuals =", nInd(x), "\n")
      cat("  No. of populations =", nPop(x), "\n")
      cat("    Minimum Observed Heterozygosity: ", round(min(df$Ho, na.rm = TRUE), 6))
      if (n.invariant > 0) {
        cat("   [Corrected:", round(min(df$Ho.adj, na.rm = TRUE), 6), "]\n")
      } else {
        cat("\n")
      }
      cat("    Maximum Observed Heterozygosity: ", round(max(df$Ho, na.rm = TRUE), 6))
      if (n.invariant > 0) {
        cat("   [Corrected:", round(max(df$Ho.adj, na.rm = TRUE), 6), "]\n")
      } else {
        cat("\n")
      }
      cat("    Average Observed Heterozygosity: ",
          round(mean(df$Ho, na.rm = TRUE), 6))
      if (n.invariant > 0) {
        cat("   [Corrected:", round(mean(df$Ho.adj, na.rm = TRUE), 6), "]\n\n")
      } else {
        cat("\n\n")
      }
      cat("    Minimum Unbiased Expected Heterozygosity: ",
          round(min(df$uHe, na.rm = TRUE), 6))
      if (n.invariant > 0) {
        cat("   [Corrected:", round(min(df$He.adj, na.rm = TRUE), 6), "]\n")
      } else {
        cat("\n")
      }
      cat("    Maximum Unbiased Expected Heterozygosity: ",
          round(max(df$uHe, na.rm = TRUE), 6))
      if (n.invariant > 0) {
        cat("   [Corrected:", round(max(df$He.adj, na.rm = TRUE), 6), "]\n")
      } else {
        cat("\n")
      }
      cat("    Average Unbiased Expected Heterozygosity: ",
          round(mean(df$uHe, na.rm = TRUE), 6))
      if (n.invariant > 0) {
        cat("   [Corrected:", round(mean(df$He.adj, na.rm = TRUE), 6), "]\n\n")
      } else {
        cat("\n\n")
      }
      
      if (n.invariant > 0) {
        cat(
          "  Average correction factor for invariant loci =",
          mean(n_loc / (n_loc + n.invariant), na.rm = TRUE),
          "\n"
        )
      } else {
        cat("  Heterozygosity estimates not corrected for uncalled invariant loci\n")
      }
    }
    
    # PRINTING OUTPUTS
    if (plot.display) {
      suppressWarnings(print(p3))
    }
    if (verbose >= 2) {
      # if (n.invariant > 0) {
      print(df)
      # } else {
      #   print(df[, c(
      #     "pop",
      #     "n.Ind",
      #     "n.Loc",
      #     "polyLoc",
      #     "monoLoc",
      #     "all_NALoc",
      #     "Ho",
      #     "HoSD",
      #     "He",
      #     "HeSD",
      #     "uHe",
      #     "uHeSD",
      #     "FIS",
      #     "FISSD"
      #   )], row.names = FALSE)
      # }
    }
  }
  
  ########### FOR METHOD BASED ON INDIVIDUAL
  
  if (method == "ind") {
    if (verbose >= 2) {
      cat(report("  Calculating observed heterozygosity for individuals\n"))
      cat(report(
        "  Note: No adjustment for invariant loci (n.invariant set to 0)\n"
      ))
    }
    # Convert to matrix
    m <- as.matrix(x)
    
    # For each individual determine counts of hets, homs and NAs
    c.na <- array(NA, nInd(x))
    c.hets <- array(NA, nInd(x))
    c.hom0 <- array(NA, nInd(x))
    c.hom2 <- array(NA, nInd(x))
    c.nloc <- array(NA, nInd(x))
    
    for (i in 1:nInd(x)) {
      c.na[i] <- sum(is.na(m[i, ]))
      c.hets[i] <-
        sum(m[i, ] == 1, na.rm = TRUE) / (nLoc(x) - c.na[i])
      c.hom0[i] <-
        sum(m[i, ] == 0, na.rm = TRUE) / (nLoc(x) - c.na[i])
      c.hom2[i] <-
        sum(m[i, ] == 2, na.rm = TRUE) / (nLoc(x) - c.na[i])
      c.nloc[i] <- (nLoc(x) - c.na[i])
    }
    
    # Join the sample sizes with the heterozygosities
    df <-
      cbind.data.frame(x@ind.names, c.hets, c.hom0, c.hom2,c.nloc)
    names(df) <-
      c("ind.name", "Ho", "f.hom.ref", "f.hom.alt","n.Loc")
    
    # Boxplot
    if (plot.display) {
      upper <- ceiling(max(df$Ho) * 10) / 10
      p1 <-
        ggplot(df, aes(y = Ho)) +
        geom_boxplot(color = plot.colors.ind[1], fill = plot.colors.ind[2]) +
        coord_flip() +
        plot.theme +
        xlim(range = c(-1, 1)) +
        ylim(0, upper) +
        ylab(" ") +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
        ggtitle("Observed Heterozygosity by Individual")
      
      # Histogram
      p2 <-
        ggplot(df, aes(x = Ho)) +
        geom_histogram(bins = 25,
                       color = plot.colors.ind[1],
                       fill = plot.colors.ind[2]) +
        coord_cartesian(xlim = c(0, upper)) +
        xlab("Observed heterozygosity") +
        ylab("Count") +
        plot.theme
    }
    
    outliers_temp <-
      ggplot_build(p1)$data[[1]]$outliers[[1]]
    outliers <-
      data.frame(ID = as.character(df$ind.name[df$Ho %in% outliers_temp]),
                 Ho = outliers_temp)
    
    # OUTPUT REPORT
    if (verbose >= 3) {
      cat("Reporting Heterozygosity by Individual\n")
      cat("No. of loci =", nLoc(x), "\n")
      cat("No. of individuals =", nInd(x), "\n")
      cat("  Minimum Observed Heterozygosity: ",
          round(min(df$Ho), 6),
          "\n")
      cat("  Maximum Observed Heterozygosity: ",
          round(max(df$Ho), 6),
          "\n")
      cat("  Average Observed Heterozygosity: ",
          round(mean(df$Ho), 6),
          "\n\n")
      cat("  Results returned as a dataframe\n\n")
      if (nrow(outliers) == 0) {
        cat("  No outliers detected\n\n")
      } else {
        cat("  Outliers detected\n")
        print(outliers)
        cat("\n")
      }
    }
    
    # PRINTING OUTPUTS
    if (plot.display) {
      p3 <- (p1 / p2) + plot_layout(heights = c(1, 4))
      print(p3)
    }
    if (verbose >= 2) {
      print(df, row.names = FALSE)
    }
  }
  
  # SAVE INTERMEDIATES TO TEMPDIR
  if (save2tmp) {
    # creating temp file names
    if (plot.display) {
      temp_plot <- tempfile(pattern = "Plot_")
    }
    temp_table <- tempfile(pattern = "Table_")
    match_call <-
      paste0(names(match.call()),
             "_",
             as.character(match.call()),
             collapse = "_")
    # saving to tempdir
    if (plot.display) {
      saveRDS(list(match_call, p3), file = temp_plot)
      if (verbose >= 2) {
        cat(report("  Saving the ggplot to session tempfile\n"))
      }
    }
    saveRDS(list(match_call, df), file = temp_table)
    if (verbose >= 2) {
      cat(report("  Saving tabulation to session tempfile\n"))
      cat(
        report(
          "  NOTE: Retrieve output files from tempdir using
                    gl.list.reports() and gl.print.reports()\n"
        )
      )
    }
  }
  
  if (verbose >= 3) {
    cat(report("  Returning a dataframe with heterozygosity values\n"))
  }
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # RETURN
  
  invisible(df)
}
