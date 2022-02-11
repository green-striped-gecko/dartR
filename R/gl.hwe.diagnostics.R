#' @name gl.hwe.diagnostics
#' @title Provides descriptive stats and plots to diagnose potential problems with H-W proportions
#' @description
#' Different causes may be responsible for lack of H-W proportions. This function helps diagnose potential problems
#' @inheritParams gl.report.hwe
#' @details
#'  This function initially runs \code{gl.report.hwe} and reports the the ternary plots. 
#'  the remaining outputs follow the recommendations from Waples (2015) paper.  
#'  These include:
#'  \itemize{
#'  \item A bar plot with the distribution of thep-values of the HWE tests. 
#'  These should be roughly uniform  across the equal-sized bins.
#'  \item A bar plot with the observed and expected number of significant HWE tests
#'   for the same locus in multiple populations (that is, the x-axis is whether 
#'   a locus results significant in 1, 2, ..., n populations. 0 indicates the 
#'   number of non-significant tests). The observed and expected number of HWE 
#'   tests should be roughly similarly distributed (if they are significant by 
#'   chance alone). 
#'   \item A table where the number of observed and expected significant HWE tests 
#'   are reported by each population, indicating whether these are due to 
#'   heterozygosity excess or deficiency. These can be used to have a clue of 
#'   potential problems (e.g. deficiency= Wahlund effect, null alleles, 
#'   non-random sampling; excess=sex linkage or different selection between sexes, 
#'   demographic changes/small Ne. See Table 1 in Wapples 2015). In the last two 
#'   columns of the table, chisquare value and its associated p-value is reported. The chisquare 
#'   is computed following Fisher's procedure for a global test. This basically 
#'   tests whether there is at least one test that is truly significant in the 
#'   series of tests conducted (De Meeûs et al 2009).
#'
#' @return A list with the plots and table.
#' @author Custodian: Carlo Pacioni -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' gl.hwe.diagnostics(x = platypus.gl)
#' @references
#' \itemize{
#' \item De Meeûs, T., Guégan, J.-F., Teriokhin, A.T., 2009. MultiTest V.1.2, a 
#' program to binomially combine independent tests and performance comparison 
#' with other related methods on proportional data. BMC Bioinformatics 10, 443-443.
#' \item Fisher, R., 1970. Statistical methods for research workers Edinburgh: 
#' Oliver and Boyd.
#' \item Waples, R. S. (2015). Testing for Hardy–Weinberg proportions: have we
#' lost the plot?. Journal of heredity, 106(1), 1-19.
#' }
#' @seealso \code{\link{gl.report.hwe}}
#' @family 
#' @export
#' @import data.table
#' @import ggplot2



gl.hwe.diagnostics <- function(x, alpha_val=0.05, brk = seq(0, 1, 1 / 20)) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jackson",
                   verbosity = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  #### Function helpers ####
  hwe.dist.plot <- function(hweout, brk) {
    
    ggplot(hweout, aes(Prob)) + geom_histogram(breaks=brk) +
      geom_hline(yintercept = nrow(hweout)/(length(brk)-1), col="red")
  }
  
  # DO THE JOB
  
  # abandoned this for now
  # FstvsFis.plot <- function(Fstats) {
  #   ggplot(Fstats, aes(Fst, Fis)) + geom_point() + geom_smooth(method = "lm") #+ facet_grid(~.id)
  # }
  # 
  nLoBynPop_plot <- function(nTimes) {
    ggplot(nTimes, aes(nPop, Freq, col=Data, fill=Data)) + 
      geom_histogram(stat="identity", position = "dodge") + 
      scale_y_log10()
  }
  
  
  # Distribution of pvlues by qual bins
  suppressWarnings(
  hweout <- gl.report.hwe(x, sig_only = F, verbose = 0)
  )
  
  hwe_dist_plot <- hwe.dist.plot(hweout, brk)
  
  # Fst vs Fis scatter plot with linear regression
  # lpops <- seppop(x)
  # lFstats <- lapply(lpops, gl.basic.stats, verbose = 0)
  # lFstats <- lapply(lFstats, "[[", "perloc")
  # Fstats <- rbindlist(l = lFstats, use.names = TRUE, idcol = TRUE)
  Fstats <- gl.basic.stats(x, verbose = 0)
  # FstatsLoc <- Fstats$perloc
  # print(FstvsFis.plot(FstatsLoc))
  
  # Number of loci out of HWE as a function of a population
  hweout.dt <- data.table(hweout)
  nTimesBypop <- hweout.dt[, .N, by=c("Locus", "Sig")]
  setkey(nTimesBypop, Sig)

  nTimesBypop.df <- as.data.frame(table(nTimesBypop["sig", N]))
  
  # Include the non-sig tests
  nTimesBypop.df <- rbind(data.frame(Var1=0, Freq=nTimesBypop["no_sig", sum(N)]), nTimesBypop.df)
  nTimesBypop.df$Data <- "Observed"

  # Generate the null distribution
  nullDist <- as.data.frame(table(rbinom(length(hweout.dt[, unique(Locus)]), 
                                         size = length(hweout.dt[, unique(Population)]), 
                                         prob = alpha_val)))
  nullDist$Data <- "Null expectation"
  
  # Compile the data for the plot
  nTimesBypop.fin <- rbind(nTimesBypop.df, nullDist)
  names(nTimesBypop.fin)[1] <- "nPop"
  nTestByPops_plot <- nLoBynPop_plot(nTimesBypop.fin)
  
  # Collate HWE tests and Fis per locus and pop
  FisPops <- data.table(Fstats$Fis, keep.rownames = TRUE)
  
  # fix the ehadings when there is only one pop
  if(length(levels(pop(x))) == 1) {
    FisPops[, dumpop := NULL]
    setnames(FisPops, "1", levels(pop(x)))
  }
  setnames(FisPops, "rn", "Locus")
  FisPopsLong <- melt(FisPops, id.vars = "Locus", variable.name = "Population", value.name = "Fis")
  FisPopsLong[, Locus := sub("^X", "", Locus)]
  hweout.dt[, Locus := gsub("-|/", replacement = ".", x = Locus)]
  hwe_Fis <- merge(hweout.dt, FisPopsLong, by = c("Locus", "Population"))
  hwe_Fis[, Deficiency := Fis > 0]
  hwe_Fis[,  Excess:= Fis < 0]
  setkey(hwe_Fis, Sig)

  hwe_summary <- hwe_Fis["sig", .(nSig=.N, nExpected = alpha_val * nLoc(x),
                                  Deficiency=sum(Deficiency, na.rm = TRUE), 
                                  Excess=sum(Excess, na.rm = TRUE),
                                  PropDeficiency=sum(Deficiency, na.rm = TRUE)/.N),
                         by=Population]
  
  chsq <- hwe_Fis[, .(ChiSquare=-2*(sum(log(Prob)))), by= Population]
  chsq[, pvalue:=pchisq(ChiSquare, 2*nLoc(x), lower.tail = FALSE)]
  hwe_summary <- merge(hwe_summary, chsq, by = "Population")
  res <- list(hwe_dist_plot, nTestByPops_plot, hwe_summary)
  print(res)
  return(res)
}

