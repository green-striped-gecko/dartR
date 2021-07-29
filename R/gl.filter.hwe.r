#' @name gl.filter.hwe
#' @title Filters loci that show significant departure from Hardy-Weinberg Equilibrium
#' @description 
#' Calculates the probabilities of agreement with H-W equilibrium based on observed
#' frequencies of reference homozygotes, heterozygotes and alternate homozygotes. 
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param subset Either a vector with population names, "each", "all" (see details)
#' [default "each"].
#' @param method_sig Method for determining statistical significance: "ChiSquare" 
#' or "Exact" [default "Exact"].
#' @param multi_comp Whether to adjust p-values for multiple comparisons [default TRUE].
#' @param multi_comp_method Method to adjust p-values for multiple comparisons: 
#' "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"
#' (see details) [default "bonferroni"]
#' @param alpha_val Level of significance for testing [default 0.05].
#' @param sample_size Minimum sample size per population required to test for 
#' Hardy-Weinberg proportions [default 5]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, 
#' progress log ; 3, progress and results summary; 5, full report 
#' [default NULL, unless specified using gl.set.verbosity]
#'
#' @details 
#' This function filters out loci departing from Hardy-Weinberg proportions in three ways:
#' \itemize{
#' \item across all populations together using: subset="all"
#' \item within each population separately using: subset="each"
#' \item within selected populations using for example: subset=c("pop1","pop2")
#' }
#' 
#' Loci are filtered out if they show HWE departure in any one population.
#' 
#' Tests for HWE are only valid if there is no population substructure (assuming
#' random mating), and the tests have sufficient power only when there is 
#' sufficient sample size (say, n individuals > 20). Note also that correction 
#' for multiple comparisons is probably required if you wish to place particular
#' importance on one or a few significant departures.
#' 
#' 
#' \code{\link[stats]{p.adjust}}
#' 
#' p.adjust.methods
#  "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#   "fdr", "none"
#
#'  \code{\link[HardyWeinberg]{HWChisqStats}}
#' Uses the exact calculations contained in the function 
#' \code{\link[HardyWeinberg]{HWExactStats}} as developed by Wigginton 
#' et al. (2005).
#' 
#' The C++ code was generously shared by Christopher Chang, and the same code is 
#' used in the program PLINK (2.0).
#' 
#' @return A genlight object with the loci departing significantly from HWE removed
#' @author Arthur Georges -- Post to \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' result <- gl.filter.hwe(x = bandicoot.gl)
#' @references 
#' \itemize{
#' \item Wigginton, J.E., Cutler, D.J., & Abecasis, G.R. (2005). A Note on Exact 
#' Tests of Hardy-Weinberg Equilibrium. American Journal of Human Genetics 76:887-893.
#' \item Graffelman, J. & Morales-Camarena, J. (2008). Graphical tests for 
#' Hardy-Weinberg equilibrium based on the ternary plot. Human Heredity 65:77-84.
#' \item Graffelman, J. (2015). Exploring Diallelic Genetic Markers: The Hardy 
#' Weinberg Package. Journal of Statistical Software 64:1-23.
#' }
#'
#' @seealso \code{\link{gl.filter.hwe}}
#' @family filters/filter reports
#'
#' @importFrom graphics polygon
#' @importFrom stats qchisq
#' @export
#'  

gl.filter.hwe <- function(x, 
                          subset = "each", 
                          multi_comp = TRUE, 
                          multi_comp_method = "bonferroni",
                          method_sig = "Exact", 
                          alpha_val = 0.05,
                          sample_size = 5,
                          verbose = NULL){
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func=funname,build="Jackson",v=verbose)
  
  # CHECK DATATYPE 
  datatype <- utils.check.datatype(x)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  # check if package is installed
  pkg <- "HardyWeinberg"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop(error("Package",pkg," needed for this function to work. Please install it."))
  }
  
  if (datatype == "SilicoDArT"){
    cat(error("  Detected Presence/Absence (SilicoDArT) data\n"))
    stop(error("Cannot calculate HWE from fragment presence/absence data. Please provide a SNP dataset.\n"))
  }
  
  if (alpha_val < 0 | alpha_val > 1){
    cat(warn("    Warning: level of significance per locus alpha must be an integer between 0 and 1, set to 0.05\n"))
    alpha_val <- 0.05
  }
  
  # DO THE JOB
  
  hold <- x
  
  #### Interpret options
  ########### for subset all
  if (subset[1] == "all") {
    if(verbose >= 2){
      cat(report("  Pooling all populations for HWE calculations\n"))
    }
    if(verbose >= 3){
      cat(warn("  Warning: Significance of tests may indicate heterogeneity among populations\n\n"))
    }
    # assigning the same population to all individuals
    pop(x) <- array("pop",dim = nInd(x))
    pop(x) <- as.factor(pop(x))
  }
  
  ########### for subset each
  if (subset[1] == "each") {
    if(verbose >= 2){
      cat(report("  Analysing each population separately\n"))
    }
  }
  
  ########### for subset selected populations
  if(subset[1] != "each" & subset[1] != "all"){
    # check whether the populations exist in the dataset
    pops_hwe_temp <- pop(x) %in% subset
    pops_hwe <- sum(pops_hwe_temp[pops_hwe_temp==TRUE])
    # if the populations are not in the dataset
    if(pops_hwe == 0){
      stop(error("Fatal Error: subset parameter must be \"each\", \"all\", or a list of populations existing in the dataset\n"))
    }
    # subsetting the populations
    x <- x[pop(x) %in% subset]
    # assigning the same population to all individuals
    pop(x) <- array("pop",dim = nInd(x))
    pop(x) <- as.factor(pop(x))
    if(verbose >= 2){
      cat(report(paste("  Pooling populations", paste(subset,collapse = " "),"together for HWE calculations\n")))
    }
    if(verbose >= 3){
      cat(warn("  Warning: Significance of tests may indicate heterogeneity among populations\n\n"))
    }
  }
  
  poplist_temp <- seppop(x)
  # filtering monomorphs
  poplist <- lapply(poplist_temp,gl.filter.monomorphs,verbose=0)
  
  # testing whether populations have heteromorphic loci 
  monomorphic_pops_temp <- unlist(lapply(poplist, nLoc))
  monomorphic_pops <- monomorphic_pops_temp[which(monomorphic_pops_temp==0)]
  
  if(length(monomorphic_pops)>0){
    if(verbose >= 2){
      cat(warn(" Warning: No heteromorphic loci in population",names(monomorphic_pops),"... skipped\n"))
      # removing pops that do not have heteromorphic loci 
      pops_to_remove <- which(names(poplist) %in% names(monomorphic_pops))
      poplist <- poplist[-pops_to_remove]
    }
  }
  
  # testing whether populations have small sample size 
  n_ind_pops_temp <- unlist(lapply(poplist, nInd))
  n_ind_pops <- n_ind_pops_temp[which(n_ind_pops_temp<=sample_size)]
  
  if(length(n_ind_pops)>0){
    if(verbose >= 2){
      cat(warn(" Warning: population",names(n_ind_pops),"has less than",sample_size,"individuals... skipped\n"))
      # removing pops that have low sample size 
      pops_to_remove_2 <- which(names(poplist) %in% names(n_ind_pops))
      poplist <- poplist[-pops_to_remove_2]
    }
  }
  
  if(length(poplist)<1){
    stop(error("No populations left after removing populations with low sample size and populations with monomorphic loci"))
  }
  
  result <- as.data.frame(matrix(nrow = 1, ncol = 10))
  colnames(result) <- c("Population","Locus", "Hom_1", "Het", "Hom_2", "N", "Prob", "Sig", "Prob.adj", "Sig.adj")
  
  for (i in poplist) {
    
    mat_HWE_temp <- t(as.matrix(i))
    mat_HWE <- matrix(nrow = nLoc(i),ncol = 3)
    colnames(mat_HWE) <- c("AA", "AB", "BB")
    mat_HWE[,"AA"] <- apply(mat_HWE_temp,1,function(y){length(y[which(y==0)])})
    mat_HWE[,"AB"] <- apply(mat_HWE_temp,1,function(y){length(y[which(y==1)])})
    mat_HWE[,"BB"] <- apply(mat_HWE_temp,1,function(y){length(y[which(y==2)])})
    
    if(method_sig == "ChiSquare"){
      p.values <- HardyWeinberg::HWChisqStats(mat_HWE)
    }
    if(method_sig == "Exact"){
      p.values <- HardyWeinberg::HWExactStats(mat_HWE)
    }
    
    total <- rowSums(mat_HWE,na.rm = T)
    
    sig2 <- rep(NA,length(p.values))
    p.values_adj <- rep(NA,length(p.values))
    bonsig2 <- rep(NA,length(p.values))
    
    # Assemble results into a dataframe
    result_temp <- cbind.data.frame(popNames(i),locNames(i),mat_HWE, total, p.values, sig2, p.values_adj , bonsig2, stringsAsFactors = FALSE)
    names(result_temp) <- c("Population","Locus", "Hom_1", "Het", "Hom_2", "N", "Prob", "Sig", "Prob.adj", "Sig.adj")
    
    result <- rbind.data.frame(result,result_temp, stringsAsFactors = FALSE)
  }
  result <- result[-1,]
  
  if(multi_comp == TRUE){
    result$Prob.adj <- stats::p.adjust(result$Prob, method = multi_comp_method)
  }
  
  if(multi_comp==F){
    result <- result[which(result$Prob <= alpha_val),]
  }
  if(multi_comp==T){
    result <- result[which(result$Prob.adj <= alpha_val),]
  }
  
  failed.loci <- as.character(unique(result$Locus))
  
  if (verbose >= 2){
    cat("  Loci examined:", nLoc(hold),"\n")
  }
  
  index <- !locNames(hold) %in% failed.loci
  hold <- hold[,index]
  hold@other$loc.metrics <- hold@other$loc.metrics[index,]

  #### Report the results
  if (verbose >= 2){
      if(multi_comp==T){
        cat("  Deleted",length(failed.loci),"loci with significant departure from HWE, after correction for multiple tests using the",multi_comp_method, "method at experiment-wide alpha =",alpha_val,"\n")
      }
      if(multi_comp==F){
        cat("  Deleted",length(failed.loci),"loci with significant departure from HWE at alpha =",alpha_val,"applied locus by locus\n")
      }
      cat("  Loci retained:",nLoc(hold),"\n\n")
      cat(important("    Adjustment of p-values for multiple comparisons vary with sample size\n"))
  }

  # ADD TO HISTORY
  nh <- length(hold@other$history)
  hold@other$history[[nh + 1]] <- match.call()
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("\nCompleted:", funname, "\n\n"))
  }
  
  # RETURN
  invisible(hold)
  
}


