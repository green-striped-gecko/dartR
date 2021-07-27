#' @name gl.report.hwe
#' @title Reports departure from Hardy-Weinberg Equilibrium
#' @description 
#' Calculates the probabilities of agreement with H-W equilibrium based on observed
#' frequencies of reference homozygotes, heterozygotes and alternate homozygotes. 
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param subset Either a vector with population names, "each", "all" (see details)
#' [default "each"].
#' @param plot.out If TRUE, will produce a Ternary Plot(s) [default FALSE].
#' @param method_sig Method for determining statistical significance: "ChiSquare" 
#' or "Exact" [default "Exact"].
#' @param multi_comp Whether to adjust p-values for multiple comparisons [default TRUE].
#' @param multi_comp_method Method to adjust p-values for multiple comparisons: 
#' "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"
#' (see details) [default "bonferroni"]]
#' @param alpha_val Level of significance for testing [default 0.05].
#' @param sample_size = 5
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, 
#' progress log ; 3, progress and results summary; 5, full report 
#' [default NULL, unless specified using gl.set.verbosity]
#'
#' @details 
#' Factors that can cause deviations from Hardy-Weinberg proportions include: 
#' mutation, finite population size, selection, population structure, age 
#' structure, assortative mating, sex linkage, nonrandom sampling and genotyping
#'  errors. Testing for Hardy-Weinberg proportions should be a process that involves
#'  a careful evaluation of the results.  
#'  
#' This function applies HWE tests to each locus in three ways:
#' \itemize{
#' \item across all populations together using: subset="all"
#' \item within each population separately using: subset="each"
#' \item within selected populations using for example: subset=c("pop1","pop2")
#' }
#' 
#' This function tests for deviations from Hardy Weinberg proportions using two 
#' different methods:
#' \itemize{
#' \item "ChiSquare"   \code{\link[HardyWeinberg]{HWChisqStats}}
#' \item "Exact"  Uses the exact calculations contained in the function 
#' \code{\link[HardyWeinberg]{HWExactStats}} as developed by Wigginton 
#' et al. (2005). The C++ code was generously shared by Christopher Chang, and the same code is 
#' used in the program PLINK (2.0).
#' }
#' 
#' This function uses correction for multiple tests using the following methods
#' \code{\link[stats]{p.adjust}}:
#' \itemize{
#' \item "holm", 
#' \item "hochberg", 
#' \item "hommel", 
#' \item "bonferroni", 
#' \item "BH", 
#' \item "BY", 
#' \item "fdr"
#' }
#'  
#' Tests for HWE are only valid if there is no population substructure (assuming
#' random mating), and the tests have sufficient power only when there is 
#' sufficient sample size (say, n individuals > 20). Note also that correction 
#' for multiple comparisons is probably required if you wish to place particular
#' importance on one or a few significant departures.
#' 
#' A Ternary Plot is optionally produced (plot.out = TRUE) -- see Graffelman et al.(2008) for 
#' further details. Implementation of the Ternary Plot is via package {HardyWeinberg} 
#' (Graffelman (2015). The plot labels loci that depart significantly from HWE 
#' as red, and those not showing significant departure as green. 
#       
#' HWTernaryPlot is a routine that draws a ternary plot for three-way genotypic compositions (AA,AB,BB), and represents the acceptance region for different tests for Hardy-Weinberg equilibrium (HWE) in the plot.
#' 
#' say the monomorphic loci are filtered
#' populations with less than 5 loci cannot be assessed tests for this in all cases
#' 
#' p.adjust {stats}
#' 
#' p.adjust.methods
# c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#   "fdr", "none")
#' colour the marker points automatically according to the result of a 
#' significance test (green markers non-significant, red markers significant).
#' #' Uses the exact calculations contained in function \code{\link{utils.prob.hwe}}
#' as developed by Wigginton et al. (2005).
#' @return returns a dataframe containing loci, counts of reference SNP homozygotes,
#' heterozygotes and alternate SNP homozygotes; probability of departure from 
#' H-W equilibrium, and per locus significance with and without Bonferroni correction.
#' @author Arthur Georges -- Post to \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' gl.report.hwe(x = bandicoot.gl,plot.out = FALSE)
#' @references 
#' \itemize{
#' \item Wigginton, J.E., Cutler, D.J., & Abecasis, G.R. (2005). A Note on Exact 
#' Tests of Hardy-Weinberg Equilibrium. American Journal of Human Genetics 76:887-893.
#' \item Graffelman, J. & Morales-Camarena, J. (2008). Graphical tests for 
#' Hardy-Weinberg equilibrium based on the ternary plot. Human Heredity 65:77-84.
#' \item Graffelman, J. (2015). Exploring Diallelic Genetic Markers: The Hardy 
#' Weinberg Package. Journal of Statistical Software 64:1-23.
#' \item Waples, R. S. (2015). Testing for Hardy–Weinberg proportions: have we 
#' lost the plot?. Journal of heredity, 106(1), 1-19.
#' }
#'
#' @seealso \code{\link{gl.filter.hwe}}
#' @family filters/filter reports
#'
#' @importFrom graphics polygon
#' @importFrom stats qchisq
#' @export
#'  

gl.report.hwe <- function(x, 
                          subset = "each", 
                          plot.out = FALSE,
                          multi_comp = TRUE, 
                          multi_comp_method = "bonferroni",
                          method_sig = "Exact", 
                          alpha_val = 0.05,
                          sample_size = 5,
                          verbose = NULL) {
  
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
      
      if(multi_comp == TRUE){
        p.values_adj <- stats::p.adjust(p.values, method = multi_comp_method)
      }
      
      # Assemble results into a dataframe
      result_temp <- cbind.data.frame(popNames(i),locNames(i),mat_HWE, total, p.values, sig2, p.values_adj , bonsig2, stringsAsFactors = FALSE)
      names(result_temp) <- c("Population","Locus", "Hom_1", "Het", "Hom_2", "N", "Prob", "Sig", "Prob.adj", "Sig.adj")

      result <- rbind.data.frame(result,result_temp, stringsAsFactors = FALSE)
    }
    result <- result[-1,]
    
    result[which(result$Prob < alpha_val),"Sig"] <- "sig" 
    result[which(result$Prob > alpha_val),"Sig"] <- "no_sig" 
    result[which(result$Prob.adj < alpha_val),"Sig.adj"] <- "sig" 
    result[which(result$Prob.adj > alpha_val),"Sig.adj"] <- "no_sig" 
    
    if(plot.out){
      # Determine the page layout for plots based on the number of populations to plot
      npops2plot <- length(poplist)
      
      if (npops2plot == 2) {
        graphics::layout(matrix(c(1,2), 1, 2, byrow = TRUE),widths=c(0.5,0.5))
        if(verbose >= 2){
          cat(report("  Plotting two ternary plots on one page\n\n"))
        }
      }
      if (npops2plot > 2) {
        graphics::layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
        if(verbose >= 2){
          cat(report("  Plotting up to four ternary plots on one page\n\n"))
        }
      } 
      
      # Plot the tertiary plots
      for (z in poplist) {
        pop_name <- popNames(z)
        xlabel <- paste0("\n\n\nPopulation: ",pop_name,"\n",method_sig," method\nalpha = ",alpha_val,"")
        result_pop <- result[which(result$Population==pop_name),]
        mat_genotypes <- as.matrix(result_pop[,c("Hom_1" ,"Het", "Hom_2" )])
        colnames(mat_genotypes) <- c("AA","AB","BB")
        
        if (method_sig == "Exact"){
          res <- HWTernaryPlot_correction(X=mat_genotypes, region = 7, vertex.cex = 1.25, alpha=alpha_val, axislab=xlabel, multi_comp = multi_comp, multi_comp_method=multi_comp_method)
        } 
        
        if (method_sig == "ChiSquare"){
          res <- HardyWeinberg::HWTernaryPlot(X=mat_genotypes, region = 2, vertex.cex = 1.25, alpha=alpha_val, axislab=xlabel, multi_comp = multi_comp, multi_comp_method=multi_comp_method)
        }
      }
    }
  
  #### Report the results
  if(multi_comp==F){
    result <- result[which(result$Prob <= alpha_val),]
  }
  if(multi_comp==T){
    result <- result[which(result$Prob.adj <= alpha_val),]
  }
  result <- result[order(result$Locus),]
  cat("    Reporting significant departures from Hardy-Weinberg Equilibrium\n")
  if (nrow(result)==0){
    cat("    No significant departures\n")
  } else {
    cat("    NB: Departures significant at the alpha level of",alpha_val,"are listed\n")
      cat(important("    Adjustment of p-values for multiple comparisons vary with sample size\n"))
      print(result, row.names=FALSE)
  }
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("\nCompleted:", funname, "\n\n"))
  }
  
  # RETURN
  invisible(result)
  
}
