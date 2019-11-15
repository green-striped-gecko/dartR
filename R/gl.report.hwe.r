#' Reports departure from Hardy-Weinberg Equilibrium
#' 
#' Calculates the probabilities of agreement with H-W equilibrium based on observed
#' frequencies of reference homozygotes, heterozygotes and alternate homozygotes. 
#' Uses the exact calculations contained in function utils.prob.hwe() as developed by Wigginton et al. (2005).
#' 
#' Tests are applied to each locus across all populations pooled (subset="all"), to each locus considered within each population treated separately
#' (subset="each") or to each locus within selected populations pooled (subset=c("pop1","pop2")). Tests for HWE are
#' only valid if there is no population substructure (assuming random mating), and the tests have sufficient power
#' only when there is sufficient sample size (say, n individuals > 20). Note also that correction for multiple comparisons is probably required
#' if you wish to place particular importance on one or a few significant departures.
#' 
#' A Ternary Plot is optionally produced -- see Graffelman et al.(2008) for further details. Implementation of the Ternary Plot is via package {HardyWeinberg} 
#' (Graffelman (2015). The plot labels loci that depart significantly from HWE as red, and those not showing significant departure as green. 
#' Two methods are used to determine significance. ChiSquare (with correction) is traditional but involves approximations; Fisher is computationally 
#' more expensive, but applies a Fisher Exact Test of departure from HWE.
#' 
#' @param x -- a genlight object containing the SNP genotypes [Required]
#' @param subset -- either, list populations to combine in the analysis | each | all [Default "each"] 
#' @param plot -- if TRUE,  will produce a Ternary Plot(s) [default FALSE]
#' @param method -- for determining the statistical signicance in the ternary plot: ChiSquare (with continuity correction) | Fisher [default "ChiSquare"]
#' @param bonf -- if TRUE, Bonferroni correction will be applied to the level of significance [default TRUE]
#' @param alpha -- level of significance for testing [default 0.05]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 3]
#' @return returns a dataframe containing loci, counts of reference SNP homozygotes, heterozygotes
#' and alternate SNP homozygotes; probability of departure from H-W equilibrium,
#' and per locus significance with and without Bonferroni Correction.
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @importFrom HardyWeinberg HWTernaryPlot
#' @export
#' @references Wigginton, J.E., Cutler, D.J., & Abecasis, G.R. (2005). A Note on Exact Tests of Hardy-Weinberg Equilibrium. American Journal of Human Genetics 76:887-893.
#' @references Graffelman, J. & Morales-Camarena, J. (2008). Graphical tests for Hardy-Weinberg equilibrium based on the ternary plot. Human Heredity 65:77-84.
#' @references Graffelman, J. (2015). Exploring Diallelic Genetic Markers: The HardyWeinberg Package. Journal of Statistical Software 64:1-23.
#' @examples
#' list <- gl.report.hwe(testset.gl,subset=c("EmmacMaclGeor", "EmmacCoopCully"),plot=TRUE,bonf=FALSE)
#' gl.report.hwe(testset.gl,subset=c("EmmacCoopCully"), plot=TRUE, verbose=3)
#' gl.report.hwe(testset.gl,subset="all", plot=TRUE, bonf=FALSE, verbose=3)
#' gl.report.hwe(testset.gl, subset="each", plot=TRUE, bonf=FALSE)

gl.report.hwe <- function(x, subset="each", plot=FALSE, method="ChiSquare", alpha=0.05, bonf=TRUE, verbose=NULL) {
  
# TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  
# SET VERBOSITY
  
  if (is.null(verbose)){ 
    if(!is.null(x@other$verbose)){ 
      verbose <- x@other$verbose
    } else { 
      verbose <- 2
    }
  } 
  
  if (verbose < 0 | verbose > 5){
    cat(paste("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n"))
    verbose <- 2
  }
  
# FLAG SCRIPT START
  
  if (verbose >= 1){
    if(verbose==5){
      cat("Starting",funname,"[ Build =",build,"]\n")
    } else {
      cat("Starting",funname,"\n")
    }
  }  
  
# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    stop("Fatal Error: genlight object required!\n")
  }
  
  if (all(x@ploidy == 1)){
    cat("  Detected Presence/Absence (SilicoDArT) data\n")
    stop("Cannot calculate HWE from fragment presence/absence data. Please provide a SNP dataset.\n")
  } else if (all(x@ploidy == 2)){
    cat("  Processing a SNP dataset\n")
  } else {
    stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)!\n")
  }

# FUNCTION SPECIFIC ERROR CHECKING

  # To be added

# DO THE JOB
  
  # Initialize a flag to indicate if populations are to be pooled or not
  flag <- 0
  
  #### Interpret options
  
  if (subset[1] == "all") {
    gl <- gl.filter.monomorphs(x,verbose=0)
    if(nPop(gl) > 1) {
      if(verbose >= 2){cat("  Pooling all populations for HWE calculations\n")}
      if(verbose >= 3){cat("  Warning: Significance of tests may indicate heterogeneity among populations\n\n")}
    } else {
      if(verbose >= 2){cat("  Calculating HWE for population",popNames(gl))}
    }  
    flag <- 1
    
  } else if (subset[1] == "each") {
    if (nPop(x) == 1){
      if(verbose >= 2){cat("  Calculating HWE for population",popNames(x))}
      gl <- gl.filter.monomorphs(x,verbose=0)
      flag <- 1
    } else {
      if(verbose >= 2){cat("  Analysing each population separately\n")}
      poplist <- seppop(x)
      flag <- 0
    }  
    
  } else if(nPop(x[pop(x) %in% subset])){
      gl <- x[pop(x) %in% subset]
      gl <- gl.filter.monomorphs(gl,verbose=0)
      if (nPop(gl) == 1){
        if(verbose >= 2){cat("  Calculating HWE for population",popNames(gl),"\n")}
        flag <- 1
      } else {
        if(verbose >= 2){cat("  Pooling populations",subset,"for HWE calculations\n")}
        if(verbose >= 3){cat("  Warning: Significance of tests may indicate heterogeneity among populations\n\n")}
        flag <- 1
      } 

  } else {
      stop("Fatal Error: subset parameter must be \"each\", \"all\", or a list of populations\n")
    }

  #### Single or Pooled populations
  
  if(flag == 1){
    result <- utils.hwe(gl, prob=alpha,verbose=0)
    mat <- array(NA,3*dim(result)[1])
    dim(mat) <- c(dim(result)[1],3)
    mat[,1] <- as.numeric(as.character(result$Hom_1)) # for God knows why
    mat[,2] <- as.numeric(as.character(result$Het))
    mat[,3] <- as.numeric(as.character(result$Hom_2))
    colnames(mat) <- c("AA", "AB", "BB")
    # Plot a single ternary plot
    if (bonf) {
      c.alpha <- alpha/length(mat[,1])
      xlabel <- paste0(method," (alpha = ",signif(c.alpha,4),", Bonferroni corrected)")
    } else {
      c.alpha <- alpha
      xlabel <- paste0(method," (alpha = ",signif(c.alpha,4),")")
    }
    if (method == "Fisher"){
      res <- HardyWeinberg::HWTernaryPlot(mat, 100, region = 7, vertex.cex = 1.25, alpha=c.alpha, signifcolour = TRUE, vbounds=TRUE, axislab=xlabel)
    } else {
      res <- HardyWeinberg::HWTernaryPlot(mat, 100, region = 2, vertex.cex = 1.25, alpha=c.alpha, signifcolour = TRUE, vbounds=TRUE, axislab=xlabel)
    }  
  }  
  
  #### Disaggregated populations
  
  if (flag == 0){

    # Generate the results for printout, count populations to plot
    count <- 0
    npops2plot <- nPop(x)
    for (i in poplist) {
      count <- count + 1
      ii <- gl.filter.monomorphs(i,verbose =0)
      if (nLoc(ii) == 0){
        if(verbose >= 2){cat("  Warning: No heteromorphic loci in population",popNames(i),"... skipped\n")}
        count <- count - 1
        npops2plot <- npops2plot - 1
        next
      }
      if (count == 1) {
        result <- utils.hwe(ii, prob=alpha,verbose=0)
        Population <- rep(names(poplist)[count],nrow(result))
        result <- cbind(Population,result)
      } else {
        r <- utils.hwe(ii, prob=alpha,verbose=0)
        Population <- rep(names(poplist)[count],nrow(r))
        r <- cbind(Population,r)
        result <- rbind(result, r)
      }
    } 
    
    # Determine the page layout for plots based on the number of populations to plot   
    layout(mat = matrix(c(1,1),1,1, byrow=FALSE))
    if (npops2plot > 1){
      if (npops2plot == 2) {
        layout(matrix(c(1,2), 1, 2, byrow = TRUE),widths=c(0.5,0.5))
        if(verbose >= 2){cat("  Plotting two ternary plots on the one page\n\n")}
      } else if (npops2plot == 3) {
        layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE),widths=c(0.5,0.5), heights=c(0.5,0.5))
        if(verbose >= 2){cat("  Plotting three ternary plots on the one page\n\n")}
      } else if (npops2plot == 4) {
        layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE),widths=c(0.5,0.5), heights=c(0.5,0.5))
        if(verbose >= 2){cat("  Plotting four ternary plots on the one page\n\n")}
      } else {
        layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE),widths=c(0.5,0.5), heights=c(0.5,0.5))
        if(verbose >= 2){cat("  Plotting four ternary plots per page\n\n")}
      }
    } else {
      if(verbose >= 2){cat("  Plotting one ternary plot\n")}
    }
 
    # Plot the tertiary plots
    for (i in poplist) {
      ii <- gl.filter.monomorphs(i,verbose=0)
      if (nLoc(ii) == 0){
        next
      }
      # Plot the tertiary plots
      a <- utils.hwe(ii, prob=alpha,verbose=0)
      mat <- array(NA,3*dim(a)[1])
      dim(mat) <- c(dim(a)[1],3)
      mat[,1] <- as.numeric(as.character(a$Hom_1)) # for God knows why
      mat[,2] <- as.numeric(as.character(a$Het))
      mat[,3] <- as.numeric(as.character(a$Hom_2))
      colnames(mat) <- c("AA", "AB", "BB")
      if (bonf) {
        c.alpha <- alpha/length(mat[,1])
        xlabel <- paste0("Population: ",popNames(i)[1])
      } else {
        c.alpha <- alpha
        xlabel <- paste0("Population: ",popNames(i)[1])
      }
      if (method == "Fisher"){
        res <- HardyWeinberg::HWTernaryPlot(mat, 100, region = 7, vertex.cex = 1.25, alpha=c.alpha, signifcolour = TRUE, vbounds=TRUE, axislab=xlabel)
      } else {
        res <- HardyWeinberg::HWTernaryPlot(mat, 100, region = 2, vertex.cex = 1.25, alpha=c.alpha, signifcolour = TRUE, vbounds=TRUE, axislab=xlabel)
      } 
    }
  }
  layout(mat = matrix(c(1,1),1,1, byrow=FALSE))
  
  #### Report the results
  
  rprob <- as.numeric(as.character(result$Prob))
  result <- result[(rprob>0 & rprob<=alpha),]
  result <- result[order(result$Locus),]
  cat("    Reporting significant departures from Hardy-Weinberg Equilibrium\n")
  if (nrow(result)==0){
    cat("    No significant departures\n")
  } else {
    cat("    NB: Departures significant at the alpha level of",alpha,"are listed\n")
    if (alpha > 0.05) {
      cat("      ns --",alpha,"< p < 0.05; * -- 0.05 < p < 0.01; ** -- 0.01 < p < 0.001; etc\n")
    } else {
      cat("      ns -- p > 0.05; * -- 0.05 < p < 0.01; ** -- 0.01 < p < 0.001; etc\n")
    }
      cat("    Critical values for significance of Bonferroni Corrected significance vary with sample size\n\n")
    print(result, row.names=FALSE)
  }  
  
  # FLAG SCRIPT END
  
  if(verbose >= 1){
    cat("Completed:",funname,"\n")
  }
  
  return(result)

}
