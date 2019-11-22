#' Reports observed and expected heterozygosity by population or by individual from SNP data. 
#' 
#' Calculates the observed and expected heterozygosities for each population (method="pop") 
#' or the observed heterozyosity for each individual (method="ind") in a genlight object.
#' 
#' Observed heterozygosity for a population takes the proportion of heterozygous
#' loci for each individual then averages over the individuals in that population. 
#' The calculations take into account missing values.
#' 
#' Expected heterozygosity for a population takes the expected proportion of
#' heterozygotes, that is, expected under Hardy-Weinberg equilibrium, for each locus, then
#' averages this across the loci for an average estimate for the population.
#' The calculations of expected heterozygosity use the unbiassed estimates of Nei, M. (1987) 
#' Molecular evolutionary genetics. New York: Columbia University Press.
#'
#' Output for method="pop" is an ordered barchart of observed heterozygosity across 
#' populations together with a table of observed and expected heterozygosity by population.
#' 
#' Observed heterozygosity for individuals is calculated as the proportion of loci that
#' are heterozygous for that individual.
#' 
#' Output for method="ind" is a histogram of heterozygosity across individuals.
#' The histogram is accompanied by a box and whisker plot presented either in standard 
#' (boxplot="standard") or adjusted for skewness (boxplot=adjusted). 
#' 
#' Refer to Tukey (1977, Exploratory Data Analysis. Addison-Wesley) for standard
#' Box and Whisker Plots and Hubert & Vandervieren (2008), An Adjusted Boxplot for Skewed
#' Distributions, Computational Statistics & Data Analysis 52:5186-5201) for adjusted
#' Box and Whisker Plots.
#' 
#' Finally, the loci that are invariant across all individuals in the dataset (that is,
#' across populations), is typically unknown. This can render estimates of heterozygosity
#' analysis specific, and so it is not valid to compare such estimates across species
#' or even across different analyses. This is a similar problem faced by microsatellites.
#' If you have an estimate of the number of invariant sequence tags (loci) in your data,
#' such as provided by gl.report.secondaries, you can specify it with the n.invariant
#' parameter to standardize your estimates of heterozygosity.
#' 
#' @param x -- a genlight object containing the SNP genotypes [Required]
#' @param method -- calculate heterozygosity by population (method='pop') or by individual (method='ind') [default 'pop']
#' @param n.invariant -- an estimate of the number of invariant sequence tags used to adjust the heterozygosity rate [default 0]
#' @param boxplot -- if 'standard', plots a standard box and whisker plot; if 'adjusted',
#' plots a boxplot adjusted for skewed distributions [default 'adjusted']
#' @param range -- specifies the range for delimiting outliers [default = 1.5 interquartile ranges]
#' @param cex.labels -- sets the size of the population labels [default 0.7]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return returns a dataframe containing population labels, heterozygosities and sample sizes
#' @export
#' @author Bernd Gruber, Arthur Georges and Renee Catullo (Post to \url{https://groups.google.com/d/forum/dartr})
#' @importFrom grDevices rainbow
#' @importFrom graphics par
#' @importFrom plyr join
#' @importFrom robustbase adjbox
#' @examples
#' out <- gl.report.heterozygosity(testset.gl,verbose=3)
#' out <- gl.report.heterozygosity(testset.gl,method='ind',verbose=3)

gl.report.heterozygosity <- function(x, 
                                     method="pop", 
                                     n.invariant=0,
                                     boxplot="adjusted",
                                     range=1.5,
                                     cex.labels=0.7,
                                     verbose=NULL) {
  
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
    stop("Fatal Error: genlight object required!")
  }
  
  if (all(x@ploidy == 1)){
    stop("  Processing  Presence/Absence (SilicoDArT) data, heterozygosity can only be calculated for SNP data\n")
    data.type <- "SilicoDArT"
  } else if (all(x@ploidy == 2)){
    if (verbose >= 2){cat("  Processing a SNP dataset\n")}
    data.type <- "SNP"
  } else {
    stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
  }

# SCRIPT SPECIFIC ERROR CHECKING

  if (!(method=="pop" | method == "ind")) {
    cat("Warning: Method must either be by population or by individual, set to method='pop'\n")
    method <- "pop"   
  }
  
  if (n.invariant < 0){
    cat("Warning: Number of invariant loci must be non-negative, set to zero\n")
    n.invariant <- 0
  }
  
  if (!(boxplot=="standard" | boxplot == "adjusted")) {
    cat("Warning: Box-whisker plots must either standard or adjusted for skewness, set to boxplot='adjusted'\n")
    boxplot <- 'adjusted'   
  }

  # Check for monomorphic loci
  
  if (x@other$loc.metrics.flags$monomorphs==FALSE) {
    if(verbose >= 1){
      cat("  Warning: genlight object contains monomorphic loci which will be factored into heterozygosity estimates\n")
    }
  }

# DO THE JOB FOR POPULATIONS
    
  if (method=="pop"){
    
  # Set a population if none is specified (such as if the genlight object has been 
    # generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  No population assignments detected, 
                             individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nInd(x))
      pop(x) <- as.factor(pop(x))
      }
    
  # Split the genlight object into a list of populations
    sgl <- seppop(x)
  
# OBSERVED HETEROZYGOSITY
  if (verbose >=2){cat("  Calculating Observed Heterozygosities, averaged across loci, for each population\n")}
    
  # Calculate heterozygosity for each population in the list
  Ho <- unlist(lapply(sgl, function(x) mean(colMeans(as.matrix(x)==1, na.rm=TRUE), na.rm=TRUE) ))
  
  # Calculate the number of loci used
  
  nl <- unlist(lapply(sgl, function(x) sum(colSums(is.na(as.matrix(x)))==0)))
  # Apply correction
  Ho.adj <- Ho*nLoc(x)/(nLoc(x)+n.invariant)
  
  # Calculate sample sizes =Number of individuals
  
  
  
  
  
  sums <- data.frame(lapply(sgl, function(x) mean(colSums(as.matrix(x)==1, na.rm=TRUE), na.rm=TRUE)))
  n <- t(sums/Ho)
  n <- cbind(row.names(n),n,nl,nLoc(x)/(nLoc(x)+n.invariant))
  
  # Join the sample sizes with the heteozygosities
  df1 <- data.frame(pop=names(Ho), Ho=as.numeric(Ho), Ho.adj=as.numeric(Ho.adj))
  df2 <- data.frame(n)
  names(df2)<- c("pop","nInd","nLoc","nLoc.adj")
  df <- plyr::join(df1,df2, by="pop")  

# EXPECTED HETEROZYGOSITY
  if (verbose >=2){cat("  Calculating Expected Heterozygosities\n")}
  
  Hexp <- array(NA,length(sgl))
  Hexp.adj <- array(NA,length(sgl))
  # For each population
  for (i in 1:length(sgl)){
    gl <- sgl[[i]]
    gl <- utils.recalc.freqhomref(gl,verbose=0)
    gl <- utils.recalc.freqhomsnp(gl,verbose=0)
    gl <- utils.recalc.freqhets(gl,verbose=0)
    p <- gl@other$loc.metrics$FreqHomRef
    q <- gl@other$loc.metrics$FreqHomSnp
    hets <- gl@other$loc.metrics$FreqHets
    p <- (2*p + hets)/2
    q <- (2*q + hets)/2
    H <- 1 - (p*p + q*q)
    Hexp[i] <- mean(H,na.rm=T)
    Hexp.adj[i] <- Hexp[i]*nLoc(x)/(nLoc(x)+n.invariant)
  }
  
  df <- data.frame(popNames(x),
                   as.numeric(table(pop(x))),
                   nl,
                   n.invariant,
                   round(df$Ho,6),
                   round(df$Ho.adj,6),
                   round(Hexp,6),
                   round(Hexp.adj,6)
                   )
  names(df) <- c("pop","nInd","nLoc","nLoc.inv","Ho","Ho.adj","He","He.adj")
  
  op <- par(mfrow=c(2,1),mai=c(1,1,1,0),pty="m")
  if(is.null(n.invariant)){
    df.ordered <- df[order(df$Ho),]
    barplot(df.ordered$Ho, names.arg=paste(df.ordered$pop, df.ordered$nInd, sep=" | "), las=2, cex.names=cex.labels, space=0, border=F, col=rainbow(nrow(df.ordered)), 
          main="Observed Heterozygosity by Population")
    df.ordered <- df[order(df$He),]
    barplot(df.ordered$He, names.arg=paste(df.ordered$pop, df.ordered$nInd, sep=" | "), las=2, cex.names=cex.labels, space=0, border=F, col=rainbow(nrow(df.ordered)), 
          main="Expected Heterozygosity by Population")
  } else {
    df.ordered <- df[order(df$Ho.adj),]
    barplot(df.ordered$Ho.adj, names.arg=paste(df.ordered$pop, df.ordered$nInd, sep=" | "), las=2, cex.names=cex.labels, space=0, border=F, col=rainbow(nrow(df.ordered)), 
          main="Observed Heterozygosity by Population")
    df.ordered <- df[order(df$He.adj),]
    barplot(df.ordered$He.adj, names.arg=paste(df.ordered$pop, df.ordered$nInd, sep=" | "), las=2, cex.names=cex.labels, space=0, border=F, col=rainbow(nrow(df.ordered)), 
          main="Expected Heterozygosity by Population")
  }
  
# OUTPUT REPORT
  if (verbose >= 3){
    cat("  Reporting Heterozygosity by Population\n")
    cat("\n  No. of loci =", nLoc(x), "\n")
    cat("  No. of individuals =", nInd(x), "\n")
    cat("  No. of populations =", nPop(x), "\n")
  
    cat("    Miniumum Observed Heterozygosity: ",round(min(df$Ho,na.rm=TRUE),6))
    if (n.invariant > 0) {
      cat("   [Corrected:",round(min(df$Ho.adj,na.rm=TRUE),6),"]\n")
    } else {
      cat("\n")
    }  
    cat("    Maximum Observed Heterozygosity: ",round(max(df$Ho,na.rm=TRUE),6))
    if (n.invariant > 0) {
      cat("   [Corrected:",round(max(df$Ho.adj,na.rm=TRUE),6),"]\n")
    } else {
      cat("\n")
    }  
    cat("    Average Observed Heterozygosity: ",round(mean(df$Ho,na.rm=TRUE),6))
    if (n.invariant > 0) {
      cat("   [Corrected:",round(mean(df$Ho.adj,na.rm=TRUE),6),"]\n\n")
    } else {
      cat("\n\n")
    }  
    
    cat("    Miniumum Expected Heterozygosity: ",round(min(df$He,na.rm=TRUE),6))
    if (n.invariant > 0) {
      cat("   [Corrected:",round(min(df$He.adj,na.rm=TRUE),6),"]\n")
    } else {
      cat("\n")
    }  
    cat("    Maximum Expected Heterozygosity: ",round(max(df$He,na.rm=TRUE),6))
    if (n.invariant > 0) {
      cat("   [Corrected:",round(max(df$He.adj,na.rm=TRUE),6),"]\n")
    } else {
      cat("\n")
    }  
    cat("    Average Expected Heterozygosity: ",round(mean(df$He,na.rm=TRUE),6))
    if (n.invariant > 0) {
      cat("   [Corrected:",round(mean(df$He.adj,na.rm=TRUE),6),"]\n\n")
    } else {
      cat("\n\n")
    }  
    
    if (n.invariant > 0){
      cat("  Average correction factor for invariant loci =",nLoc(x)/(nLoc(x)+n.invariant),"\n")
    } else {
      cat("  Heterozygosity estimates not corrected for uncalled invariant loci\n")
    }
  
    if(verbose >= 3){
      if (n.invariant > 0 ) {
        print(df)
      } else {
        print(df[,c("pop","nInd","nLoc","Ho","He")])
      } 
    }
    
  }
  
  }
  
  # DO THE JOB FOR INDIVIDUALS
  
  if (method=="ind"){
    if(verbose >= 2){
      cat("  Calculating observed heterozygosity for individuals\n")
      cat("  Note: No adjustment for invariant loci (n.invariant set to 0)\n")
    }
    # Convert to matrix
    m <- as.matrix(x)
    
    # For each individual determine counts of hets, homs and NAs
    c.na <- array(NA, nInd(x))
    c.hets <- array(NA, nInd(x))
    c.hom0 <- array(NA, nInd(x))
    c.hom2 <- array(NA, nInd(x))
    for (i in 1:nInd(x)){
      c.na[i] <- sum(is.na(m[i,]))
      c.hets[i] <- sum(m[i,]==1,na.rm=TRUE)/(nLoc(x)-c.na[i])
      c.hom0[i] <- sum(m[i,]==0,na.rm=TRUE)/(nLoc(x)-c.na[i])
      c.hom2[i] <- sum(m[i,]==2,na.rm=TRUE)/(nLoc(x)-c.na[i])
    }
    
    # Join the sample sizes with the heteozygosities
    df <- cbind.data.frame(x@ind.names, c.hets, c.hom0, c.hom2)
    names(df)<- c("ind.name", "Ho", "f.hom.ref", "f.hom.alt")
    
    # Prepare for plotting
    # Save the prior settings for mfrow, oma, mai and pty, and reassign
    op <- par(mfrow = c(2, 1),  mai=c(0.5,0.5,0.5,0.5),pty="m")
    # Set margins for first plot
    #par(mai=c(1,0.5,0.5,0.5))
    # Plot Box-Whisker plot
    if (boxplot == "standard"){
      whisker <- boxplot(df$Ho, horizontal=TRUE, col='red', range=range, main = "Heterozygosity by Individual")
      if (length(whisker$out)==0){
        if(verbose >= 1){cat("  Standard boxplot, no adjustment for skewness\n")}
      } else {
        outliers <- data.frame(spacer="     ",ID=as.character(df$ind.name[df$Ho %in% whisker$out]),
                          Ho=whisker$out
        )
        if(verbose >= 1){cat("  Standard boxplot, no adjustment for skewness\n")}
      }
      
    } else {
      whisker <- robustbase::adjbox(data=df, x=as.numeric(df$Ho),
           horizontal = TRUE,
           col='red',
           range=range,
           main = "Heterozygosity by Individual")
      if (length(whisker$out)==0){
        if(verbose >= 1){cat("  Boxplot adjusted to account for skewness\n")}
      } else {
        outliers <- data.frame(ID=as.character(df$ind.name[df$Ho %in% whisker$out]),
          Ho=whisker$out
          )
        if(verbose >= 1){cat("  Boxplot adjusted to account for skewness\n")}
      }
    }  
    # Set margins for second plot
    #par(mai=c(0.5,0.5,0,0.5))  
    # Plot Histogram
      hist(c.hets, col='red', main=NULL, breaks=100)
      
    # OUTPUT REPORT
      if (verbose >= 3){
        cat("Reporting Heterozygosity by Individual\n")
        cat("No. of loci =", nLoc(x), "\n")
        cat("No. of individuals =", nInd(x), "\n")

        cat("  Miniumum Observed Heterozygosity: ",round(min(df$Ho),6),"\n")
        cat("  Maximum Observed Heterozygosity: ",round(max(df$Ho),6),"\n")
        cat("  Average Observed Heterozygosity: ",round(mean(df$Ho),6),"\n\n")
        cat("  Results returned as a dataframe\n\n")
        if (length(whisker$out)==0){
          cat("  No outliers detected\n")
        } else {  
          cat("  Outliers detected -- \n")
          print(outliers)
        }  
      }   
  }  

  # Reset the par options    
    par(op)  
    
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

  return(df) 

}
