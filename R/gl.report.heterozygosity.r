#' Reports observed and expected hetrozygosity by population
#' 
#' Calculates the observed and expected heterozygisities by population from a genlight object
#' and plots as an ordered barchart of observed heterozygosity. Also calculates the estimated variance
#' for expected heterozygosity.
#' 
#' Expected heterozygosity is calculated according to Nei, M. (1987) Molecular evolutionary genetics. New York: Columbia University Press, using the Pegas R package.
#' 
#' @param x -- a genlight object containing the SNP genotypes [Required]
#' @param plot -- if TRUE, plots a barchart of observed heterozygosity [default FALSE]
#' @param cex.labels -- sets the size of the population labels [default 0.7]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return a dataframe containing population labels, heterozygosities and sample sizes
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @importFrom grDevices rainbow
#' @importFrom graphics par
#' @importFrom plyr join
#' @importFrom pegas heterozygosity
#' @examples
#' tmp <- gl.report.heterozygosity(testset.gl,plot=TRUE,v=3)

gl.report.heterozygosity <- function(x, plot=FALSE, cex.labels=0.7, v=2) {
  
# ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required for gl.drop.pop.r!\n"); stop("Execution terminated\n")
  }
  if (v > 0) {
    cat("Starting gl.report.heterozygosity: Obs and Exp\n")
  }
  if (v < 0 | v > 5){
    cat("Warning: Verbosity must take on an integer value between 0 and 5, set to 3\n")
    v <- 3
  }
  
# DO THE DEED  

# Split the genlight object into a list of populations
  sgl <- seppop(x)
  
# OBSERVED HETEROZYGOSITY  
  if (v >=2){cat("  Calculating Observed Heterozygosities by population\n")}
  
  # Calculate heterozygosity for each population in the list
    Ho <- data.frame(lapply(sgl, function(x) mean(colMeans(as.matrix(x, na.rm=TRUE)==1), na.rm=TRUE) ))
    
  # Calculate sample sizes
    sums <- data.frame(lapply(sgl, function(x) mean(colSums(as.matrix(x, na.rm=TRUE)==1), na.rm=TRUE) ))
    n <- t(sums/Ho)
    n <- cbind(row.names(n),n)
    
  # Join the sample sizes with the heteozygosities
    df1 <- data.frame(pop=names(Ho), Ho=as.numeric(Ho[1,]))
    df2 <- data.frame(n)
    names(df2)<- c("pop","nInd")
    df <- plyr::join(df1,df2, by="pop")
    
  # Plot the results
  if (plot){
    par(mfrow=c(1,1),mai=c(2.5,1,0.5,0.2),pty="m")
    df.ordered <- df[order(df$Ho),]
    barplot(df.ordered$Ho, names.arg=paste(df.ordered$pop, df.ordered$nInd, sep=" | "), las=2, cex.names=cex.labels, space=0, border=F, col=rainbow(nrow(df.ordered)), main="Observed Heterozygosity")
  }
  
# EXPECTED HETEROZYGOSITY
    if (v >=2){cat("  Calculating Expected Heterozygosities by population\n")}
  
  freq <- array(NA,2)
  He <- array(NA,length(sgl))
  Hv <- array(NA,length(sgl))
  ##Hef <- array(NA,length(sgl))
  ##He.n <- array(NA,length(sgl))

  for (i in 1:length(sgl)){
    gl <- sgl[[i]]
    
    # Calculate the number of heterozygous and homozygous individuals
    tbl <- table(as.matrix(gl)[,i],useNA = "always")
    if (is.na(tbl["0"])) {tbl["0"] <- 0}
    if (is.na(tbl["1"])) {tbl["1"] <- 0}
    if (is.na(tbl["2"])) {tbl["2"] <- 0}
    
    # Calculate the allele frequencies
    freq[1] <- (2*tbl[1] + tbl[2])
    freq[2] <- (tbl[2] + 2*tbl[3])
    
    # Calculate the expected heterozygosity and its variance
    tmp <- pegas::heterozygosity(freq, variance=TRUE)
    He[i] <- tmp[1]
    Hv[i] <- tmp[2]
    
    # Sample size
    ##He.n[i] <- nInd(gl)

  }
  ##df <- data.frame(popNames(x),as.numeric(table(pop(x))),Ho,round(Ho.n,1),He, Hv, He.n)
  ##names(df) <- c("pop","nInd","Ho","Ho.n","He", "He.var","He.n")
  df <- data.frame(popNames(x),as.numeric(table(pop(x))),df$Ho,He,Hv)
  names(df) <- c("pop","nInd","Ho","He","He.var")

# PLOT OBSERVED AGAINST EXPECTED 
  
#  if (plot){
#    par(mai=c(1,1,1,1),pty="s")
#    p <- ggplot(df, aes(x=Ho, y=He))
#    p <- p + geom_point(cex=2,col="red") 
#    p <- p + geom_abline(intercept = 0, slope = 1)
#    p <- p + xlab("Observed Heterozygosity") + ylab("Expected Heterozygosity")
#    p <- p + theme(panel.background = element_rect(fill = "white"),
#                 plot.margin = margin(2, 2, 2, 2, "cm"),
#                 plot.background = element_rect(fill = "grey90",colour = "black",size = 1),
#                 axis.title=element_text(size=15),
#                 axis.text=element_text(size=15)
#    )
#    plot(p)
#  }  

# OUTPUT REPORT
  if (v >= 3){
    cat("Reporting Heterozygosity by Population\n")
    cat("No. of loci =", nLoc(x), "\n")
    cat("No. of individuals =", nInd(x), "\n")
    cat("No. of populations =", nPop(x), "\n")
  
    cat("  Miniumum Observed Heterozygosity: ",round(min(df$Ho),4),"\n")
    cat("  Maximum Observed Heterozygosity: ",round(max(df$Ho),4),"\n")
    cat("  Average Observed Heterozygosity: ",round(mean(df$Ho),4),"\n\n")
  
    print(df)
  }  
  
  
# FLAG SCRIPT END
  if (v >= 1) {
    cat("\ngl.report.heterozygosity Completed\n")
  }
  
  # Return the result
  return(df) 
}