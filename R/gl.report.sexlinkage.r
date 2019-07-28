#' Identify loci that are sex linked in specimens in a genlight \{adegenet\} object
#'
#' Alleles unique to the Y or W chromosome and monomorphic on the X chromosomes will appear in the SNP dataset as 
#' genotypes that are heterozygotic in all individuals of the heterogametic sex and homozygous in all individuals 
#' of the homogametic sex.
#' 
#' This script will identify loci with alleles that behave in this way, as putative sex specific SNP markers.
#' 
#' Sex of the individuals for which sex is known with certainty is to be held in the variable x@other$ind.metrics$sex, 
#' as M for male, F for female, NA otherwise. The script abbreviates the entries here to the first character. So coding of "Female" and "Male" works as well. Character are also converted to upper cases.
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param as.sex -- specify the individual metric variable that spedifies known sexes (M or F) [default 'sex']
#' @param t.het -- tolerance, that is tm=0.05 means that 5% of the heterogametic sex can be homozygous and still 
#' be regarded as consistent with a sex specific marker [default 0]
#' @param t.hom -- tolerance, that is tf=0.05 means that 5% of the homogametic sex can be heterozygous and still
#' be regarded as consistent with a sex specific marker [default 0]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return The list of sex specific loci
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' result <- gl.report.sexlinkage(testset.gl)

# Last amended 3-Feb-19

gl.report.sexlinkage <- function(x, as.sex="sex", t.het=0, t.hom=0, verbose=2) {

# TIDY UP FILE SPECS

  funname <- match.call()[[1]]

# FLAG SCRIPT START

  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }

  if (verbose > 0) {
    cat("Starting",funname,"\n")
  }

# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }

  # Work around a bug in adegenet if genlight object is created by subsetting
    x@other$loc.metrics <- x@other$loc.metrics[1:nLoc(x),]

  # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nLoc(x))
      pop(x) <- as.factor(pop(x))
    }

  # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x, verbose=0)
    if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {cat("  Warning: genlight object contains monomorphic loci\n")}

# SCRIPT SPECIFIC ERROR CHECKS
    
    If(!(as.sex %in% names(x@other$ind.metrics))) {
      cat("  Fatal Error: 'sex' or the individual metric specified as the sex variable does not exist\n")
      stop()
    }
    
# DO THE JOB

# Extract the sex variable from whereever it may be -- might need alteration    
  sex <- toupper(substr(x@other$ind.metrics$sex,1,1))

# Extract the data for the females
  matf <- as.matrix(x[sex=="F"])
  # For each individual
    f <- array(data=NA, dim=c(ncol(matf),3))
    for (i in 1:ncol(matf)) {
      for (j in 1:3) {
        f[i,j] <- length(which(matf[,i]==(j-1)))
      }  
    }
    dff <- data.frame(f)
    row.names(dff) <- locNames(x)
    colnames(dff) <- c("F0","F1","F2")

# Extract the data for the males
  matm <- as.matrix(x[sex=="M"])
# For each individual
  m <- array(data=NA, dim=c(ncol(matm),3))
  for (i in 1:ncol(matm)) {
    for (j in 1:3) {
      m[i,j] <- length(which(matm[,i]==(j-1)))
    }  
  }
  dfm <- data.frame(m)
  row.names(dfm) <- locNames(x)
  colnames(dfm) <- c("M0","M1","M2")
  
  
  # Save the prior settings for mfrow, oma, mai and pty, and reassign
  op <- par(mfrow = c(1, 2), oma=c(1,1,1,1), mai=c(0.5,0.5,0.5,0.5),pty="m")
  # Set margins for first plot
  par(mai=c(1.5,1.5,1.5,1.5))
  f <- dff$F1/(dff$F0+dff$F1+dff$F2)
  m <- dfm$M1/(dfm$M0+dfm$M1+dfm$M2)

  # Plot Box-Whisker plot
  plot(x=f,y=m,xlab="Female Heterozygosity",ylab="Male Heterozygosity",col='red')
  p1 <- ggplot(df, aes(x=f, y=m-f)) + 
    geom_point(colour='red') +
    labs(y='Ho Difference', x='Female Heterozygosity') +
    xlim(0, 1) + ylim(-1, 1) +
    geom_hline(yintercept=c(0,1,-1))
  z <- m-f
  p2 <- adjbox(c(0,m),
         horizontal = F,
         col='red',
         range=range,
         main = "Box and Whisker Plot")
  
  grid.arrange(p1, p2, nrow = 1)
  par(op)

  

# Combine the two files
  Trimmed_Sequence <- x@other$loc.metrics$TrimmedSequence
  df <- cbind(dff,dfm,Trimmed_Sequence)
  a <- strsplit(row.names(df), split="-")
  a <- do.call(rbind,a)
  a <- strsplit(a[,1], split="\\|")
  a <- do.call(rbind,a)
  a <- as.numeric(a[,2])
  
  df$Trimmed_Sequence <- as.character(df$Trimmed_Sequence)
  b <- substr(df$Trimmed_Sequence,1,a)
  c <- substr(df$Trimmed_Sequence,a+1,a+1)
  c <- tolower(c)
  d <- substr(df$Trimmed_Sequence,a+2,nchar(df$Trimmed_Sequence))
  
  df$Trimmed_Sequence <- paste0(b,c,d)
  df$AvgCountRef <- x@other$loc.metrics$AvgCountRef
  df$AvgCountSnp <- x@other$loc.metrics$AvgCountSnp
  
# Check for hets in all males, homs in all females (XY); ditto for ZW
  sumf <- df$F0+df$F1+df$F2
  summ <- df$M0+df$M1+df$M2
  # Pull loci that are 100% homozygous for females and 100% heterozygous for males
  index <- ((df$F0/(sumf)>=(1-t.hom) | df$F2/(sumf)>=(1-t.hom)) & df$M1/(summ)>=(1-t.het))
  zw <- df[index,]
  # Pull loci that are 100% homozygous for males and 100% heterozygous for females
  index <- ((df$M0/(summ)>=(1-t.hom) | df$M2/(summ)>=(1-t.hom)) & df$F1/(sumf)>=(1-t.het))
  xy <- df[index,]
  
if (nrow(zw) == 0){
  cat("  No sex linked markers consistent with female heterogamety (ZZ/ZW)\n")
} else {
  cat("\n  Sex linked loci consistent with female heterogamety (ZZ/ZW)\n")
  cat(paste("    Threshold proportion for homozygotes in the heterozygotic sex (ZW)",t.hom,";\n"))
  cat(paste("    for heterozygotes in the homozygotic sex (ZZ)",t.het,"\n"))
  cat("    0 = homozygous reference; 1 = heterozygous; 2 = homozygous alternate\n")
  print(zw)
  cat("  Note: Snp location in Trimmed Sequence indexed from 0 not 1, SNP position in lower case\n")
  cat("  Note: The most reliable putative markers will have AvgCount for Ref or Snp 10 or more, approx one half of that for the other\n")
}
  
if (nrow(xy) == 0){
  cat("  No sex linked markers consistent with male heterogamety (XX/XY)\n")
} else {
  cat("\n  Sex linked loci consistent with male heterogamety (XX/XY)\n")
  cat(paste("    Threshold proportion for homozygotes in the heterozygotic sex (XY)",t.hom,"\n"))
  cat(paste("    for heterozygotes in the homozygotic sex (XX)",t.het,"\n"))
  cat("    0 = homozygous reference; 1 = heterozygous; 2 = homozygous alternate\n")
  print(xy)
  cat("  Note: Snp location in Trimmed Sequence indexed from 0 not 1, SNP position in lower case\n")
  cat("  Note: The most reliable putative markers will have AvgCount for Ref or Snp 10 or more, one ca half the other\n")
  }

l <- list(zw,xy)

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

return(l)

}
