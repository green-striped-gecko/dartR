#' Identify loci that are sex linked in specimens in a genlight \code{adegenet} object
#'
#' Alleles unique to the Y or W chromosome and monomorphic on the X chromosomes will appear in the SNP dataset as 
#' genotypes that are heterozygotic in all individuals of the heterogametic sex and homozygous in all individuals 
#' of the homogametic sex.
#' 
#' This script will identify loci with alleles that behave in this way, as putative sex specific SNP markers.
#' 
#' Sex of the individuals for which sex is known with certainty can be provided via a factor (equal to the length of the number of individuals) or to be held in the variable \code{x@other$ind.metrics$sex}.
#' Coding is: M for male, F for female, U or NA for unknown/missing. The script abbreviates the entries here to the first character. So coding of "Female" and "Male" works as well. Character are also converted to upper cases.
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param sex -- factor that defines the sex of individuals. See explanation above.
#' @param t.het -- tolerance, that is t.het=0.05 means that 5\% of the heterogametic sex can be homozygous and still be regarded as consistent with a sex specific marker [default 0]
#' @param t.hom -- tolerance, that is t.hom=0.05 means that 5\% of the homogametic sex can be heterozygous and still be regarded as consistent with a sex specific marker [default 0]
#' @param t.pres -- tolerance, that is t.pres=0.05 means that a silicodart marker can be present in either of the sexes and still be regarded as a sex-linked marker. [default 0]
#' @param plot -- creates a plot that shows the heterozygosity of males and females at each loci for SNP data or percentage of present/absent in the case of silicodart data.
#' be regarded as consistent with a sex specific marker [default 0]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return two list of sex specific loci, for XX/XY and ZZ/ZW systems.
#' @export
#' @author Arthur Georges, Bernd Gruber & Floriaan Devloo-Delvan (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' cat("does not work yet")
#' #result <- gl.sexlinkage(testset.gl)

# Last amended 3-Feb-19

gl.report.sexlinked <- function(x,sex=NULL, t.het=0, t.hom=0,t.pres=0, plot=TRUE,verbose=NULL) {

# TIDY UP FILE SPECS

  funname <- match.call()[[1]]

# FLAG SCRIPT START
  # set verbosity
  if (is.null(verbose) & !is.null(x@other$verbose)) verbose=x@other$verbose
  if (is.null(verbose)) verbose=2
 

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
  if (all(x@ploidy ==1)){
    if (verbose>0) ("Processing Presence/Absence (SilicoDArT) data [ploidy=1]), use gs.report.sexlinkage for this kind of data\n")
    data.type <- "SilicoDArT"
  } else if (all(x@ploidy == 2)){
    if (verbose >= 2){cat("  Processing a SNP dataset [ploidy=2]\n")}
    data.type <- "SNP"
  } else {
    stop("Fatal Error: Ploidy must be universally 2 (SNP data)")
  }
  
  
  
  
  # Work around a bug in adegenet if genlight object is created by subsetting
      if (nLoc(x)!=nrow(x@other$loc.metrics)) { stop("The number of rows in the loc.metrics table does not match the number of loci in your genlight object!")  }

  # Set a population if none is specified (such as if the genlight object has been generated manually)
  #  if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
  #    if (verbose >= 2){ cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")}
  #    pop(x) <- factor(rep("pop1",nInd(x)))
  #  }

  

  # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x, verbose=0)
    if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {cat("  Warning: genlight object contains monomorphic loci\n")}

# DO THE JOB
  
# sex should be provided as it is not a default setting, if not provided it will be searched here: reproducibility
  if (is.null(sex)) sex <- x@other$ind.metrics$sex
  
  if (is.null(sex)) stop("No definition for the sex of individuals is provided. If not provided via the function call it needs to be at gl@other$ind.metrics$sex. Please refer to the help pages how you can define sex of individuals.")

  
  if (length(sex) != nInd(x)) stop("The number of individuals and the number of entries defining the sex do not match. Check your genlight object and your sex defining column.")
  
  sex <- as.character(sex)
  UP <- toupper(sex)
  UP <- substring(UP, 1, 1)
  sex[UP == 'F'] <- 'F'
  sex[UP == 'M'] <- 'M'
  sex <- ifelse(sex=='F' | sex=='M', sex, 'U')
  sex[is.na(sex)] <- 'U'  
  
  
  if (data.type=="SNP") {  #for SNP data
# Extract the data for the females
  matf <- as.matrix(x[sex=="F",])
  # For each individual
    f <- array(data=NA, dim=c(ncol(matf),3))
    for (i in 1:ncol(matf)) {
      for (j in 1:3) {
        dummy <- sum(matf[,i]==(j-1), na.rm=T)
        if (is.na(dummy)) dummy <- 0
        f[i,j] <- dummy
      }  
    }
    dff <- data.frame(f)
    row.names(dff) <- locNames(x)
    colnames(dff) <- c("F0","F1","F2")

# Extract the data for the males
  matm <- as.matrix(x[sex=="M",])
# For each individual
  m <- array(data=NA, dim=c(ncol(matm),3))
  for (i in 1:ncol(matm)) {
    for (j in 1:3) {
      dummy <- sum(matm[,i]==(j-1), na.rm=T)
      if (is.na(dummy)) dummy <- 0
      m[i,j] <- dummy
    }  
  }
  dfm <- data.frame(m)
  row.names(dfm) <- locNames(x)
  colnames(dfm) <- c("M0","M1","M2")

# Combine the two files
  
  df <- cbind(dff,dfm)
  
  df$read.depth <- x@other$loc.metrics$rdepth
  
  
# Check for hets in all males, homs in all females (XY); ditto for ZW
  sumf <- df$F0+df$F1+df$F2
  summ <- df$M0+df$M1+df$M2
  # Pull loci that are 100% homozygous for females and 100% heterozygous for males
  indexxy <- ((df$F0/(sumf)>=(1-t.hom) | df$F2/(sumf)>=(1-t.hom)) & df$M1/(summ)>=(1-t.het))
  xy <- cbind(locnr =which(indexxy==TRUE), df[indexxy,])
 
  # Pull loci that are 100% homozygous for males and 100% heterozygous for females
  indexzw <- ((df$M0/(summ)>=(1-t.hom) | df$M2/(summ)>=(1-t.hom)) & df$F1/(sumf)>=(1-t.het))
  zw <- cbind(locnr =which(indexzw==TRUE),df[indexzw,])
 
  if(verbose>0) {
    
    cat("Number of females:",sum(sex=="F"),"\n")
    cat("Number of males:",sum(sex=="M"),"\n")
    cat("Sexratio females:(males+females):",round(sum(sex=="F")/(sum(sex=="F")+sum(sex=="M")),2),"\n")
  }
    if (nrow(zw) == 0){
  if(verbose>0) cat("  No sex linked markers consistent with female heterogamety (ZZ/ZW)\n")
} else {
  if(verbose>0) cat("\n  Sex linked loci consistent with female heterogamety (ZZ/ZW)\n")
  if(verbose>0) cat(paste("    Threshold proportion for homozygotes in the heterozygotic sex (ZW)",t.hom,";\n"))
  if(verbose>0) cat(paste("    for heterozygotes in the homozygotic sex (ZZ)",t.het,"\n"))
  if(verbose>0) cat("    0 = homozygous reference; 1 = heterozygous; 2 = homozygous alternate\n")
  if(verbose>0) print(zw)
  if(verbose>0) cat("  Note: The most reliable putative markers will have read depth of 10 or more.\n")
}
if (nrow(xy) == 0){
  if(verbose>0) cat("  No sex linked markers consistent with male heterogamety (XX/XY)\n")
} else {
  if(verbose>0) cat("\n  Sex linked loci consistent with male heterogamety (XX/XY)\n")
  if(verbose>0) cat(paste("    Threshold proportion for homozygotes in the heterozygotic sex (XY)",t.hom,"\n"))
  if(verbose>0) cat(paste("    for heterozygotes in the homozygotic sex (XX)",t.het,"\n"))
  if(verbose>0) cat("    0 = homozygous reference; 1 = heterozygous; 2 = homozygous alternate\n")
  if(verbose>0) print(xy)
  if(verbose>0) cat("  Note: The most reliable putative markers will have read depth for Ref or Snp 10 or more, one ca half the other\n")
  }

  if (plot)  {
      # Set margins for first plot
    df$fhet <- dff$F1/(dff$F0+dff$F1+dff$F2)
    df$mhet <- dfm$M1/(dfm$M0+dfm$M1+dfm$M2)
    gg <- ggplot(df, aes(x=df$fhet, y=df$mhet))+geom_rect(xmin=0, xmax=t.het, ymin=1-t.het, ymax=1, fill="darkgrey")+geom_text(x=0, y=1.03, label="XX/XY")+geom_rect(xmin=1, xmax=1-t.het, ymin=0, ymax=t.het, fill="darkgrey")+geom_text(x=1, y=-0.02, label="ZZ/ZW")  +geom_point(color=indexxy+indexzw+1,  alpha = 0.3, size=2)+xlab("Female Heterozygosity")+ ylab("Male Heterozygosity")+xlim(0,1)+ylim(0,1)
    print(gg)
  }
l <- list(xxxy=xy, zzzw=zw, plot=gg)
} #end if data.type='SNP'
  
if (data.type=="SilicoDArT")
{
  
  matf <- as.matrix(x[sex=="F"])
  # For each individual
  f <- array(data=NA, dim=c(ncol(matf),2))
  for (i in 1:ncol(matf)) {
    for (j in 1:2) {
      dummy <- sum(matf[,i]==(j-1), na.rm=T)
      if (is.na(dummy)) dummy <- 0
      f[i,j] <- dummy
    }  
  }
  dff <- data.frame(f)
  row.names(dff) <- locNames(x)
  colnames(dff) <- c("F0","F1")
  
  # Extract the data for the males
  matm <- as.matrix(x[sex=="M"])
  # For each individual
  m <- array(data=NA, dim=c(ncol(matm),2))
  for (i in 1:ncol(matm)) {
    for (j in 1:2) {
      dummy <- sum(matm[,i]==(j-1), na.rm=T)
      if (is.na(dummy)) dummy <- 0
      m[i,j] <- dummy
    }  
  }
  dfm <- data.frame(m)
  row.names(dfm) <- locNames(x)
  colnames(dfm) <- c("M0","M1")
  
  # Combine the two files
  
  df <- cbind(dff,dfm)
  
  df$read.depth <- x@other$loc.metrics$AvgReadDepth
  
  # Check for hets in all males, homs in all females (XY); ditto for ZW
  sumf <- df$F0+df$F1
  summ <- df$M0+df$M1
  # Pull loci that are 100% present in  females and 0% in males
  indexzw <- (df$F1/(sumf)>=(1-t.pres)  & df$M0/(summ)>=(1-t.pres))
  zw <- cbind(locnr =which(indexzw==TRUE), df[indexzw,])
  # Pull loci that are 100% present in males and 0% in females
  indexxy <- (df$M1/(summ)>=(1-t.pres)  & df$F0/(sumf)>=(1-t.pres))
  xy <- cbind(locnr =which(indexxy==TRUE),df[indexxy,])
  if (nrow(zw) == 0){
    if(verbose>0) cat("  No sex linked markers consistent with female heterogamety (ZZ/ZW)\n")
  } else {
    if(verbose>0) cat("\n  Sex linked loci consistent with female heterogamety (ZZ/ZW)\n")
    if(verbose>0) cat(paste("    Threshold proportion for homozygotes in the heterozygotic sex (ZW)",t.hom,";\n"))
    if(verbose>0) cat(paste("    for heterozygotes in the homozygotic sex (ZZ)",t.het,"\n"))
    if(verbose>0) cat("    0 = homozygous reference; 1 = heterozygous; 2 = homozygous alternate\n")
    if(verbose>0) print(zw)
    if(verbose>0) cat("  Note: The most reliable putative markers will have read depth of 10 or more.\n")
  }
  if (nrow(xy) == 0){
    if(verbose>0) cat("  No sex linked markers consistent with male heterogamety (XX/XY)\n")
  } else {
    if(verbose>0) cat("\n  Sex linked loci consistent with male heterogamety (XX/XY)\n")
    if(verbose>0) cat(paste("    Threshold proportion for homozygotes in the heterozygotic sex (XY)",t.hom,"\n"))
    if(verbose>0) cat(paste("    for heterozygotes in the homozygotic sex (XX)",t.het,"\n"))
    if(verbose>0) cat("    0 = homozygous reference; 1 = heterozygous; 2 = homozygous alternate\n")
    if(verbose>0) print(xy)

    if(verbose>0) cat("  Note: The most reliable putative markers will have read depth for Ref or Snp 10 or more, one ca half the other\n")
  }
  
  if (plot)  {
    # Set margins for first plot
    df$fhet <- dff$F1/(dff$F0+dff$F1)
    df$mhet <- dfm$M1/(dfm$M0+dfm$M1)
    gg <- ggplot(df, aes(x=df$fhet, y=df$mhet))+geom_rect(xmin=0, xmax=t.pres, ymin=1-t.pres, ymax=1, fill="darkgrey")+geom_text(x=0, y=1.03, label="XX/XY")+geom_rect(xmin=1, xmax=1-t.pres, ymin=0, ymax=t.pres, fill="darkgrey")+geom_text(x=1, y=-0.02, label="ZZ/ZW")  +geom_point(color=indexxy+indexzw+1,  alpha = 0.3, size=2)+xlab("% present in females")+ ylab("% present in males")+xlim(0,1)+ylim(0,1)
    print(gg)
  }
  l <- list(xxxy=xy, zzzw=zw, plot=gg) 
}
  
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
return(l)
}


