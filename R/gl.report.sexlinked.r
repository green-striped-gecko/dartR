#' @name gl.report.sexlinked
#' 
#' @title Identify loci that are sex linked
#'
#' @description 
#' Alleles unique to the Y or W chromosome and monomorphic on the X chromosomes 
#' will appear in the SNP dataset as genotypes that are heterozygotic in all 
#' individuals of the heterogametic sex and homozygous in all individuals of the 
#' homogametic sex. This function identifies loci with alleles that behave in 
#' this way, as putative sex specific SNP markers.
#'
#' @param x Name of the genlight object containing the SNP or presence/absence
#'  (SilicoDArT) data [required].
#' @param sex Factor that defines the sex of individuals. See explanation in details [default NULL].
#' @param t.het Tolerance in the heterogametic sex, that is t.het=0.05 means that 5\% of the heterogametic sex can be homozygous and still be regarded as consistent with a sex specific marker [default 0].
#' @param t.hom Tolerance in the homogametic sex, that is t.hom=0.05 means that 5\% of the homogametic sex can be heterozygous and still be regarded as consistent with a sex specific marker [default 0].
#' @param t.pres Tolerance in presence, that is t.pres=0.05 means that a silicodart marker can be present in either of the sexes and still be regarded as a sex-linked marker [default 0].
#' @param plot Whether produce a plot of the results [default TRUE].
#' @param plot_theme Theme for the plot. See Details for options [default theme_dartR()].
#' @param plot_colours List of three color names for the not sex-linked loci, for the sex-linked loci and for the area in which sex-linked loci appear [default three_colors].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default NULL, unless specified using gl.set.verbosity].
#'
#' @details 
#' Sex of the individuals for which sex is known with certainty can be provided 
#' via a factor (equal to the length of the number of individuals) or to be held 
#' in the variable \code{x@other$ind.metrics$sex}.
#' Coding is: M for male, F for female, U or NA for unknown/missing. 
#' The script abbreviates the entries here to the first character. So, coding of
#' "Female" and "Male" works as well. Character are also converted to upper cases.
#' 
#''\strong{ Function's output }
#'
#' This function creates a plot that shows the heterozygosity of males and females at each loci 
#' for SNP data or percentage of present/absent in the case of SilicoDArT data.
#'
#'  Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'
#' @return Two lists of sex-linked loci, one for XX/XY and one for ZZ/ZW systems and a plot.
#'
#' @author Arthur Georges, Bernd Gruber & Floriaan Devloo-Delva (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' out <- gl.report.sexlinked(testset.gl)
#'
#' @family reporting functions
#'
#' @export
#'  

gl.report.sexlinked <- function(x,
                                sex = NULL, 
                                t.het = 0, 
                                t.hom = 0,
                                t.pres = 0,
                                plot = TRUE,
                                plot_theme = theme_dartR(), 
                                plot_colours = three_colors, 
                                verbose = NULL) {
  
  # TRAP COMMAND
  
  funname <- match.call()[[1]]
  
  # SET VERBOSITY
  
  verbose <- gl.check.verbosity(verbose)
  
  # CHECKS DATATYPE 
  
  datatype <- utils.check.datatype(x)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
  # Check for monomorphic loci
  tmp <- gl.filter.monomorphs(x, verbose=0)
  if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {
    cat(warn("  Warning: genlight object contains monomorphic loci\n"))
  }
  
  # FLAG SCRIPT START
  
  if (verbose >= 1) {
    if (verbose == 5) {
      cat(report("\n\nStarting", funname, "[ Build =", 
                 build, "]\n\n"))
    } else {
      cat(report("\n\nStarting", funname, "\n\n"))
    }
  }

# DO THE JOB
  
# sex should be provided as it is not a default setting, if not provided it will be searched here:
  if (is.null(sex)){
    sex <- x@other$ind.metrics$sex
    }
  
  if (is.null(sex)){
    stop(error("No definition for the sex of individuals is provided. If not 
               provided via the function call it needs to be at gl@other$ind.metrics$sex. 
               Please refer to the help pages how you can define sex of individuals."))
  }

  if (length(sex) != nInd(x)){
    stop(error("The number of individuals and the number of entries defining the 
         sex do not match. Check your genlight object and your sex defining column."))
  }
  
  sex <- as.character(sex)
  UP <- toupper(sex)
  UP <- substring(UP, 1, 1)
  sex[UP == 'F'] <- 'F'
  sex[UP == 'M'] <- 'M'
  sex <- ifelse(sex=='F' | sex=='M', sex, 'U')
  sex[is.na(sex)] <- 'U'  
  
  ########### FOR SNP data
  
  if (datatype=="SNP"){ 
    
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
    cat("Sex ratio females:(males+females):",round(sum(sex=="F")/(sum(sex=="F")+sum(sex=="M")),2),"\n")
  }
    if (nrow(zw) == 0 & verbose>0){
    cat(important("  No sex linked markers consistent with female heterogamety (ZZ/ZW)\n"))
      } else {
  if(verbose>0){
    cat("\n  Sex linked loci consistent with female heterogamety (ZZ/ZW)\n")
    cat(paste("    Threshold proportion for homozygotes in the heterozygotic sex (ZW)",t.hom,";\n"))
    cat(paste("    for heterozygotes in the homozygotic sex (ZZ)",t.het,"\n"))
    cat("    0 = homozygous reference; 1 = heterozygous; 2 = homozygous alternative\n\n")
    print(zw)
    cat(important("  Note: The most reliable putative markers will have read depth of 10 or more.\n\n"))
  } 
}
  
if (nrow(xy) == 0 & verbose>0 ){
  cat(important("  No sex linked markers consistent with male heterogamety (XX/XY)\n"))
} else {
  if(verbose>0){
    cat("\n  Sex linked loci consistent with male heterogamety (XX/XY)\n")
    cat(paste("    Threshold proportion for homozygotes in the heterozygotic sex (XY)",t.hom,"\n"))
    cat(paste("    for heterozygotes in the homozygotic sex (XX)",t.het,"\n"))
    cat("    0 = homozygous reference; 1 = heterozygous; 2 = homozygous alternative\n\n")
    print(xy)
    cat(important("  Note: The most reliable putative markers will have read depth for Ref or Snp 10 or more, one ca half the other\n\n"))
  }
  }

  if (plot){

    df$fhet <- dff$F1/(dff$F0+dff$F1+dff$F2)
    df$mhet <- dfm$M1/(dfm$M0+dfm$M1+dfm$M2)
    df$xy <- indexxy
    df$zw <- indexzw
    df$test <-  df$xy + df$zw 
    df_sex_linked <- df[which(df$test==1),]
    df_no_sex_linked <- df[which(df$test==0),]
    
    
    gg <- ggplot()+
      geom_rect(aes(xmin=0, xmax=t.het, ymin=1-t.het, ymax=1), fill=three_colors[3],alpha=1/2)+
      geom_text(x=0, y=1.03, aes(label="XX/XY"))+
      geom_rect(aes(xmin=1, xmax=1-t.het, ymin=0, ymax=t.het), fill=three_colors[3],alpha=1/2)+
      geom_text(x=1, y=-0.02, aes(label="ZZ/ZW"))  +
      geom_point(data =  df_no_sex_linked, aes(x=fhet, y=mhet),alpha=1/3, size=2,color=three_colors[1] )+
      geom_point(data = df_sex_linked, aes(x=fhet, y=mhet),alpha=1/2, size=3,color=three_colors[2] )+
      xlab("Female Heterozygosity")+ 
      ylab("Male Heterozygosity")+
      xlim(0,1)+
      ylim(0,1)+
      plot_theme
    
    print(gg)
  }
  
l <- list(xxxy=xy, zzzw=zw, plot=gg)

}
  
  ########### FOR SilicoDArT data
  
if (datatype=="SilicoDArT"){
  
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
  
  if (nrow(zw) == 0 & verbose>0){
    cat(important("  No sex linked markers consistent with female heterogamety (ZZ/ZW)\n"))
  } else {
    if(verbose>0){
      cat("\n  Sex linked loci consistent with female heterogamety (ZZ/ZW)\n")
      cat(paste("    Threshold proportion for homozygotes in the heterozygotic sex (ZW)",t.hom,";\n"))
      cat(paste("    for heterozygotes in the homozygotic sex (ZZ)",t.het,"\n"))
      cat("    0 = homozygous reference; 1 = heterozygous; 2 = homozygous alternative\n\n")
      print(zw)
      cat(important("  Note: The most reliable putative markers will have read depth of 10 or more.\n\n"))
    } 
  }
  if (nrow(xy) == 0){
    if(verbose>0) cat(important("  No sex linked markers consistent with male heterogamety (XX/XY)\n"))
  } else {
    if(verbose>0){
      cat("\n  Sex linked loci consistent with male heterogamety (XX/XY)\n")
      cat(paste("    Threshold proportion for homozygotes in the heterozygotic sex (XY)",t.hom,"\n"))
      cat(paste("    for heterozygotes in the homozygotic sex (XX)",t.het,"\n"))
      cat("    0 = homozygous reference; 1 = heterozygous; 2 = homozygous alternative\n\n")
      print(xy)
      cat(important( "  Note: The most reliable putative markers will have read depth for Ref or Snp 10 or more, one ca half the other\n\n"))
    } 
  }
  
  if (plot)  {
    fhet <- mhet <- NULL
    df$fhet <- dff$F1/(dff$F0+dff$F1)
    df$mhet <- dfm$M1/(dfm$M0+dfm$M1)
    df$xy <- indexxy
    df$zw <- indexzw
    df$test <-  df$xy + df$zw 
    df_sex_linked <- df[which(df$test==1),]
    df_no_sex_linked <- df[which(df$test==0),]
    
    gg <- ggplot()+
      geom_rect(aes(xmin=0, xmax=t.pres, ymin=1-t.pres, ymax=1), fill=three_colors[3],alpha=1/2)+
      geom_text(x=0, y=1.03, aes(label="XX/XY"))+
      geom_rect(aes(xmin=1, xmax=1-t.pres, ymin=0, ymax=t.pres), fill=three_colors[3],alpha=1/2)+
      geom_text(x=1, y=-0.02, aes(label="ZZ/ZW"))  +
      geom_point(data =  df_no_sex_linked, aes(x=fhet, y=mhet),alpha=1/3, size=2,color=three_colors[1] )+
      geom_point(data = df_sex_linked, aes(x=fhet, y=mhet),alpha=1/2, size=3,color=three_colors[2] )+
      xlab("% present in females")+ 
      ylab("% present in males")+
      xlim(0,1)+
      ylim(0,1)+
      plot_theme
    
    print(gg)
  }
  
  l <- list(xxxy=xy, zzzw=zw, plot=gg) 
}
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("\n\nCompleted:", funname, "\n\n"))
  }
  
  # RETURN
  
  invisible(l)
  
}


