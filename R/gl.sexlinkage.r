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
#' @param t.het -- tolerance, that is tm=0.05 means that 5% of the heterogametic sex can be homozygous and still 
#' be regarded as consistent with a sex specific marker [default 0]
#' @param t.hom -- tolerance, that is tf=0.05 means that 5% of the homogametic sex can be heterozygous and still
#' be regarded as consistent with a sex specific marker [default 0]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return The list of sex specific loci
#' @export
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' result <- gl.sexlinkage(testset.gl)

gl.sexlinkage <- function(x, t.het=0, t.hom=0, v=2) {

  if (v > 0) {
    cat("Starting gl.sexlinkage: Identifying sex linked loci\n")
  }
  
  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required for gl.sexlinkage!\n"); stop("Execution terminated\n")
  }

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

# Combine the two files
  Trimmed_Sequence <- x@other$loc.metrics$TrimmedSequence
  df <- cbind(dff,dfm,Trimmed_Sequence)
  a <- strsplit(row.names(df), split="-")
  a <- do.call(rbind,a)
  a <- as.numeric(a[,1])
  
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
  df_zero_sum <- df[(sumf == 0) | (summ == 0), ]
  df <- df[(sumf > 0) & (summ > 0), ]
  zw <- df[df$F1/(sumf[(sumf > 0) & (summ > 0)]) >= (1 - t.hom) &
             df$M1/(summ[(sumf > 0) & (summ > 0)]) <= (0 + t.het), ]
  xy <- df[df$F1/(sumf[(sumf > 0) & (summ > 0)]) <= (0 + t.het) &
             df$M1/(summ[(sumf > 0) & (summ > 0)]) >= (1 - t.hom), ]
  
  if (nrow(zw) == 0){
    cat("No sex linked markers consistent with female heterogamety (ZZ/ZW)\n")
  } else {
    cat("\nSex linked loci consistent with female heterogamety (ZZ/ZW)\n")
    cat(paste("  Threshold proportion for homozygotes in the heterozygotic sex (ZW)",t.hom,";\n")) 
    cat(paste("  for heterozygotes in the homozygotic sex (ZZ)",t.het,"\n"))
    cat("  0 = homozygous reference; 1 = heterozygous; 2 = homozygous alternate\n")
    print(zw)
    cat("Note: Snp location in Trimmed Sequence indexed from 0 not 1, SNP position in lower case\n")
    cat("Note: The most reliable putative markers will have AvgCount for Ref or Snp 10 or more, one ca half the other\n")
  }
  if (nrow(xy) == 0){
    cat("No sex linked markers consistent with male heterogamety (XX/XY)\n")
  } else {
    cat("\nSex linked loci consistent with male heterogamety (XX/XY)\n")
    cat(paste("  Threshold proportion for homozygotes in the heterozygotic sex (XY)",t.hom,";\n")) 
    cat(paste("  for heterozygotes in the homozygotic sex (XX)",t.het,"\n"))
    cat("  0 = homozygous reference; 1 = heterozygous; 2 = homozygous alternate\n")
    print(xy)
    cat("Note: Snp location in Trimmed Sequence indexed from 0 not 1, SNP position in lower case\n")
    cat("Note: The most reliable putative markers will have AvgCount for Ref or Snp 10 or more, one ca half the other\n")
  }
  if (nrow(df_zero_sum) != 0){
    cat("\nFound loci with zero alleles for one sex\n")
    print(df_zero_sum)
  } 
  
  l <- list(zw,xy, df_zero_sum)
  
  l
  
}
