#' Calculates departure from Hardy-Weinberg Equilibrium. Utility script not meant 
#' for end users.
#' 
#' Calculates the probabilities of agreement with H-W equilibrium based on observed
#' frequencies of reference homozygotes, heterozygotes and alternate homozygotes. 
#' Uses the exact calculations contained in function utils.prob.hwe() as developed by
#' Wigginton, JE, Cutler, DJ, and Abecasis, GR.
#' 
#' @param x -- a genlight object containing the SNP profiles for a population [Required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @param prob -- level of significance [Default 0.05]
#' @return Locus, Hom_1, Het, Hom_2, N, Prob, Sig, BonSig)
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' 

utils.hwe <- function (x, 
                       prob=0.05, 
                       verbose=NULL) {
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func=funname,build="Jackson",v=verbose)
  
  # CHECK DATATYPE 
  datatype <- utils.check.datatype(x)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  # Check monomorphs have been removed up to date
  if (x@other$loc.metrics.flags$monomorphs == FALSE){
    if (verbose >= 2){
      cat(warn("  Warning: Dataset contains monomorphic loci which will be included in the ",funname," calculations\n"))
    }  
  }

# DO THE JOB
  
  m <- as.matrix(x)  
  # Initialize arrays
  hom.ref <- array(data=NA, dim=ncol(m))
  hom.snp <- hom.ref
  het <- hom.ref
  total <- hom.ref 
  p.values <- hom.ref
  boncrit <- hom.ref
  sig <- hom.ref
  sig2 <- hom.ref
  bonsig <- hom.ref
  bonsig2 <- hom.ref
  columns <- colnames(m)
  
  # For each locus, calculate counts of homozygotes and heterozygotes
  # then probability of departure from H-W equilibrium.
  for (i in 1:ncol(m)) {
    hom.ref[i] <- length(which(m[,i]==0))
    hom.snp[i] <- length(which(m[,i]==2))
    het[i] <- length(which(m[,i]==1)) 
    # Significance
    p.values[i] <- utils.prob.hwe(het[i], hom.ref[i], hom.snp[i])
  }
  
  
    
  if (p.values[i] < 0) {
      sig2[i] <- NA
    } else if (p.values[i] > 0.05) {
      sig2[i] <- "ns"
    } else if (p.values[i] <=0.05 && p.values[i] > 0.01){
      sig2[i] <- "*"
    } else if (p.values[i] <= 0.01 && p.values[i] > 0.001){
      sig2[i] <- "**"
    } else if (p.values[i] <= 0.001) {
      sig2[i] <- "***"
    }
    # Sample size taking into account NAs
    total[i] <- hom.ref[i]+het[i]+hom.snp[i]
    # Significance adjusted for the number of tests
    boncrit[i] <- prob/total[i]
    bonsig[i] <- ((p.values[i] >= 0) && (p.values[i] < boncrit[i]))
    if (bonsig[i] == FALSE){
      bonsig2[i] <- "ns"
    } else if (bonsig[i] == TRUE){
      bonsig2[i] <- "*"
    }
  

  # Assemble results into a dataframe
  result <- data.frame(cbind(hom.ref, het, hom.snp, total, p.values, sig2, bonsig2))
  result <- cbind(columns, result)
  names(result) <- c("Locus", "Hom_1", "Het", "Hom_2", "N", "Prob", "Sig", "BonSig")
  
  # Report the results
  pc <- sum(sig,na.rm=TRUE)*100/length(sig)
  pc <- round(pc, digits=1)
  pc2 <- sum(bonsig,na.rm=TRUE)*100/length(bonsig)
  pc2 <- round(pc2, digits=1)
  #cat("\nPopulations combined:",nPop(x),",",nInd(x),"individuals -- ")
  #cat(levels(pop(x)),"\n")
  #cat(paste("  Number of loci examined:",ncol(as.matrix(x)),"\n"))
  #cat(paste("  Number of polymorphic loci examined:",ncol(x),"\n"))
  #cat(paste("  Polymorphic loci that depart significantly from HWe:",pc,"%\n"))
  #cat(paste("  Polymoprhic loci that depart significantly from HWe (Bonferroni Corrected):",pc2,"%\n"))

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

  return(result)
  
}
