#' Calculates the genomic relatedness matrix
#' 
#' The G matrix is calculated by centering the allele frequency matrix  of the second allele by substracting 2 times the allefrequency
#'@param x -- a genlight object 
#'@param plotheatmap -- a switch if a heatmap should be shown [Default:TRUE] 
#'@param return.imputed switch if loci with imputed data should be returned (see ?A.mat in package rrBLUP)
#'@param ... parameters passed to function A.mat from package rrBLUP
#'@param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#'@return a genomic relatedness matrix 
#'@importFrom stats heatmap cov var
#'@importFrom rrBLUP A.mat
#'@export
#'  
#'@examples
#'gl.grm(bandicoot.gl[1:5,1:10],plotheatmap=TRUE)  


gl.grm <- function(x, plotheatmap=TRUE, return.imputed=FALSE, verbose=NULL, ...){
  
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
      stop("  Detected Tag Presence/Absence (SilicoDArT) data. Please provide a SNP dataset\n")
    } else if (all(x@ploidy == 2)){
      cat("  Processing a SNP dataset\n")
    } else {
      stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)!\n")
    }

# DO THE JOB    
  
G <- A.mat(as.matrix(x)-1,return.imputed = return.imputed)
#if (plotheatmap & return.imputed==FALSE) heatmap(G) else heatmap(G$A) ####   BERND, G$A THROWS AND ERROR
if (plotheatmap & return.imputed==FALSE) heatmap(G) else heatmap(G)

# ff <- as.matrix(gl)
# alf <- colMeans(ff, na.rm = T)/2
# pjm <- matrix(rep(alf,nInd(gl)), ncol=nLoc(gl), nrow=nInd(gl), byrow=T)
# W <- ff  - (2*pjm)
# het <- 2*sum(alf *(1-alf))
#             
# G <- (W %*% t(W) )/het
# GG <-G
# ii<- !(colMeans(is.na(G))==1)
# GG <- GG[, ii ]
# ii<- !(rowMeans(is.na(G))==1)
# GG <- GG[ii,  ]
# 
# if (plotheatmap & nrow(GG)>0) heatmap(GG)

# FLAG SCRIPT END

  if(verbose >= 1){
    cat("Completed:",funname,"\n")
  }  

  return (G)
}