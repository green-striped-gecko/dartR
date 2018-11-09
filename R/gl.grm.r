#' Calculates the genomic relatedness matrix
#' 
#' The G matrix is calculated by centering the allele frequency matrix  of the second allele by substracting 2 times the allefrequency
#'@param gl -- a genlight object 
#'@param plotheatmap -- a switch if a heatmap should be shown [Default:TRUE] 
#'@param return.imputed switch if loci with imputed data should be returned (see ?A.mat in package rrBLUP)
#'@param ... parameters passed to function A.mat from package rrBLUP
#'@return a genomic relatedness matrix 
#'@importFrom stats heatmap cov var
#'@importFrom rrBLUP A.mat
#'@export
#'  
#'@examples
#'gl.grm(foxes.gl[1:5,1:10])  


gl.grm <- function(gl, plotheatmap=TRUE, return.imputed=FALSE, ...)
{
G <- A.mat(as.matrix(gl)-1,return.imputed = return.imputed, ...)
if (plotheatmap & return.imputed==FALSE) heatmap(G) else heatmap(G$A)

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
return (G)
}