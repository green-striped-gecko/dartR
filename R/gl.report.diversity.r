#' Calculate diversity indices for SNPs
#'
#' !!Just an intro placeholder!! This script takes a genlight object and calculates alpha and beta diversity for q=0:2. Formulas are taken from Sherwin et al. 2017. The paper describes nicely the relationship between the different q levels and how they relate to population genetic processes such as dispersal and selection. For all indices the entropies (H) and corrosponding effective numbers Hill numbers (D), which reflect the amount of entities that are needed to get the observed valuea are calculated. In a nutshell the alpha indices between the different q-values should be similar if there are no deviation from expected allele frequencies and occurrences (e.g. all loci in HWE & equilibrium). If there is a deviation of an index this links to a process causing it such as dispersal, selection or strong drift. For a detailed explanation of all the indices, we recommend to resort to the literature provided below.
#'
#' @param gl genlight object containing the SNP genotypes  [required]
#' @param spectrumplot switch to provide a plot [default TRUE]
#' @param confiplot switch if confidence intervals (1 sd) should be drawn [default FALSE]
#' @param silent -- if FALSE, function returns an object, otherwise NULL [default TRUE]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @param pbar report on progress. Silent if set to FALSE. [default TRUE]
#' @param table prints a tabular output to the console either 'D'=D values, or 'H'=H values or 'DH','HD'=both or 'N'=no table.
#' 
#' @return if silent==TRUE, returns NULL; otherwise returns a list of entropy indices for each level of q and equivalent numbers for alpha and beta diversity.
#' @export
#' @importFrom graphics arrows
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr}), Contributors: William B. Sherwin, Alexander Sentinella 
#' @references 
#' Sherwin, W.B., Chao, A., Johst, L., Smouse, P.E. (2017). Information Theory Broadens the Spectrum of Molecular Ecology and Evolution. TREE 32(12) 948-963. doi:10.1016/j.tree.2017.09.12
#' 
#' Chao et al. 2014
#' 

### To be done:
# adjust calculation of betas for population sizes (switch)
# check that it works on SilicoDArT datasets


gl.report.diversity <- function(gl, spectrumplot=TRUE, confiplot=FALSE, pbar=TRUE, table="DH", silent=TRUE, verbose=NULL) {
  
# TIDY UP FILE SPECS
  
  build ='Jacob'
  funname <- match.call()[[1]]
  # Note does not draw upon or modify the loc.metrics.flags
  
# FLAG SCRIPT START
  # set verbosity
  if (is.null(verbose) & !is.null(x@other$verbose)) verbose=x@other$verbose
  if (is.null(verbose)) verbose=2
 
  
  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }
  
  cat("Starting",funname,"\n")
  
# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    stop("Fatal Error: genlight object required!\n")
  }
  
  if (all(x@ploidy == 1)){
    cat("  Detected Presence/Absence (SilicoDArT) data\n")
    stop("Cannot calculate minor allele frequences for Tag presence/absence data. Please provide a SNP dataset.\n")
  } else if (all(x@ploidy == 2)){
    cat("  Processing a SNP dataset\n")
  } else {
    stop("Fatal Error: Ploidy must be universally 1 (Tag P/A data) or 2 (SNP data)!\n")
  }
  
# FUNCTION SPECIFIC ERROR CHECKING
  
# DO THE JOB

if  (is.null(pop(gl)))  pop(gl)<- factor(rep("all", nInd(gl)))

#split in pops
pops <- seppop(gl)

if (pbar){
pb <- txtProgressBar(0,8, style = 3, width = 20) 
if(verbose >=2){cat(" Counting missing loci...           ")}
}
#number of missing loci

nlocpop <- lapply(pops, function(x) sum(!is.na(colMeans(as.matrix(x), na.rm = T))))
#nlocpop$totalloc <- nLoc(gl)

if (pbar){
setTxtProgressBar(pb, 1)
  if(verbose >=2){cat(" Calculating zero_H/D_alpha ...           ")}
##0Halpha (average number of alleles, ignoring missing values and sd)
}
zero_H_alpha_es <- lapply(pops, function(x) {
  dummys <- ((colMeans(as.matrix(x), na.rm=T) %% 2)>0)+1-1 
  return( list(estH=mean(dummys,  na.rm = T), sdH=sd(dummys, na.rm=T),estD =mean(dummys, na.rm = T)+1, sdD=sd(dummys, na.rm = T))) 
})
zero_H_alpha <- unlist(lapply(zero_H_alpha_es, function(x) x[[1]]))
zero_H_alpha_sd <- unlist(lapply(zero_H_alpha_es, function(x) x[[2]]))
zero_D_alpha <- unlist(lapply(zero_H_alpha_es, function(x) x[[3]]))
zero_D_alpha_sd <- unlist(lapply(zero_H_alpha_es, function(x) x[[4]]))

if (pbar){
setTxtProgressBar(pb, 2)
  if(verbose >=2){cat(" Calculating one_H/D_alpha ...          ")}
}
### one_H_alpha
one_H_alpha_es <- lapply(pops, function(x) {
  p<- colMeans(as.matrix(x), na.rm = T)/2
  p <- p[!is.na(p)] #ignore loci with just missing data
  logp <- ifelse(!is.finite(log(p)), 0, log(p))
  log1_p <- ifelse(!is.finite(log(1-p)), 0, log(1-p))
  
  dummys <- -(p*logp + (1-p)*log1_p)
  
  return(list(estH=mean(dummys), sdH=sd(dummys), estD=mean(exp(dummys)), sdD=sd(exp(dummys)), dummys=dummys))
})
one_H_alpha <- unlist(lapply(one_H_alpha_es, function(x) x[[1]]))
one_H_alpha_sd <- unlist(lapply(one_H_alpha_es, function(x) x[[2]]))
one_D_alpha <- unlist(lapply(one_H_alpha_es, function(x) x[[3]]))
one_D_alpha_sd <- unlist(lapply(one_H_alpha_es, function(x) x[[4]]))
if (pbar){
setTxtProgressBar(pb, 3)
  if(verbose >=2){cat(" Calculating two_H/D_alpha ...           ")}
}
#two_H_alpha
two_H_alpha_es <- lapply(pops, function(x) {
  p<- colMeans(as.matrix(x), na.rm = T)/2

  p <- p[!is.na(p)] #ignore loci with just missing data
  dummys <- (1-(p*p+(1-p)*(1-p)))
  
  return(list(estH=mean(dummys), sdH=sd(dummys), estD=mean(1/(1-dummys)), sdD=sd(1/(1-dummys)) , dummys=dummys ))
})
two_H_alpha <- unlist(lapply(two_H_alpha_es, function(x) x[[1]]))
two_H_alpha_sd <- unlist(lapply(two_H_alpha_es, function(x) x[[2]]))
two_D_alpha <- unlist(lapply(two_H_alpha_es, function(x) x[[3]]))
two_D_alpha_sd <- unlist(lapply(two_H_alpha_es, function(x) x[[4]]))



#initiallize betas as NA 
mat_zero_H_beta <- NA
mat_one_H_beta <- NA
mat_two_H_beta <- NA
npops <- length(pops)

if (npops>1)
{
if(pbar){
  setTxtProgressBar(pb, 4)
  if(verbose >=2){cat(" Counting pairwise missing loci...")}
}
pairs <- t(combn(npops,2))
### pairwise missing loci 
nlocpairpop <- apply(pairs,1, function(x)  {
  pop1 <- pops[[x[1]]]
  pop2 <- pops[[x[2]]]
  pp1 <- colMeans(as.matrix(pop1), na.rm = T)/2
  pp2 <- colMeans(as.matrix(pop2), na.rm = T)/2
  index <- !is.na(pp1) & !is.na(pp2)
  return(sum(index))
} )
mat_nloc_pops <- matrix(NA, nrow = npops, ncol = npops)
mat_nloc_pops[lower.tri(mat_nloc_pops)] <- nlocpairpop
colnames(mat_nloc_pops) <- rownames(mat_nloc_pops) <- names(pops)

if (pbar){
setTxtProgressBar(pb, 5)
  if(verbose >=2){cat(" Calculating zero_H/D_beta ...         ")}
}
#zero_H_beta
zero_H_beta_es <- apply(pairs,1, function(x)  {
  
  pop1 <- pops[[x[1]]]
  pop2 <- pops[[x[2]]]
  pp1 <- colMeans(as.matrix(pop1), na.rm = T)/2
  
  
  pp1 <- ifelse(pp1>0 & pp1<1, 0.5 , pp1)
  pp2 <- colMeans(as.matrix(pop2), na.rm = T)/2
  pp2 <- ifelse(pp2>0 & pp2<1, 0.5 , pp2)
  
  index <- !is.na(pp1) & !is.na(pp2)
  pp1 <- pp1[index]
  pp2 <- pp2[index]
  
  dummys <- abs(pp1-pp2)
  #  mat_zero_D_beta <- mat_zero_H_beta+1
  return(list(estH=mean(dummys), sdH=sd(dummys),estD= mean(dummys)+1,sdD =sd(dummys)))
  
} )

zero_H_beta <- unlist(lapply(zero_H_beta_es, function(x) x[[1]]))
zero_H_beta_sd<-unlist(lapply(zero_H_beta_es, function(x) x[[2]]))
zero_D_beta <- unlist(lapply(zero_H_beta_es, function(x) x[[3]]))
zero_D_beta_sd<-unlist(lapply(zero_H_beta_es, function(x) x[[4]]))

mat_zero_H_beta <- matrix(NA, nrow = npops, ncol = npops)
mat_zero_H_beta[lower.tri(mat_zero_H_beta)] <- zero_H_beta
mat_zero_H_beta[pairs] <- zero_H_beta_sd
colnames(mat_zero_H_beta) <- rownames(mat_zero_H_beta) <- names(pops)

mat_zero_D_beta <- matrix(NA, nrow = npops, ncol = npops)
mat_zero_D_beta[lower.tri(mat_zero_D_beta)] <- zero_D_beta
mat_zero_D_beta[pairs] <- zero_D_beta_sd
colnames(mat_zero_D_beta) <- rownames(mat_zero_D_beta) <- names(pops)

if (pbar){
  setTxtProgressBar(pb, 6)
  if(verbose >=2){cat(" Calculating one_H/D_beta ...    ")}
}

#one_H_beta
# calculate one_H_alpha_all for combined pops
p<- colMeans(as.matrix(gl), na.rm = T)/2
i0 <- which(!is.na(p))#ignore loci with just missing data
#p <- p[!is.na(p)] #ignore loci with just missing data
logp <- ifelse(!is.finite(log(p)), 0, log(p))
log1_p <- ifelse(!is.finite(log(1-p)), 0, log(1-p))
one_H_alpha_all <- -(p*logp + (1-p)*log1_p)

one_H_beta_es <- apply(pairs,1, function(x) { 
  i1 <- which(!is.na(colMeans(as.matrix(pops[[x[1]]]), na.rm = T)/2))
  i2 <- which(!is.na(colMeans(as.matrix(pops[[x[2]]]), na.rm = T)/2))
  tt <- table(c(i0,i1,i2))
  index <-as.numeric(names(tt)[tt==3])
  dummys <- one_H_alpha_all[i0 %in% index]-(one_H_alpha_es[[x[1]]]$dummys[i1 %in% index]+ one_H_alpha_es[[x[2]]]$dummys[i2 %in% index])/2
  return(list(estH=mean(dummys), sdH=sd(dummys), estD=mean(exp(dummys)), sdD=sd(exp(dummys))))
})

one_H_beta <- unlist(lapply(one_H_beta_es, function(x) x[[1]]))
one_H_beta_sd<-unlist(lapply(one_H_beta_es, function(x) x[[2]]))
one_D_beta <- unlist(lapply(one_H_beta_es, function(x) x[[3]]))
one_D_beta_sd<-unlist(lapply(one_H_beta_es, function(x) x[[4]]))

mat_one_H_beta <- matrix(NA, nrow = npops, ncol = npops)
mat_one_H_beta[lower.tri(mat_one_H_beta)] <- one_H_beta
mat_one_H_beta[pairs] <- one_H_beta_sd
colnames(mat_one_H_beta) <- rownames(mat_one_H_beta) <- names(pops)

mat_one_D_beta <- matrix(NA, nrow = npops, ncol = npops)
mat_one_D_beta[lower.tri(mat_one_D_beta)] <- one_D_beta
mat_one_D_beta[pairs] <- one_D_beta_sd
colnames(mat_one_D_beta) <- rownames(mat_one_D_beta) <- names(pops)

if (pbar){
  setTxtProgressBar(pb, 7)
  cat(" Calculating two_H/D_beta...    ")
}


#two_H_beta
# calculate two_H_alpha_all for combined pops
p<- colMeans(as.matrix(gl), na.rm = T)/2
i0 <- which(!is.na(p))#ignore loci with just missing data
#p <- p[!is.na(p)] #ignore loci with just missing data
two_H_alpha_all <- (1-(p*p+(1-p)*(1-p)))  

two_H_beta_es <- apply(pairs,1, function(x) {
  
  i1 <- which(!is.na(colMeans(as.matrix(pops[[x[1]]]), na.rm = T)/2))
  i2 <- which(!is.na(colMeans(as.matrix(pops[[x[2]]]), na.rm = T)/2))
  tt <- table(c(i0,i1,i2))
  index <-as.numeric(names(tt)[tt==3])
  
  m2Ha <- (two_H_alpha_es[[x[1]]]$dummys[i1 %in% index]+ two_H_alpha_es[[x[2]]]$dummys[i2 %in% index])/2
  dummys <- ((two_H_alpha_all[i0 %in% index]-m2Ha)/(1-m2Ha))*(npops/(npops-1))
#Johst-D
    return(list(estH=mean(dummys), sdH=sd(dummys), estD=mean(exp(dummys)), sdD=sd(exp(dummys))) )
 # return((two_H_beta_all-m2Ha)/two_H_beta_all  )
})

two_H_beta <- unlist(lapply(two_H_beta_es, function(x) x[[1]]))
two_H_beta_sd<-unlist(lapply(two_H_beta_es, function(x) x[[2]]))
two_D_beta <- unlist(lapply(two_H_beta_es, function(x) x[[3]]))
two_D_beta_sd<-unlist(lapply(two_H_beta_es, function(x) x[[4]]))

mat_two_H_beta <- matrix(NA, nrow = npops, ncol = npops)
mat_two_H_beta[lower.tri(mat_two_H_beta)] <- two_H_beta
mat_two_H_beta[pairs] <- two_H_beta_sd
colnames(mat_two_H_beta) <- rownames(mat_two_H_beta) <- names(pops)

mat_two_D_beta <- matrix(NA, nrow = npops, ncol = npops)
mat_two_D_beta[lower.tri(mat_two_D_beta)] <- two_D_beta
mat_two_D_beta[pairs] <- two_D_beta_sd
colnames(mat_two_D_beta) <- rownames(mat_two_D_beta) <- names(pops)

} # npops>1  


if (pbar){
  setTxtProgressBar(pb, 8)
  cat(" Done.                               ")
}


if (spectrumplot)
{

fs <- cbind(zero_D_alpha, one_D_alpha, two_D_alpha) 
cx <- max(1-(max(-12+nrow(fs),0)*0.025),0.5)
bb<-  barplot(fs, 
              beside = T, 
              names.arg = rep(rownames(fs),3), 
              ylim=c(1,2.15), 
              main="q-profile", 
              col=rainbow(npops), 
              las=2, 
              xpd=FALSE, 
              cex.names=cx)
 text(colMeans(bb), rep(2.1,3), labels = c("q=0", "q=1", "q=2")) 

 sds <- cbind(zero_D_alpha_sd, one_D_alpha_sd, two_D_alpha_sd)
 up <- fs+sds
 low <- fs-sds
if (confiplot)
  for (i in 1:ncol(bb))
 {
   for (ii in 1:nrow(bb))
   {
     arrows(x0 = bb[ii,i], y0 = up[ii,i], x1=bb[ii,i], y1=low[ii,i] ,angle = 90, code = 3, length = 0)
   }
 }

} #end spectrumplot

if (!is.na(match(table, c("H", "DH","HD"))))
{
tt <- data.frame(nloci=unlist(nlocpop),"m_0Ha"= zero_H_alpha, "sd_0Ha"= zero_H_alpha_sd,"m_1Ha"= one_H_alpha, "sd_1Ha"= one_H_alpha_sd, "m_2Ha"= two_H_alpha, "sd_2Ha"= two_H_alpha_sd)
print(knitr::kable(tt, digits = 3))
if (npops>1)
{
cat("\n\npairwise non-missing loci")  
print(knitr::kable(mat_nloc_pops, digits = 3))

cat("\n\n0_H_beta")
print(knitr::kable(mat_zero_H_beta, digits = 3))
cat("\n\n1_H_beta")
print(knitr::kable(mat_one_H_beta, digits = 3))
cat("\n\n2_H_beta")
print(knitr::kable(mat_two_H_beta, digits = 3))
}
}

if (!is.na(match(table, c("D", "DH","HD"))))
{
  tt <- data.frame(nloci=unlist(nlocpop),"m_0Da"= zero_D_alpha, "sd_0Da"= zero_D_alpha_sd,"m_1Da"= one_D_alpha, "sd_1Da"= one_D_alpha_sd, "m_2Da"= two_D_alpha, "sd_2Da"= two_D_alpha_sd)
  print(knitr::kable(tt, digits = 3))
  if (npops>1)
  {
  cat("\n\npairwise non-missing loci")  
  print(knitr::kable(mat_nloc_pops, digits = 3))
    cat("\n\n0_D_beta")  
  print(knitr::kable(mat_zero_D_beta, digits = 3))
  cat("\n\n1_D_beta")
  print(knitr::kable(mat_one_D_beta, digits = 3))
  cat("\n\n2_D_beta")
  print(knitr::kable(mat_two_D_beta, digits = 3))
  }
}
if (npops>1) out <-
  list(  nlocpop=unlist(nlocpop), 
         nlocpairpop = mat_nloc_pops,
         
         zero_H_alpha=    zero_H_alpha,
         zero_H_alpha_sd= zero_H_alpha_sd,
         one_H_alpha=     one_H_alpha,
         one_H_alpha_sd=  one_H_alpha_sd,
         two_H_alpha=     two_H_alpha,
         two_H_alpha_sd=  two_H_alpha_sd,
         
         zero_D_alpha   = zero_D_alpha, 
         zero_D_alpha_sd= zero_D_alpha_sd,
         one_D_alpha    = one_D_alpha,
         one_D_alpha_sd = one_D_alpha_sd,
         two_D_alpha    = two_D_alpha,
         two_D_alpha_sd = two_D_alpha_sd,
         
         zero_H_beta=mat_zero_H_beta, 
         one_H_beta=mat_one_H_beta, 
         two_H_beta=mat_two_H_beta, 
         
         zero_D_beta = mat_zero_D_beta, 
         one_D_beta = mat_one_D_beta, 
         two_D_beta = mat_two_D_beta) else out <- 
  list(  nlocpop=unlist(nlocpop), 
         
         zero_H_alpha=    zero_H_alpha,
         zero_H_alpha_sd= zero_H_alpha_sd,
         one_H_alpha=     one_H_alpha,
         one_H_alpha_sd=  one_H_alpha_sd,
         two_H_alpha=     two_H_alpha,
         two_H_alpha_sd=  two_H_alpha_sd,
         
         zero_D_alpha   = zero_D_alpha, 
         zero_D_alpha_sd= zero_D_alpha_sd,
         one_D_alpha    = one_D_alpha,
         one_D_alpha_sd = one_D_alpha_sd,
         two_D_alpha    = two_D_alpha,
         two_D_alpha_sd = two_D_alpha_sd)

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

  if(silent==TRUE){
    return(NULL)
  } else{
    return(out)
  } 

}




