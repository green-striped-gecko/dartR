#' Calculates an identity by descent matrix
#' 
#' This function calculates the mean probability of identity by descent (IBD) across loci that would result from all the possible crosses of the individuals analyzed. IBD is calculated by an additive relationship matrix approach developed by Endelman and Jannink (2012) as implemented in the function \link[rrBLUP]{A.mat} (package rrBLUP). Two or more alleles are identical by descent (IBD) if they are identical copies of the same ancestral allele in a base population. The additive relationship matrix is a theoretical framework for estimating a relationship matrix that is consistent with an approach to estimate the probability that the alleles at a random locus are identical in state (IBS).
#' This function also plots a heatmap, and a dendrogram, of IBD values where each diagonal element has a mean that equals 1+f, where f is the inbreeding coefficient (i.e. the probability that the two alleles at a randomly chosen locus are IBD from the base population). As this probability lies between 0 and 1, the diagonal elements range from 1 to 2. Because the inbreeding coefficients are expressed relative to the current population, the mean of the off-diagonal elements is -(1+f)/n, where n is the number of loci. Individual names are shown in the margins of the heatmap and colors represent different populations.
#'@param x -- a genlight object 
#'@param plotheatmap -- a switch if a heatmap should be shown [Default:TRUE] 
#'@param ... parameters passed to function A.mat from package rrBLUP
#'@param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#'@return an identity by descent matrix and a heatmap plot
#'@export
#'@references 
#' Endelman, J. B. (2011). Ridge regression and other kernels for genomic selection with r package rrblup. The Plant Genome 4, 250.
#' Endelman, J. B. , Jannink, J.-L. (2012). Shrinkage estimation of the realized relationship matrix. G3: Genes, Genomics, Genetics 2, 1405.
#'@examples
#'gl.grm(bandicoot.gl[1:20,])  

gl.grm <- function(x, plotheatmap=TRUE, verbose=NULL, ...){

# CHECK IF PACKAGES ARE INSTALLED
  pkg <- "rrBLUP"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop("Package ",pkg," needed for this function to work. Please install it.") }   
  pkg <- "gplots"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop("Package ",pkg," needed for this function to work. Please install it.") }   
  
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

     # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nInd(x))
      pop(x) <- as.factor(pop(x))
    }
    
# DO THE JOB    
    
# function to replicate defaults colors of ggplot
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# assigning colors to populations
    colors_pop_temp <- gg_color_hue(length(levels(pop(x))))
    names(colors_pop_temp) <- as.character(levels(x$pop))
    cols_pops <- as.data.frame(cbind(names(colors_pop_temp),colors_pop_temp))
    colnames(cols_pops) <- c("pop","colour")
    df_pops <- as.data.frame(x$pop)
    colnames(df_pops) <- c("pop")
    cols_pops <- merge(df_pops,cols_pops,by="pop")

# calculating the realized additive relationship matrix
  
G <- rrBLUP::A.mat(as.matrix(x)-1, ...)

# creating color palette for probability of identity by descent
cool = rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
cols = c(rev(cool), rev(warm))
mypalette <- colorRampPalette(cols)(255)

if(plotheatmap==T){
# plotting heatmap 
gplots::heatmap.2(G,col=mypalette,dendrogram="column",ColSideColors=as.character(cols_pops$colour),RowSideColors=as.character(cols_pops$colour), trace = "none", density.info = "none",scale="none",main="Probability of identity by descent")
legend(0,0.8,legend=unique(cols_pops$pop),fill=unique(cols_pops$colour),cex=0.75,title="Populations")
}

# FLAG SCRIPT END

  if(verbose >= 1){
    cat("Completed:",funname,"\n")
  }  

  return (G)
}