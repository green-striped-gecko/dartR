#' Generate percentage allele frequencies by locus and population
#'
#' This is a support script, to take SNP data or SilocoDArT presence/absence data grouped into populations in a genelight object \{adegenet\}
#' and generate a table of allele frequencies for each population and locus
#'
#' @param x -- name of the genlight object containing the SNP or Tag P/A (SilicoDArT) data [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return A matrix with allele (SNP data) or presence/absence frequencies (Tag P/A data) broken down by population and locus
#' @export
# #' @importFrom plyr
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' m <-  gl.percent.freq(testset.gl)
#' m

gl.percent.freq<- function(x, verbose=NULL) {

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
    if (verbose >= 2){cat("  Processing Presence/Absence (SilicoDArT) data\n")}
    data.type <- "SilicoDArT"
  } else if (all(x@ploidy == 2)){
    if (verbose >= 2){cat("  Processing a SNP dataset\n")}
    data.type <- "SNP"
  } else {
    stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
  }
    
# SCRIPT SPECIFIC ERROR CHECKING
    
  # Checking for and removing monomorphic loci
    if(!(x@other$loc.metrics.flags$monomorphs==TRUE)){
      if (verbose >= 1) {cat("Warning: Monomorphic loci retained, used in calculations\n")}
    } 
# DO THE JOB
  x2 <- seppop(x)
  x2_list <- lapply(x2,as.matrix)

  if(data.type=="SilicoDArT"){

    if (verbose >= 2) {
      cat("Starting gl.percent.freq: Calculating Tag P/A frequencies for populations\n")
    }
    # Treat SilicoDArT as biallelic, no heterozygotes
    x2_list <- lapply(x2_list,function(x){
      x[x==1] <- 2
      return(x)
    })

  } else {

    if (verbose >= 1) {
      cat("Starting gl.percent.freq: Calculating allele frequencies for populations\n")
    }
    if (verbose >= 3){
      cat("  This may take some time -- be patient\n")
    }
  }

  loc_names <- lapply(x2_list,colnames)
  
  nmissing_temp <- lapply(x2_list, is.na)
  nmissing <- lapply(nmissing_temp, colSums)
  
  n_temp <- lapply(x2_list,nrow)
  n <- lapply(n_temp,rep,nLoc(x))
  
  nobs_temp <- lapply(nmissing,unname)
  nobs <- Map("-",n,nobs_temp)
  
  sum_res <- lapply(x2_list, colSums,na.rm=T)

  f <- lapply(x2_list,colMeans, na.rm = T)
  f <- lapply(f, "/", 2)
  f <- lapply(f, "*", 100)
  f <- lapply(f, "round", 2)

  m <- Map(cbind,names(sum_res),loc_names,sum_res,nobs,nmissing,f,n)
  m <- lapply(m,cbind,1:nLoc(x))
  m <- lapply(m,as.data.frame)
  m <- rbind.fill(m)

  colnames(m) <- c("popn","locus","sum","nobs","nmissing","frequency","n","loc_order")
  
  m$popn <- as.factor(m$popn)
  m$locus <- as.factor(m$locus)
  m$sum <- as.numeric(m$sum)
  m$nobs <- as.numeric(m$nobs)
  m$nmissing <- as.numeric(m$nmissing)
  m$frequency <- as.numeric(m$frequency)
  m$n <- as.numeric(m$n)
  m$loc_order <- as.numeric(m$loc_order)
  
  m <- m[order(m$loc_order,m$popn),]
  m <- m[,-ncol(m)]
  
  rownames(m) <- NULL

# FLAG SCRIPT END

  if (verbose >= 1) {
    cat("Completed:",funname,"\n")
  }
  
  return(m)
  
  }

