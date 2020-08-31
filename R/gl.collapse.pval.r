#' Collapse a fixed distance matrix by amalgamating populations for which pairwise fixed differences are not significant
#'
#' This script takes the output from gl.collapse and further collapses the fixed difference matrix
#' based on the pvalue associated with each comparison. The results are subsets of populations (OTUs) for which
#' diagnosability is demonstrated in the sample set, but non-significant. 
#' 
#' @param fd -- name of the list containing the collapsed gl object and associated distance matricies output by gl.collapse run with test=TRUE [required]
#' @param recode.table -- name of the new recode.table to receive the new population reassignments 
#' arising from the amalgamation of populations [tmp.csv]
#' @param outpath -- path where to save the output file [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath="." when calling this function to direct output files to your working directory.
#' @param delta -- threshold for the level of difference between two populations that will be regarded as operationally fixed [Default 0.02]
#' @param reps number of repetitions in the simulations to estimate false positives. [Default 1000].
#' @param alpha -- significance level for test of false positives [default 0.05]
#' @param plot -- if TRUE, plot a PCoA with the new groupings [default FALSE]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return A list containing the gl object with the new collapsed populations and the following square matricies
#'         [[1]] $gl -- the input genlight object;
#'         [[2]] $fd -- raw fixed differences;
#'         [[3]] $pcfd -- percent fixed differences;
#'         [[4]] $nobs -- mean no. of individuals used in each comparison;
#'         [[5]] $nloc -- total number of loci used in each comparison;
#'         [[6]] $expobs -- the expected count of false positives for each comparison [by simulation], otherwise NAs
#'         [[7]] $prob -- the significance of the count of fixed differences [by simulation]. These should all be significant (< alpha)
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})

gl.collapse.pval <- function(fd, 
                             recode.table="tmp.csv",
                             outpath=tempdir(),
                             delta=0.02, 
                             reps=1000, 
                             alpha=0.05, 
                             plot=FALSE,
                             verbose=NULL) {
  
# TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  
# SET VERBOSITY
  
  if (is.null(verbose)){ 
    
      verbose <- 2
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
  
  stop("This script has been discontinued. Requires subjective judgement.\n  Amalgamate non-significant pairs manually using gl.fixed.diff with test=TRUE and gl.merge.pop\n")

  if ( verbose >= 2){
    cat("  Amalgamating populations based on non-significant fixed differences\n")
  }
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
  # Checking count of populations
  if(nPop(fd$gl) < 2) {
    stop("Fatal Error: Distance calculation requires at least two populations, one or none present\n")
  }

  if (!("fd" %in% names(fd)) ||
      !("pcfd" %in% names(fd)) ||
      !("nobs" %in% names(fd)) ||
      !("nloc" %in% names(fd)) ||
      !("pval" %in% names(fd)) ||
      !("expobs" %in% names(fd))) {
    stop("Fatal Error: fd must be a list produced by gl.collapse\n")
  }

  # DO THE JOB      
  
  if (is.na(fd$pval[1,1])) {
    cat("\nWarning: gl.collapse needed to be run with test set to TRUE, running now\n")
    fd <- gl.fixed.diff(fd$gl,tloc=0,test=TRUE,delta=delta,reps=reps,mono.rm=TRUE,pb=FALSE, verbose=verbose)
    cat("\n")
  } 
  
  # Print out the results of the prior aggregation, and the pvalue matrix 
  if(verbose >= 2) {
    
      cat("\n\nInitial Populations\n",variable.names(fd$pval),"\n")
    
      cat("\nRaw fixed differences for the collapsed matrix\n")
      print(fd$fd)
      cat("\n")
      
      cat("\nP-values for the collapsed matrix\n")
      print(fd$pval)
      cat("\n")
      
      cat("Sample sizes")
      print(table(pop(fd$gl)))
      cat("\n")
      
  }
  
  # Check for non-significant pairwise differences
  
  if (all(fd$pval <= alpha)){
    
    if (verbose >= 2){
      cat("No amalgamation of populations at alpha >=",alpha,"-- all differences significant\n")
      cat("  Analysis complete\n")
    }
    
    l <- list(gl=fd$gl,fd=fd$fd,pcfd=fd$pcfd,nobs=fd$nobs,nloc=fd$nloc,expobs=fd$expobs,pval=fd$pval)

  } else {
  # Collapse the matrix
    
  # Store the number of populations in the matrix
  npops <- dim(fd$pval)[1]
  
  # Extract the column names
  pops <- variable.names(fd$pval)
  
  # Initialize a list to hold the populations that do not differ significantly (>= alpha)
  zero.list <- list()

  # For each pair of populations
  if (verbose >= 2){
    cat("Identifying and amalgamating pairs of populations that are not significantly different\n")
  }  
    for(i in 1:npops){
      zero.list[[i]] <- c(rownames(fd$pval)[i])
      for (j in 1:npops) {
        if (fd$pval[i,j] >= alpha) {
          zero.list[[i]] <- c(zero.list[[i]],rownames(fd$pval)[i],rownames(fd$pval)[j])
          zero.list[[i]] <- unique(zero.list[[i]])
        }
      }
      zero.list[[i]] <- sort(zero.list[[i]])
    }
    # Pull out the unique aggregations 
      zero.list <- unique(zero.list)
  
  # Amalgamate populations
 
  if (length(zero.list) >= 2) {
    for (i in 1:(length(zero.list)-1)) {
      for (j in 2:length(zero.list)) {
        if (length(intersect(zero.list[[i]],zero.list[[j]])) > 0 ) {
          zero.list[[i]] <- union(zero.list[[i]],zero.list[[j]])
          zero.list[[j]] <- union(zero.list[[i]],zero.list[[j]])
        }
      }
    }
    for (i in 1:length(zero.list)) {
      zero.list <- unique(zero.list)
    }
  }  
  
  for (i in 1:length(zero.list)) {
    # Create a group label
    if (length(zero.list[[i]])==1) {
      replacement <- zero.list[[i]][1]
    } else {
      replacement <- paste0(zero.list[[i]][1],"+")
    }
    if(verbose >= 2) {
      cat(paste0("Group:",replacement,"\n"))
      print(as.character(zero.list[[i]]))
      cat("\n")
    }
    # Create a dataframe with the pop names and their new group names  
    if (i==1) {
      df <- rbind(data.frame(zero.list[[i]],replacement, stringsAsFactors = FALSE))
    } else {
      df <- rbind(df,data.frame(zero.list[[i]],replacement, stringsAsFactors = FALSE))
    }
  }
  
  # Create a recode table corresponding to the aggregations
    recode.table <- "tmp.csv"
    outfilespec <- file.path(outpath,recode.table)
    write.table(df, file=outfilespec, sep=",", row.names=FALSE, col.names=FALSE)
  
  # Recode the data file (genlight object)
    x2 <- gl.recode.pop(fd$gl, pop.recode=outfilespec, verbose=0)
    if (verbose >= 2){
      cat("Recalculating fixed difference matrix for collapsed populations\n\n")
    }
    fd2 <- gl.fixed.diff(x2,tloc=0,test=TRUE,delta=delta,reps=reps,mono.rm=TRUE,pb=FALSE, verbose=verbose)
  
  # Display the fd matrix
    if (verbose >= 3) {
      
      cat("\n\nRaw fixed differences for the collapsed matrix\n")
      print(fd2$fd)
      cat("\n")

      cat("\nP-values for the collapsed matrix\n")
      print(fd2$pval)
      cat("\n")
      
      cat("Sample sizes")
      print(table(pop(x2)))
      cat("\n")
      
      if (any(fd2$pval > alpha)){
        cat("  Warning: Some resulting pairwise differences between populations are now non-signicant\n")
        cat("  Consider running gl.collapse.pval again.\n\n")
      }
    }
    l <- list(gl=x2,fd=fd2$fd,pcfd=fd2$pcfd,nobs=fd2$nobs,nloc=fd2$nloc,expobs=fd2$expobs,pval=fd2$pval)
  }
    
    # Plot the results  
    if (plot){
      pcoa <- gl.pcoa(x2,verbose=verbose)
      tmp <- gl.pcoa.plot(pcoa,x2)
      show(tmp)
    }  
    
    # Return the matricies
    if (verbose >= 4) {
      cat("Returning a list containing the new genlight object and square matricies, as follows:\n",
          "         [[1]] $gl -- input genlight object;\n",
          "         [[2]] $fd -- raw fixed differences;\n",
          "         [[3]] $pcfd -- percent fixed differences;\n",
          "         [[4]] $nobs -- mean no. of individuals used in each comparison;\n",
          "         [[5]] $nloc -- total number of loci used in each comparison;\n",
          "         [[6]] $expobs -- if test=TRUE, the expected count of false positives for each comparison [by simulation]\n",
          "         [[7]] $prob -- if test=TRUE, the significance of the count of fixed differences [by simulation]\n")
    }

  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat("Completed:",funname,"\n")
  }  
  
    return(l)
}