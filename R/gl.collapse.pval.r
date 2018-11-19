#' Collapse a fixed distance matrix by amalgamating populations for which fixed differences are not significant
#'
#' This script takes the output from gl.collapse recursive and further collapses the fixed difference matrix
#' based on the pvalue associated with each comparison. The results are subsets of populations (OTUs) for which
#' diagnosability is non-significant. A recode table is generated applied to the genlight object to reflect 
#' the resultant OTUs. 
#' 
#' The recode table and final distance matrix are stored to disk as csv files. 
#'
#' @param fd -- name of the list of matricies output by gl.collapse.recursive run with test=TRUE [required]
#' @param prefix -- a string to be used as a prefix in generating the matrix of fixed differences (stored to disk) and the recode
#' table (also stored to disk) [default "fd_sig"]
#' @param delta Default 0.02
#' @param reps number of repetition. Default 1000.
#' @param alpha -- significance level for test of false positives [default 0.05]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A list containing the gl object x and the following square matricies
#'         [[1]] $gl -- the input genlight object;
#'         [[2]] $fd -- raw fixed differences;
#'         [[3]] $pcfd -- percent fixed differences;
#'         [[4]] $nobs -- mean no. of individuals used in each comparison;
#'         [[5]] $nloc -- total number of loci used in each comparison;
#'         [[6]] $expobs -- if test=TRUE, the expected count of false positives for each comparison [by simulation], otherwise NAs
#'         [[7]] $prob -- if test=TRUE, the significance of the count of fixed differences [by simulation], otherwise NAs
#' @export
#' @author Arthur Georges (bugs? Post to https://groups.google.com/d/forum/dartr)

gl.collapse.pval <- function(fd, prefix="fd.sig", delta=0.02, reps=1000, alpha=0.05, v=2) {
  
  if (v > 0) {
    cat("Starting gl.collapse.pval: Amalgamating populations with no statistically significant fixed differences\n")
  }
  if (!("fd" %in% names(fd)) ||
      !("pcfd" %in% names(fd)) ||
      !("nobs" %in% names(fd)) ||
      !("nloc" %in% names(fd)) ||
      !("pval" %in% names(fd)) ||
      !("expobs" %in% names(fd))) {
    cat("Fatal Error: fd must be a list produced by gl.collapse.recursive\n"); stop("Execution terminated\n")
  }
  if (fd$pval[1,1] == -1) {
    cat("Fatal Error: gl.collapse.recursive needs to be run with test set to TRUE\n"); stop("Execution terminated\n")
  }
  # Store the number of populations in the matrix
  npops <- dim(fd$pval)[1]
  
  # Extract the column names
  pops <- variable.names(fd$pval)
  
  # Initialize a list to hold the populations that differ by <= tpop
  zero.list <- list()

  # For each pair of populations
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
  if (length(zero.list) > 1) {
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
  
  # Print out the results of the aggregations 
  if(v > 1) {cat("\n\nPOPULATION GROUPINGS\n")}
  
  for (i in 1:length(zero.list)) {
    # Create a group label
    if (length(zero.list[[i]])==1) {
      replacement <- zero.list[[i]][1]
    } else {
      replacement <- paste0(zero.list[[i]][1],"+")
    }
    if(v > 1) {
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
    recode.table <- paste0(prefix,"_recode_sig.csv")
    write.table(df, file=recode.table, sep=",", row.names=FALSE, col.names=FALSE)
  
  # Recode the data file (genlight object)
  x2 <- gl.recode.pop(fd$gl, pop.recode=recode.table, v=v)
  fd2 <- gl.fixed.diff(x2,test=TRUE,delta=delta,reps=reps,v=v)
  
  # Return the matricies
  if (v > 1) {
    cat("Returning a list containing the following square matricies:\n",
        "         [[1]] $gl -- input genlight object;\n",
        "         [[2]] $fd -- raw fixed differences;\n",
        "         [[3]] $pcfd -- percent fixed differences;\n",
        "         [[4]] $nobs -- mean no. of individuals used in each comparison;\n",
        "         [[5]] $nloc -- total number of loci used in each comparison;\n",
        "         [[6]] $expobs -- if test=TRUE, the expected count of false positives for each comparison [by simulation]\n",
        "         [[7]] $prob -- if test=TRUE, the significance of the count of fixed differences [by simulation]\n")
  }

  if(setequal(levels(pop(x2)),levels(pop(fd$gl)))) { 
    if (v > 1) {cat(paste("\nPOPULATION GROUPINGS\n     No populations collapsed at alpha >=", alpha,"\n"))}
    if (v > 0) {
      cat("Completed gl.collapse.pval\n")
    }
    l <- list(gl=fd$gl,fd=fd$fd,pcfd=fd$pcfd,nobs=fd$nobs,nloc=fd$nloc,expobs=fd$expobs,pval=fd$pval)
    return(l)
  } else {
    if (v > 1) {
      cat("\nPOPULATION GROUPINGS\n")
      print(table(pop(x2)))
    }
    if (v > 0) {
      cat("Completed gl.collapse.pval\n\n")
    }
    l <- list(gl=x2,fd=fd2$fd,pcfd=fd2$pcfd,nobs=fd2$nobs,nloc=fd2$nloc,expobs=fd2$expobs,pval=fd2$pval)
    return(l)
  }
}
  