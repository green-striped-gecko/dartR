#' @name gl.report.replicates
#' @title Identify replicated individuals 
#' @description
#' Identify replicated individuals 
#' @param x Name of the genlight object containing the SNP data [required].
#' @param loc_threshold Minimum number of loci required to asses that two 
#' individuals are replicates [default 100].
#' @param perc_geno Mimimum percentage of genotypes in which two individuals 
#' should be the same [default 0.99]. 
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot_theme User specified theme [default theme_dartR()].
#' @param plot_colors Vector with two color names for the borders and fill
#' [default two_colors].
#' @param bins Number of bins to display in histograms [default 100].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @details
#' This function uses an C++ implementation, so package Rcpp needs to be 
#' installed and it is therefore fast (once it has compiled the function after 
#' the first run).
#' 
#' Ideally, in a large dataset with related and unrelated individuals and 
#' several replicated individuals, such as in a capture/mark/recapture study, 
#' the first histogram should have four "peaks". The first peak should represent
#'  unrelated individuals, the second peak should correspond to second-degree 
#'  relationships (such as cousins), the third peak should represent 
#'  first-degree relationships (like parent/offspring and full siblings), and
#'   the fourth peak should represent replicated individuals. 
#'   
#' In order to ensure that replicated individuals are properly identified, it's
#'  important to have a clear separation between the third and fourth peaks in 
#'  the second histogram. This means that there should be bins with zero counts 
#'  between these two peaks.
#' @return A list with three elements:
#'\itemize{
#'\item table.rep: A dataframe with pairwise results of percentage of same 
#'genotypes between two individuals, the number of loci used in the comparison 
#'and the missing data for each individual.
#'\item ind.list.drop: A vector of replicated individuals to be dropped.
#' Replicated individual with the least missing data is reported.
#'\item ind.list.rep: A list of of each individual that has replicates in the 
#'dataset, the name of the replicates and the percentage of the same genotype.
#'  }
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' res_rep <- gl.report.replicates(platypus.gl, loc_threshold = 500, 
#' perc_geno = 0.85)
#' @family report functions
#' @export

gl.report.replicates <- function(x,
                            loc_threshold = 100,
                            perc_geno = 0.99,
                            plot.out = TRUE,
                            plot_theme = theme_dartR(),
                            plot_colors = two_colors,
                            bins = 100,
                            verbose = NULL
){
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jody",
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # DO THE JOB
  xx <- as.matrix(x)
  
  # number of same genotypes pairwise
  #to hack package checking...
  SameGeno <- function() {}  
  Rcpp::cppFunction(
    "NumericMatrix SameGeno(NumericMatrix x) {
  int nrow = x.nrow();
  NumericMatrix out(nrow,nrow);
  for (int i=0; i<(nrow-1); i++) {
     for (int j=(i+1); j<nrow; j++) {
     out(j,i) = sum(na_omit(x(i,_)==x(j,_)) );
     }
  }
  return out;
}"
  )
  
  SameGeno_tmp <- SameGeno(xx)
  # number of no NAs in either or both genotypes pairwise
  #to hack package checking...
  SameGenoNA <- function() {}  
  Rcpp::cppFunction(
    "NumericMatrix SameGenoNA(NumericMatrix x) {
  int nrow = x.nrow();
  NumericMatrix out(nrow,nrow);
  for (int i=0; i<(nrow-1); i++) {
     for (int j=(i+1); j<nrow; j++) {
     out(j,i) = sum( !is_na( x(i,_) == x(j,_) ));
     }
  }
  return out;
}"
  )
  
  SameGenoNA_tmp <- SameGenoNA(xx)
  colnames(SameGenoNA_tmp) <- indNames(x)
  rownames(SameGenoNA_tmp) <- indNames(x)
  
  mat_same <- SameGeno_tmp/SameGenoNA_tmp
  
  colnames(mat_same) <- indNames(x)
  rownames(mat_same) <- indNames(x)
  
  col_same <- as.data.frame(as.table(mat_same))
  colnames(col_same) <- c("ind1","ind2","perc")
  col_noNas <- as.data.frame(as.table(SameGenoNA_tmp))
  
  col_same$nloc <- col_noNas$Freq
  
  col_same <- col_same[complete.cases(col_same$perc),]
  # holding complete dataset for plotting 
  col_same_hold <- col_same
  
  col_same <- col_same[which(col_same$nloc > loc_threshold),]
  col_same <- col_same[which(col_same$perc > perc_geno),]
  
  col_same <- col_same[order(col_same$perc,decreasing = TRUE),]
  
  unique_ind <- col_same
  unique_ind$ind1 <- as.character(unique_ind$ind1)
  unique_ind$ind2 <- as.character(unique_ind$ind2)
  
  ind_NA <- 1 - rowSums(is.na(as.matrix(x))) / nLoc(x)
  ind1_NA <- data.frame(ind1 = indNames(x), ind1_NA = ind_NA)
  ind2_NA <- data.frame(ind2 = indNames(x), ind2_NA = ind_NA)

  unique_ind <- merge(unique_ind,ind1_NA,by="ind1")
  unique_ind <- merge(unique_ind,ind2_NA,by="ind2")
  unique_ind <- cbind(ind1=unique_ind[,2],unique_ind[,-2])
  
  # getting vector of individuals (replicates) to drop based on missing data
  
  unique_ind_tmp <- unique_ind
  unique_ind_tmp$test_stat <- unique_ind_tmp$ind1_NA > unique_ind_tmp$ind2_NA
  
  ind_list <- NULL
  
    for (y in 1:nrow(unique_ind_tmp)) {
      if (unique_ind_tmp[y, "test_stat"] == TRUE) {
        ind_tmp <- unique_ind_tmp[y, "ind1"]
      } else{
        ind_tmp <- unique_ind_tmp[y, "ind2"]
      }
      
      if (ind_tmp %in% ind_list) {
        next
      } else{
        ind_list <- c(ind_list, ind_tmp)
      }
      
    }
  
  # getting list of replicated individuals
  
  ind_mat <- as.matrix(reshape2::acast(unique_ind, ind1~ind2, value.var="perc"))
  
  ind_list_rep <- apply(ind_mat,2,function(x){
    x[which(!is.na(x))]
  })
  
  # Histograms
  p1_col_same_hold <- col_same_hold
  p1 <-
    ggplot(p1_col_same_hold, aes(x = perc)) + 
    geom_histogram(bins = bins, color = plot_colors[1],fill = plot_colors[2]) +
    xlab(" ") +
    ylab("Count") +
    plot_theme
  
  p2_col_same_hold <- col_same_hold[which(col_same_hold$perc>0.8),]
  p2 <-
    ggplot(p2_col_same_hold, aes(x = perc)) + 
    geom_histogram(bins = bins, color = plot_colors[1],fill = plot_colors[2]) +
    xlab("Percentage of same genotype") +
    ylab("Count") +
    plot_theme
  
  # PRINTING OUTPUTS
  if (plot.out) {
    # using package patchwork
    p3 <- (p1 / p2) 
    print(p3)
  }
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # RETURN
  
  return(list(table.rep = unique_ind,
              ind.list.drop = ind_list,
              ind.list.rep = ind_list_rep))
  
}