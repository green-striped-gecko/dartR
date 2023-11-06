
gl.ID.replicate <- function(x,
                            loc_threshold = 100,
                            perc_geno = 0.99
){
  
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
  col_same <- col_same[which(col_same$nloc > loc_threshold),]
  col_same <- col_same[which(col_same$perc > perc_geno),]
  
  col_same <- col_same[order(col_same$perc,decreasing = TRUE),]
  
  return(col_same)
  
}