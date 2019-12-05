#' Create a dataframe suitable for input to package \{Demerelate\} from a genlight \{adegenet\} object
#'
#' @param gl -- name of the genlight object containing the SNP data [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A dataframe suitable as input to package \{Demerelate\}
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' df <- gl2demerelate(testset.gl)

gl2demerelate <- function(gl, verbose=2) {

# ERROR CHECKING
  
  if(!is(gl, "genlight")) {
    cat("Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }

  if (verbose < 0 | verbose > 5){
    cat("    Warning: Parameter verbose must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }

# FLAG SCRIPT START

  if (verbose >= 1) {
    cat("Starting gl2demelerate: Outputing data to demelerate format\n")
  }
  
# DO THE JOB
  

# Convert the genlight data to a form expected by demerelate
  # Strategy is to create two identical genlight objects, converted to matricies, then
  # to have one hold one allelic state for each locus and the other to hold the alternate
  # allelic state for each locus. The locus names are relabled _1 and _2 in each matrix
  # and then the two matricies are concatenated side by side. The resultant matrix is ordered
  # on locus to bring the two allelic states for a locus back to adjacency. Format tidying up
  # and bob's your uncle.

  x1 <- as.matrix(gl) # gl is a genlight object
  x2 <- as.matrix(gl)

  x1[x1 == 2] <- 2 # homozygote alternate
  x1[x1 == 1] <- 1 # heterozygote
  x1[x1 == 0] <- 1 # homozygote reference
  colnames(x1) <- gsub(" ","",paste(colnames(x1),"_1"))

  x2[x2 == 2] <- 2 # homozygote alternate
  x2[x2 == 1] <- 2 # heterozygote
  x2[x2 == 0] <- 1 # homozygote reference
  colnames(x2) <- gsub(" ","",paste(colnames(x2),"_2"))

  x <- cbind(x1,x2)
  x <- x[,order(colnames(x))]
  #x[is.na(x)] <- 0 # Related uses zero as missing

# Tidy up the locus names
  colnames(x) <- gsub("-","_",colnames(x), fixed=TRUE)
  colnames(x) <- gsub("/","",colnames(x), fixed=TRUE)
  colnames(x) <- gsub("|","_",colnames(x), fixed=TRUE)

#Add the individual names
  
  df <- cbind.data.frame(row.names(x), factor(pop(gl)), x, stringsAsFactors=FALSE)
  df[,1] <- as.character(df[,1])

# Convert to a dataframe suitable for input to package demerelate
  #df <- data.frame(x, stringsAsFactors = FALSE)
  names(df)[1] <- "Sample-ID"
  names(df)[2] <- "Population"
#  df[,2] <- factor(df[,2])

  if (verbose >= 1) {
    cat("Completed: gl2demelerate\n")
  }

  return(df)
}

