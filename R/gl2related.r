#' Create a dataframe suitable for input to package \{related\}, from a genlight (SNP) \{adegenet\} object
#'
#' @param gl -- name of the genlight object containing the SNP data or a genind object containing presence absence data [required]
#' @return A dataframe suitable as input to package \{related\}
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' \dontrun{
#' #install package via
#' #install.packages("related", repos="http://R-Forge.R-project.org")
#' df <- gl2related(gl)
#' }

gl2related <- function(gl) {

  if(class(gl)!="genlight") {
    cat("Fatal Error: genlight object required for gl2related.pop.r!\n"); stop()
  }

# Convert the genlight data to a form expected by related
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
  x[is.na(x)] <- 0 # Related uses zero as missing

# Tidy up the locus names
  colnames(x) <- gsub("-","_",colnames(x), fixed=TRUE)
  colnames(x) <- gsub("/","",colnames(x), fixed=TRUE)
  colnames(x) <- gsub("|","_",colnames(x), fixed=TRUE)

#Add the individual names
  x <- cbind(as.character(row.names(x)), x)

# Convert to a dataframe suitable for input to package related
  df <- data.frame(x, stringsAsFactors = FALSE)

}

