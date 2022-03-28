#' structure util functions
#' 
#' These functions were copied from package strataG, which is no longer on CRAN (maintained by Eric Archer)
#' @export
#' @author Bernd Gruber (bugs? Post to
#'  \url{https://groups.google.com/d/forum/dartr}); original implementation of
#'   Eric Archer \url{https://github.com/EricArcher/strataG}
#' @param x a genind object
#' @return a gtypes object 


utils.structure.genind2gtypes <- 
function (x) 
{
##################################################
  df2gtypes <- function (x, ploidy, id.col = 1, strata.col = 2, loc.col = 3, 
            sequences = NULL, schemes = NULL, description = NULL, other = NULL) 
  {
    if (!(is.matrix(x) | is.data.frame(x))) {
      stop("'x' must be a matrix or data.frame")
    }
    x <- as.data.frame(x)
    if (!is.null(id.col)) {
      if (is.character(id.col)) 
        id.col <- match(id.col, colnames(x))
      id.col <- as.numeric(id.col)
    }
    if (!is.null(strata.col)) {
      if (is.character(strata.col)) 
        strata.col <- match(strata.col, colnames(x))
      strata.col <- as.numeric(strata.col)
    }
    loc.col <- as.numeric(loc.col)
    if (loc.col < max(id.col, strata.col, loc.col)) {
      stop("'loc.col' must be greater than 'id.col' and 'strata.col'")
    }
    ind.names <- if (is.null(id.col)) {
      if (is.null(rownames(x))) 
        1:nrow(x)
      else rownames(x)
    }
    else {
      x[, id.col]
    }
    ind.names <- as.character(ind.names)
    strata <- if (is.null(strata.col)) 
      NULL
    else x[, strata.col]
    gen.data <- x[, loc.col:ncol(x), drop = FALSE]
    methods::new("gtypes", gen.data = gen.data, ploidy = ploidy, 
                 ind.names = ind.names, strata = strata, schemes = schemes, 
                 sequences = sequences, description = description, other = if (is.null(other)) 
                   list()
                 else other)
  }
##################################################  
  
  
  
  gen.mat <- adegenet::genind2df(x, usepop = TRUE, oneColPerAll = TRUE)
  gen.mat[gen.mat == "NA"] <- NA
  has.pop <- !is.null(x@pop)
  df2gtypes(x = gen.mat, ploidy = x@ploidy[1], id.col = NULL, 
            strata.col = if (has.pop) 
              1
            else NULL, loc.col = if (has.pop) 
              2
            else 1, schemes = x@strata, other = list(genind = adegenet::other(x)))
}
