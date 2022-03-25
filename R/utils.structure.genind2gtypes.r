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
