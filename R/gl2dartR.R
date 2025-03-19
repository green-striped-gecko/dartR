#' Convert a genlight object to a dartR object
#'
#' This function converts a `genlight` object into a `dartR` object by 
#' changing its class attribute.
#' It is used to convert legacy data sets to the new dartR format. 
#'
#' @param x An object of class `genlight` to be converted.
#' @param save Logical. If `TRUE`, the converted object is saved as an file.
#' 
#' @return The input object with class changed to `"dartR"` and its package attribute set to `"dartR.base"`.
#' @export
#' @examples
#' gl <- glSim(10, 100, ploidy = 2, ind.names=1:10, loc.names=1:100)  # Simulating a genlight object
#' gl <- gl2dartR(gl)
#' pop(gl)<- rep("A",10)

gl2dartR <- function(x,  filename=NULL, file.path=getwd()) {
  if (!is(x, "genlight")) {
    stop("Input must be a genlight object.")
  }
  
  if (is(x,"genlight"))
  {
    class(x) <- "dartR"
    attr(class(x), "package") <- "dartR.base"
  if (!is.null(filename)) {
    save_path <- file.path(file.path, filename)
    cat(report("Saving the converted object as a file under: ", save_path ))
    saveRDS(x, file=save_path)
  }
    return(x)
  }
}

