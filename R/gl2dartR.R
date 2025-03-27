#' Convert a genlight object to a dartR object
#'
#' This function converts a `genlight` object into a `dartR` object by 
#' changing its class attribute.
#' It is used to convert legacy data sets to the new dartR format. 
#'
#' @param x An object of class `genlight` to be converted.
#' @param filename A character string specifying the name of the file to save 
#' the converted object. [default is gl.rds]
#' @param file.path A character string specifying the path to save the file.
#' @return The input object with class changed to `"dartR"` and its package attribute set to `"dartR.base"`.
#' @export
#' @examples
#' simgl <- glSim(10, 100, ploidy = 2, indnames=1:10, locnames=1:100)  # Simulating a genlight object
#' simgl <- gl2dartR(simgl)
#' pop(simgl)<- rep("A",10)
#' indNames(simgl) <- paste0("ind",1:10)
#' gl.smearplot(simgl, verbose=0)

gl2dartR <- function(x,  filename=NULL, file.path=tempdir()) {
  if (!is(x, "genlight")) {
    stop("Input must be a genlight object.")
  }
  
  if (is(x,"genlight"))
  {
    class(x) <- "dartR"
    attr(class(x), "package") <- NULL
  if (!is.null(filename)) {
    save_path <- file.path(file.path, filename)
    cat(report("Saving the converted object as a file under: ", save_path ))
    saveRDS(x, file=save_path)
  }
    return(x)
  }
}

