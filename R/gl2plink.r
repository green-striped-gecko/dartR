#' Converts a genlight object to plink file format
#'
#' This function exports a genlight object into plink format and save it into a file
#' @param x -- genlight object
#' @param outfile -- name (path) of the output plink file [default plink.csv]
#' @param outpath -- path of the output file. Default is to tempdir(). If to be saved in the current working directory change to "."
#' @export
#' @author Bernd Guber (glbugs@@aerg.canberra.edu.au)
#' @examples
#' gl2plink(testset.gl)


gl2plink <- function(x, outfile="plink.csv", outpath=tempdir()) {
tt <- as.matrix(x)
write.csv(tt, file=file.path(outpath, outfile), row.names = TRUE, na = "-9")
}












