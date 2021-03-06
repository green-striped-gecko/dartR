% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl.filter.monomorphs.r
\name{gl.filter.monomorphs}
\alias{gl.filter.monomorphs}
\title{Remove monomorphic loci, including those with all NAs}
\usage{
gl.filter.monomorphs(x, verbose = NULL)
}
\arguments{
\item{x}{-- name of the input genlight object [required]}

\item{verbose}{-- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2, unless specified using gl.set.verbosity]}
}
\value{
A genlight object with monomorphic ( and all NA) loci removed
}
\description{
This script deletes monomorphic loci from a genlight \{adegenet\} object
}
\details{
A DArT dataset will not have monomorphic loci, but they can arise, along with loci that are scored all NA, when populations or individuals are deleted.
Retaining monomorphic loci unnecessarily increases the size of the dataset and will affect some calculations.

Note that for SNP data, NAs likely represent null alleles; in tag presence/absence data, NAs represent missing values (presence/absence could not 
be reliably scored)
}
\examples{
# SNP data
  result <- gl.filter.monomorphs(testset.gl, verbose=3)
# Tag P/A data
  result <- gl.filter.monomorphs(testset.gs, verbose=3)
}
\author{
Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
}
