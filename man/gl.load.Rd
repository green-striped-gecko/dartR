% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl.load.r
\name{gl.load}
\alias{gl.load}
\title{Loads an object from compressed binary format produced by gl.save()}
\usage{
gl.load(file, verbose = NULL)
}
\arguments{
\item{file}{Name of the file to receive the binary version of the object
[required].}

\item{verbose}{Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
progress log; 3, progress and results summary; 5, full report
[default 2 or as specified using gl.set.verbosity].}
}
\value{
The loaded object
}
\description{
This is a wrapper for readRDS()
}
\details{
The script loads the object from the current workspace and returns the
 gl object.
}
\examples{
gl.save(testset.gl,file.path(tempdir(),'testset.rds'))
gl <- gl.load(file.path(tempdir(),'testset.rds'))
}
\seealso{
\code{\link{gl.save}}
}
\author{
Custodian: Arthur Georges -- Post to
\url{https://groups.google.com/d/forum/dartr}
}
