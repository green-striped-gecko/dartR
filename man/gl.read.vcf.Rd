% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl.read.vcf.r
\name{gl.read.vcf}
\alias{gl.read.vcf}
\title{Converts a vcf file into a genlight object}
\usage{
gl.read.vcf(vcffile, verbose = NULL)
}
\arguments{
\item{vcffile}{A vcf file (works only for diploid data) [required].}

\item{verbose}{Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
progress log; 3, progress and results summary; 5, full report
[default 2, unless specified using gl.set.verbosity].}
}
\value{
A genlight object.
}
\description{
This function needs package vcfR, please install it. The converted genlight
object does not have individual metrics. You need to add them 'manually' to
the other$ind.metrics slot.
}
\examples{
\dontrun{
obj <- gl.read.vcf(system.file('extdata/test.vcf', package='dartR'))
}
}
\author{
Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
}
