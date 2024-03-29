% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl2snapp.r
\name{gl2snapp}
\alias{gl2snapp}
\title{Converts a genlight object to nexus format suitable for phylogenetic analysis
 by SNAPP (via BEAUti)}
\usage{
gl2snapp(x, outfile = "snapp.nex", outpath = tempdir(), verbose = NULL)
}
\arguments{
\item{x}{Name of the genlight object containing the SNP data [required].}

\item{outfile}{File name of the output file (including extension)
[default "snapp.nex"].}

\item{outpath}{Path where to save the output file
[default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath='.'
when calling this function to direct output files to your working directory.}

\item{verbose}{Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
progress log; 3, progress and results summary; 5, full report
[default 2 or as specified using gl.set.verbosity].}
}
\value{
returns no value (i.e. NULL)
}
\description{
The output nexus file contains the SNP data and relevant PAUP command lines
suitable for BEAUti.
}
\examples{
gl2snapp(testset.gl)
}
\references{
Bryant, D., Bouckaert, R., Felsenstein, J., Rosenberg, N.A. and
RoyChoudhury, A. (2012). Inferring species trees directly from biallelic
genetic markers: bypassing gene trees in a full coalescent analysis.
 Molecular Biology and Evolution 29:1917-1932.
}
\author{
Custodian: Arthur Georges (Post to
\url{https://groups.google.com/d/forum/dartr})
}
