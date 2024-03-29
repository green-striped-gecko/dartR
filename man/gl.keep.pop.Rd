% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl.keep.pop.r
\name{gl.keep.pop}
\alias{gl.keep.pop}
\title{Removes all but the specified populations from a genlight object}
\usage{
gl.keep.pop(
  x,
  pop.list,
  as.pop = NULL,
  recalc = FALSE,
  mono.rm = FALSE,
  verbose = NULL
)
}
\arguments{
\item{x}{Name of the genlight object containing the SNP or presence/absence
(SilicoDArT) data [required].}

\item{pop.list}{A list of populations to be kept [required].}

\item{as.pop}{Assign another metric to represent population [default NULL].}

\item{recalc}{Recalculate the locus metadata statistics [default FALSE].}

\item{mono.rm}{Remove monomorphic loci [default FALSE].}

\item{verbose}{Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
progress log; 3, progress and results summary; 5, full report
[default 2 or as specified using gl.set.verbosity]}
}
\value{
A genlight object with the reduced data
}
\description{
Individuals are assigned to populations based on the specimen metadata data
file (csv) used with gl.read.dart().

The script, having deleted the specified populations, optionally identifies
resultant monomorphic loci or loci with all values missing and deletes them
(using gl.filter.monomorphs.r). The script also optionally recalculates
statistics made redundant by the deletion of individuals from the dataset.

The script returns a genlight object with the new population assignments and
the recalculated locus metadata.

#' See more about data manipulation in the [tutorial](http://georges.biomatix.org/storage/app/media/uploaded-files/tutorial4dartrdatamanipulation22-dec-21-3.pdf).
}
\examples{
 # SNP data
   gl2 <- gl.keep.pop(testset.gl, pop.list=c('EmsubRopeMata', 'EmvicVictJasp'))
   gl2 <- gl.keep.pop(testset.gl, pop.list=c('EmsubRopeMata', 'EmvicVictJasp'),
   mono.rm=TRUE,recalc=TRUE)
   gl2 <- gl.keep.pop(testset.gl, pop.list=c('Female'),as.pop='sex')
 # Tag P/A data
   gs2 <- gl.keep.pop(testset.gs, pop.list=c('EmsubRopeMata','EmvicVictJasp'))

}
\seealso{
\code{\link{gl.drop.pop}} to drop rather than keep specified populations
}
\author{
Custodian: Arthur Georges -- Post to
\url{https://groups.google.com/d/forum/dartr}
}
