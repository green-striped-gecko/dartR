% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl.drop.pop.r
\name{gl.drop.pop}
\alias{gl.drop.pop}
\title{Removes specified populations from a genlight object}
\usage{
gl.drop.pop(
  x,
  pop.list,
  as.pop = NULL,
  recalc = FALSE,
  mono.rm = FALSE,
  verbose = NULL
)
}
\arguments{
\item{x}{Name of the genlight object containing SNP genotypes or Tag P/A data
(SilicoDArT) [required].}

\item{pop.list}{A list of populations to be removed [required].}

\item{as.pop}{Temporarily assign another metric to represent population for
the purposes of deletions [default NULL].}

\item{recalc}{Recalculate the locus metadata statistics [default FALSE].}

\item{mono.rm}{Remove monomorphic loci [default FALSE].}

\item{verbose}{Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
progress log; 3, progress and results summary; 5, full report
 [default 2 or as specified using gl.set.verbosity].}
}
\value{
A genlight object with the reduced data
}
\description{
Individuals are assigned to populations based on the specimen metadata
file (csv) used with \code{\link{gl.read.dart}}.
The script, having deleted populations, optionally identifies resultant
monomorphic loci or loci with all values missing and deletes them
(using gl.filter.monomorphs.r). The script also optionally
recalculates statistics made redundant by the deletion of individuals from
the dataset.

The script returns a genlight object with the new population assignments and
the recalculated locus metadata.

#' See more about data manipulation in the [tutorial](http://georges.biomatix.org/storage/app/media/uploaded-files/tutorial4dartrdatamanipulation22-dec-21-3.pdf).
}
\examples{
 # SNP data
   gl2 <- gl.drop.pop(testset.gl,
   pop.list=c('EmsubRopeMata','EmvicVictJasp'),verbose=3)
   gl2 <- gl.drop.pop(testset.gl, pop.list=c('EmsubRopeMata','EmvicVictJasp'),
   mono.rm=TRUE,recalc=TRUE)
   gl2 <- gl.drop.pop(testset.gl, pop.list=c('Male','Unknown'),as.pop='sex',verbose=3)
 # Tag P/A data
   gs2 <- gl.drop.pop(testset.gs, pop.list=c('EmsubRopeMata','EmvicVictJasp'))

}
\seealso{
\code{\link{gl.keep.pop}} to keep rather than drop specified
 populations
}
\author{
Custodian: Arthur Georges -- Post to
\url{https://groups.google.com/d/forum/dartr}
}
