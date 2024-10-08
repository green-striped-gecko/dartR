% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl.report.reproducibility.r
\name{gl.report.reproducibility}
\alias{gl.report.reproducibility}
\title{Reports summary of RepAvg (repeatability averaged over both alleles 
for
each locus) or reproducibility (repeatability of the scores for fragment
presence/absence)}
\usage{
gl.report.reproducibility(
  x,
  plot.out = TRUE,
  plot_theme = theme_dartR(),
  plot_colors = two_colors,
  save2tmp = FALSE,
  verbose = NULL
)
}
\arguments{
\item{x}{Name of the genlight object containing the SNP or presence/absence
(SilicoDArT) data [required].}

\item{plot.out}{If TRUE, displays a plot to guide the decision on a filter
threshold [default TRUE].}

\item{plot_theme}{Theme for the plot. See Details for options
[default theme_dartR()].}

\item{plot_colors}{List of two color names for the borders and fill of the
plots [default two_colors].}

\item{save2tmp}{If TRUE, saves any ggplots and listings to the session
temporary directory (tempdir) [default FALSE].}

\item{verbose}{Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
progress log; 3, progress and results summary; 5, full report
[default 2, unless specified using gl.set.verbosity].}
}
\value{
An unaltered genlight object
}
\description{
SNP datasets generated by DArT have an index, RepAvg, generated by
reproducing the data independently for 30% of loci. RepAvg is the proportion
of alleles that give a repeatable result, averaged over both alleles for each
locus.

In the case of fragment presence/absence data (SilicoDArT), repeatability is
the percentage of scores that are repeated in the technical replicate
 dataset.
}
\details{
The function displays a table of minimum, maximum, mean and quantiles for
 repeatbility against possible thresholds that might subsequently be
 specified in \code{\link{gl.filter.reproducibility}}.

 If plot.out=TRUE, display also includes a boxplot and a histogram to guide
 in the selection of a threshold for filtering on repeatability.

 If save2tmp=TRUE, ggplots and relevant tabulations are saved to the
 session's temp directory (tempdir)

 For examples of themes, see:
   \itemize{
 \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and 
 \item
  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
 }
}
\examples{
 \donttest{
# SNP data
  out <- gl.report.reproducibility(testset.gl)
  }
# Tag P/A data
  out <- gl.report.reproducibility(testset.gs)

}
\seealso{
\code{\link{gl.filter.reproducibility}}

Other report functions: 
\code{\link{gl.report.bases}()},
\code{\link{gl.report.callrate}()},
\code{\link{gl.report.diversity}()},
\code{\link{gl.report.hamming}()},
\code{\link{gl.report.heterozygosity}()},
\code{\link{gl.report.hwe}()},
\code{\link{gl.report.ld.map}()},
\code{\link{gl.report.locmetric}()},
\code{\link{gl.report.maf}()},
\code{\link{gl.report.monomorphs}()},
\code{\link{gl.report.overshoot}()},
\code{\link{gl.report.pa}()},
\code{\link{gl.report.parent.offspring}()},
\code{\link{gl.report.rdepth}()},
\code{\link{gl.report.secondaries}()},
\code{\link{gl.report.sexlinked}()},
\code{\link{gl.report.taglength}()}
}
\author{
Custodian: Arthur Georges -- Post to
\url{https://groups.google.com/d/forum/dartr}
}
\concept{report functions}
