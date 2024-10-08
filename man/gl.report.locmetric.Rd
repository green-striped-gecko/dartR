% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl.report.locmetric.r
\name{gl.report.locmetric}
\alias{gl.report.locmetric}
\title{Reports summary of the slot $other$loc.metrics}
\usage{
gl.report.locmetric(
  x,
  metric,
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

\item{metric}{Name of the metric to be used for filtering [required].}

\item{plot.out}{Specify if plot is to be produced [default TRUE].}

\item{plot_theme}{Theme for the plot. See Details for options
[default theme_dartR()].}

\item{plot_colors}{List of two color names for the borders and fill of the
plots [default two_colors].}

\item{save2tmp}{If TRUE, saves any ggplots and listings to the session
temporary directory (tempdir) [default FALSE].}

\item{verbose}{Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
progress log; 3, progress and results summary; 5, full report
[default NULL, unless specified using gl.set.verbosity].}
}
\value{
An unaltered genlight object.
}
\description{
This script uses any field with numeric values stored in $other$loc.metrics
to produce summary statistics (mean, minimum, average, quantiles), histograms
and boxplots to assist the decision of choosing thresholds for the filter
function \code{\link{gl.filter.locmetric}}.
}
\details{
The function \code{\link{gl.filter.locmetric}} will filter out the
 loci with a locmetric value below a specified threshold.

The fields that are included in dartR, and a short description, are found
below. Optionally, the user can also set his/her own field by adding a vector
into $other$loc.metrics as shown in the example. You can check the names of
all available loc.metrics via: names(gl$other$loc.metrics).

\itemize{
\item SnpPosition - position (zero is position 1) in the sequence tag of the
defined SNP variant base.
\item CallRate - proportion of samples for which the genotype call is
non-missing (that is, not '-' ).
\item OneRatioRef - proportion of samples for which the genotype score is 0.
\item OneRatioSnp - proportion of samples for which the genotype score is 2.
\item FreqHomRef - proportion of samples homozygous for the Reference allele.
\item FreqHomSnp - proportion of samples homozygous for the Alternate (SNP)
allele.
\item FreqHets - proportion of samples which score as heterozygous, that is,
scored as 1.
\item PICRef - polymorphism information content (PIC) for the Reference allele.
\item PICSnp - polymorphism information content (PIC) for the SNP.
\item AvgPIC - average of the polymorphism information content (PIC) of the
reference and SNP alleles.
\item AvgCountRef - sum of the tag read counts for all samples, divided by the
number of samples with non-zero tag read counts, for the Reference allele row.
\item AvgCountSnp - sum of the tag read counts for all samples, divided by the
number of samples with non-zero tag read counts, for the Alternate (SNP) allele
row.
\item RepAvg - proportion of technical replicate assay pairs for which the
marker score is consistent.
\item rdepth - read depth.
}

\strong{ Function's output }

The minimum, maximum, mean and a tabulation of quantiles of the locmetric
values against thresholds rate are provided. Output also includes a boxplot
and a histogram.

Quantiles are partitions of a finite set of values into q subsets of (nearly)
equal sizes. In this function q = 20. Quantiles are useful measures because
they are less susceptible to long-tailed distributions and outliers.

 Plots and table were saved to the temporal directory (tempdir) and can be
 accessed with the function \code{\link{gl.print.reports}} and listed with
 the function \code{\link{gl.list.reports}}. Note that they can be accessed
 only in the current R session because tempdir is cleared each time that the
 R session is closed.

 Examples of other themes that can be used can be consulted in:
  \itemize{
 \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
 \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
 }
}
\examples{
# adding dummy data
test <- testset.gl
test$other$loc.metrics$test <- 1:nLoc(test)
# SNP data
out <- gl.report.locmetric(test,metric='test')

# adding dummy data
test.gs <- testset.gs
test.gs$other$loc.metrics$test <- 1:nLoc(test.gs)
# Tag P/A data
out <- gl.report.locmetric(test.gs,metric='test')

}
\seealso{
\code{\link{gl.filter.locmetric}}, \code{\link{gl.list.reports}},
 \code{\link{gl.print.reports}}

Other report functions: 
\code{\link{gl.report.bases}()},
\code{\link{gl.report.callrate}()},
\code{\link{gl.report.diversity}()},
\code{\link{gl.report.hamming}()},
\code{\link{gl.report.heterozygosity}()},
\code{\link{gl.report.hwe}()},
\code{\link{gl.report.ld.map}()},
\code{\link{gl.report.maf}()},
\code{\link{gl.report.monomorphs}()},
\code{\link{gl.report.overshoot}()},
\code{\link{gl.report.pa}()},
\code{\link{gl.report.parent.offspring}()},
\code{\link{gl.report.rdepth}()},
\code{\link{gl.report.reproducibility}()},
\code{\link{gl.report.secondaries}()},
\code{\link{gl.report.sexlinked}()},
\code{\link{gl.report.taglength}()}
}
\author{
Luis Mijangos (Post to \url{https://groups.google.com/d/forum/dartr})
}
\concept{report functions}
