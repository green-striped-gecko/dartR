% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl.sim.offspring.r
\name{gl.sim.offspring}
\alias{gl.sim.offspring}
\title{Simulates a specified number of offspring based on alleles provided by
 potential father(s) and mother(s)}
\usage{
gl.sim.offspring(fathers, mothers, noffpermother, sexratio = 0.5)
}
\arguments{
\item{fathers}{Genlight object of potential fathers [required].}

\item{mothers}{Genlight object of potential mothers simulated [required].}

\item{noffpermother}{Number of offspring per mother [required].}

\item{sexratio}{The sex ratio of simulated offspring 
(females / females +males, 1 equals 100 percent females) [default 0.5.].}
}
\value{
A genlight object with n individuals.
}
\description{
This takes a population (or a single individual) of fathers (provided as a
genlight object) and mother(s) and simulates offspring based on 'random'
 mating. It can be used to simulate population dynamics and check the effect
 of those dynamics and allele frequencies, number of alleles. Another
 application is to simulate relatedness of siblings and compare it to actual
 relatedness found in the population to determine kinship.
}
\examples{
#Simulate 10 potential fathers
gl.fathers <- glSim(10, 20, ploidy=2)
#Simulate 10 potential mothers
gl.mothers <- glSim(10, 20, ploidy=2)
gl.sim.offspring(gl.fathers, gl.mothers, 2, sexratio=0.5)
}
\author{
Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
}
