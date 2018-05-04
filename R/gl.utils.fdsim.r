#' Estimate the rate of false positives in a fixed difference analysis
#'
#' This is a support script, called by gl.collapse.recursive.
#' 
#' The script takes two populations and generates allele frequency profiles for them. It then samples
#’ an allele frequency for each, at random, and estimates a sampling distribution for those two allele
#’ frequencies. Drawing two samples from those sampling distributions, it calculates whether or not they
#’ represent a fixed difference. This is applied to all loci, and the number of fixed differences so
#’ generated are counted, as an expectation. The script distinguished between true fixed differences (with
#’ a tolerance of delta), and false positives. The simulation is repeated a given number of times 
#’ (default=1000)to provide an expectation of the number of false positives, given the observed allele 
#’ frequency profiles and the sample sizes. The probability of the observed count of fixed differences is
#’ greater than the expected number of false positives is calculated.
#'
#' @param gl -- name of the genlight containing the SNP genotypes [required]
#' @param poppair -- labels of two populations for comparison in the form c(popA,popB) [required]
#' @param obs -- observed number of fixed differences between the two populations [required]
#' @param reps -- number of replications to undertake in the simulation [default 1000]
#' @param delta -- the threshold value for the minor allele frequency to regard the difference between two populations to be fixed [default 0.02]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A list containing the following square matricies
#'         [[1]] observed fixed differences;
#'         [[2]] mean expected number of false positives for each comparison;
#'         [[3]] standard deviation of the no. of false positives for each comparison;
#'         [[4]] probability the observed fixed differences arose by chance for each comparison;
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)

gl.utils.fdsim <- function(gl, poppair, obs=NULL, reps=1000, delta=0.02, v=2) {
  
  if (v > 0) {
    cat("Starting gl.utils.fdsim: Using simulation to estimate frequency of false positives\n")
  }
  # Determine data type
  if(!(class(gl)=="genlight")){
    cat("Fatal Error: Input data must be a genlight object\n")
    stop("Execution terminated\n")
  }
  
   if (!(poppair[1] %in% levels(pop(gl)))){
    cat("Fatal Error: Population A mislabelled\n")
     stop("Execution terminated\n")
  }
  if (!(poppair[2] %in% levels(pop(gl)))){
    cat("Fatal Error: Population B mislabelled\n")
    stop("Execution terminated\n")
  }
  
  # Extract the data for the two nominated populations
  if (length(poppair) == 2) {
    pair <- gl.keep.pop(gl, poppair, recalc=FALSE, mono.rm = TRUE, v=v)
  } else {
    cat("Fatal Error: Must specify two populations labels from the genlight object, e.g. poppair=c(popA,popB)\n")
    stop("Execution terminated\n")
  }  

  if (v > 1) {cat("Comparing ",poppair[1],"and",poppair[2],"\n")}
  if (v > 2) {
    cat("Sample sizes:",table(pop(pair)),"\n")
    cat("No. of loci:",nLoc(pair),"\n")
    cat("No. of observed fixed differences:",obs,"\n")
  }

  # Calculate the percentage frequencies
  rf <- gl.percent.freq(pair, v=v)
  
  # Disaggregate the data for the two populations
  rfA <- rf[rf$popn==poppair[1],]
  rfB <- rf[rf$popn==poppair[2],]
  # hist(rfA$frequency)
  # hist(rfB$frequency)
  
  # Initialize storage vectors
  simA <- array(data=NA,dim=nLoc(pair))
  simB <- array(data=NA,dim=nLoc(pair))
  fd <- array(data=NA,dim=nLoc(pair))
  falsepos <- array(data=NA,dim=reps)
  
  # Caclulate the observed fixed differences
  if (is.null(obs)) {
    fdmat <- gl.fixed.diff(pair)
    obs <- fdmat[2,1]
  }

  # Calculate the rate of false positives, given delta and actual sample sizes
  if (v > 1) {
    cat("Calculating false positive rate with",reps,"replications. Please be patient\n")
  }
  for (j in 1:reps){ # Repeat reps times
    for (i in 1:nLoc(pair)) {
      # Eliminate from consideration loci for which either frequency is missing
      if (is.na(rfA$frequency[i]) || is.na(rfB$frequency[i])) {
        simA[i]<-NA
        simB[i]<-NA
      # Eliminate from consideration fixed differences or near fixed differences, for given delta
      } else if ((((rfA$frequency[i]/100) < delta) && (1-(rfB$frequency[i]/100)) < delta)
           || (((rfB$frequency[i]/100) < delta) && (1-(rfA$frequency[i]/100)) < delta)) {
        simA[i]<-NA
        simB[i]<-NA
      # Calculate the number of fixed differences to arise by chance
      } else {
        szA <- rfA$nobs[i]*2
        szB <- rfB$nobs[i]*2
        pbA <- rfA$frequency[i]/100
        pbB <- rfB$frequency[i]/100
        # Sample an allele frequency from each population, expressed as a proportion
        simA[i] <- (rbinom(n=1, size=szA, prob=pbA))/szA
        simB[i] <- (rbinom(n=1, size=szB, prob=pbB))/szB
        # Apply Equation 5
        fd[i] <- ((1-simA[i])**szA)*((simB[i])**szB)+((simA[i])**szA)*((1-simB[i])**szB)
      }
    }
    # Count the number of false positives
    falsepos[j] <- sum(fd, na.rm=TRUE)
  }
  # Calculate the mean and standard deviation for the false positive rate
  mn <- mean(falsepos)
  sdev <- sd(falsepos)
  # Calculate the probability of the observed FD count, given the simulated result
  nprob <- pnorm(obs, mean=mn, sd=sdev, lower.tail=FALSE)
  
  if (v > 2) {
    cat("Threshold minor allele frequency for generating a false positive:",delta,"\n")
    cat("Estimated mean count of false positives:",round(mean(falsepos),1),"\n")
    cat("Estimated SD of the count of false positives:",round(sd(falsepos),2),"\n")
    cat("Prob that observed count of",obs,"are false positives:",nprob,"\n")
  }
  
  l <- list(observed=obs,mnexpected=mn,sdexpected=sdev,prob=nprob)
  
  if (v > 0) {
    cat("Completed gl.utils.fdsim\n\n")
  }
  return(l)
}
