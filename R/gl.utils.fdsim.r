#' Estimate the rate of false positives in a fixed difference analysis
#'
#' This is a support script, called by gl.collapse.recursive.
#' 
#' The script takes two populations, generates allele frequency profiles for them, then takes two samples
#' from those profiles to estimate the rate by which fixed differences are generated for a given pair
#' of sample sizes. It then calculates the rate of fixed differences for effectively samples of
#' infinite size, and subtracts the two to get an expected rate of false positives. The expected
#' rate of false positives is used to compute a threshold above which the two populations can be
#' considered distinct.
#'
#' @param gl -- name of the genlight containing the SNP genotypes [required]
#' @param poppair -- labels of two populations for comparison in the form c(popA,popB) [required]
#' @param obs -- observed number of fixed differences between the two populations [required]
#' @param reps -- number of replications to undertake in the simulation [default 1000]
#' @param delta -- the threshold value for the minor allele frequency to regard the difference between two populations to be fixed [default 0.02]
#' @param v -- verbosity = 0, silent; 1, brief; 2, verbose [default 1]
#' @return Expected count of fixed differences arising from loci that are not considered to represent a fixed difference between the two focal populations.
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples

gl.utils.fdsim<- function(gl, poppair, obs=NULL, reps=1000, delta=0.02, v=1) {
  
  # Determine data type
  if(class(gl)=="genlight"){
    cat("Using SNP data from a genlight object\n")
  } else {
    cat("Fatal Error: Input data must be a genlight object\n")
    stop()
  }
  
  if (is.null(obs)) {
    cat("Warning: Number of observed fixed differences not supplied, set to NA.\n")
    obs <- NA
  }
  
  if (!(poppair[1] %in% levels(pop(gl)))){
    cat("Fatal Error: Population A mislabelled\n")
    stop()
  }
  if (!(poppair[2] %in% levels(pop(gl)))){
    cat("Fatal Error: Population B mislabelled\n")
    stop()
  }
  
  # Extract the data for the two nominated populations
  if (length(poppair) == 2) {
    pair <- gl.keep.pop(gl, poppair, recalc=FALSE, mono.rm = TRUE)
  } else {
    cat("Fatal Error: Must specify two populations labels from the genlight object, e.g. poppair=c(popA,popB)\n")
    stop()
  }  

  if (v==2) {
    cat("Comparing two populations",poppair[1],"and",poppair[2],"\n")
    cat("Sample sizes:",table(pop(pair)),"\n")
    cat("No. of loci:",nLoc(pair),"\n")
    cat("No. of observed fixed differences:",obs,"\n")
  }

  # Calculate the percentage frequencies
  rf <- gl.percent.freq(pair)
  
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
  
  # Calculate the rate of false positives, given delta and actual sample sizes
  cat("Calculating false positve rate with",reps,"replications for",poppair[1],"vs",poppair[2],". Please be patient\n\n")
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
  
  if (v==2) {
    cat("Threshold minor allele frequency for generating a false positive:",delta,"\n")
    cat("Estimated mean count of false positives:",round(mean(falsepos),1),"\n")
    cat("Estimated SD of the count of false positives:",round(sd(falsepos),2),"\n")
    cat("Prob that observed count of",obs,"are false positives:",nprob,"\n")
  }
  
  return(mn)
  
}
