####### Likelihood functions needed by OutFLANK


EffectiveNumberSamplesMLE=function(FstVect, Fstbar, NumberOfSamples, SmallestFstInTrimmedList, LargestFstInTrimmedList){
  #This function should find the maximum likelihood value 
  #of the effective number of samples, for a given list of
  #Fst values.
  
  #The FstVect should already have been purged of NaN values and of loci with 
  #too low heterozygosity or MAF. 
  
  sortedFst=FstVect[order(FstVect)]  
  
  #The Minimum Fst considered in the trimmed data is the larger of the amount
  #specified by the user or the mean FSt over 100. This is to prevent extremely
  #small Fsts from causing estiamtion errors (Especially when R interprets a
  #small Fst as FSt=0.)
  LowTrimPoint=max(Fstbar/100,SmallestFstInTrimmedList)
  
  trimmedFstVect =FstVect[which((FstVect>=LowTrimPoint)&(FstVect<=LargestFstInTrimmedList))]
  
  trimmedFstArray=as.array(trimmedFstVect)
  
  localNLLAllData=function(dfInferred){
    localNLLOneLocus=function(Fst){
      negLLdfFstTrim(Fst,dfInferred,Fstbar,LowTrimPoint,LargestFstInTrimmedList)
    }
    sum(localNLLOneLocus(trimmedFstVect))
  }
  
  optim(NumberOfSamples, localNLLAllData, lower=2, method="L-BFGS-B")$par
}

IncompleteGammaFunction=function(a, z) {
  #equivalence to Mathematica Gamma[a,z] according to 
  #   http://r.789695.n4.nabble.com/Incomplete-Gamma-function-td833545.html
  pgamma(z,a,lower.tail=FALSE)*gamma(a)
}

negLLdfFstTrim=function(Fst, dfInferred, Fstbar, LowTrimPoint, HighTrimPoint){
  #Fst is the Fst from a locus, and dfInferred is the candidate value for the
  #degrees of freedom for the chi-squared distribution of neutral Fst, and
  #Fstbar is the mean Fst of all neutral loci (sequentially inferred from
  #non-outlier loci) LowTrimPoint and HighTrimPoint are the values of the lowest
  #and highest Fst values allowed to be included in the Fst list.
  #
  #Finds contribution to the negative log likelihood of a given locus' Fst for a
  #given dfInferred #CHECKED AGAINST MATHEMATICA DERIVATION##
  
  df=dfInferred
  
  1/(2*Fstbar)*(df * Fst +df * Fstbar * log(2) - df * Fstbar *log(df)-(df-2)*Fstbar * log(Fst)+df * Fstbar * log(Fstbar) + 2*Fstbar * log(-IncompleteGammaFunction(df/2,df*HighTrimPoint/(2*Fstbar))+IncompleteGammaFunction(df/2,df*LowTrimPoint/(2*Fstbar))))
}






