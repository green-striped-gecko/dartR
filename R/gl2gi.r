#' Converts a genlight object to genind object
#' 
#' @param gl -- a genlight object
#' @param v -- level of verbosity. v=0 is silent, v=1 returns more detailed output during conversion.
#' @return A genind object, with all slots filled.
#' @export
#' @author Bernd Gruber (bugs? Post to \url{https://groups.google.com/d/forum/dartr})
#' @details this function uses a faster version of df2genind (from the adgegenet package)

gl2gi <- function(gl, v=1) {

if (v==1) {
  cat("Start conversion....\n")
  ptm <- proc.time()[3]
  cat("Please note conversion of bigger data sets will take some time!\n" )
  cat("Once finished, we recommend to save the object using >save(object, file=\"object.rdata\")\n")
}
#convert to genind....
x <- as.matrix(gl[,])
if (v==1) pb <- txtProgressBar(min=0, max=1, style=3, initial=NA)

for (i in 1:nrow(x))
  {
  for (ii in 1:ncol(x))
    {

    inp <- x[i,ii]
    if (!is.na(inp))
      {
      if (inp==0) x[i,ii] <- "A/A" else if (inp==1) x[i,ii] <- "A/B" else if (inp==2) x[i,ii] <- "B/B"
      }
    }
if (v==1)   setTxtProgressBar(pb, i/nrow(x))
  }
  
if (v==1) {
  cat("\nMatrix converted.. Prepare genind object...\n")
  close(pb)
}
gen<-df2genind(x[,], sep="/", ncode=1, ind.names=gl@ind.names, pop = gl@pop, ploidy=2)#, probar=probar)
gen@other <- gl@other

if (v==1)cat(paste("Finished! Took", round(proc.time()[3]-ptm),"seconds.\n") )
gen
}


