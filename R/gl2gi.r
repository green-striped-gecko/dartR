#' Converts a genlight object to genind object
#' 
#' @param gl -- a genlight object
#' @param probar -- switch on/off progress bar
#' @return A genind object, with all slots filled.
#' @export
#' @author Bernd Gruber (glbugs@@aerg.canberra.edu.au)
#' @details this function uses a faster version of df2genind (from the adgegenet package)

gl2gi <- function(gl, probar=TRUE) {

cat("Start conversion....\n")
ptm <- proc.time()[3]
cat("Please note conversion of bigger data sets will take some time!\n" )
cat("Once finished, we recommend to save the object using >save(object, file=\"object.rdata\")\n")
#convert to genind....
x <- as.matrix(gl[,])
if (probar) pb <- txtProgressBar(min=0, max=1, style=3, initial=NA)

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
if (probar)   setTxtProgressBar(pb, i/nrow(x))
  }
  
cat("\nMatrix converted.. Prepare genind object...\n")
if (probar) close(pb)
gen<-df2genind(x[,], sep="/", ncode=1, ind.names=gl@ind.names, pop = gl@pop, ploidy=2)#, probar=probar)
gen@other <- gl@other

cat(paste("Finished! Took", round(proc.time()[3]-ptm),"seconds.\n") )
gen
}


