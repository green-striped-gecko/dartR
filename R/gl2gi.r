#' Converts a genlight object to genind object
#' 
#' @param x -- a genlight object
#' @param probar -- if TRUE, a progress bar will be displayed for long loops [default = TRUE]
#' @param v -- level of verbosity. v=0 is silent, v=1 returns more detailed output during conversion.
#' @return A genind object, with all slots filled.
#' @export
#' @author Bernd Gruber (bugs? Post to \url{https://groups.google.com/d/forum/dartr})
#' @details this function uses a faster version of df2genind (from the adgegenet package)

gl2gi <- function(gl, probar=TRUE, v=2) {

if (v==1) {
  cat("Start conversion....\n")
  ptm <- proc.time()[3]
  cat("Please note conversion of bigger data sets will take some time!\n" )
  cat("Once finished, we recommend to save the object using >save(object, file=\"object.rdata\")\n")
}
#convert to genind....
x <- as.matrix(x[,])

if (probar) {pb <- txtProgressBar(min=0, max=1, style=3, initial=NA)}

  if (is.null(gl@loc.all))  {
    gl@loc.all <- rep("A/T", nLoc(gl))
    gl@loc.all[1]<- "C/G"
  }


homs1 <- paste(substr(gl@loc.all,1,1),"/",substr(gl@loc.all,1,1), sep = "")
hets <-  gl@loc.all
homs2 <- paste(substr(gl@loc.all,3,3),"/",substr(gl@loc.all,3,3), sep = "")
xx <- matrix(NA, ncol=ncol(x), nrow=nrow(x))
for (i in 1:nrow(x))
  {
  for (ii in 1:ncol(x))
    {

    inp <- x[i,ii]
    if (!is.na(inp))
      {
      
      if (inp==0) xx[i,ii] <- homs1[ii] else if (inp==1) xx[i,ii] <- hets[ii] else if (inp==2) xx[i,ii] <- homs2[ii]
      } else xx[i,ii]="-/-"
    }
  if (probar) {setTxtProgressBar(pb, i/nrow(x))}

if (v==1) {
  cat("\nMatrix converted.. Prepare genind object...\n")
}
  if (probar) {close(pb)}

gen<-df2genind(xx[,], sep="/", ncode=1, ind.names=gl@ind.names, pop = gl@pop, ploidy=2,  NA.char = "-")#, probar=probar)
gen@other <- gl@other
locNames(gen)<- locNames(gl)

if (v==1)cat(paste("Finished! Took", round(proc.time()[3]-ptm),"seconds.\n") )
gen
}


