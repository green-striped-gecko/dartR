#' Export DArT genlight object \{adegenet\} to faststructure format (to run faststructure elsewhere)
#'
#' Recodes in the quite specific faststructure format (e.g first six columns need to be there, but are ignored...check faststructure documentation (if you find any :-( )))
#'

#' The script writes out the a file in faststructure format. 
#' @param gl -- genlight object
#' @param outfile -- filename of the output fastA file [default genlight.fasta]
#' @param probar switch to show/hide progress bar
#' @return NULL
#' @export
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @import adegenet
#' @importFrom utils getTxtProgressBar setTxtProgressBar txtProgressBar 

gl2faststructure <- function(gl, outfile="gl.str", probar=TRUE) 
{
  x <- as.matrix(gl)
  #add six dummy colums
  nc <- ncol(x)+6
  if (probar) pb <- txtProgressBar(min=0, 1, style=3, initial=NA)
  zz <- file(outfile, "w")
  for (i in 1:nrow(x))
  {
    dummy <- rbind(x[i,], x[i,])
    index <- colSums(dummy, na.rm=T)==2
    dummy[, index] <- c(0,2)
    dummy <- ifelse(is.na(dummy), -9, dummy)
    dummy <- ifelse(dummy==0,1, dummy)
    dummy <- cbind(i,i,i,i,i,i,dummy)
    write(t(dummy), file=outfile, sep="\t", ncolumns = nc, append=TRUE)  
    if (probar)   setTxtProgressBar(pb, i/nrow(x))  
  }
  close(zz)
  if (probar) close(pb)
  #gi <-read.structure("project_data.str", ask=F, n.ind=50, n.loc=100, col.pop=0, col.lab = 1, onerowperind = TRUE, sep="\t")
  cat(paste0("Saved faststructure file: ",getwd(),"/", outfile, "\n") )
  cat(paste("Consists of", nrow(x), "individuals and ", ncol(x), "loci."))
  return(NULL)
}
