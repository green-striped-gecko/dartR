#' Filter for CloneID to select only unique SNPs
#'
#' @param gl a genlight object created via read.dart (needs to have a cloneID as provided by dart)
#' @return filtered genlight object, with unique cloneIDs
#' @export
#' @examples{
#' }



gl.filter.cloneid <- function(gl)

{  
  
  # Work around a bug in adegenet if genlight object is created by subsetting

  if (nLoc(gl)!=nrow(x@other$loc.metrics)) { stop("The number of rows in the loc.metrics table does not match the number of loci in your genlight object!")  }

  
  
nclones <- table(gl@other$loc.metrics$clone)


print(table(nclones)  )



ind <- NA

cl <- gl@other$loc.metrics$clone
for (i in 1:length(nclones))
{
  
  sel <- which(cl== names(nclones)[i])
  
  ind[i] <- sample(as.character(sel), 1)
  
}

ind <- as.numeric(ind)
  

glf <- gl[, ind]
glf@other$loc.metrics <- gl@other$loc.metrics[ind,]

#add to history
nh <- length(glf@other$history)
glf@other$history[[nh + 1]] <- match.call()
return(glf)
}


