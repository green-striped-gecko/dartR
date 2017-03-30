#' Filter for CloneID to select only unique SNPs
#'
#' @param gl a genlight object created via read.dart (needs to have a cloneID as provided by dart)
#' @return filtered genlight object, with unique cloneIDs
#' @export
#' @examples{
#' }



gl.filter.cloneid <- function(gl)

{  
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


return(glf)
}


