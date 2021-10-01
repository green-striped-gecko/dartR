#' Simulate emigration between populations
#' 
#' A function that allows to exchange individuals of populations within a genlight object (=simulate emigration between populations). There are two ways to specify emigration. If an emi.table is provided (a square matrix of dimension of the populations that specifies the emigration from column x to row y), then emigration is deterministic in terms of numbers of individuals as specified in the table. If perc.mig and emi.m are provided, then emigration is probabilistic. The number of emigrants is determined by the population size times the perc.mig and then the population where to migrate to is taken from the relative probability in the columns of the emi.m table. Be aware if the diagonal is non zero then migration can occur into the same patch. So most often you want to set the diagonal of the emi.m matrix to zero. Which individuals is moved is random, but the order is in the order of populations. It is possible that an individual moves twice within an emigration call(as there is no check, so an individual moved from population 1 to 2 can move again from population 2 to 3). 
#' @param x genlight or list of genlight objects
#' @param perc.mig percentage of individuals that migrate. [emigrates = nInd times perc.mig]
#' @param emi.m probabilistic emigration matrix. emigrate from=column to=row
#' @param emi.table if presented emi.m matrix is ignored. Deterministic emigration as specified in the matrix [a square matrix of dimenstion of the number of populations]. e.g. an entry in the 'emi.table[2,1]<- 5' means that five individuals emigrate from population 1 to population 2 [from=columns and to=row].
#' @return a list or a single [depends on the input] genlight object, where emigration between population has happened
#' @author Custodian: Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' x <- possums.gl
#' #one individual moves from every population to 
#' #every other population
#' emi.tab <- matrix(1, nrow=nPop(x), ncol=nPop(x))
#' diag(emi.tab)<- 0
#' np <- gl.sim.emigration(x, emi.table=emi.tab)
#' np
#' @export

gl.sim.emigration <- function (x, perc.mig=NULL, emi.m=NULL, emi.table = NULL) 
{
  if (length(x)==1) p <- seppop(x) else p <- x
  pn <- names(p)
  n.pops = length(p)
  migs <- matrix(0, nrow = n.pops, ncol = n.pops)
  migrants <- NA
 if(is.null(emi.table)){
   #convert disdis into prob
    emi.m = emi.m /rep(colSums(emi.m), each = ncol(emi.m))  
    for (i in 1:n.pops) {
      pop.size <- nrow(p[[i]])
      if (!is.null(pop.size)) {
        migrants[i] <- sum(ifelse(runif(pop.size) < perc.mig, 
                                  1, 0))
        fromto <- table(sample(1:n.pops, migrants[i], 
                               replace = T, prob = emi.m[, i]))
        to <- as.numeric(names(fromto))
        migs[i,to] <- migs[i,to] + fromto 
      }
    }
  }  else migs <- emi.table
  
  diag(migs) <- 0 #do not care about staying
  for (from in 1:n.pops) {
    for (to in 1:n.pops) {
      psize <- nrow(p[[from]])
      m <- migs[to, from]
      if (m > 0) {
        ind.from <- sample(1:psize, m, replace = F)
        p[[to]] <- rbind(p[[to]], p[[from]][ind.from,])
        p[[from]] <- p[[from]][-ind.from, ]
      }
    }
  }
  for (i in 1:n.pops)   pop(p[[i]])<- rep(pn[i],nrow(p[[i]]))
  #return list or single population (depending on the input)
  if (length(x)==1)  xout <- do.call("rbind", p) else xout <- p
  return(xout)
}

