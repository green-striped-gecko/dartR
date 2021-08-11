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
  for (i in 1:n.pop)   pop(p[[i]])<- rep(pn[i],nrow(p[[i]]))
  #return list or single population (depending on the input)
  if (length(x)==1)  xout <- do.call("rbind", p) else xout <- p
  return(xout)
}

