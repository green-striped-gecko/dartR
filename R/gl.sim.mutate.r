#' Simulate mutations within a genlight object
#'
#' This script is intended to be used within the simulation framework of dartR. 
#' It adds the ability to add a constant mutation rate across all loci. Only 
#' works currently for biallelic data sets (SNPs). Mutation rate is checking for 
#' all alleles position and mutations at loci with missing values are ignored 
#' and in principle 'double mutations' at the same loci can occur, but should be 
#' rare.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param mut.rate Constant mutation rate over nInd*nLoc*2 possible locations
#'  [default 1e-6]
#' @return Returns a genlight object with the applied mutations
#' @export
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' b2 <- gl.sim.mutate(bandicoot.gl,mut.rate=1e-4 )
#' #check the mutations that have occurred
#' table(as.matrix(bandicoot.gl), as.matrix(b2))

gl.sim.mutate <- function(x,
                          mut.rate = 1e-06) {
    nm <- rbinom(1, nInd(x) * nLoc(x) * 2, mut.rate)
    for (ii in 1:nm) {
        ri <- sample(1:nInd(x), 1)
        rl <- sample(1:nLoc(x), 1)
        cs <- as.matrix(x)[ri, rl]
        if (!is.na(cs))
        {
            xx <- as.matrix(x[ri,])
            
            if (!cs %% 2)
                nv <- 1
            else
                nv = sample(c(0, 2), 1)
            xx[rl] <- nv
            x@gen[[ri]] <- new("SNPbin", xx)
        }  #end missing
    }
    return(x)
}
