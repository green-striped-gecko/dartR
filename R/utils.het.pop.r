#' Calculates expected mean expected heterozygosity per population
#' @param x A genlight object containing the SNP genotypes [required].
#' @param t_het A string specifying the type of expected heterozygosity to be
#' calculated. Options are "He" for expected heterozygosity and "Ho" for observed
#' @return A vector with the mean expected heterozygosity for each population
#' @export
#' @author Bernd Gruber & Luis Mijangos (bugs? Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' out <- utils.het.pop(testset.gl,t_het="He")

utils.het.pop <- function(x,t_het) {
    # Split the genlight object into a list of populations
    sgl <- seppop(x)
    res_het <- array(NA, length(sgl))
    # For each population
    for (i in 1:length(sgl)) {
        gl <- sgl[[i]]
        t <- as.matrix(gl)
        p <- colMeans(t == 0, na.rm = TRUE)
        q <- colMeans(t == 2, na.rm = TRUE)
        hets <- colMeans(t == 1, na.rm = TRUE)
        if(t_het=="He"){
          p <- (2 * p + hets) / 2
          q <- (2 * q + hets) / 2
          H <- 1 - (p * p + q * q)
          res_het[i] <- round(mean(H, na.rm = TRUE), 6)
        }
        if(t_het=="Ho"){
          res_het[i] <- round(mean(hets, na.rm = TRUE), 6)
        }

    }
    invisible(res_het)
}

ind.count <- function(x) {
    # the loci that are completely missing
    loci.na <-
        which(colSums(is.na(as.matrix(x))) == nrow(as.matrix(x)))
    # the number of samples in the matrix the number of non-genotyped
    # samples remove the loci that are completely missing
    if (length(loci.na) > 0) {
        nind <-
            mean(nrow(as.matrix(x)) - colSums(is.na(as.matrix(x)))[-loci.na])
        # the number of samples in the matrix the number of
        # non-genotyped samples
    } else {
        nind <- mean(nrow(as.matrix(x)) - colSums(is.na(as.matrix(x))))
    }
    
    return(nind)
}

