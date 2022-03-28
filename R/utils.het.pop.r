#' Calculates expected mean expected heterozygosity per population
#' @param x A genlight object containing the SNP genotypes [required].
#' @return A vector with the mean expected heterozygosity for each population
#' @export
#' @author Bernd Gruber & Luis Mijangos (bugs? Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' out <- utils.het.pop(testset.gl)

utils.het.pop <- function(x) {
    # Split the genlight object into a list of populations
    sgl <- seppop(x)
    Hexp <- array(NA, length(sgl))
    # For each population
    for (i in 1:length(sgl)) {
        gl <- sgl[[i]]
        t <- as.matrix(gl)
        p <- colMeans(t == 0, na.rm = T)
        q <- colMeans(t == 2, na.rm = T)
        hets <- colMeans(t == 1, na.rm = T)
        p <- (2 * p + hets) / 2
        q <- (2 * q + hets) / 2
        H <- 1 - (p * p + q * q)
        Hexp[i] <- round(mean(H, na.rm = T), 6)
    }
    invisible(Hexp)
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

