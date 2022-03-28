################################################################################
############## utils functions for gl.report.hwe ###############################
################################################################################

####################################### GenerateSamples #######################
GenerateSamples <- function(n = 5) {
    # generates all possible samples of size n.
    Res <- NULL
    for (i in 0:n) {
        AA <- i
        for (j in 0:(n - i)) {
            AB <- j
            BB <- (n - (AA + AB))
            sam <- c(AA, AB, BB)
            Res <- rbind(Res, sam)
        }
    }
    rownames(Res) <- 1:nrow(Res)
    colnames(Res) <- c("AA", "AB", "BB")
    return(Res)
}

############################################## CritSam #######################
CritSam <- function(n, Dpos, alphalimit, pvaluetype) {
    X <- GenerateSamples(n)
    ncomp <- nrow(X)
    Res <- NULL
    Ds <- NULL
    pval <- NULL
    fA <- NULL
    for (i in 1:nrow(X)) {
        fA <- c(fA, (2 * X[i, 1] + X[i, 2]) / (2 * n))
        Ds <-
            c(Ds, HardyWeinberg::HWChisq(X[i,], verbose = FALSE)$D)
        pval <-
            c(
                pval,
                HardyWeinberg::HWExact(
                    X[i,],
                    alternative = "two.sided",
                    pvaluetype = pvaluetype,
                    verbose = F
                )$pval
            )
    }
    
    Y <- data.frame(X[, 1], X[, 2], X[, 3], fA, Ds, pval)
    colnames(Y) <- c("AA", "AB", "BB", "fA", "Ds", "pval")
    if (Dpos)
        Y <- Y[Y$Ds > 0,]
    else
        Y <- Y[Y$Ds < 0,]
    Y <- Y[Y$pval < alphalimit,]
    fre <- unique(fA)
    for (i in 1:length(fre)) {
        Ys <- Y[Y$fA == fre[i],]
        if (nrow(Ys) > 0) {
            indi <- which.max(Ys$pval)
            Ys <- Ys[indi,]
            Res <- rbind(Res, c(Ys$AA, Ys$AB, Ys$BB))
        }
    }
    Xn <- Res / n
    return(list(Xn = Xn, Ds = Ds, fA = fA))
}

############################################ CritSam_Chi #######################
CritSam_Chi <- function(n, Dpos, alphalimit, cc) {
    X <- GenerateSamples(n)
    ncomp <- nrow(X)
    Res <- NULL
    Ds <- NULL
    pval <- NULL
    fA <- NULL
    for (i in 1:nrow(X)) {
        fA <- c(fA, (2 * X[i, 1] + X[i, 2]) / (2 * n))
        Ds <-
            c(Ds,
              HardyWeinberg::HWChisq(X[i,], cc = cc, verbose = FALSE)$D)
        pval <-
            c(pval,
              HardyWeinberg::HWChisq(X[i,], cc = cc, verbose = FALSE)$pval)
    }
    
    Y <- data.frame(X[, 1], X[, 2], X[, 3], fA, Ds, pval)
    colnames(Y) <- c("AA", "AB", "BB", "fA", "Ds", "pval")
    if (Dpos)
        Y <- Y[Y$Ds > 0,]
    else
        Y <- Y[Y$Ds < 0,]
    Y <- Y[Y$pval < alphalimit,]
    fre <- unique(fA)
    for (i in 1:length(fre)) {
        Ys <- Y[Y$fA == fre[i],]
        if (nrow(Ys) > 0) {
            indi <- which.max(Ys$pval)
            Ys <- Ys[indi,]
            Res <- rbind(Res, c(Ys$AA, Ys$AB, Ys$BB))
        }
    }
    Xn <- Res / n
    return(list(Xn = Xn, Ds = Ds, fA = fA))
}
