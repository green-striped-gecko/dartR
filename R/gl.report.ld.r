##################################### pairwise LD function across all loci
#' Calculates pairwise population based Linkage Disequilibrium across all loci
#' using the specified number of cores
#' @description
#' This function is implemented in a parallel fashion to speed up the process.
#' There is also the ability to restart the function if crashed by specifying
#' the chunk file names or restarting the function exactly in the same way as in
#' the first run. This is implemented as sometimes due to connectivity loss
#' between cores the function my crash half way. Also remove loci with have only
#'  missing value before running the function.
#'
#' @param x A genlight or genind object created (genlight objects are internally
#'  converted via \code{\link{gl2gi}} to genind) [required].
#' @param name Character string for rdata file. If not given genind object name
#' is used [default NULL].
#' @param save Switch if results are saved in a file [default TRUE].
#' @param outpath Folder where chunks and results are saved (if save=TRUE)
#' [default tempdir()].
#' @param nchunks How many subchunks will be used (the less the faster, but if
#' the routine crashes more bits are lost) [default 2].
#' @param ncores How many cores should be used [default 1].
#' @param chunkname The name of the chunks for saving [default NULL].
#' @param probar if TRUE, a progress bar is displayed for long loops
#' [default = TRUE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @return Returns calculation of pairwise LD across all loci between
#' subpopulations. This functions uses if specified many cores on your computer
#' to speed up. And if save is used can restart (if save=TRUE is used) with the
#' same command starting where it crashed. The final output is a data frame that
#'  holds all statistics of pairwise LD between loci. (See ?LD in package
#'  genetics for details).
#' @export
#' @import foreach
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})

gl.report.ld <- function(x,
                         name = NULL,
                         save = TRUE,
                         outpath = tempdir(),
                         nchunks = 2,
                         ncores = 1,
                         chunkname = NULL,
                         probar = FALSE,
                         verbose = NULL) {
    if (!(requireNamespace("doParallel", quietly = TRUE))) {
        stop(error(
            "Package doParallel needed for this function to work. Please install it."
        ))
    }
    
    if (!(requireNamespace("data.table", quietly = TRUE))) {
        stop(error(
            "Package data.table needed for this function to work. Please install it."
        ))
    }
    
    if (!(requireNamespace("foreach", quietly = TRUE))) {
        stop(error(
            "Package foreach needed for this function to work. Please install it."
        ))
    } else {
        # SET VERBOSITY
        verbose <- gl.check.verbosity(verbose)
        
        # FLAG SCRIPT START
        funname <- match.call()[[1]]
        utils.flag.start(func = funname,
                         build = "Jackson",
                         verbosity = verbose)
        
        # CHECK DATATYPE
        datatype <-
            utils.check.datatype(x, verbose = verbose)
        
        # convert genlight to genind
        if (class(x) == "genlight") {
            gi <- gl2gi(x)
        }
        
        # library(doParallel) library(adegenet) library(data.table)
        if (verbose >= 2) {
            cat(report(paste(
                "  Calculating LD for all pairs of loci...\n"
            )))
            cat(report(
                paste("  Using", ncores, "cores in", nchunks, " chunks.\n")
            ))
        }
        if (verbose >= 3) {
            cat(report(
                "  Depending on the number of loci this may take a while...\n"
            ))
            cat(
                report(
                    "  nchunks specifies the number of steps in the progress bar and the number of intermediate saves, but slows the computation a bit. nchunks = 1 is fastest.\n"
                )
            )
        }
        
        if (verbose >= 2) {
            cat(report(paste(
                "  Separating all", length(locNames(gi)), "loci...\n"
            )))
        }
        flush.console()
        
        # convert into list of
        slg <- seploc(gi)
        for (i in 1:length(slg)) {
            slg[[i]] <- slg[[i]]@tab
        }
        
        if (verbose >= 2) {
            cat(report(
                paste(
                    "  Generating all possible pairs:",
                    length(slg) * (length(slg) - 1) / 2,
                    "...\n"
                )
            ))
        }
        flush.console()
        allp <- combn(length(slg), 2)
        resnames <-
            c("loc1", "loc2", "D", "Dprime", "r", "R2", "n", "X2", "p")
        lddone <- NULL
        chunknr <- 0  #to make sure old chunks are not overridden
        if (!is.null(chunkname)) {
            if (verbose >= 2) {
                cat(report("  You specified results from a previous run ...\n"))
                cat(report(
                    paste(
                        "  Loooking for LD_chunks_",
                        chunkname,
                        "files.\n"
                    )
                ))
            }
            chunkfiles <-
                list.files(pattern = paste0("LD_chunks_", chunkname))
            if (length(chunkfiles > 0)) {
                if (verbose >= 2) {
                    cat(report(paste(
                        "  Found",
                        length(chunkfiles),
                        "file(s).\n"
                    )))
                }
                for (i in 1:length(chunkfiles)) {
                    load(file.path(
                        outpath,
                        paste0("LD_chunks_", chunkname, "_", i, ".rdata")
                    ))
                    lddone[[i]] <- ldc
                }
                chunknr <- length(chunkfiles)
                lddone <- data.table::rbindlist(lddone)
                
                if (verbose >= 2) {
                    cat(report(paste(
                        "  Found", nrow(lddone), "pairs...\n"
                    )))
                }
                data.table::setnames(lddone, resnames)
                done <- nrow(lddone)
                if (done == ncol(allp)) {
                    if (verbose >= 2) {
                        cat(
                            important(
                                "  Already everyting is calculated. If you want to recalculate please delete al LD_chunk files or specify a different chunkname.\n Aborting function...\n"
                            )
                        )
                        return(lddone)
                    }
                }
                allp <- allp[,-c(1:done)]
                if (verbose >= 2) {
                    cat(report(
                        paste(
                            "  ...only",
                            ncol(allp),
                            "pairs left to be done...\n"
                        )
                    ))
                }
            } else {
                if (verbose >= 2) {
                    cat(important(
                        paste(
                            "  No chunkfiles with LD_chunks_",
                            chunkname,
                            "_x.rdata found. \nRestart to calculation of all pairs.\n"
                        )
                    ))
                }
            }
        }
        if (verbose >= 2) {
            cat(report(paste(
                "  Calculate LD for all pairs...\n"
            )))
        }
        flush.console()
        n <- ncol(allp)
        ptm <- proc.time()[3]
        # ld functions here
        
        LD.fast <- function(g1, g2, ...) {
            prop.A <- colMeans(g1, na.rm = T) / 2
            if (length(prop.A) != 2)
                return(list(NA, NA, NA, NA, NA, NA, NA, NA))
            names(prop.A) <- c("A", "B")
            prop.B <- colMeans(g2, na.rm = T) / 2
            if (length(prop.B) != 2)
                return(list(NA, NA, NA, NA, NA, NA, NA, NA))
            names(prop.B) <- c("A", "B")
            
            major.A <- names(prop.A)[which.max(prop.A)]
            major.B <- names(prop.B)[which.max(prop.B)]
            pA <- max(prop.A, na.rm = TRUE)
            pB <- max(prop.B, na.rm = TRUE)
            pa <- 1 - pA
            pb <- 1 - pB
            
            Dmin <- max(-pA * pB,-pa * pb)
            pmin <- pA * pB + Dmin
            
            Dmax <- min(pA * pb, pB * pa)
            pmax <- pA * pB + Dmax
            
            mja <- which.max(colSums(g1, na.rm = T))
            mjb <- which.max(colSums(g2, na.rm = T))
            counts <- table(g1[, mja], g2[, mjb], useNA = "no")
            
            # counts <- counts[c('1','2'),c('1','2')]
            
            # counts <- table( allele.count(g1, major.A), allele.count(g2, major.B) ) rownames(counts)<-1:2 colnames(counts)<-1:2
            
            n3x3 <- matrix(0, nrow = 3, ncol = 3)
            colnames(n3x3) <- rownames(n3x3) <- 0:2
            
            # ensure the matrix is 3x3, with highest frequency values in upper left
            for (i in rownames(counts)) {
                for (j in colnames(counts)) {
                    n3x3[3 - as.numeric(i), 3 - as.numeric(j)] <-
                        counts[i, j]
                }
            }
            
            loglik <- function(pAB, ...) {
                (2 * n3x3[1, 1] + n3x3[1, 2] + n3x3[2, 1]) * log(pAB) + (2 * n3x3[1, 3] + n3x3[1, 2] + n3x3[2, 3]) * log(pA - pAB) + (2 * n3x3[3,
                                                                                                                                               1] + n3x3[2, 1] + n3x3[3, 2]) * log(pB - pAB) + (2 * n3x3[3, 3] + n3x3[3, 2] + n3x3[2, 3]) * log(1 - pA - pB + pAB) + n3x3[2,
                                                                                                                                                                                                                                                                          2] * log(pAB * (1 - pA - pB + pAB) + (pA - pAB) * (pB - pAB))
            }
            
            # SAS code uses: s <- seq(pmin+0.0001,pmax-0.0001,by=0.0001) lldmx <- loglik(s) maxi <- which.max(lldmx) pAB <- s[maxi]
            # but this should be faster:
            solution <-
                optimize(
                    loglik,
                    lower = pmin + .Machine$double.eps,
                    upper = pmax - .Machine$double.eps,
                    maximum = TRUE
                )
            pAB <- solution$maximum
            
            estD <- pAB - pA * pB
            if (estD > 0) {
                estDp <- estD / Dmax
            } else{
                estDp <- estD / Dmin
            }
            
            n <- sum(n3x3)
            
            corr <- estD / sqrt(pA * pB * pa * pb)
            
            dchi <- (2 * n * estD ^ 2) / (pA * pa * pB * pb)
            dpval <- 1 - pchisq(dchi, 1)
            
            retval <-
                list(
                    call = match.call(),
                    D = estD,
                    `D'` = estDp,
                    r = corr,
                    `R^2` = corr ^ 2,
                    n = n,
                    `X^2` = dchi,
                    `P-value` = dpval
                )
            class(retval) <- "LD"
            retval
        }
        
        # split into nchunks steps for the progress bar...
        runs <- 1:n
        if (nchunks > length(runs)) {
            nchunks <- length(runs)
        }
        chunks <-
            function(x, n)
                split(x, cut(seq_along(x), n, labels = FALSE))
        if (nchunks < 2) {
            splitruns <-
                list(runs)
        }  else{
            splitruns <- chunks(runs, nchunks)
        }
        
        ldchunks <- list()
        cl <-
            parallel::makeCluster(ncores)#adjust the number of cores of your computer!!!!
        doParallel::registerDoParallel(cl)
        if (probar) {
            pbar <-
                txtProgressBar(
                    min = 0,
                    max = nchunks,
                    style = 3,
                    initial = NA
                )
        }
        for (i in 1:nchunks) {
            iter <- splitruns[[i]]
            if (ncores < 2) {
                it2 <- list(iter)
            } else{
                it2 <- chunks(iter, ncores)
            }
            ip <- NA
            ll <-
                foreach::foreach(ip = 1:length(it2), .combine = rbind) %dopar% {
                    res <- matrix(NA, ncol = 9, nrow = length(it2[[ip]]))
                    for (ii in 1:length(it2[[ip]])) {
                        l1 <- allp[1, it2[[ip]][ii]]
                        l2 <- allp[2, it2[[ip]][ii]]
                        s1 <- slg[[l1]]
                        s2 <- slg[[l2]]
                        if ((length(colSums(s1, na.rm = T)) + length(colSums(s2, na.rm = T))) == 4) {
                            r <- LD.fast(s1, s2)
                            res[ii,] <-
                                c(l1, l2, do.call(cbind, r[2:8]))
                        } else {
                            res[ii,] <- c(l1, l2, rep(NA, 7))
                        }
                    }
                    res
                }
            ldchunks[[i]] <- as.data.frame(ll)
            if (probar) {
                setTxtProgressBar(pbar, i)
            }
            ldc <- ldchunks[[i]]
            save(ldc, file = file.path(
                outpath,
                paste0("LD_chunks_", chunkname, "_", i + chunknr, ".rdata")
            ))
        }
        parallel::stopCluster(cl)
        LDres2 <- data.table::rbindlist(ldchunks)
        # LDres2 <- t(do.call(cbind, ldchunks))
        data.table::setnames(LDres2, resnames)
        if (!is.null(lddone)) {
            LDres2 <-
                data.table::rbindlist(list(lddone, LDres2))
        }
        # colnames(LDres2)<- resnames
        
        if (verbose >= 2) {
            cat(report(
                paste(
                    "\n  No. of Simulations:",
                    n,
                    ". Took",
                    round(proc.time()[3] - ptm),
                    "seconds.\n"
                )
            ))
        }
        if (save) {
            if (!is.null(name)) {
                nobj <- name
                filename <- paste0("LD_", name, ".rdata")
            } else {
                filename = "LDallp.rdata"
                nobj <- "LDallp"
            }
            
            if (verbose >= 2) {
                cat(report(
                    paste0(
                        "\n  Results are saved as object ",
                        nobj,
                        " under ",
                        filename,
                        ".\n"
                    )
                ))
                (cat(report(
                    paste(
                        "  Once you have checked you can delete your LD_chunks_",
                        chunkname,
                        "files.\n"
                    )
                )))
            }
            assign(nobj, LDres2)
            save(list = nobj,
                 file = file.path(outpath, filename))
        }
        
        # FLAG SCRIPT END
        
        if (verbose >= 1) {
            cat(report("Completed:", funname, "\n"))
        }
        
        return(LDres2)
        
    }
}
