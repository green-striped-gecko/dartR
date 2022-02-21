#' @name gl2geno
#' @title Converts a genlight object to geno format from package LEA
#' @description
#' The function converts a genlight object (SNP or presence/absence
#'  i.e. SilicoDArT data) into a file in the 'geno' and the 'lfmm' formats from 
#'  (package LEA).
#' @param x Name of the genlight object containing the SNP or presence/absence
#'  (SilicoDArT) data [required].
#' @param outfile File name of the output file [default 'gl_geno'].
#' @param outpath Path where to save the output file
#' [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath='.'
#'  when calling this function to direct output files to your working directory.
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @return NULL
#' @author Custodian: Luis Mijangos (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' # SNP data
#' gl2geno(testset.gl)
#' # Tag P/A data
#' gl2geno(testset.gs)
#' @export

gl2geno <- function(x,
                    outfile = "gl_geno",
                    outpath = tempdir(),
                    verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbosity = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # FUNCTION SPECIFIC ERROR CHECKING
    hold <- x
    
    x <- gl.filter.allna(x, verbose = 0)
    
    if (nLoc(hold) > nLoc(x) & verbose > 0) {
        cat(warn(
            "  Loci with all missing data has been removed for conversion. \n"
        ))
    }
    
    # DO THE JOB
    
    outfilespec <- file.path(outpath, outfile)
    
    dat <- as.matrix(x)
    n = nInd(x)
    L = nLoc(x)
    
    # Convert allelic data into absence/presence data at each locus Results are stored in the 'dat.binary' object
    
    if (datatype == "SilicoDArT") {
        dat <- as.matrix(dat)
        dat.binary <- NULL
        for (j in 1:L) {
            allele <- sort(unique(dat[, j]))
            for (i in allele[allele >= 0]) {
                dat.binary <- cbind(dat.binary, dat[, j] == i)
            }
            LL <- ncol(dat.binary)
            ind <- which(dat[, j] < 0)
            if (length(ind) != 0) {
                dat.binary[ind, (LL - length(allele) + 2):LL] <- -9
            }
        }
    }
    
    if (datatype == "SNP") {
        dummy_final <- NULL
        for (r in 1:nrow(dat)) {
            dummy <- rbind(dat[r, ], dat[r, ])
            index <- colSums(dummy, na.rm = T) == 2
            dummy[, index] <- c(0, 2)
            dummy <- ifelse(is.na(dummy), -9, dummy)
            dummy <- ifelse(dummy == 0, 1, dummy)
            dummy_final <- rbind(dummy_final, dummy)
        }
        
        dat.binary <- NULL
        for (j in 1:L) {
            allele = sort(unique(dummy_final[, j]))
            for (i in allele[allele >= 0]) {
                dat.binary <- cbind(dat.binary, dummy_final[, j] == i)
            }
            LL <- dim(dat.binary)[2]
            ind <- which(dummy_final[, j] < 0)
            if (length(ind) != 0) {
                dat.binary[ind, (LL - length(allele) + 2):LL] <- -9
            }
        }
    }
    
    # Compute a genotype count for each allele (0,1,2 or 9 for a missing value) The results are stored in 'genotype'
    
    n <- nrow(dat.binary)
    L <- ncol(dat)
    LL <- ncol(dat.binary)
    
    if (datatype == "SNP") {
        n = n / 2
        genotype <- matrix(NA, nrow = n, ncol = LL)
        for (i in 1:n) {
            genotype[i, ] <- dat.binary[2 * i - 1, ] + dat.binary[2 * i, ]
            genotype[i, (genotype[i, ] < 0)] <- NA
        }
        if (LL == 2 * L) {
            genotype <- genotype[, seq(2, LL, by = 2)]
        }
    }
    
    if (datatype == "SilicoDArT") {
        genotype <- dat.binary
        for (i in 1:n) {
            genotype[i, (genotype[i, ] < 0)] <- NA
        }
        if (LL == 2 * L) {
            genotype <- genotype[, seq(2, LL, by = 2)]
        }
    }
    
    genotype[is.na(genotype)] <- 9
    lst.monomorphic <- apply(
        genotype,
        2,
        FUN = function(x) {
            length(unique(x[x != 9]))
        }
    )
    
    if (sum(lst.monomorphic == 1) > 0) {
        if (verbose > 0) {
            cat(warn(
                "  Monomorphic alleles generated during conversion were removed. \n"
            ))
        }
        genotype <- genotype[, lst.monomorphic > 1]
    }
    
    # write lfmm
    write.table(
        genotype,
        paste(outfilespec, ".lfmm", sep = ""),
        col.names = FALSE,
        row.names = FALSE,
        sep = " "
    )
    # write geno
    write.table(
        t(genotype),
        paste(outfilespec, ".geno", sep = ""),
        col.names = FALSE,
        row.names = FALSE,
        sep = ""
    )
    
    if (verbose > 0) {
        cat(report("  Output files:", paste(
            paste0(outfile, ".geno", ".lfmm."), sep = ""
        ), "\n"))
    }
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    invisible(NULL)
}
