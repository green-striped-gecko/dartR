
#'Converts a genlight object into a sfs input file
#'
#'The output of this function is suitable for analysis in fastsimcoal2 or dada.
#'
#'It saves a derived sfs, assuming that the reference allele is the ancestral,
#'and a MAF sfs.
#'
#'At this stage this function caters only for diploid organisms, for samples
#'from one population only, and for genotypes without missing data. Note that
#'sfs uses frequencies considered \bold{independent}, data are assumed to be
#'from independent (i.e. not linked) loci. This means that only one site per tag
#'should be considered 9i.e. secondaries should be removed). If no monomorphic
#'site estimates is provided (with \code{n.invariant.tags}), the sfs will only
#'include the number of monomorphic sites in the data (but this will be a biased
#'estimates as it doesn't take into account the invariant tags that have not
#'been included. This will affect parameter estimates in the analyses). Note
#'that the number of invariant tags can be estimated with
#'\code{gl.report.secondaries}. In a limited number of cases, ascertainment bias
#'can be explicitly modelled in fastsimcoal2. See fastsimcoal2 manual for
#'details.
#'
#'
#'It expects a dartR formatted genlight object, but it should also  work with
#'other genlight objects.
#'
#'@param outfile_root The root of the name of the output file
#'@inheritParams gl.report.heterozygosity
#'@inheritParams gl2vcf
#'@return A list with two elements: the DAF and MAF.
#'@author Custodian: Carlo Pacioni (Post to
#'  \url{https://groups.google.com/d/forum/dartr})
#'@seealso \code{\link{gl.report.heterozygosity}},
#'  \code{\link{gl.report.secondaries}}, \code{\link{utils.n.var.invariant}}
#'@export
#'@references Excoffier L., Dupanloup I., Huerta-Sánchez E., Sousa V. C. and
#'  Foll M. (2013) Robust demographic inference from genomic and SNP data. PLoS
#'  genetics 9(10)

gl2sfs <- function(x,
                   n.invariant.tags = 0,
                   outfile_root = "gl2sfs",
                   outpath = tempdir(),
                   verbose = NULL) {
    #---------- Helper ---------------#
    count.freq.ref <- function(x) {
        if (x == 0)
            return(2)
        else if (x == 1)
            return(1)
        else
            return(0)
    }
    #---------------------------------#
    
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbosity = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    if (!all(x$ploidy == 2))
        stop(error("This function currently caters only for diploid organisms"))
    if (length(unique(x$pop)) > 1)
        stop(error("This function currently caters only for samples from one population only"))
    glm <- as.matrix(x)
    if (sum(is.na(glm)) > 0)
        stop(error("This function  currently caters only for genotypes without missing data"))
    
    daf <- table(colSums(glm))
    
    names(daf) <-
        paste("d0", seq(0, length(daf) - 1), sep = "_")
    
    daf[1] <- daf[1] + n.invariant.tags
    writeLines(c(
        " 1 observations",
        paste(paste("d0", seq(
            0, length(daf) - 1
        ), sep = "_"), collapse = " "),
        paste(as.character(daf), collapse = " ")
    ),
    con = file.path(outpath, paste0(outfile_root, "_DAFpop0.obs")))
    
    freq.ref <-
        colSums(apply(glm, FUN = count.freq.ref, MARGIN = c(1, 2)))
    freq.alt <- colSums(glm)
    maf <- table(c(freq.ref, freq.alt))
    maf <- maf[which(as.numeric(names(maf)) <= nrow(glm))]
    maf[1] <- maf[1] + n.invariant.tags
    if (length(maf) == nrow(glm))
        maf[length(maf)] <- maf[length(maf)] / 2
    writeLines(c(
        " 1 observations",
        paste(paste("d0", seq(
            0, length(maf) - 1
        ), sep = "_"), collapse = " "),
        paste(as.character(maf), collapse = " ")
    ),
    con = file.path(outpath, paste0(outfile_root, "_MAFpop0.obs")))
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(list(DAF = daf, MAF = maf))
}

