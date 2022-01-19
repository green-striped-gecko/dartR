
#'Converts a genlight object into a sfs input file
#'
#'The output of this function is suitable for analysis in fastsimcoal2 or dada.
#'
#'It saves a derived sfs, assuming that the reference allele is the ancestral,
#'and a MAF sfs.
#'
#'At this stage this function caters only for diploid organisms, for samples
#'from one population only, and for genotypes without missing data.
#'
#'If no invariant site estimates is provided (with \code{n.invariant}), it will
#'estimate the number of invariant site from the sequenced tags (but this will
#'be a biased estimates as it doesn't take into account the invariant tags).
#'Note that the invariant sites can be estimated with
#'\code{gl.report.secondaries}.
#'
#'
#'It expects a dartR formatted genlight object, but it should also  work with
#'other genlight objects. If \code{n.invariant == 0} the genlight object needs
#'to have a \code{data.frame} in \code{other} called loc.metrics with a column
#'named \code{TrimmedSequence} that contains a character string whose length
#'corresponds to the number of sites of the allele. This information is used to
#'obtained the number of non-polymorphic sites.
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
                   n.invariant = 0,
                   outfile_root = "gl2sfs",
                   outpath = tempdir()) {
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
    if (n.invariant == 0) {
        if (isFALSE("TrimmedSequence" %in% names(x$other$loc.metrics))) {
            stop(
                error(
                    "The column 'TrimmedSequence' is not present in x$other$loc.metrics, but it is needed for this function"
                )
            )
        }
        inv <-
            sum(nchar(as.character(
                x$other$loc.metrics$TrimmedSequence
            ))) - x@n.loc
    } else {
        inv <- n.invariant
    }
    
    daf[1] <- daf[1] + inv
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
    maf[1] <- maf[1] + inv
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

