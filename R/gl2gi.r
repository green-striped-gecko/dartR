#' Converts a genlight object to genind object
#' @param x A genlight object [required].
#' @param probar If TRUE, a progress bar will be displayed for long loops
#' [default TRUE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @return A genind object, with all slots filled.
#' @export
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @details This function uses a faster version of df2genind (from the adegenet
#'  package)

gl2gi <- function(x,
                  probar = FALSE,
                  verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # convert to genind....
    x2 <- as.matrix(x[,])
    
    if (probar) {
        pb <- txtProgressBar(
            min = 0,
            max = 1,
            style = 3,
            initial = NA
        )
    }
    
    if (is.null(x@loc.all)) {
        x@loc.all <- rep("A/T", nLoc(x))
        x@loc.all[1] <- "C/G"
    }
    
    
    homs1 <-
        paste(substr(x@loc.all, 1, 1), "/", substr(x@loc.all, 1, 1), sep = "")
    hets <- x@loc.all
    homs2 <-
        paste(substr(x@loc.all, 3, 3), "/", substr(x@loc.all, 3, 3), sep = "")
    xx <- matrix(NA, ncol = ncol(x2), nrow = nrow(x2))
    for (i in 1:nrow(x2)) {
        for (ii in 1:ncol(x2)) {
            inp <- x2[i, ii]
            if (!is.na(inp)) {
                if (inp == 0)
                    xx[i, ii] <- homs1[ii]
                else if (inp == 1)
                    xx[i, ii] <- hets[ii]
                else if (inp == 2)
                    xx[i, ii] <- homs2[ii]
            } else
                xx[i, ii] <-"-/-"
        }
        if (probar) {
            setTxtProgressBar(pb, i / nrow(x2))
        }
    }
    if (probar) {
        close(pb)
    }
    
    if (verbose >= 1) {
        cat(report("Matrix converted.. Prepare genind object...\n"))
    }
    
    gen <-
        df2genind(
            xx[,],
            sep = "/",
            ncode = 1,
            ind.names = x@ind.names,
            pop = x@pop,
            ploidy = 2,
            NA.char = "-"
        )  #, probar=probar)
    gen@other <- x@other
    locNames(gen) <- locNames(x)
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    # RETURN
    
    return(gen)
    
}
