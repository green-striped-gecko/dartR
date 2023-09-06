#' @name utils.dart2genlight
#' Converts DarT to genlight.
#' Internal function called by gl.read.dart()
#' @description Converts a DArT file (read via \code{read.dart}) into an
#' genlight object \code{\link{adegenet}}. 
#' @param dart A dart object created via read.dart [required].
#' @param ind.metafile Optional file in csv format with metadata for each
#' individual (see details for explanation) [default NULL].
#' @param covfilename Depreciated, use parameter ind.metafile.
#' @param probar Show progress bar [default TRUE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report [default NULL].
#' @return A genlight object. Including all available slots are filled.
#' loc.names, ind.names, pop, lat, lon (if provided via the ind.metadata file)
#' @details
#' The ind.metadata file needs to have very specific headings. First a heading
#' called id. Here the ids have to match the ids in the dart object
#' \code{colnames(dart[[4]])}. The following column headings are optional.
#' pop: specifies the population membership of each individual. lat and lon
#' specify spatial coordinates (in decimal degrees WGS1984 format). Additional
#' columns with individual metadata can be imported (e.g. age, gender).
#' 
#'@author Custodian: Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})

utils.dart2genlight <- function(dart,
                                ind.metafile = NULL,
                                covfilename = NULL,
                                probar = TRUE,
                                verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbose = verbose)
    
    # DO THE JOB
    
    if (is.null(ind.metafile)) {
        ind.metafile <- covfilename
    }
    
    #### out contains the dart data
    nind <- dart[["nind"]]
    nsnp <- dart[["nsnp"]]
    sraw <- dart[["covmetrics"]]
    nrows <- dart[["nrows"]]  #check if nrows are provided...
    service <- as.character(unlist(unname(dart[["service"]])))
    plate_location <-
        as.character(unlist(unname(dart[["plate_location"]])))
    
    if (is.null(nrows)) {
        cat(report(
            "nrows not provided. Trying to guess if one row or two row format...\n"
        ))
        gnrows <-3 - max(dart$gendata, na.rm = TRUE)
        
        if (gnrows == 1 | gnrows == 2) {
            nrows <- gnrows
            cat(report(
                paste(
                    "Should be ",
                    nrows,
                    " row(s) format. Please check if this is the case. Trying to proceed...\n"
                )
            ))
        } else {
            stop(
                error(
                    "Cannot be guessed. The dart format must be either one row or two row format and needs to be provided via nrows=1 or 2.\n"
                )
            )
        }
    }
    
    if (sum(c("SNP", "SnpPosition") %in% names(sraw)) != 2) {
        stop(error(
            "Could not find SNP or SnpPosition in Dart file. Check you headers!!!"
        ))
    }
    
    if (verbose >= 2) {
        cat(report("Starting conversion....\n"))
        cat(report(paste0("Format is ", nrows, " rows.\n")))
        cat(report(
            "Please note conversion of bigger data sets will take some time!\n"
        ))
        cat(
            report(
                "Once finished, we recommend to save the object using save(object, file=\"object.rdata\")\n"
            )
        )
    }
    
    if (probar) {
        pb <- txtProgressBar(
            min = 0,
            max = 1,
            style = 3,
            initial = NA
        )
    }
    
    sdata <- dart[["gendata"]]
    # every second line only....
    esl <-seq(nrows, nrow(sdata), nrows)
    
    pos <- sraw$SnpPosition[esl]
    alleles <- as.character(sraw$SNP)[esl]
    a1 <- substr(alleles, nchar(alleles) - 2, nchar(alleles))
    a2 <- sub(">", "/", a1)
    locname <- paste(sraw$uid[esl], a2, sep = "-")
    geninddata <- matrix(NA, nrow = nsnp, ncol = nind)
    
    if (nrows == 2) {
        for (i in 1:nind) {
            isnp <-paste(sdata[esl - 1, i], sdata[esl, i], sep = "/")
            g <- isnp
            g <- gsub("0/1", 2, g)
            g <- gsub("1/0", 0, g)
            g <- gsub("1/1", 1, g)
            g <- gsub("NA/NA", NA, g)
            geninddata[, i] <- as.numeric(g)
            if (probar) {
                setTxtProgressBar(pb, i / nind)
            }
        }
    } else {
        for (i in 1:nind) {
            isnp <-sdata[esl, i]
            g <- isnp
            g <- 3 - g
            g <- ifelse(g == 3, 0, g)
            geninddata[, i] <- g
            if (probar) {
                setTxtProgressBar(pb, i / nind)
            }
        }
    }
    gout <-
        new(
            "dartR",
            gen = t(geninddata),
            ploidy = 2,
            ind.names = colnames(sdata),
            loc.names = locname,
            loc.all = a2,
            position = pos,
            parallel = FALSE
        )
    
    if (probar) {
        close(pb)
    }
    
    # refactor data.frame
    df <-
        as.data.frame(lapply(sraw[esl, ], function(x)
            if (is.factor(x))
                factor(x)
            else
                x))
    
    gout@other$loc.metrics <- df
    df_ind.metrics <-
        as.data.frame(matrix(nrow = nInd(gout), ncol = 2))
    colnames(df_ind.metrics) <- c("service", "plate_location")
    gout@other$ind.metrics <- df_ind.metrics
    
    # additional metadata and long lat to the data file are stored in other
    
    if (!is.null(ind.metafile)) {
        if (verbose >= 2) {
            cat(report(
                paste("Adding individual metrics:", ind.metafile, ".\n")
            ))
        }
        ###### population and individual file to link AAnumbers to populations...
        ind.cov <-
            read.csv(ind.metafile,
                     header = T,
                     stringsAsFactors = T)
        # is there an entry for every individual
        
        id.col <-match("id", names(ind.cov))
        
        if (is.na(id.col)) {
            stop(error("Fatal Error: There is no id column\n"))
        } else {
            ind.cov[, id.col] <-
                trimws(ind.cov[, id.col], which = "both")  #trim spaces
            
            if (length(ind.cov[, id.col]) != length(unique(ind.cov[, id.col]))) {
                cat(error(
                    "Individual names are not unique. You need to change them!\n"
                ))
                stop()
            }
            
            # reorder
            if (length(ind.cov[, id.col]) != length(names(sdata))) {
                cat(
                    warn(
                        "Ids for individual metadata does not match the number of ids in the SNP data file. Maybe this is fine if a subset matches.\n"
                    )
                )
                nam.indmeta <- ind.cov[, id.col]
                nam.dart <- names(sdata)
                
                nm.indmeta <-
                    nam.indmeta[!nam.indmeta %in% nam.dart]
                nm.inddart <-
                    nam.dart[!nam.dart %in% nam.indmeta]
                if (length(nm.indmeta) > 0) {
                    cat(warn("ind.metafile ids not matched were:\n"))
                    print(nm.indmeta)
                }
                if (length(nm.inddart) > 0) {
                    cat(warn("DArT file ids not matched were:\n"))
                    print(nm.inddart)
                }
            }
            
            ord <- match(names(sdata), ind.cov[, id.col])
            ord <- ord[!is.na(ord)]
            
            if (length(ord) > 1 & length(ord) <= nind) {
                if (verbose >= 2) {
                    cat(report(
                        paste(
                            "  Ids for individual metadata (at least a subset of) are matching!\n"
                        )
                    ))
                    cat(report(
                        paste(
                            "  Found ",
                            length(ord == nind),
                            "matching ids out of",
                            nrow(ind.cov),
                            "ids provided in the ind.metadata file.\n "
                        )
                    ))
                }
                ord2 <-
                    match(ind.cov[ord, id.col], indNames(gout))
                gout <- gout[ord2, ]
            } else {
                stop(error(
                    "Fatal Error: Individual ids are not matching!!!!\n"
                ))
            }
        }
        
        pop.col <-match("pop", names(ind.cov))
        
        if (is.na(pop.col)) {
            if (verbose >= 1) {
                cat(
                    warn(
                        "Warning: There is no pop column, created one with all pop1 as default for all individuals\n"
                    )
                )
            }
            pop(gout) <- factor(rep("pop1", nInd(gout)))
        } else {
            pop(gout) <- as.factor(ind.cov[ord, pop.col])
            if (verbose >= 2) {
                cat(report(" Added population assignments.\n"))
            }
        }
        
        lat.col <-match("lat", names(ind.cov))
        lon.col <-match("lon", names(ind.cov))
        if (verbose >= 2) {
            if (is.na(lat.col)) {
                cat(
                    warn(
                        "Warning: Individual metrics do not include a latitude (lat) column\n"
                    )
                )
            }
            if (is.na(lon.col)) {
                cat(
                    warn(
                        "Warning: Individual metrics do not include a longitude (lon) column\n"
                    )
                )
            }
        }
        if (!is.na(lat.col) & !is.na(lon.col)) {
            gout@other$latlon <- ind.cov[ord, c(lat.col, lon.col)]
            rownames(gout@other$latlon) <- ind.cov[ord, id.col]
            if (verbose >= 2) {
                cat(report("  Added latlon data.\n"))
            }
        }
        
        # known.col <- names( ind.cov) %in% c('id','pop', 'lat', 'lon') known.col <- ifelse(is.na(known.col), , known.col) other.col <-
        # names(ind.cov)[!known.col]
        other.col <- names(ind.cov)
        if (length(other.col) > 0) {
            gout@other$ind.metrics <- ind.cov[ord, other.col, drop = FALSE]
            rownames(gout@other$ind.metrics) <-
                ind.cov[ord, id.col]
            if (verbose >= 2) {
                cat(report(
                    paste(
                        " Added ",
                        other.col,
                        " to the other$ind.metrics slot.\n"
                    )
                ))
            }
        }
    }
    
    ord3 <- match(indNames(gout), names(sdata))
    
    gout@other$ind.metrics$service <- service[ord3]
    gout@other$ind.metrics$plate_location <-
        plate_location[ord3]
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report(paste("Completed:", funname, "\n")))
    }
    
    return(gout)
    
}

