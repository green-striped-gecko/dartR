#' Reads SNP data from a csv file into a genlight object
#'
#' This script takes SNP genotypes from a csv file, combines them with
#' individual and locus metrics and creates a genlight object.
#'
#' The SNP data need to be in one of two forms. SNPs can be coded 0 for
#' homozygous reference, 2 for homozygous alternate, 1 for heterozygous, and NA 
#' for missing values; or the SNP data can be coded A/A, A/C, C/T, G/A etc,
#' and -/- for missing data. In this format, the reference allele is the most 
#' frequent allele, as used by DArT. Other formats will throw an error.
#'
#' The SNP data need to be individuals as rows, labeled, and loci as columns,
#' also labeled. If the orientation is individuals as columns and loci by rows,
#'  then set transpose=TRUE.
#'
#' The individual metrics need to be in a csv file, with headings, with a
#'  mandatory id column corresponding exactly to the individual identity labels
#'  provided with the SNP data and in the same order.
#'
#' The locus metadata needs to be in a csv file with headings, with a mandatory
#' column headed AlleleID corresponding exactly to the locus identity labels
#' provided with the SNP data and in the same order.
#'
#' Note that the locus metadata will be complemented by calculable statistics
#' corresponding to those that would be provided by Diversity Arrays Technology
#' (e.g. CallRate).
#'
#' @param filename Name of the csv file containing the SNP genotypes [required].
#' @param transpose If TRUE, rows are loci and columns are individuals
#' [default FALSE].
#' @param ind.metafile Name of the csv file containing the metrics for
#' individuals [optional].
#' @param loc.metafile Name of the csv file containing the metrics for
#' loci [optional].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @return A genlight object with the SNP data and associated metadata included.
#' @export
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' csv_file <- system.file('extdata','platy_test.csv', package='dartR')
#' ind_metadata <- system.file('extdata','platy_ind.csv', package='dartR')
#' gl  <- gl.read.csv(filename = csv_file, ind.metafile = ind_metadata)

gl.read.csv <- function(filename,
                        transpose = FALSE,
                        ind.metafile = NULL,
                        loc.metafile = NULL,
                        verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jackson",
                     verbosity = verbose)
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    if (is.null(loc.metafile) & verbose > 0) {
        cat(
            warn(
                "Warning: Locus metafile not provided, locus metrics will be
        calculated where this is possible\n"
            )
        )
    }
    
    if (is.null(ind.metafile) & verbose > 0) {
        cat(
            warn(
                "Warning: Individual metafile not provided, pop set to 'A' for all individuals\n"
            )
        )
    }
    
    # DO THE JOB
    
    # FIRST THE SNP DATA
    
    # Create the SNP data matrix, indNames and LocNames
    
    df0 <-
        read.csv(file = filename,
                 header = FALSE,
                 stringsAsFactors = TRUE)
    
    if (transpose) {
        df0 <- t(df0)
    }
    
    numrows <- dim(df0)[1]  # Individuals plus labels if any
    numcols <- dim(df0)[2]  # Loci plus labels if any
    
    if (verbose > 0) {
        cat(
            report(
                "Input data should be a csv file with individuals as rows, loci as columns\n"
            )
        )
        cat("  ",
            numcols - 1,
            "loci, confirming first 5:",
            as.matrix(df0[1, 2:6]),
            "\n")
        cat("  ",
            numrows - 1,
            "individuals, confirming first 5:",
            as.matrix(df0[2:6, 1]),
            "\n")
        cat(important(
            "    If these are reversed, re-run the script with transpose=TRUE\n"
        ))
    }
    data <- as.matrix(df0[2:numrows, 2:numcols])
    
    loci <- df0[1, 2:numcols]
    loci <- as.character(as.matrix(loci))
    individuals <- df0[2:numrows, 1]
    individuals <- as.character(individuals)
    
    if (length(unique(individuals)) != length(individuals)) {
        cat(
            error(
                "Fatal Error: Individual labels are not unique, check and edit your input file\n"
            )
        )
        stop()
    }
    if (length(unique(loci)) != length(loci)) {
        cat(error(
            "Fatal Error: AlleleID not unique, check and edit your input file\n"
        ))
        stop()
    }
    
    # Validate and convert the SNP data
    
    test <- paste0(data, collapse = "")
    test <- gsub("NA", "9", test)
    test <- gsub(" ", "", test)
    if (nchar(test) > nrow(data) * ncol(data)) {
        if (verbose >= 2) {
            cat(
                report(
                    "Character data detected, assume genotypes are of the
            form C/C, A/T, C/G, -/- etc\n"
                )
            )
        }
        # Check that this is true
        s1 <- paste(data, collapse = " ")
        s1 <- gsub("/", " ", s1)
        s1 <- toupper(s1)
        s2 <- unlist(strsplit(s1, " "))
        tmp <- table(s2)
        if (all(names(tmp) %in% c("A", "C", "G", "T", "-")) == F) {
            cat(
                error(
                    "Fatal Error: Genotypes must be defined by the letters A, C, G, T or missing -\n"
                )
            )
            stop()
        }
        # Check that the data are bi-allelic
        for (i in 1:dim(data)[2]) {
            v1 <- data[, i]
            v1 <- paste(v1, collapse = " ")
            v1 <- gsub("/", " ", v1)
            v1 <- gsub("- ", "", v1)
            v1 <- toupper(v1)
            v1 <- unlist(strsplit(v1, " "))
            tmp <- table(v1)
            tmp <- tmp[order(as.numeric(tmp),decreasing = T)]
            if (length(names(tmp)) > 2) {
                cat(error("Fatal Error: Loci are not bi-allelic\n"))
                stop()
            }
            
            # Step through and convert data to 0, 1, 2, NA
            homRef <- paste0(names(tmp)[1],"/",names(tmp)[1])
            homAlt <- paste0(names(tmp)[2],"/",names(tmp)[2])
            het1 <- paste0(names(tmp)[1],"/",names(tmp)[2])
            het2 <- paste0(names(tmp)[2],"/",names(tmp)[1])
            missing <- "-/-"
  
                data[,i] <- gsub(homRef,"0",data[,i])
                data[,i] <- gsub(homAlt,"2",data[,i])
                data[,i] <- gsub(het1,"1",data[,i])
                data[,i] <- gsub(het2,"1",data[,i])
                data[,i] <- gsub(missing,NA,data[,i])
        
        }
        if (verbose >= 2) {
            cat(report("  Data confirmed as biallelic\n"))
        }
        
        if (verbose >= 2) {
            cat(report("  SNP coding converted to 0, 1, 2 and NA\n"))
        }
        
        data <- apply(data, 2, as.numeric)
        
    } else {
        if (verbose >= 2) {
            cat(
                report(
                    "  Numeric data detected, assume genotypes are 0 = homozygous reference, 1 = heterozygous, 2 = homozygous alternate\n"
                )
            )
        }
        # Check that this is true
        data <- apply(data, 2, as.numeric)
        s1 <- paste(data, collapse = " ")
        s2 <- unlist(strsplit(s1, " "))
        tmp <- table(s2)
        if (!(names(tmp) == "0" ||
              names(tmp) == "1" ||
              names(tmp) == "2" || names(tmp) == "NA")) {
            cat(
                error(
                    "Fatal Error: Genotypes must be defined by the numbers 0, 1, 2 or missing NA\n"
                )
            )
            stop()
        }
    }
    
    # Create a genlight object
    
    gl <-
        new(
            "genlight",
            data,
            ploidy = 2,
            loc.names = loci,
            ind.names = individuals
        )
    
    pop(gl) <- array("A", nInd(gl))
    gl <- gl.compliance.check(gl, verbose = verbose)
    # gl@other$loc.metrics <- data.frame(CloneID = locNames(gl), AlleleID = locNames(gl))
    gl@other$ind.metrics <-
        data.frame(id <-
                       indNames(gl), pop = array("A", nInd(gl)))
    
    # NOW THE LOCUS METADATA
    
    if (!is.null(loc.metafile)) {
        loc.metrics <-
            read.csv(
                file = loc.metafile,
                header = TRUE,
                stringsAsFactors = TRUE
            )
        if (!("AlleleID" %in% names(loc.metrics))) {
            cat(
                error(
                    "Fatal Error: mandatory AlleleID column absent
                                               from the locus metrics file\n"
                )
            )
        }
        for (i in 1:nLoc(gl)) {
            if (loc.metrics[i, 1] != gl@other$loc.metrics$AlleleID[i]) {
                stop(
                    error(
                        "Fatal Error: AlleleID in the locus metrics file does not correspond with",
                        "AlleleID in the input data file, or they are not in the same order\n"
                    )
                )
            }
        }
        gl@other$loc.metrics <- loc.metrics
        
    }
    gl <- gl.recalc.metrics(gl, verbose = 0)
    if (verbose >= 2) {
        cat(report(
            paste(
                " Added or updated ",
                names(gl@other$loc.metrics),
                "to the other$ind.metrics slot.\n"
            )
        ))
    }
    
    # NOW THE INDIVIDUAL METADATA
    
    if (!is.null(ind.metafile)) {
        ind.metrics <-
            read.csv(
                file = ind.metafile,
                header = TRUE,
                stringsAsFactors = TRUE,
                fileEncoding = "UTF-8-BOM"
            )
        if (!("id" %in% names(ind.metrics))) {
            cat(
                error(
                    "Fatal Error: mandatory id column absent from the individual metadata file\n"
                )
            )
            stop()
        }
        for (i in 1:nInd(gl)) {
            if (ind.metrics[i, 1] != gl@other$ind.metrics$id[i]) {
                cat(
                    error(
                        "Fatal Error: id in the individual metrics file does not correspond with",
                        "id in the input data file, or they are not in the same order\n"
                    )
                )
                stop()
            }
        }
        if (!("pop" %in% names(ind.metrics))) {
            cat(
                warn(
                    "  Warning: pop column absent from the individual metadata file, setting to 'A'\n"
                )
            )
            
            gl@other$ind.metrics <- ind.metrics
            gl@other$ind.metrics$id <- individuals
            gl@other$ind.metrics$pop <- array("A", nInd(gl))
            pop(gl) <- gl@other$ind.metrics$pop
        } else {
            gl@other$ind.metrics <- ind.metrics
            gl@other$ind.metrics$id <- individuals
            gl@other$ind.metrics$pop <- ind.metrics$pop
            pop(gl) <- gl@other$ind.metrics$pop
        }
        if (verbose >= 2) {
            cat(report(
                paste(
                    " Added ",
                    names(gl@other$ind.metrics),
                    " to the other$ind.metrics slot.\n"
                )
            ))
        }
    }
    
    # MAKE COMPLIANT
    gl <- gl.compliance.check(gl, verbose = verbose)
    
    # ADD TO HISTORY (add the first entry)
    gl@other$history <- list()
    gl@other$history[[1]] <- match.call()
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(gl)
    
}
