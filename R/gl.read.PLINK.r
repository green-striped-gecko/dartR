#' @name gl.read.PLINK
#' @title Reads PLINK data file into a genlight object
#' @description This function imports PLINK data into a genlight object and append available metadata. 
#' 
#' @details  Additional metadata can be included
#' passing .csv files. These will be appended to the existing metadata present in the PLINK files. This 
#' function handles .ped or .bed file (with the associate files - e.g. .fam, .bim). However, if a .ped 
#' file is provided, PLINK needs to be installed and it is used to convert 
#' the .ped into a .bed, which is then converted into a genlight.
#' 
#' The locus metadata needs to be in a csv file with headings, with a mandatory
#' column headed AlleleID corresponding exactly to the locus identity labels
#' provided with the SNP data
#'
#' @param filename Fully qualified path to PLINK input file (without including the extension)
#' @param ind.metafile Name of the csv file containing the metrics for
#' individuals [optional].
#' @param loc.metafile Name of the csv file containing the metrics for
#' loci [optional].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @inheritParams utils.plink.run
#' @return A genlight object with the SNP data and associated metadata included.
#' @export
#' @author Custodian: Carlo Pacioni -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @rawNamespace import(data.table, except = c(melt,dcast))
#' @examples

gl.read.PLINK <- function(filename,
                        ind.metafile = NULL,
                        loc.metafile = NULL,
                        plink.cmd="plink", 
                        plink.path="path",
                        verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jackson",
                     verbosity = verbose)
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    # check if packages are installed
    pkg <- "snpStats"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      stop(error(
        "Package",
        pkg,
        " needed for this function to work. Please install it."
      ))
    }
    
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
    dir.out <- tempdir() # use a tmp dir to handle transformation
    plink.fns <- list.files(path=dirname(filename), pattern=paste0("^", basename(filename)), full.names=TRUE, recursive=TRUE)
    plink.fn <- plink.fns[grep(".bed$", plink.fns)]
    if(length(plink.fn)>1) { 
      stop(
      error(
        "Found more than one .bed file and don't know which one to use\n"
       )
     )
    } else {
      if(length(plink.fn) == 0) {
        plink.fn <- plink.fns[grep(".ped$", plink.fns)]
        if(length(plink.fn) == 0 | length(plink.fn)>1) {
          stop(
            error(paste(
              "Found no .bed files and", length(plink.fn), ".ped file(s). This function needs either one .bed or .ped file\n"
            )
            )
          )
      } else {
        if(length(plink.fn) == 1) {
          utils.plink.run(dir.in = dirname(filename), plink.cmd = plink.cmd, plink.path = plink.path, out = basename(filename), 
                          syntax = paste("--file", basename(filename), "--make-bed"))
          plink.fn <- list.files(path=dirname(filename), pattern=".bed$", full.names=TRUE, recursive=TRUE)
        } else {
          stop(
            error(
              "Couldn't find any .bed or .ped file\n"
            )
          )
        }
        
      }
    }
    
    snpMatrix <- snpStats::read.plink(bed=sub(x=plink.fn, pattern=".bed$", ""))
    gen <- snpMatrix[[1]]
    fam <- snpMatrix[[2]]
    map <- snpMatrix[[3]]
    row.names(gen) <- fam$member
    write.SnpMatrix(gen, file.path(dir.out, "snpMatrixComb.txt"))
    
    suppressWarnings(
      genCombdt <- data.table::fread(file=file.path(dir.out, "snpMatrixComb.txt"))
    )
    setnames(genCombdt, "V1", "id")
    
    genCombNA <- as.matrix(genCombdt[, -1, with=FALSE])
    row.names(genCombNA) <- fam$member
    gl <- new("genlight", gen=genCombNA,
               ind.names=fam$member,
               loc.names=map$snp.name,
               chromosome=map$chromosome,
               position=map$position, # Need to confirm that we want here the chr position and that the SNP position on the read goes only in the loc.metrics
               other=list(ind.metrics=cbind(id=fam$member, 
                                            snpStats::row.summary(gen), # This adds Call.rate, Certain.calls, Heterozygosity
                                            Family=fam$pedigree,
                                            Father=fam$father,
                                            Mother=fam$mother,
                                            Sex=fam$sex),
                          loc.metrics=map)
    )
   
    if (length(unique(fam$member)) != length(fam$member)) {
        cat(
            error(
                "Fatal Error: Individual labels are not unique, check and edit your input file\n"
            )
        )
        stop()
    }
    if (length(unique(map$snp.name)) != length(map$snp.name)) {
        cat(error(
            "Fatal Error: AlleleID not unique, check and edit your input file\n"
        ))
        stop()
    }
    

    pop(gl) <- array("A", nInd(gl))
    names(gl@other$loc.metrics)[2] <- "AlleleID"
    gl@other$loc.metrics <- gl@other$loc.metrics[, c("AlleleID", names(gl@other$loc.metrics)[-2])] # Fix order of cols
    gl <- gl.compliance.check(gl, verbose = verbose)
    
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
        
        if(nrow(loc.metrics) != nLoc(gl)) stop(
          error(
            "Fatal Error: the locus metrics file does not have the same number of loci of the input data file\n"
          )
        )
        
        if(!all(loc.metrics[, "AlleleID"] %in% gl@other$loc.metrics$AlleleID))
                stop(
                    error(
                        "Fatal Error: AlleleID in the locus metrics file does not correspond with",
                        "AlleleID in the input data file\n"
                    )
                )
            }
        }
        row.names(loc.metrics) <- loc.metrics$AlleleID
        gl@other$loc.metrics <- cbind(gl@other$loc.metrics, loc.metrics[map$snp.name,])
        
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
        fam$member %in% ind.metrics[, "id"] 
        if(sum(fam$member %in% ind.metrics[, "id"]) < length(fam$member)) stop(
          error(
            paste(
            "Fatal Error: there are", sum(fam$member %in% ind.metrics[, "id"]), 
            "individuals id that match the ones in the data, but", length(fam$member),
            "individuals genotyped\n")
          )
        )
          
        
        if (!("pop" %in% names(ind.metrics))) {
            cat(
                warn(
                    "  Warning: pop column absent from the individual metadata file, setting to 'A'\n"
                )
            )
            ind.metrics$pop <- array("A", nInd(gl))
        }
        
        gl@other$ind.metrics <- cbind(gl@other$ind.metrics, ind.metrics)
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
