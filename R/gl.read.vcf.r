#' Converts a vcf file into a genlight object
#'
#' This function needs package vcfR, please install it. 
#' @param vcffile A vcf file (works only for diploid data) [required].
#' @param ind.metafile Optional file in csv format with metadata for each
#' individual (see details for explanation) [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @details
#' The ind.metadata file needs to have very specific headings. First a heading
#' called id. Here the ids have to match the ids in the dartR object. 
#' The following column headings are optional.
#' pop: specifies the population membership of each individual. lat and lon
#' specify spatial coordinates (in decimal degrees WGS1984 format). Additional
#' columns with individual metadata can be imported (e.g. age, gender).
#' @return A genlight object.
#' @export
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' \dontrun{
#' obj <- gl.read.vcf(system.file('extdata/test.vcf', package='dartR'))
#' }

gl.read.vcf <- function(vcffile,
                        ind.metafile = NULL,
                        verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jackson",
                     verbose = verbose)
    
    x <- NULL
    
    pkg <- "vcfR"
      if (!(requireNamespace(pkg, quietly = TRUE))) {
        cat(error(
          "Package",
          pkg,
          " needed for this function to work. Please install it.\n"
        ))
        return(-1)
      } 
    
        vcf <- vcfR::read.vcfR(file = vcffile, verbose = verbose)
        myRef <- vcfR::getREF(vcf)
        myAlt <- vcfR::getALT(vcf)
        chrom <- vcfR::getCHROM(vcf)
        pos <- vcfR::getPOS(vcf) 
        loc.all <- paste0(myRef,"/",myAlt)
        x <- vcfR::vcfR2genlight(vcf)

        # adding SNP information from VCF
        info_tmp_1 <- vcf@fix[,6:7]
        info_tmp_2 <- vcfR::getINFO(vcf)
        if(any(is.na(info_tmp_2[1]) | is.na(info_tmp_1[1]))==TRUE){
          info <- info_tmp_1
          colnames(info) <- c("QUAL","FILTER")
        }else{
          info_tmp_2 <- as.data.frame(do.call(rbind,stringr::str_split(info_tmp_2,pattern = "=|;")))
          info <- info_tmp_2[,seq(2,ncol(info_tmp_2),2)]
          info <- cbind(info_tmp_1,info)
          col.names.info <- c("QUAL","FILTER",unname(unlist(info_tmp_2[1,seq(1,ncol(info_tmp_2),2)])))
          if(length(col.names.info)!=  length(colnames(info))){
            message(warn(
              "  Locus information is not formatted correctly. One reason could be that a field could have missing values."))
            info <- info_tmp_1
            colnames(info) <- c("QUAL","FILTER")
          }else{
          colnames(info) <- col.names.info
          }
        }
      
        # identify which SNPs have more than 2 alleles
        if("AC" %in% colnames(info)){
          more_alleles <- grep(",",info$AC)
          if(length(more_alleles)!=0){
            info <- info[-more_alleles,]
            info[] <- lapply(info, as.numeric)
            x@loc.all <- loc.all[-more_alleles]
            x@chromosome <- as.factor(chrom[-more_alleles])
            x@position <- pos[-more_alleles]
          }else{
            x@loc.all <- loc.all
            x@chromosome <- as.factor(chrom)
            x@position <- pos
          }
        }else{
        ALT <- vcfR::getALT(vcf)
        more_alleles <- grep(pattern = ",",ALT)
        if(length(more_alleles)>0){
          info <- info[-more_alleles,]
          x@loc.all <- loc.all[-more_alleles]
          x@chromosome <- as.factor(chrom[-more_alleles])
          x@position <- pos[-more_alleles]
        }else{
          x@loc.all <- loc.all
          x@chromosome <- as.factor(chrom)
          x@position <- pos
        }
        }
        
        ploidy(x) <- 2
        x <- gl.compliance.check(x)
        
        x$other$loc.metrics <- cbind(x$other$loc.metrics,info)
        x$other$loc.metrics$QUAL <- as.numeric(x$other$loc.metrics$QUAL)
        
        # additional metadata and long lat to the data file are stored in other
        
        if (!is.null(ind.metafile)) {
          if (verbose >= 2) {
            cat(report(
              paste("Adding individual metrics:", ind.metafile, ".\n")
            ))
          }
          ###### population and individual file to link numbers to populations...
          ind.cov <- read.csv(ind.metafile,  
                              header = TRUE, 
                              stringsAsFactors = TRUE)
          # is there an entry for every individual
          
          id.col <- match("id", names(ind.cov))
          
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
            if (length(ind.cov[, id.col]) != length(indNames(x))) {
              cat(
                warn(
                  "Ids for individual metadata does not match the number of ids in the SNP data file. Maybe this is fine if a subset matches.\n"
                )
              )
              nam.indmeta <- ind.cov[, id.col]
              nam.dart <- indNames(x)
              
              nm.indmeta <- nam.indmeta[!nam.indmeta %in% nam.dart]
              nm.inddart <- nam.dart[!nam.dart %in% nam.indmeta]
              if (length(nm.indmeta) > 0) {
                cat(warn("ind.metafile ids not matched were:\n"))
                print(nm.indmeta)
              }
              if (length(nm.inddart) > 0) {
                cat(warn("DArT file ids not matched were:\n"))
                print(nm.inddart)
              }
            }
            
            ord <- match(indNames(x), ind.cov[, id.col])
            ord <- ord[!is.na(ord)]
            
            if (length(ord) > 1 & length(ord) <= nInd(x)) {
              if (verbose >= 2) {
                cat(report(
                  paste(
                    "  Ids for individual metadata (at least a subset of) are matching!\n"
                  )
                ))
                cat(report(
                  paste(
                    "  Found ",
                    length(ord == nInd(x)),
                    "matching ids out of",
                    nrow(ind.cov),
                    "ids provided in the ind.metadata file.\n "
                  )
                ))
              }
              ord2 <- match(ind.cov[ord, id.col], indNames(x))
              x <- x[ord2, ]
            } else {
              stop(error(
                "Fatal Error: Individual ids are not matching!!!!\n"
              ))
            }
          }
          
          pop.col <- match("pop", names(ind.cov))
          
          if (is.na(pop.col)) {
            if (verbose >= 1) {
              cat(
                warn(
                  "  Warning: There is no pop column, created one with all pop1 as default for all individuals\n"
                )
              )
            }
            pop(x) <- factor(rep("pop1", nInd(x)))
          } else {
            pop(x) <- as.factor(ind.cov[ord, pop.col])
            if (verbose >= 2) {
              cat(report(" Added population assignments.\n"))
            }
          }
          
          lat.col <- match("lat", names(ind.cov))
          lon.col <- match("lon", names(ind.cov))
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
            x@other$latlon <- ind.cov[ord, c(lat.col, lon.col)]
            rownames(x@other$latlon) <- ind.cov[ord, id.col]
            if (verbose >= 2) {
              cat(report("  Added latlon data.\n"))
            }
          }
          
          other.col <- names(ind.cov)
          if (length(other.col) > 0) {
            x@other$ind.metrics <- ind.cov[ord, other.col, drop = FALSE]
            rownames(x@other$ind.metrics) <- ind.cov[ord, id.col]
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
        
        # add history
        x@other$history <- list(match.call())
        x <- gl.recalc.metrics(x)
        
        if (verbose > 2) {
            cat(
                important(
                    "Genlight object does not have individual metrics. You need to add them 'manually' to the @other$ind.metrics slot.\n"
                )
            )
        }
        return(x)
    
}
