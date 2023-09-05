#' Converts a vcf file into a genlight object
#'
#' This function needs package vcfR, please install it. The converted genlight
#' object does not have individual metrics. You need to add them 'manually' to
#' the other$ind.metrics slot.
#' @param vcffile A vcf file (works only for diploid data) [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @return A genlight object.
#' @export
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' \dontrun{
#' obj <- gl.read.vcf(system.file('extdata/test.vcf', package='dartR'))
#' }

gl.read.vcf <- function(vcffile,
                        verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jackson",
                     verbosity = verbose)
    
    x <- NULL
    
    pkg <- "vcfR"
      if (!(requireNamespace(pkg, quietly = TRUE))) {
        cat(error(
          "Package",
          pkg,
          " needed for this function to work. Please install it.\n"
        ))
        return(-1)
      } else {
        vcf <- vcfR::read.vcfR(file = vcffile, verbose = verbose)
        myRef <- vcfR::getREF(vcf)
        myAlt <- vcfR::getALT(vcf)
        loc.all <- paste0(myRef,"/",myAlt)
        x <- vcfR::vcfR2genlight(vcf)
        x@loc.all <- loc.all
        
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
          colnames(info) <- c("QUAL","FILTER",unname(unlist(info_tmp_2[1,seq(1,ncol(info_tmp_2),2)])))
        }
      
        # identify which SNPs have more than 2 alleles
        if("AC" %in% colnames(info)){
          more_alleles <- grep(",",info$AC)
          info <- info[-more_alleles,]
          info[] <- lapply(info, as.numeric)
        }else{
        ALT <- vcfR::getALT(vcf)
        more_alleles <- which(stringr::str_length(ALT) >1)
        info <- info[-more_alleles,]
        }
        
        ploidy(x) <- 2
        x <- gl.compliance.check(x)
        
        x$other$loc.metrics <- cbind(x$other$loc.metrics,info)
        
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
    
}
