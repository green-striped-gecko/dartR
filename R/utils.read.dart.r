#' Imports DarT data to R
#'
#' Internal function called by gl.read.dart
#' @param filename Path to file (csv file only currently) [required].
#' @param nas A character specifying NAs [default '-'].
#' @param topskip A number specifying the number of rows to be skipped. If not
#' provided the number of rows to be skipped are 'guessed' by the number of rows
#' with '*' at the beginning [default NULL].
#' @param service_row The row number in which the information of the DArT
#' service is contained [default 1].
#' @param plate_row The row number in which the information of the plate
#' location is contained [default 3].
#' @param lastmetric Specifies the last non genetic column [default 'RepAvg'].
#' Be sure to check if that is true, otherwise the number of individuals will
#' not match. You can also specify the last column by a number.
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log ; 3, progress and results summary; 5, full report [default NULL].
#' @return A list of length 5. #dart format (one or two rows) #individuals,
#' #snps, #non genetic metrics, #genetic data (still two line format, rows=snps,
#'  columns=individuals)

utils.read.dart <- function(filename,
                            nas = "-",
                             topskip = NULL,
                             lastmetric = "RepAvg",
                            service_row = 1,
                            plate_row = 3,
                            verbose = NULL) {
  
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbosity = verbose)
    
    # DO THE JOB
    
    if (is.null(topskip)) {
        if (verbose >= 2) {
            cat(report("  Topskip not provided.\n "))
        }
        tdummy <-
            read.csv(
                filename,
                na.strings = nas,
                check.names = FALSE,
                nrows = 20,
                header = FALSE,
                stringsAsFactors = TRUE
            )
        
        nskip <- sum(tdummy[, 1] == "*")
        if (nskip > 0) {
            topskip <- nskip
            if (verbose >= 2) {
                cat(report(paste(
                    "Setting topskip to", nskip, ".\n"
                )))
            }
        } else {
            stop(
                error(
                    "Could not determine the number of rows that need to be skipped. Please provide it manually by setting the topskip parameter.\n"
                )
            )
        }
    }
    
    if (verbose >= 2) {
        cat(report("  Reading in the SNP data\n"))
    }
    snpraw <-
        read.csv(
            filename,
            na.strings = nas,
            skip = topskip,
            check.names = FALSE,
            stringsAsFactors = TRUE
        )
    
    if (is.character(lastmetric)) {
        lmet <- which(lastmetric == names(snpraw))
        if (length(lmet) == 0) {
            stop(error(
                paste(
                    "Could not determine number of data columns based on",
                    lastmetric,
                    "!\n"
                )
            ))
        }
    } else {
        lmet <- lastmetric
    }
    
    service <- NA 
    plate_location <- NA 
    if(exists("tdummy")){
    # extracting service information
    service <- tdummy[service_row, (lmet + 1):ncol(tdummy)]
    # extracting plate information
    plate <-
        unlist(unname(tdummy[plate_row, (lmet + 1):ncol(tdummy)]))
    plate_row_res <-
        unlist(unname(tdummy[(plate_row + 1), (lmet + 1):ncol(tdummy)]))
    plate_col_res <-
        unlist(unname(tdummy[(plate_row + 2), (lmet + 1):ncol(tdummy)]))
    plate_location <-
        paste0(plate, "-", plate_row_res, plate_col_res)
    }
    
    ind.names <- colnames(snpraw)[(lmet + 1):ncol(snpraw)]
    ind.names <- trimws(ind.names, which = "both")  #trim for spaces
    
    if (length(ind.names) != length(unique(ind.names))) {
        cat(
            warn(
                "Warning: Individual names are not unique, adding '_n' to replicates (but not the first instance) to render them unique.\n"
            )
        )
        ind.names <-
            make.unique(as.character(ind.names), sep = "_")
    }
    
    stdmetricscols <- 1:lmet
    # if (length(stdmetricscols) != length(stdmetrics)) { cat(paste('\nCould not find all standard metrics.\n',stdmetrics[!(stdmetrics
    # %in% names(snpraw) )] ,' is missing.\n Carefully check the spelling of your headers!\n')) stop() } if (!is.null(addmetrics)) {
    # addmetricscols <- which( names(snpraw) %in% addmetrics ) if (length(addmetricscols) != length(addmetrics)) { cat(paste('\nCould not
    # find all additional metrics.\n',addmetrics[!(addmetrics %in% names(snpraw) )] ,' is missing.\n Carefully check the spelling of your
    # headers! or set addmetrics to NULL\n')) stop() } stdmetricscols <- c(stdmetricscols, addmetricscols) }
    
    covmetrics <- snpraw[, stdmetricscols]
    
    ##### Various checks (first are there two rows per allele?  we do not need cloneid any more....  covmetrics$CloneID =
    ##### as.character(covmetrics$CloneID) check that there are two lines per locus... covmetrics = separate(covmetrics, CloneID, into =
    ##### c('clid','clrest'),sep = '\\|', extra='merge')
    
    # covmetrics$AlleleID = as.character(covmetrics$AlleleID)
    
    # check that there are two lines per locus... covmetrics = separate(covmetrics, AlleleID, into = c('allid','alrest'),sep = '\\|',
    # extra='merge')
    
    if("MarkerName" %in% colnames(covmetrics)){
      covmetrics$clone <- sub("\\|.*", "", covmetrics$MarkerName, perl = TRUE)
      spp <- sub(".+-+(\\d{1,3}):.+", "\\1", covmetrics$MarkerName)
    }else{
      covmetrics$clone <- sub("\\|.*", "", covmetrics$AlleleID, perl = TRUE)
      spp <- sub(".+-+(\\d{1,3}):.+", "\\1", covmetrics$AlleleID)
    }

    #### find uid within allelid
    covmetrics$uid <- paste(covmetrics$clone, spp, sep = "-")
    ### there should be only twos (and maybe fours)
    tt <- table(table(covmetrics$uid))

    # Testing for SNPs with the same ID
    # Jason Carling e-mail
    # the marker discovery and calling pipeline which we use (called ds14) has a 
    # mechanism to prevent this situation from occurring. When the clusters are 
    # parsed, the software will not allow more than a single marker with the same 
    # CloneID, SNP position, and SNP variant. However, there are a number of 
    # scenarios I can think of which could occur after the initial marker outputs
    # are produced which could violate this rule, resulting in the situation you 
    # have described. One example could be the combination of markers from multiple
    # software outputs into a single new report file. This can occur either 
    # manually, or using SNP re-calling software which calls the markers based on an 
    # input definition file, rather than running the marker discovery de novo. 
    # Secondly, marker renaming is sometimes undertaken, when newly discovered 
    # sequences are added to our CloneID database, and these IDs are retrospectively
    # incorporated into the marker report to replace temporary IDs. I have not yet 
    # determined where the marker name clash was introduced in this case, but I can 
    # answer your question by saying that unfortunately, we cannot completely 
    # prevent this sort of problem for the reasons indicated above due to the
    # potential for processes which modify names or combine data sets after the 
    # initial marker discovery and output.
    
    if(length(tt)>1){
      # if format is two rows
      if(as.numeric(names(tt[2]))==4){
        
        t1 <- as.data.frame(cbind(1:nrow(covmetrics),covmetrics$uid))
        t2 <- table(t1$V2)
        t3 <- names(t2[t2>2])
        t4 <- unlist(lapply(t3,function(x){
          which(t1$V2==x)
        }))
        t5  <- t1[t4,]
        # removing SNPs with the same ID
        
        t5$V1 <- as.numeric(t5$V1)
        covmetrics <- covmetrics[-t5$V1,]
        
        snpraw <- snpraw[-t5$V1,]
          
        cat(warn(
          "  There were",tt[2],"SNPs with the same ID. These SNPs will be removed. SNPs names are:",paste(unique(t5$V2),collapse = " "),"\n"
          
        ))
        
      }else{
        
        t1 <- as.data.frame(cbind(1:nrow(covmetrics),covmetrics$uid))
        t2 <- table(t1$V2)
        t3 <- names(t2[t2>1])
        t4 <- unlist(lapply(t3,function(x){
          which(t1$V2==x)
        }))
        t5  <- t1[t4,]
        t5$V1 <- as.numeric(t5$V1)
        # removing SNPs with the same ID
        covmetrics <- covmetrics[-t5$V1,]
        snpraw <- snpraw[-t5$V1,]
        
        cat(warn(
          "  There were",tt[2],"SNPs with the same ID. These SNPs will be removed. SNPs names are:",paste(unique(t5$V2),collapse = " "),"\n"
          
        ))
        
        
      }
      
      
    }
      
    datas <- snpraw[, (lmet + 1):ncol(snpraw)]
    
    nrows <-NULL
    if (is.null(nrows)) {
        gnrows <-3 - max(datas, na.rm = TRUE)  #if max(datas==1) then two row format, if two then one row format
        
        if (gnrows == 1 | gnrows == 2) {
            nrows <- gnrows
            if (verbose >= 2) {
                cat(report(paste(
                    "  Detected", nrows, "row format.\n"
                )))
            }
        } else {
            stop(
                error(
                    "The DArT format must be either 1row or 2row. This does not seem to be the case here.\n"
                )
            )
        }
        
    }
    
    if (nrows != as.numeric(names(tt[1]))) {
      cat(
        warn(
          "  Warning: The no. rows per Clone does not fit with nrow format. Most likely your data are not read in correctly!\n"
        )
      )
    }
    
    if (verbose >= 2) {
      cat(report(
        paste(
          "Number of rows per clone (should be only ",
          nrows,
          "s):",
          names(tt),
          "\n "
        )
      ))
    }
    
    if (verbose >= 2) {
        cat(report("Added the following locus metrics:\n"))
        cat(report(paste(
            paste(names(snpraw)[stdmetricscols], collapse = " "), ".\n"
        )))
    }
    
    nind <- ncol(datas)
    nsnp <- nrow(covmetrics) / nrows
    
    if (verbose >= 2) {
        cat(important(
            paste(
                "Recognised:",
                nind,
                "individuals and",
                nsnp,
                " SNPs in a",
                nrows,
                "row format using",
                filename,
                "\n"
            )
        ))
    }
    
    out <-
        list(
            nrows = nrows,
            nind = nind,
            nsnp = nsnp,
            covmetrics = covmetrics,
            gendata = datas,
            service = service,
            plate_location = plate_location
        )
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(out)
    
}
