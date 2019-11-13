#' Import DarT data to R
#'
#' Internal function called by gl.read.dart
#' @param filename path to file (csv file only currently)
#' @param nas a character specifying NAs (default is "-")
#' @param topskip a number specifying the number of rows to be skipped. If not provided the number of rows to be skipped are "guessed" by the number of rows with "*" at the beginning.
#' @param lastmetric specifies the last non genetic column (Default is "RepAvg"). Be sure to check if that is true, otherwise the number of individuals will not match. You can also specify the last column by a number.
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return a list of length 5. #dart format (one or two rows) #individuals, #snps, #non genetic metrics, #genetic data (still two line format, rows=snps, columns=individuals)

utils.read.dart <- function(filename, nas = "-", topskip=NULL,  lastmetric ="RepAvg",verbose=2){
  
# TIDY UP FILE SPECS
  
  build <- "Jacob"
  funname <- match.call()[[1]]

# FLAG SCRIPT START

  
  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }
  
  if (verbose >= 1){
    cat(paste("Starting",funname,"\n"))
  }
  
# DO THE JOB
  

  if (is.null(topskip)) {
    if (verbose >= 2){
      cat("  Topskip not provided. ") 
    }
    tdummy <- read.csv(filename,   na.strings=nas,  check.names=FALSE, nrows = 20, header=FALSE)
  
    nskip <- sum(tdummy[,1] == "*"  )
    if (nskip > 0) { 
      topskip <- nskip
      if (verbose >= 2){
        cat(paste("Setting topskip to",nskip,".\n"))
      }  
    } else {
      stop("Could not determine the number of rows that need to be skipped. Please provide it manually by setting the topskip parameter.\n") 
    }
  }

  if (verbose >= 2){
    cat("  Reading in the SNP data\n")
  }
  snpraw <- read.csv(filename, na.strings=nas, skip = topskip, check.names=FALSE)

  if (is.character(lastmetric)) {
    lmet <- which(lastmetric==names(snpraw))
    if (length(lmet)==0)  {
      stop (paste("Could not determine number of data columns based on", lastmetric,"!\n"))
    }  
  } else {
    lmet  <- lastmetric
  }  
  
  ind.names <- colnames(snpraw)[(lmet+1):ncol(snpraw) ]
  ind.names <- trimws(ind.names, which = "both") #trim for spaces
  if (length(ind.names)!= length(unique(ind.names))) {
    cat("Warning: Individual names are not unique, adding '_n' to replicates (but not the first instance) to render them unique.\n")
    ind.names <- make.unique(as.character(ind.names), sep = "_")
  }  
  
  datas <- snpraw[, (lmet+1):ncol(snpraw)]
  
  nrows = NULL
  if (is.null(nrows)) {
    gnrows = 3-max(datas, na.rm = TRUE)  #if max(datas==1) then two row format, if two then one row format
    
    if (gnrows==1 | gnrows==2)  {
      nrows <-gnrows
      if (verbose >= 2){
        cat(paste("  Detected",nrows,"row format.\n"))
      }  
    } else {
      stop("The DArT format must be either 1row or 2row. This does not seem to be the case here.\n")
    }
    
  } 
  
  stdmetricscols <- 1:lmet
  # 
  # if (length(stdmetricscols) != length(stdmetrics))
  # { cat(paste("\nCould not find all standard metrics.\n",stdmetrics[!(stdmetrics %in% names(snpraw)   )]
  #             ," is missing.\n Carefully check the spelling of your headers!\n"))
  #   stop()
  # }
  # 
  # if (!is.null(addmetrics)) 
  # {
  #   addmetricscols <- which(  names(snpraw)   %in% addmetrics )
  #   if (length(addmetricscols) != length(addmetrics))
  #   { cat(paste("\nCould not find all additional metrics.\n",addmetrics[!(addmetrics %in% names(snpraw)   )]
  #               ," is missing.\n Carefully check the spelling of your headers! or set addmetrics to NULL\n"))
  #     stop()
  #   }
  #   stdmetricscols <- c(stdmetricscols, addmetricscols)
  # } 
  
  if (verbose >= 2){
    cat ("Added the following locus metrics:\n")
    cat (paste(paste(names(snpraw)[stdmetricscols], collapse=" "),".\n"))
  }
  covmetrics <-  snpraw[,stdmetricscols]
  
  #####Various checks (first are there two rows per allele?
  # we do not need cloneid any more....  
  #covmetrics$CloneID = as.character(covmetrics$CloneID)
  #check that there are two lines per locus...
  #covmetrics = separate(covmetrics, CloneID, into  = c("clid","clrest"),sep = "\\|", extra="merge")
  
  #covmetrics$AlleleID = as.character(covmetrics$AlleleID)
  
  #check that there are two lines per locus...
  #covmetrics = separate(covmetrics, AlleleID, into  = c("allid","alrest"),sep = "\\|", extra="merge")
  covmetrics$clone <- (sub("\\|.*","",covmetrics$AlleleID, perl=T))
  spp <- ( sub(".+-+(\\d{1,3}):.+","\\1",covmetrics$AlleleID))
  
  
  #### find uid within allelid 
  covmetrics$uid <- paste(covmetrics$clone, spp,sep="-")
  ### there should be only twos (and maybe fours)
  tt <- table(table(covmetrics$uid) )
  if (verbose >= 2){
    cat(paste("Number of rows per clone (should be only ", nrows,"s):", names(tt),"\n "))
  }
  if (nrows!=as.numeric(names(tt))) {
    cat("  Warning: The no. rows per Clone does not fit with nrow format. Most likely your data are not read in correctly!\n") 
  }  
  nind <- ncol(datas)
  nsnp <- nrow(covmetrics)/nrows
  
  if (verbose >= 2){
    cat(paste("Recognised:", nind, "individuals and",nsnp," SNPs in a",nrows,"row format using", filename,"\n"))
  }
  
  out <- list(nrows=nrows, nind=nind, nsnp=nsnp, covmetrics= covmetrics, gendata =datas)
  
# FLAG SCRIPT END
  
  if (verbose > 0) {
    cat(paste("Completed:",funname,"\n"))
  }
  
  return(out)
  
}
