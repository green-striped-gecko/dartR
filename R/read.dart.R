#' Import DarT data to R
#'
#' @param filename path to file (csv file only currently)
#' @param nas a character specifying NAs (default is "-")
#' @param topskip a number specifying the number of rows to be skipped. If not provided the number of rows to be skipped are "guessed" by the number of rows with "*" at the beginning.
#' @param stdmetrics a vector of column headings that are extracted. AlleleID and its format is compulsory, the rest are needed for filtering.
#' @param addmetrics add additional headers/columns by name
#' @param lastmetric specifies the last non genetic column (Default is "RepAvg"). Be sure to check if that is true, otherwise the number of individuals will not match. You can also specify the last column by a number.
#' @return a list of length 5. #dart format (one or two rows) #individuals, #snps, #non genetic metrics, #genetic data (still two line format, rows=snps, columns=individuals)
#' @export
#' @examples{
#' dartfile <- system.file("extdata","testset_SNPs_2Row.csv", package="dartR")
#' dart <-read.dart(dartfile)
#' }



read.dart <- function(filename, nas = "-", topskip=NULL, stdmetrics =c("AlleleID", "SNP","SnpPosition","RepAvg","CallRate", "AvgCountRef", "AvgCountSnp", "FreqHomRef", "FreqHomSnp", "FreqHets","OneRatioSnp"), addmetrics=NULL, lastmetric ="RepAvg")
{
  
  if (is.null(topskip))
  {
  cat("Topskip not provided. Try to guess topskip...\n")    
  tdummy <- read.csv(filename,   na.strings=nas,  check.names=FALSE, nrows = 20, header=FALSE)
  
  nskip <- sum(tdummy[,1]=="*"  )
  if (nskip>0) { topskip <- nskip; cat(paste("Set topskip to ", nskip,". Trying to proceed...\n"))} else {
    stop("Could not determine topskip (the number of rows that need to be skipped. Please provide it manually.\n") 
    
  }
  }

  
  

  
  snpraw <- read.csv(filename,   na.strings=nas, skip = topskip, check.names=FALSE)


  if (is.character(lastmetric))
  {
    lmet <- which(lastmetric==names(snpraw))
    if (length(lmet)==0)  stop (paste("Could not determine data columns based on", lastmetric,"!\n"))
  } else lmet  <- lastmetric
  
  ind.names <- colnames(snpraw)[(lmet+1):ncol(snpraw) ]
  ind.names <- trimws(ind.names, which = "both") #trim for spaces
  if (length(ind.names)!= length(unique(ind.names))) stop("Individual names are not unique. You need to change them!\n")
  
  datas <- snpraw[, (lmet+1):ncol(snpraw)]
  
  nrows = NULL
  if (is.null(nrows)) 
  {
    cat("Trying to determine if one row or two row format...\n")
    gnrows = 3-max(datas, na.rm = TRUE)  #if max(datas==1) then two row format, if two then one row format
    
    if (gnrows==1 | gnrows==2)  {nrows <-gnrows;  cat(paste("Found ", nrows , " row(s) format. Proceed...\n"))} else stop("The dart format either one row or two row format. This does not seem to be the case here.\n")
  } 
  
  
  
  stdmetricscols <- which(  names(snpraw)   %in% stdmetrics )
  
  if (length(stdmetricscols) != length(stdmetrics))
  { cat(paste("\nCould not find all standard metrics.\n",stdmetrics[!(stdmetrics %in% names(snpraw)   )]
              ," is missing.\n Carefully check the spelling of your headers!\n"))
    stop()
  }
  
  if (!is.null(addmetrics)) 
  {
    addmetricscols <- which(  names(snpraw)   %in% addmetrics )
    if (length(addmetricscols) != length(addmetrics))
    { cat(paste("\nCould not find all additional metrics.\n",addmetrics[!(addmetrics %in% names(snpraw)   )]
                ," is missing.\n Carefully check the spelling of your headers! or set addmetrics to NULL\n"))
      stop()
    }
    stdmetricscols <- c(stdmetricscols, addmetricscols)
  } 
  cat ("Added the following covmetrics:\n")
  cat (paste(paste(names(snpraw)[stdmetricscols], collapse=" "),".\n"))
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
  cat(paste("Number of rows per Clone. Should be only ", nrows,"s:", names(tt),"\n "))
  if (nrows!=as.numeric(names(tt))) cat("!!!!!Number of rows per clone does not fit with nrow format. Most likely your data are not read in correctly.!!!!!\n") 
  nind <- ncol(datas)
  nsnp <- nrow(covmetrics)/nrows
  
  cat(paste("Recognised:", nind, "individuals and",nsnp," SNPs in a",nrows,"row format using", filename,"\n"))
  
  out <- list(nrows=nrows, nind=nind, nsnp=nsnp, covmetrics= covmetrics, gendata =datas)
  
  out
  
  
}
