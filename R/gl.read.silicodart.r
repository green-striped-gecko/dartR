#' Import presence/absence data from SilicoDArT to genlight \{agegenet\} format (ploidy=1)
#'
#' DaRT provide the data as a matrix of entities (individual animals) across the top and
#' attributes (P/A of sequenced fragment) down the side in a format that is unique to DArT. This program
#' reads the data in to adegenet format for consistency with
#' other programming activity. The script may require modification as DArT modify their
#' data formats from time to time.
#'
#' gl.read.silicodart() opens the data file (csv comma delimited) and skips the first n=topskip lines. The script assumes
#' that the next line contains the entity labels (specimen ids) followed immediately by the SNP data for the first locus.
#' It reads the presence/absence data into a matrix of 1s and 0s,
#' and inputs the locus metadata and specimen metadata. The locus metadata comprises a series of columns of values for
#' each locus including the essential columns of CloneID and the desirable variables Reproducibility and PIC.
#' Refer to documentation provide by DArT for an explanation of these columns.
#'
#' The specimen metadata provides the opportunity
#' to reassign specimens to populations, and to add other data relevant to the specimen. The key variables are id (specimen identity
#' which must be the same and in the same order as the SilicoDArT file, each unique), pop (population assignment), lat (latitude, optional)
#' and lon (longitude, optional). id, pop, lat, lon are the column headers in the csv file. Other optional columns can be added.
#'
#' The data matrix, locus names (forced to be unique), locus metadata, specimen names,
#' specimen metadata are combined into a genInd object. Refer to the documentation for \{adegenet\} for further details.
#'
#' @param filename -- name of csv file containing the SilicoDArT data [required]
#' @param ind.metafile -- name of csv file containing metadata assigned to each entity (individual) [default NULL]
#' @param nas -- missing data character [default "-"]
#' @param topskip -- number of rows to skip before the header row (containing the specimen identities) [optional]
#' @param lastmetric -- specifies the last non genetic column (Default is "Reproducibility"). Be sure to check if that is true, otherwise the number of individuals will not match. You can also specify the last column by a number. [default Reproducibility]
#' @param probar show progress bar
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2, or as set by gl.set.verbose()]
#' @return An object of class \code{genlight} with ploidy set to 1, containing the presence/absence data, and locus and individual metadata
#' @export
#' @author Bernd Gruber and Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' silicodartfile <- system.file("extdata","testset_SilicoDArT.csv", package="dartR")
#' metadata <- system.file("extdata",ind.metafile ="testset_metadata_silicodart.csv", package="dartR")
#' testset.gs <- gl.read.silicodart(filename = silicodartfile, ind.metafile = metadata)
 



gl.read.silicodart <- function(filename, 
                               ind.metafile=NULL,  
                               nas="-",  
                               topskip=NULL, 
                               lastmetric="Reproducibility",
                               probar=TRUE,
                               verbose=NULL) {

# TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  
# SET VERBOSITY
  
  if (is.null(verbose)){ 
    verbose <- 2
  } 
  
  if (verbose < 0 | verbose > 5){
    cat(paste("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n"))
    verbose <- 2
  }
  
  if(verbose == 0 & probar){
    probar=FALSE
    #cat("  Note: Progress bar set to FALSE for verbosity==0\n")
  }
  
# FLAG SCRIPT START
  
  if (verbose >= 1){
    if(verbose==5){
      cat("Starting",funname,"[ Build =",build,"]\n")
    } else {
      cat("Starting",funname,"\n")
    }
  }
  
# DO THE JOB  
  
  if (verbose >= 2){
    cat("Reading data from file:", filename,"\n")
    cat("  This may take some time, please wait!\n")
  }
  
  if (is.null(topskip)) {
    if (verbose >= 2){cat("Topskip not provided. Guessing topskip...\n")}   
    tdummy <- read.csv(filename,   na.strings=nas,  check.names=FALSE, nrows = 20, header=FALSE)
    
    nskip <- sum(tdummy[,1] == "*"  )
    if (nskip > 0) { 
      topskip <- nskip; cat(paste("Set topskip to ", nskip,". Proceeding ...\n"))
    } else {
      stop("Could not determine topskip (the number of rows that need to be skipped. Please provide it manually.\n") 
    }
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
  if (length(ind.names) != length(unique(ind.names))) {
    cat("The following labels for individuals are not unique:\n")
    cat(ind.names[duplicated(ind.names)])
    cat("\n")
    cat("Warning: Individual names are not unique. Rendering locus names unique with sequential suffix _1, _2 for duplicates.\n")
    ind.names <- make.unique(ind.names)
  }  
  
  datas <- snpraw[, (lmet+1):ncol(snpraw)]
  nrows=1  #there is no two row SilicoFormat??
  stdmetricscols <- 1:lmet
  
  if (verbose >= 2){ 
    cat ("Added the following locus metrics:\n")
    cat (paste(paste(names(snpraw)[stdmetricscols], collapse=" "),".\n"))
  }  
  covmetrics <-  snpraw[,stdmetricscols]
  
  nind <- ncol(datas)
  #nsnp <- nrow(covmetrics)/nrows
  
  if (verbose >= 2){cat(paste("Recognised:", nind, "individuals and",nsnp,"Sequence Tags using", filename,"\n"))}
 
 
  if(max(datas,na.rm=TRUE)!=1 || min(datas,na.rm=TRUE)!=0) {
    stop("Fatal Error: Tag P/A data must be 0 or 1!\n")
  }

  if (verbose >= 2){
    cat("\nStarting conversion to a genlight object ....\n")
    cat("  Please note conversion of bigger data sets will take some time!\n")
    cat("  Once finished, we recommend you save the object using saveRDS(object, file=\"object.rdata\")\n")
  }


  #create unique locnames based on cloneID
    index <- unique(covmetrics$CloneID[which(duplicated(covmetrics$CloneID))])
    if (length(index>0)){
      cat("Warning: Locus names [CloneIDs] are not unique!\n")
      cat("         Rendering locus names unique with sequential suffix _1, _2 for duplicates.\n")
      for (i in 1:length(index)) {
        loc <- index[i]
        i2 <- which(covmetrics$CloneID %in% loc)  
        covmetrics$CloneID[i2] <- paste0(covmetrics$CloneID[i2],"_",1:length(i2))
      }
    }

    glout <- new("genlight", gen=t(datas),ploidy=1, ind.names=ind.names, loc.names=covmetrics$CloneID )

  #add loc.metrics

    glout@other$loc.metrics <- covmetrics

    if (!is.null(ind.metafile)){
      cat(paste("  Adding individual metadata:", ind.metafile,".\n"))
      ind.cov <- read.csv(ind.metafile, header=T, stringsAsFactors=T)
      # is there an entry for every individual
      id.col = match( "id", names(ind.cov))
  
      if (is.na(id.col)) {
        stop("Fatal Error: There is no id column\n")
      } else {
        ind.cov[,id.col]<- trimws(ind.cov[,id.col], which="both")  #trim spaces
        if (length(ind.cov[,id.col])!= length(unique(ind.cov[,id.col]))) {
          stop("Fatal Error: Individual names are not unique. You need to change them!\n")
        }  
      #reorder
        if (length(ind.cov[,id.col]) !=length(names(datas))){
          cat ("Warning: Ids for individual metadata does not match in number the ids in the SNP data file. Maybe this is fine if a subset matches.\n") 
        } 
        ord <- match(names(datas), ind.cov[,id.col])
        ord <- ord[!is.na(ord)]
    
        if (length(ord)>1 & length(ord)<=nind ){
          cat (paste("Ids for individual metadata (at least a subset of) are matching!\nFound ", length( ord ==nind),"matching ids out of" , nrow(ind.cov), "ids provided in the ind.metadata file. Subsetting loci now!.\n "))
          ord2 <- match(ind.cov[ord,id.col], indNames(glout))
          glout <- glout[ord2,]
        } else {
          stop("Fatal Error: Ids are not matching!!!!\n")
        }
      }
  
      pop.col = match( "pop", names(ind.cov))
  
      if (is.na(pop.col)) {
        cat ("Please note: there is no pop column\n") 
        pop(out) <- array(NA,nInd(glout))
        cat("Created pop column with NAs\n")
      } else {
        pop(glout) <- as.factor(ind.cov[ord,pop.col])
        cat("Added pop factor.\n")
      }
  
      lat.col = match( "lat", names(ind.cov))
      lon.col = match( "lon", names(ind.cov))
  
      if (is.na(lat.col)) {cat ("Please note: there is no lat column\n") }
      if (is.na(lon.col)) {cat ("Please note: there is no lon column\n") }
      if (!is.na(lat.col) & !is.na(lon.col)){
        glout@other$latlong <- ind.cov[ord,c(lat.col, lon.col)]
        rownames(glout@other$latlong)  <-  ind.cov[ord,id.col]
        cat(" Added latlon data.\n" )
      }
  
      # known.col <- names( ind.cov) %in% c("id","pop", "lat", "lon")
      # known.col <- ifelse(is.na(known.col), , known.col)
      # other.col <- names(ind.cov)[!known.col]
      other.col <- names(ind.cov)
      if (length(other.col>0) ){
        glout@other$ind.metrics<-ind.cov[ord,other.col, drop=FALSE]
        rownames(glout@other$ind.metrics) <- ind.cov[ord,id.col]
        cat(paste(" Added ",other.col," to the other$ind.metrics slot.\n"))
      }
    }

    if (is.null(glout@other$history)) {
      glout@other$history <- list(match.call())
    }
    #add recalc flags (TRUE=up-to-date, FALSE=no longer valid)
    #all potential headers that can be relculated
    recalc.flags <-  c( "AvgPIC", "OneRatioRef","OneRatioSnp", "PICRef", "PICSnp", "CallRate",  "maf", "FreqHets" ,"FreqHomRef" , "FreqHomSnp", 
                    "monomorphs", "OneRatio", "PIC")
    glout@other$loc.metrics.flags <-  data.frame(matrix(TRUE, nrow=1, ncol=length(recalc.flags)))
    names(glout@other$loc.metrics.flags) <- recalc.flags
    glout@other$verbose <- 2

# FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(paste("Completed:",funname,"\n"))
    cat("Genlight object with ploidy=1 created to hold Tag P/A data\n")
  }
  
  return(glout)

}
