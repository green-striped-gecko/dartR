#' Import DarT data into R and conver it to a genlight object
#' 
#' This function is a wrapper function that allows you to convert you dart file into a genlight object in one step. In previous versions you had to use read.dart and then dart2genlight. In case you have individual metadata for each individual/sample you can specify as before in the dart2genlight command the file that combines the data.
#'
#' @param filename path to file (csv file only currently)
#' @param nas a character specifying NAs (default is "-")
#' @param topskip a number specifying the number of rows to be skipped. If not provided the number of rows to be skipped are "guessed" by the number of rows with "*" at the beginning.
#' @param stdmetrics a vector of column headings that are extracted. AlleleID and its format is compulsory, the rest are needed for filtering.
#' @param addmetrics add additional headers/columns by name
#' @param lastmetric specifies the last non genetic column (Default is "RepAvg"). Be sure to check if that is true, otherwise the number of individuals will not match. You can also specify the last column by a number.
#' @param covfilename the name of the file that has entails additional information on individuals. For the require format check 
#' @param probar show progress bar
#' @return a dart genlight object that contains individuals [if data were provided] and loci meta data [from a DArT report]. The dart genlight object can then be fed into a number of initial screening, export and export functions provided by the package. For some of the function it is necessary to have the metadata that was provided from DArT. Please check the vignette for more information. Additional information can also be found in the help documents for  \code{\link{read.dart}} and \code{\link{dart2genlight}}. 
#' @export
#' @examples{
#' dartfile <- system.file("extdata","testset_SNPs_2Row.csv", package="dart")
#' covfilename <- system.file("extdata","testset_metadata.csv", package="dart")
#' gl <- read.dart2genlight(dartfile, covfilename = covfilename, probar=TRUE)
#' }



read.dart2genlight <- function(filename, covfilename=NULL, nas = "-", topskip=NULL, stdmetrics =c("CloneID", "SNP","SnpPosition","RepAvg","CallRate", "AvgCountRef", "AvgCountSnp", "FreqHomRef", "FreqHomSnp", "FreqHets","OneRatioSnp"), addmetrics=NULL, lastmetric ="RepAvg", probar=TRUE)
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
  covmetrics$clone <- (sub("\\|.*","",covmetrics$CloneID, perl=T))
  spp <- ( sub(".+-+(\\d{1,3}):.+","\\1",covmetrics$CloneID))
  
  
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
  
  dart <- out
  #### out contains the dart data
  nind <- dart[["nind"]]
  nsnp <- dart[["nsnp"]]
  sraw <- dart[["covmetrics"]]
  nrows <- dart[["nrows"]] #check if nrows are provided...
  
  
  if (is.null(nrows))
  {
    cat("nrows not provided. Trying to guess if one row or two row format...\n")
    gnrows = 3 - max(dart$gendata, na.rm = TRUE) 
    
    if (gnrows==1 | gnrows==2)  {nrows <-gnrows;  cat(paste("Should be ", nrows , " row(s) format. Please check if this is the case. Trying to proceed...\n"))} else stop("Cannot be guessed. The dart format must be either one row or two row format and needs to be provided via nrows=1 or 2.\n")
  }
  
  
  
  
  if (sum(c("SNP", "SnpPosition") %in% names(sraw))!=2)
  {
    cat("Could not find SNP or SnpPosition in Dart file. Check you headers!!!")
    stop()
  }
  
  cat("Start conversion....\n")
  cat(paste0("Format is ", nrows," rows.\n"))
  cat("Please note conversion of bigger data sets will take some time!\n")
  cat("Once finished, we recommend to save the object using save(object, file=\"object.rdata\")\n")
  
  if(probar) pb <- txtProgressBar(min=0, max=1, style=3, initial=NA)
  
  
  
  
  
  sdata <- dart[["gendata"]]
  #every second line only....
  esl = seq(nrows,nrow(sdata),nrows)
  
  pos <- sraw$SnpPosition[esl]
  alleles <- as.character(sraw$SNP)[esl]
  a1 <- substr(alleles,nchar(alleles)-2,nchar(alleles))
  a2 <-  sub(">","/", a1)
  locname <- paste(sraw$uid[esl], a2,sep="-")
  geninddata <- matrix(NA, nrow=nsnp, ncol=nind)
  
  
  if (nrows==2)
  {
    for (i in 1:nind)
    {
      isnp = paste(sdata[esl-1,i],sdata[esl,i], sep="/")
      g <- isnp
      g <- gsub("0/1",2,g)
      g <- gsub("1/0",0,g)
      g <- gsub("1/1",1,g)
      g <- gsub("NA/NA",NA,g)
      geninddata[,i] <- as.numeric(g)
      if (probar) setTxtProgressBar(pb, i/nind)
    }
  } else 
    
  {
    for (i in 1:nind)
    {
      isnp = sdata[esl,i]
      g <- isnp
      g <- 3-g
      g <- ifelse(g==3, 0, g)
      geninddata[,i] <- g
      if (probar) setTxtProgressBar(pb, i/nind)
    }  
  }
  gout <- new("genlight", gen=t(geninddata), ploidy=2, ind.names=colnames(sdata), loc.names=locname ,loc.all=a2, position=pos, parallel=F)
  
  if (probar) close(pb)
  
  
  
  #refactor data.frame
  df <- as.data.frame(lapply( sraw[esl,], function (x) if (is.factor(x)) factor(x) else x)) 
  
  
  gout@other$loc.metrics <- df
  
  ####
  #additional covariates and long lat to the data file are stored in other
  
  if (!is.null(covfilename))
  {
    cat(paste("Try to add covariate file:", covfilename,".\n"))
    ###### population and individual file to link AAnumbers to populations...
    ind.cov <- read.csv(covfilename, header=T, stringsAsFactors=T)
    # is there an entry for every individual
    
    id.col = match( "id", names(ind.cov))
    
    
    
    if (is.na(id.col)) {cat ("There is no id column\n") ;stop()} else {
      
      if (length(ind.cov[,id.col])!= length(unique(ind.cov[,id.col]))) {cat("Individual names are not unique. You need to change them!\n"); stop()}  
      
      
      #reorder
      if (length(ind.cov[,id.col]) !=length(names(sdata)))  {cat ("Ids of covariate file does not match the number of ids in the genetic file. Maybe this is fine if a subset matches.\n") } 
      
      ord <- match(names(sdata), ind.cov[,id.col])
      ord <- ord[!is.na(ord)]
      
      
      if (length(ord)>1 & length(ord)<=nind ) 
      {cat (paste("Ids of covariate file (at least a subset of) are matching!\nFound ", length( ord ==nind),"matching ids out of" , nrow(ind.cov), "ids provided in the covariate file. Subsetting snps now!.\n "))
        ord2 <- match(ind.cov[ord,id.col], indNames(gout))
        gout <- gout[ord2,]
      }else {cat("Ids are not matching!!!!\n");stop()}
    }
    
    
    pop.col = match( "pop", names(ind.cov))
    
    if (is.na(pop.col)) {cat ("Please note:there is no pop column\n") }  else {
      
      pop(gout) <- as.factor(ind.cov[ord,pop.col])
      cat("Added pop factor.\n")
    }
    
    lat.col = match( "lat", names(ind.cov))
    lon.col = match( "lon", names(ind.cov))
    
    if (is.na(lat.col)) {cat ("Please note:there is no lat column\n") }
    if (is.na(lon.col)) {cat ("Please note:there is no lon column\n") }
    if (!is.na(lat.col) & !is.na(lon.col))
    {
      gout@other$latlong <- ind.cov[ord,c(lat.col, lon.col)]
      rownames(gout@other$latlong)  <-  ind.cov[ord,id.col]
      cat("Added latlon data.\n" )
    }
    
    # known.col <- names( ind.cov) %in% c("id","pop", "lat", "lon")
    # known.col <- ifelse(is.na(known.col), , known.col)
    # other.col <- names(ind.cov)[!known.col]
    other.col <- names(ind.cov)
    if (length(other.col>0) )
    {
      gout@other$ind.metrics<-ind.cov[ord,other.col, drop=FALSE]
      rownames(gout@other$ind.metrics) <- ind.cov[ord,id.col]
      cat(paste("Added ",other.col," to the other$ind.metrics slot.\n"))
    }
  }
  gout 
  
}
