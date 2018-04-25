#' Import SNP data from DArT and convert to  genlight \{agegenet\} format (gl)
#'
#' DaRT provide the data as a matrix of entities (individual turtles) across the top and
#' attributes (SNP loci) down the side in a format that is unique to DArT. This program
#' reads the data in to adegenet format (genlight) for consistency with
#' other programming activity. The script or the data may require modification as DArT modify their
#' data formats from time to time.
#'
#' gl.read.dart() opens the data file (csv comma delimited) and skips the first n=topskip lines. The script assumes
#' that the next line contains the entity labels (specimen ids) followed immediately by the SNP data for the first locus.
#' It reads the SNP data into a matrix of 1s and 0s,
#' and inputs the locus metadata and specimen metadata. The locus metadata comprises a series of columns of values for
#' each locus including the essential columns of AlleleID, SNP, SnpPostion and the desirable variables REpAvg and AvgPIC.
#' Refer to documentation provide by DArT for an explanation of these columns.
#'
#' The specimen metadata provides the opportunity
#' to reassign specimens to populations, and to add other data relevant to the specimen. The key variables are id (specimen identity
#' which must be the same and in the same order as the DArTSeq file, each unique), pop (population assignment), lat (latitude, optional)
#' and lon (longitude, optional). id, pop, lat, lon are the column headers in the csv file. Other optional columns can be added.
#'
#' The SNP matrix, locus names (constructed from the AlleleID, SNP and SnpPosition to be unique), locus metadata, specimen names,
#' specimen metadata are combined into a genlight object. Refer to the genlight documentation (Package adegenet) for further details.
#'
#' @param datafile -- name of csv file containing the DartSeq data in 2-row format (csv) [required]
#' @param topskip -- number of rows to skip before the header row (containing the specimen identities [required]
#' @param nmetavar -- number of columns containing the locus metadata (e.g. AlleleID, RepAvg) [required]
#' @param nas -- missing data character [default "-"]
#' @param ind.metafile -- name of csv file containing metadata assigned to each entity (individual) [default NULL]
#' @param pbar -- display progress bar [FALSE]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return An object of class ("genlight") containing the SNP data, and locus and individual metadata
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @export
#' @examples
#' \dontrun{
#' gl <- gl.read.dart.2row(datafile="SNP_DFwt15-1908_scores_2Row.csv", topskip=6, 
#' nmetavar=16, nas="-", ind.metafile="metadata.csv" )
#' }

# Last edit:25-Apr-18

gl.read.dart.2row <- function(datafile,topskip,nmetavar, nas="-",ind.metafile=NULL,pbar=TRUE,v=2) {

# INPUT THE DATA TO PRELIMINARY STORAGE
  
  if ( v > 0) {cat("Starting gl.read.dart.2row: Reading DArT csv file\n")}
  if (v > 1){
    cat("  Reading data from file:", datafile,"\n")
    cat("    This may take some time, please wait!\n")
  }  
  x <- read.csv(datafile, na.strings=nas, skip = topskip, check.names=FALSE)
  if (v > 1){
    cat("  The following locus metadata was identified: ", names(x[1:nmetavar]),"\n")
  }  
# Error checks
  if(any(names(x) == "AlleleID")) {
    if (v > 1){cat("    includes key variable AllelID\n")}
  } else {
      cat("Fatal Error: Dataset does not include key variable AlleleID!\n"); stop()
  }
  if(any(names(x) == "SNP")) {
    if (v > 1){cat("  includes key variable SNP\n")}
  } else {
    cat("Fatal Error: Dataset does not include key variable SNP!\n"); stop()
  }
  if(any(names(x) == "SnpPosition")) {
    if (v > 1){cat("  includes key variable SnpPosition\n")}
  } else {
    cat("Fatal Error: Dataset does not include key variable SnpPosition!\n"); stop()
  }
    if(any(names(x) == "AvgPIC")) {
      if (v > 1){cat("  includes key variable AvgPIC\n")}
  } else {
    cat("  Warning: Dataset does not include variable AvgPIC which may limit your options in later analyses!\n")
  }
  if(any(names(x) == "TrimmedSequence")) {
    if (v > 1){cat("  includes key variable TrimmedSequence\n")}
  } else {
    cat("  Warning: Dataset does not include variable TrimmedSequence which may limit your options in later analyses!\n")
  }
  if(any(names(x) == "RepAvg")) {
    if (v > 1){cat("  includes key variable RepAvg\n")}
  } else {
    cat("  Warning: Dataset does not include variable RepAvg which may limit your filtering options in later analyses!\n")
  }

# Extract names of the entities (individuals)
  ind.names <- colnames(x)[(nmetavar+1):ncol(x)]
  if (v > 1){cat("  Data identified for ",ncol(x)-nmetavar, "individuals, ", nrow(x)/2, "loci")}

  # More error checks
  if (length(ind.names)!= length(unique(ind.names))) {
    cat("  Warning: Specimen names are not unique!\n")
    cat("         Duplicated names:\n")
    noccur <- table(ind.names)
    cat(paste("              ",names(noccur[noccur>1])),"\n")
    cat("         Rendering specimen names unique with sequential suffix _1, _2 etc\n")
    ind.names <- make.unique(ind.names, sep="_")
  }
# Extract the SNP data
  snpdata <- x[, (nmetavar+1):ncol(x)]
  
# Extract the standard metadata for loci
 locus.metadata <- x[, 1:nmetavar]
 
# More error checks
  if(max(snpdata,na.rm=TRUE)!=1 || min(snpdata,na.rm=TRUE)!=0) {
    cat("Fatal Error: SNP data must be 0 or 1!\n"); stop()
  }
 
# Calculate number of entities (individuals) and attributes (loci)
  nind <- ncol(snpdata)
  nloci <- nrow(locus.metadata)/2

# CONVERT TO GENLIGHT FORMAT

  if (v > 1){
    cat("\n  Starting conversion to genlight object ....\n")
    cat("    Please note conversion of bigger data sets will take some time!\n")
  }  
  if(pbar){
    pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
    getTxtProgressBar(pb)
  }
  
# 2 Row format -- Consider every second line only
  seqby2 = seq(2,nrow(snpdata),2)
  
# Extract the SNP position
  pos <- locus.metadata$SnpPosition[seqby2]
  
# Extract the SNP transitions (metavariable SNP)
  state <- as.character(locus.metadata$SNP)[seqby2]
  state2 <- substr(state,nchar(state)-2,nchar(state))
  state <- sub(":","-", state)
  state <- sub(">","/", state)
  state2 <- sub(">","/", state2)
  uid <- as.character(locus.metadata$AlleleID)[seqby2]
  a.list <- strsplit(uid, "F")
  uid <- sapply(a.list, "[", 1)
  
# Create locus names
  locname <- paste(uid, state, sep="")
  if (length(locname)!= length(unique(locname))) {
    cat("Warning: Locus names are not unique!\n")
    cat("         Duplicated names:\n")
    noccur <- table(locname)
    cat(paste("              ",names(noccur[noccur>1])),"\n")
    cat("         Rendering locus names unique with sequential suffix _1, _2 etc\n")
    locname <- make.unique(locname, sep="_")
  }
  
# Initialize the data matrix
  x <- matrix(NA, nrow=nloci, ncol=nind)
  
# For each individual, convert the 2 line data to required one line format (0, homozygous reference; 1, heterozygous; 2 homozygous mutant)
  for (i in 1:nind) {
    isnp = paste(snpdata[seqby2-1,i],snpdata[seqby2,i], sep="/")
    g <- isnp
    g <- gsub("0/1",2,g)
    g <- gsub("1/0",0,g)
    g <- gsub("1/1",1,g)
    g <- gsub("NA/NA",NA,g)
    x[,i] <- as.numeric(g)
    if (pbar){setTxtProgressBar(pb, i/nind)}
  }
# Create the genlight object
  gl <- new("genlight", gen=t(x), ploidy=2, ind.names=colnames(snpdata), loc.names=locname ,loc.all=state2, position=pos, parallel=F)

  close(pb)

# Add in the standard metadata
  gl@other$loc.metrics <- locus.metadata[seqby2,]

# Add in extra metadata -- population assignments
if (!is.null(ind.metafile)) {
  if (v > 1){
    cat("Adding population assignments and additional individual metadata from file :", ind.metafile,"\n")
  }
  ind.metadata <- read.csv(ind.metafile, header=T, stringsAsFactors=T)
  
# Remove leading and trailing spaces that could lead to a spurious mismatch
  ind.metadata$id <- gsub("^\\s+|\\s+$", "", ind.metadata$id)
  
# Check for an entry for every individual
  id.col = match( "id", names(ind.metadata))
  if (is.na(id.col)) {
    cat ("Fatal Error: No id column present!\n") ;stop()
    } else {
    if (sum(ind.metadata[,id.col] == names(snpdata)) == nind ) {
      if (v > 1) {cat ("  Ids of individual metadata file match!\n")}
    }else {
      cat("Fatal Error: Ids in files ",datafile,"and ",ind.metafile," do not match\n     or not in the same order!\n\n");stop()
    }
  }
  pop.col = match( "pop", names(ind.metadata))
  
# Check for population assignment
  if (is.na(pop.col)) {
    if (v > 1) {cat ("  Warning: No pop column present\n")}
  } else {
    pop(gl) <- as.factor(ind.metadata[,pop.col])
    if (v > 1){cat("  Populations assigned to individuals\n")}
  }
  
# Check for latitude and longitude data
  lat.col = match( "lat", names(ind.metadata))
  lon.col = match( "lon", names(ind.metadata))
  if (is.na(lat.col)) {
    if (v > 1) {cat ("Warning: No lat column present\n")}
  }
  if (is.na(lon.col)) {
    if (v > 1) {cat ("Warning: No lon column present\n")}
  }
  if (!is.na(lat.col) & !is.na(lon.col))  {
    gl@other$latlong <- ind.metadata[,c(lat.col, lon.col)]
    if (v > 1) {cat("  Added latlon data\n" )}
  }
# Check for other metadata
  other.col <- names(ind.metadata)
  if (length(other.col>0) ) {
    gl@other$ind.metrics<-ind.metadata[,other.col]
    if (v > 1) {cat("Added ",other.col," to the other$ind.metrics slot\n")}
  }
}
# Report
  if (v > 0) {cat("gl.read.dart.2row completed: Genlight object created\n")}

  return <- gl

}
