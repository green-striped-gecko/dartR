#' Import presence/absence data from SilicoDArT and convert to  genind \{agegenet\} format
#'
#' DaRT provide the data as a matrix of entities (individual animals) across the top and
#' attributes (P/A of sequenced fragment) down the side in a format that is unique to DArT. This program
#' reads the data in to adegenet format (genind) for consistency with
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
#' @param datafile -- name of csv file containing the SilicoDArT data [required]
#' @param topskip -- number of rows to skip before the header row (containing the specimen identities) [required]
#' @param nmetavar -- number of columns containing the locus metadata (e.g. CloneID, Reproducibility) [required]
#' @param nas -- missing data character [default "-"]
#' @param ind.metafile -- name of csv file containing metadata assigned to each entity (individual) [default NULL]
#' @return An object of class ("genInd") containing the presence/absence data, and locus and individual metadata
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' \donttest{
#' glind <- gl.read.silicodart(datafile="SNP_DFwt15-1908_scores_2Row.csv", topskip=6, 
#' nmetavar=16, nas="-", ind.metafile="metadata.csv" )
#' }

# Debug
#   datafile <- "silicodart.header.fixed.csv"
#   topskip <- 4
#   nmetavar <- 19
#   ind.metafile <- "metadata.csv"
#   nas="4"
#

gl.read.silicodart <- function(datafile, topskip, nmetavar, nas="-", ind.metafile=NULL) {

# INPUT THE DATA TO PRELIMINARY STORAGE

  cat("Reading data from file:", datafile,"\n")
  cat("  This may take some time, please wait!\n")
  x <- read.csv(datafile, na.strings=nas, skip = topskip, check.names=FALSE)
  cat("The following metadata for loci was identified: ", names(x[1:nmetavar]),"\n")
# Error checks
  if(any(names(x) == "CloneID")) {
    cat("  includes key variable CloneID\n")
  } else {
      cat("Fatal Error: Dataset does not include key variable CloneID!\n"); stop()
  }
  if(any(names(x) == "Reproducibility")) {
    cat("  includes key variable Reproducibility\n")
  } else {
    cat("Warning: Dataset does not include variable Reproducibility which may limit your filtering options in later analyses!\n")
  }

# Extract names of the entities (individuals)
  ind.names <- colnames(x)[(nmetavar+1):ncol(x)]
  cat("Data identified for ",ncol(x)-nmetavar, "individuals, ", nrow(x), "loci")
# More error checks
  if (length(ind.names)!= length(unique(ind.names))) {
    cat("Warning: Specimen names are not unique!\n")
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
  nloci <- nrow(locus.metadata)

# CONVERT TO GENIND FORMAT

cat("\nStarting conversion to genInd object ....\n")
cat("Please note conversion of bigger data sets will take some time!\n")

# Create locus names

  locname <- x$CloneID
  if (length(locname)!= length(unique(locname))) {
    cat("Warning: Locus names are not unique!\n")
    cat("         Duplicated names:\n")
    noccur <- table(locname)
    cat(paste("              ",names(noccur[noccur>1])),"\n")
    cat("         Rendering locus names unique with sequential suffix _1, _2 etc\n")
    locname <- make.unique(x$CloneID, sep="_")
  }

# Create the genInd object
  glind <- new("genind", tab=t(snpdata), ploidy=2, ind.names=colnames(snpdata), loc.names=locname, parallel=F)
  colnames(glind@tab) <- locname # Seems to ignore the loc.names=locname in line above

# Add in the standard metadata
  glind@other$loc.metrics <- locus.metadata

# Add in extra metadata -- population assignments
if (!is.null(ind.metafile)) {
  cat("Adding population assignments and other additional individual metadata from file :", ind.metafile,"\n")
  ind.metadata <- read.csv(ind.metafile, header=T, stringsAsFactors=T)
  # Remove leading and trailing spaces that could lead to a spurious mismatch
  ind.metadata$id <- gsub("^\\s+|\\s+$", "", ind.metadata$id)
  # Check that the number of individuals in the metafile is the same as in the dataset
#  if(length(ind.metadata[1])!=nInd(gl)) {
#    cat("Fatal Error: Number of individuals in metadata file does not equal number of individuals in the dataset\n"); stop()
#  }
  # Check for an entry for every individual
  id.col = match( "id", names(ind.metadata))
  if (is.na(id.col)) {
    cat ("Fatal Error: No id column present!\n") ;stop()
    } else {
    if (sum(ind.metadata[,id.col] == names(snpdata))== nind ) {
      cat ("Ids of individual metadata file match!\n")
    }else {
        cat("Fatal Error: Ids in files ",datafile,"and ",ind.metafile," do not match!\n\n");stop()
    }
  }
  pop.col = match( "pop", names(ind.metadata))
  # Check for population assignment
  if (is.na(pop.col)) {
    cat ("Warning: No pop column present\n")
  } else {
    pop(glind) <- as.factor(ind.metadata[,pop.col])
    cat("Populations assigned to individuals\n")
  }
  # Check for latitude and longitude data
  lat.col = match( "lat", names(ind.metadata))
  lon.col = match( "lon", names(ind.metadata))
  if (is.na(lat.col)) {
   cat ("Warning: No lat column present\n")
  }
  if (is.na(lon.col)) {
    cat ("Warning: No lon column present\n")
  }
  if (!is.na(lat.col) & !is.na(lon.col))  {
    glind@other$latlong <- ind.metadata[,c(lat.col, lon.col)]
    cat("Added latlon data\n" )
  }
  # Check for other metadata
  known.col <- names( ind.metadata) %in% c("id","pop", "lat", "lon")
  # known.col <- ifelse(is.na(known.col), , known.col)
  other.col <- names(ind.metadata)[!known.col]
  if (length(other.col>0) )
  {
    glind@other$ind.metrics<-ind.metadata[,other.col]
    cat("Added ",other.col," to the other$ind.metrics slot\n")
  }
}
# Report
  cat("GenInd object created")

  return <- glind

}
