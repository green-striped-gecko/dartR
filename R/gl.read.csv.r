#' Read SNP data from a csv file into a genlight object 
#'
#' This script takes SNP genotypes from a csv file, combines them with individual and locus metadata and creates a genlight object.
#' 
#' The SNP data need to be in one of two forms. SNPs can be coded 0 for homozygous reference, 2 for homozygous alternate, and 1 for heterozygous with NA for missing values; or
#' it can be coded A/A, A/C, C/T, G/A etc, with -/- as missing. Other formats will thow an error.
#' 
#' The individual metadata needs to be in a csv file, with headings, with a mandatory id column corresponding exactly to the identities provided with the SNP data. The locus metadata needs
#' to be in a csv file with headings. Note that the locus metadata will be complemented by calculable statistics corresponding to those that would be provided by Diversity Arrays
#' Technology (e.g. CallRate)
#'
#' @param filename -- name of the csv file containing the SNP genotypes [required]
#' @param ind.labels -- if TRUE, then the first column of the csv file is the individual/specimen/sample labels [default TRUE]
#' @param loc.labels -- if TRUE, then the first row of the csv file is the locus labels [default TRUE]
#' @param ind.metafile -- name of the csv file containing the metadata for individuals [optional]
#' @param loc.metafile -- name of the csv file containing the metadata for loci [optional]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return a genlight object with the SNP data and associated metadata included.
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})

# Last amended 3-Feb-19

gl.read.csv <- function(filename, ind.labels=TRUE, loc.labels=TRUE, ind.metafile=NULL, loc.metafile=NULL, verbose=2){
  
# TIDY UP FILE SPECS

  funname <- match.call()[[1]]

# FLAG SCRIPT START

  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }

  if (verbose > 0) {
    cat("Starting",funname,"\n")
  }

# STANDARD ERROR CHECKING
  
# FUNCTION SPECIFIC ERROR CHECKING

# DO THE JOB

  #FIRST THE SNP DATA

  # Create the SNP data matrix, indNames and LocNames
  df0 <- read.csv(file=filename,header=FALSE)
  
  rowstart <- 1  
  if (loc.labels) {
    rowstart <- 2
  }
  colstart <- 1
  if (ind.labels){
    colstart <- 2
  }
  numrows <- dim(df0)[1]
  numcols <- dim(df0)[2]

  df <- df0[rowstart:numrows,colstart:numcols]
    data <- as.matrix(df)
  loci <- df0[1,colstart:numcols]
    loci <- as.character(as.matrix(loci))
  individuals <- df0[rowstart:numrows,1]
    individuals <- as.character(individuals)
    
  test <- paste(data,collapse="")
  test <- gsub("NA","9",test)
  if(nchar(test) > numrows*numcols) {
      if(verbose >=2){
        cat("Character data detected, assume genotypes are of the form C/C, A/T, C/G, -/- etc\n")
      }
      # Check that this is true
      s1 <- paste(data, collapse=" ")
      s1 <- gsub("/"," ",s1)
      s1 <- toupper(s1)
      s2 <- unlist(strsplit(s1, " "))
      tmp <- table(s2)
      names(tmp)
      if (!(names(tmp) == "A" || names(tmp) == "C" || names(tmp) == "G" || names(tmp) == "T" || names(tmp) == "-")){
        cat("Fatal Error: Genotypes must be defined by the letters A, C, G, T or missing -\n")
        stop()
      }
    # Check that the data are biallelic
      for (i in 1:dim(data)[2]){
        v1 <- data[,i]
        v1 <- paste(v1, collapse=" ")
        v1 <- gsub("/"," ",v1)
        v1 <- gsub("- ","",v1)
        v1 <- toupper(v1)
        v1 <- unlist(strsplit(v1, " "))
        tmp <- table(v1)
        if (length(names(tmp)) > 2) {
          cat("Fatal Error: Loci are not biallelic\n")
          stop()
        } 
      }
      if ( verbose >= 2) {cat ("  Data confirmed as biallelic\n")}
      
    # Step through and convert data to 0, 1, 2, NA
      homRef <- paste0(names(tmp)[1],"/",names(tmp)[1])
      homAlt <- paste0(names(tmp)[2],"/",names(tmp)[2])
      het1 <- paste0(names(tmp)[1],"/",names(tmp)[2])
      het2 <- paste0(names(tmp)[2],"/",names(tmp)[1])
      missing <- "-/-"
      for (i in 1:dim(data)[2]){
        data[,i] <- toupper(data[,i])
        data[,i] <- gsub(homRef,"0",data[,i])
        data[,i] <- gsub(homAlt,"2",data[,i])
        data[,i] <- gsub(het1,"1",data[,i])
        data[,i] <- gsub(het2,"1",data[,i])
        data[,i] <- gsub(missing,NA,data[,i])
      }
      cat( "  SNP coding converted to 0,1,2,NA\n")

      data <- apply(data,2,as.numeric)
  
  } else {
      if (verbose >= 2) {cat("Numeric data detected, assume genotypes are 0 = homozygous reference, 1 = heterozygous, 2 = homozygous alternate\n")}
      # Check that this is true
      s1 <- paste(data, collapse=" ")
      s2 <- unlist(strsplit(s1, " "))
      tmp <- table(s2)
      names(tmp)
      if (!(names(tmp) == "0" || names(tmp) == "1" || names(tmp) == "2" || names(tmp) == "NA")){
        cat("Fatal Error: Genotypes must be defined by the numbers 0, 1, 2 or missing NA\n")
        stop()
      }
      data <- apply(data,2,as.numeric)
  }    
      
  gl <- as.genlight(data)
  locNames(gl) <- loci
  indNames(gl) <- individuals

  # SECOND THE LOCUS METADATA
  if(!is.null(loc.metafile)){
    loc.metadata <- read.csv(file=loc.metafile,header=TRUE)
    gl@other$loc.metrics <- loc.metadata
    gl <- gl.recalc.metrics(gl)
  } else {
    gl <- gl.recalc.metrics(gl)
  }
  if (verbose >= 2){cat(paste(" Added or updated ",names(gl@other$loc.metrics)," to the other$ind.metrics slot.\n"))}
  
  # THIRD THE INDIVIDUAL METADATA
  if(!is.null(ind.metafile)){
    ind.metadata <- read.csv(file=ind.metafile,header=TRUE)
    if (!("id" %in% names(ind.metadata))) {cat("Fatal Error: mandatory id column absent from the individual metadata file\n")}
    if (!("pop" %in% names(ind.metadata))) {
      cat("  Warning: pop column absent from the individual metadata file, setting to NA\n")
      gl@other$ind.metrics$pop <- array(NA,nInd(gl))
    }
    gl@other$ind.metrics <- ind.metadata
  } else {
    gl@other$ind.metrics$id <- individuals
    gl@other$ind.metrics$pop <- array(NA,nInd(gl))
  }
  if (verbose >= 2){cat(paste(" Added ",names(gl@other$ind.metrics)," to the other$ind.metrics slot.\n"))}
 
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
      
      return(gl)
      
}    


  
