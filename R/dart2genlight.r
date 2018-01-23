#' Convert DarT to genlight
#' 
#' @description converts a dart file (read via \code{read.dart}) into an genlight object \code{\link{adegenet}}. Internal function called by gl.read.dart
#' @param dart a dart object created via read.dart
#' @param ind.metafile optional file in csv format with covariates for each individual (see details for explanation)
#' @param covfilename depreciated, use parameter ind.metafile
#' @param probar show progress bar
#' @return a genlight object is returned. Including all available slots are filled. loc.names, ind.names, pop, lat, lon (if provided via the covariate file)
#' @details the covariate file needs to have very specific headings. First an heading called id. Here the ids have to match the ids in the dart object \code{colnames(dart[[4]])}. The following column headings are optional. pop: specifies the population membership of each individual. lat and lon specify spatial coordinates (perferable in decimal degrees WGS1984 format). Additional columns with individual covariates can be imported (e.g. age, gender).
#' @examples
#' \dontrun{
#' dgl <- dart2genlight(dart, "covariates.csv")
#' }


dart2genlight <- function(dart, ind.metafile=NULL, covfilename=NULL, probar = TRUE)
 {

if (is.null(ind.metafile)) {ind.metafile <- covfilename}
 
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

if (!is.null(ind.metafile))
{
cat(paste("Try to add covariate file:", ind.metafile,".\n"))
###### population and individual file to link AAnumbers to populations...
ind.cov <- read.csv(ind.metafile, header=T, stringsAsFactors=T)
# is there an entry for every individual

id.col = match( "id", names(ind.cov))



if (is.na(id.col)) {cat ("There is no id column\n") ;stop()} else {
  ind.cov[,id.col]<- trimws(ind.cov[,id.col], which="both")  #trim spaces

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

