#' Converts genlight objects to STRUCTURE formated files
#'
#' This function exports genlight objects to STRUCTURE formatted files (be aware there is a gl2faststruture version as well). It is based on the code provided by Lindsay Clark (originally for genind objects: [Lindsay happy to link to your github repo if you want]) and this function is basically a wrapper around his genind2structure function.
#' @param gl -- genlight containing lat longs  [required]
#' @param pops -- switch if population column should be added
#' @param outfile -- name (path) of the output shape file
#' @param outpath -- path of the output file. Default is to tempdir(). If to be saved in the current working directory change to "."
#' @param v -- verbosity: if v=0 no output, v=1 reports name and path of output file. default 1
#' @export
#' @author Lindsay V. Clark [your email], wrapper by Bernd Gruber 
#' @examples
#' \dontrun{
#' gl2structure(testset.gl)
#'}



gl2structure <- function(gl, pops=FALSE, outfile="gl.str",  outpath=getwd(), v=1){
  if(!"genlight" %in% class(gl)){
    warning("Function was designed for genlight objects.")
  }
  # internally convert to genind
  gi <- gl2gi(gl, v = 0)
  # get the max ploidy of the dataset
  pl <- max(gi@ploidy)
  # get the number of individuals
  S <- nInd(gi)
  # column of individual names to write; set up data.frame
  tab <- data.frame(ind=rep(indNames(gi), each=pl))
  # column of pop ids to write
  if(pops){
    popnums <- 1:nPop(gi)
    names(popnums) <- as.character(unique(pop(gi)))
    popcol <- rep(popnums[as.character(pop(gi))], each=pl)
    tab <- cbind(tab, data.frame(pop=popcol))
  }
  loci <- locNames(gi)
  # add columns for genotypes
  tab <- cbind(tab, matrix(-9, nrow=dim(tab)[1], ncol=nLoc(gi), dimnames=list(NULL,loci)))
  # begin going through loci
  for(L in loci){
    thesegen <- gi@tab[,grep(paste(L, ".", sep=""), dimnames(gi@tab)[[2]], fixed=TRUE), drop=FALSE] # genotypes by locus
    al <- 1:dim(thesegen)[2] # numbered alleles
    for(s in 1:S){
      if(all(!is.na(thesegen[s,]))){
        tabrows <- (1:dim(tab)[1])[tab[[1]] == indNames(gi)[s]] # index of rows in output to write to
        tabrows <- tabrows[1:sum(thesegen[s,])] # subset if this is lower ploidy than max ploidy
        tab[tabrows,L] <- rep(al, times = thesegen[s,])
      }
    }
  }
  # export table
  write.table(tab, file=file.path(outpath, outfile), sep="\t", quote=FALSE, row.names=FALSE)
  if (v==1)  cat(paste("Structure file saved as:", outfile,"\nin folder:",outpath))
}