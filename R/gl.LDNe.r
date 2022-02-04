#' @name gl.LDNe
#' @title Estimate effective population size using the Linkage Disequilibrium method based on NeEstimator (V2)
#' @description
#' This function is basically a convinience function that runs the LD Ne estimator using Neestimator2 (\url{http://www.molecularfisherieslaboratory.com.au/neestimator-software/}) within R using the provided genlight object. To be able to do so, the software has to be downloaded from their website and the appropriate executable Ne2-1 has to be copied into the path as specified in the function (see example below).  
#' @param x Name of the genlight object containing the SNP data [required].
#' @param outfile File name of the output file [default 'genepopLD.txt'] with all results from Neestimator 2.
#' @param outpath Path where to save the output file
#' [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath='.'
#'  when calling this function to direct output files to your working directory.
#'  @param neest.path path to the folder of the   NE2-1 file. [Please note there are 3 different executables depending on your OS: Ne2-1.exe (=Windows), Ne2-1M (=Mac), Ne2-L (=Linux)]. You only need to point to the folder (the function will recognise which OS you are running) 
#' @param critical (vector of) Critical values that are used to remove alleles based on their minor allele frequency. This can be done before using the gl.filter.maf function, therefore the default is set to 0 (no loci are removed). To run for MAF 0 and MAF 0.05 at the same time specify: critical = c(0,0.05)
#' @param singleton.rm use this if you want to remove singleton alleles (=TRUE) [default set to TRUE].
#' @param mating use formula for Random mating='random' or monogamy= 'monogamy'. Default setting is 'random'.
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @return invisible results as table
#' @author Custodian: Bernd Gruber (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' \dontrun{
#' # SNP data (use two populations and only the first 100 SNPs)
#' pops <- possums.gl[1:60,1:100]
#' nes <- gl.LDNe(pops, outfile="popsLD.txt", outpath=tempdir(),
#' neest.path = "./path_to Ne-21",
#' critical=c(0,0.05), singleton.rm=TRUE, mating='random')
#' nes
#' }
#' @export




gl.LDNe <- function(  x, 
                      outfile = "genepopLD.txt", 
                      outpath=tempdir(),
                      neest.path ="./",
                      critical = 0,
                      singleton.rm=TRUE,
                      mating='random',
                      verbose=NULL)
{
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jody",
                   verbosity = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  #works only with SNP data
  if (datatype!="SNP") cat(error(
    "Only SNPs (diploid data can be transformed into genepop format!\n" ))
  
  # DO THE JOB
  xx <- gl2genepop(x, outfile="dummy.gen", outpath=tempdir())
  
  
  if (singleton.rm==TRUE) critical[length(critical)+1]<- 1
  #copy info file to tempdir
  info <- NA
  info[1] <- "1"
  info[2] <- ""  #path of input file
  info[3] <- "dummy.gen"  #input file
  info[4] <- 2  #Genepop format
  if (substr(outpath, nchar(outpath),500)!="/") op <- paste0(outpath,"/") else op<- outpath
  info[5] <- paste0(op) #path of output file
  info[6] <- outfile #output file
  info[7] <- length(critical)
  info[8] <- paste(critical, collapse = " ")
  mm <- pmatch(mating, c("random","mono"))-1
  if (mm==0 | mm==1)  info[9]<- mm else cat(error("mating is not either 'random' or 'monogamy'. Please check"))
  con <- file(file.path(tempdir(),"infodummy"),"w")
  writeLines(info, con)
  close(con)
  
  
  if (Sys.info()["sysname"] == "Windows") 
  {
   prog <- "Ne2-1.exe"
   cmd <- "Ne2-1.exe i:infodummy"
    
  }  
  if (Sys.info()["sysname"] == "Linux") 
  {
    prog <- "Ne2-1L"
    cmd <- "./Ne2-L i:infodummy"
  }  
  if (Sys.info()["sysname"] == "Darwin") 
  {
    prog <- "Ne2-1M"
    cmd <- "./Ne2-1M i:infodummy"
  }  
  
  #check if file program can be found
  if(file.exists(file.path(neest.path,prog)))
    file.copy(file.path(neest.path,prog), to = tempdir(), overwrite = TRUE ) else cat(error(paste0("Cannot find ", prog," in the specified folder given by neest.path: ",neest.path))) 
  
    #change into temddir (run it there)
    old.path=getwd()
    setwd(tempdir())
    ret <- system(cmd)
    setwd(old.path)
    res <- readLines(file.path(outpath, outfile))
    cat(message(paste0("The results are saved in: ", file.path(outpath, outfile))))
  return(res)
}
  

#nes <- gl.LDNe(pops, outfile="popsLD2.txt", outpath="d:/temp",neest.path = "d:/programms/Neestimator/",critical=c(0), singleton.rm=FALSE, mating='mon')
#nes

