#' @name gl2genepop
#' @title Converts a genlight object to genepop format (and file)
#' @description
#' The genepop format is used by several external applications (for example
#' Neestimator2
#' (\url{http://www.molecularfisherieslaboratory.com.au/neestimator-software/}).
#' So the main idea is to create the genepop file and then run the other
#' software externally. As a feature the genepop file also returned as an
#' invisible data.frame by the function.
#' @param x Name of the genlight object containing the SNP data [required].
#' @param outfile File name of the output file [default 'genepop.gen'].
#' @param outpath Path where to save the output file
#' [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath='.'
#'  when calling this function to direct output files to your working directory.
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @return invisible data frame in genepop format
#' @author Custodian: Bernd Gruber (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' \dontrun{
#' # SNP data
#' geno <- gl2genepop(testset.gl[1:3,1:9])
#' head(geno)
#' }
#' @export

gl2genepop <- function (x,
                        outfile = "genepop.gen",
                        outpath = tempdir(),
                        verbose = NULL) {

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
  if (datatype != "SNP") {
    cat(error(
      "  Only SNPs (diploid data can be transformed into genepop format!\n"
    ))
  }
                           
      if (is.null(pop(x))) {
    cat(important("Your genlight object does not have a population definition. Therefore the function assumes you the whole genlight object to be one population!\n"))
    pop(x) <- rep("Pop1", nInd(x))
    }
  
  # DO THE JOB
  
  #convert to genind
  x<- x[order(pop(x)),]
  x <- gl2gi(x, verbose = 0,probar = FALSE)
  data <- as.matrix(x)
  pop_names <- x@pop
  
  if (all(unlist(unique(x@all.names)) %in% c("A", "T", 
                                             "C", "G"))) {
    m_type <- "snp"
    #message("Your dataset is treated as a SNP dataset.\n            Alleles initially coded A, T, C, G were respectively coded\n            01, 02, 03 and 04")
    colnames(data) <- gsub(colnames(data), pattern = "\\.A", 
                           replacement = ".01")
    colnames(data) <- gsub(colnames(data), pattern = "\\.T", 
                           replacement = ".02")
    colnames(data) <- gsub(colnames(data), pattern = "\\.C", 
                           replacement = ".03")
    colnames(data) <- gsub(colnames(data), pattern = "\\.G", 
                           replacement = ".04")
  }
  
  loci_names_l <- x@loc.fac
  loc_all <- data.frame(col = colnames(data))
  if (all(stringr::str_count(colnames(data), "\\.") == 
          1) != TRUE) {
    stop("The columns' names of x@tab must be of form 'locus.allele' with only 1\n         '.' between locus and allele")
  }
  loc_all <- tidyr::separate(loc_all, col = 1, sep = "\\.", 
                             into = c("locus", "allele"))
  loci_names <- as.character(loci_names_l[-which(duplicated(loci_names_l))])
  n.loci <- length(loci_names_l[-which(duplicated(loci_names_l))])
  data_gpop <- data.frame(id = paste(pop_names, "_", 
                                     row.names(data), ",", sep = ""))
  for (i in 1:n.loci) {
    loc <- loci_names[i]
    a <- c()
    for (j in 1:nrow(data)) {
      col_loc <- which(loc_all[, "locus"] == loc)
      hom <- which(data[j, col_loc] == 2)
      het <- which(data[j, col_loc] == 1)
      if (length(hom) != 0) {
        a[j] <- paste(loc_all[col_loc[hom], "allele"], 
                      loc_all[col_loc[hom], "allele"], sep = "")
      }
      else if (length(het) != 0) {
        if (as.character(loc_all[col_loc[het[1]], "allele"]) < 
            as.character(loc_all[col_loc[het[2]], "allele"])) {
          a[j] <- paste(loc_all[col_loc[het[1]], "allele"], 
                        loc_all[col_loc[het[2]], "allele"], 
                        sep = "")
        }
        else {
          a[j] <- paste(loc_all[col_loc[het[2]], "allele"], 
                        loc_all[col_loc[het[1]], "allele"], 
                        sep = "")
        }
      }
      else {
        if (nchar(loc_all[1, "allele"]) == 3) {
          a[j] <- "000000"
        }
        else {
          a[j] <- "0000"
        }
      }
    }
    data_gpop <- cbind(data_gpop, a)
  }
  colnames(data_gpop) <- c("ID", as.character(loci_names))
  data_gpop[, ] <- apply(data_gpop,c(1,2), as.character)
  
  
  dummy<- paste("Genepop output. Loci:", nLoc(x), "Populations:",nPop(x))
  dummy[2] <- paste(locNames(x),collapse=",")
  cs <- c(cumsum(table(pop(x))))
  from=c(1,(cs[-length(cs)]+1))
  to= cs
  for (i in 1:nPop(x))
  {
    da <- apply(data_gpop[from[i]:to[i],], 1, function(y) paste(y, collapse = " "))  
    dummy <- c(dummy,"Pop",da)
  }
  
    data_gpop2<- dummy
    utils::write.table(data_gpop2, file = file.path(outpath,outfile), quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat(report(
    "The genepop file is saved as: ", file.path(outpath, outfile,"\n")
  ))
                
                  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
                
                  # RETURN
                
    invisible(data.frame(lines=data_gpop2))
    
}