#' @name gl2migraten
#' @title Converts a genlight object into a migrate-n SNP input file
#' @description
#' This function exports a genlight object into a migrate-n SNP input file
#' 
#' @details The input file generated is compatible with migrate-n (Beerli ) 
#' v3.7.2 onward. For the (limited) testing I have done with migrate-n and SNP data, 
#' migrate-n seems to be able to handle up to 20k SNPs (v3.7.2) and 3-5k SNPs in v4
#' and later. That seems to depend on some memory issues that Peter Beerli may be able
#' to resolve (or may have resolved by the time you use this function).
#' 
#' This function needs PLINK installed, which is used for intermediate data manipulation.
#' 
#' migrate-n and its manual can be downloaded from:
#' \url{https://peterbeerli.com/migrate-html5/index.html}
#' 
#' 
#' @param x Name of the genlight object containing the SNP data [required].
#' @param outfile File name of the output file [default 'migrateHapMap'].
#' @param outpath Path where to save the output file
#' [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath='.'
#'  when calling this function to direct output files to your working directory.
#'  @inheritParams utils.plink.run
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @return NULL
#' @references
#' Beerli, P. (2004). Effect of unsampled populations on the estimation of 
#' population sizes and migration rates between sampled populations. 
#' Molecular Ecology 13, 827-836.
#' 
#' Beerli, P. and Felsenstein, J. (2001). Maximum likelihood estimation of a 
#' migration matrix and effective population sizes in n subpopulations by using 
#' a coalescent approach. Proceedings of the National Academy of Sciences 98, 4563-4568. 
#' doi: 10.1073/pnas.081068098.
#' 
#' Beerli, P. and Palczewski, M. (2010). Unified framework to evaluate panmixia 
#' and migration direction among multiple sampling locations. Genetics 185, 313-326.


#' @rawNamespace import(data.table, except = c(melt,dcast))
#' @export
#' @author Custodian: Carlo Pacioni (Post to
#'  \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' \donttest{

#' }
gl2migraten <- function(x, 
                        plink.cmd="plink", 
                        plink.path="path",
                        outfile = "migrateHapMap",
                        outpath = tempdir(),
                        verbose=NULL) {
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jody",
                   verbosity = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # DO THE JOB
  
  #----------------- Helper FUN -------------------------#
  ifNA <- function(x) {
    return(ifelse(is.na(x), 0, x))
  }
  #------------------------------------------------------#
  # get the data into plink's tped format which is easier to convert
  tmp.out <- tempdir()
  gl2plink(x, plink_path = plink.path, bed_file = FALSE, outfile = "plink.out", outpath = tmp.out)
  
  use.plink(tmp.out, plink.cmd=plink.cmd, plink.path, out="tplink.out", 
            # rm bfile because creating a .bed doesn't work for me on W
            # syntax="--bfile plink.out --recode --transpose")
            syntax="--file plink.out --recode --transpose")
  
  tped.dt <- data.table::fread(file.path(plink.out, "tplink.out"))
  popIDsdip <- rep(pop(x), each=2)
  setnames(tped.dt, new = c("Chr", "snp.name", "gen.pos", "phys.pos", popIDsdip))
  firstLine <- paste("H", length(unique(popIDs)), nrow(tped.dt))
  
  write(x = firstLine, file = file.path(outpath, outfile))
  
  for(p in unique(popIDs)) {
    if (verbose > 2) {
      cat(report("Start converting data for population", p, "\n"))
    }
    
    popLine <- paste(sum(popIDsdip == p), p)
    write(popLine, file = file.path(outpath, outfile), append = TRUE)
    for(i in seq_len(nrow(tped.dt))) {
      # get the allele in the dataset across all pops
      variants <- names(table(as.matrix(tped.dt[i, -c(1:4)])))
      
      # Freq within each pops
      pop_p <- table(as.matrix(tped.dt[i, names(tped.dt) == p, with=FALSE]))
      
      snps <- paste(tped.dt[i, snp.name], variants[1], ifNA(pop_p[variants[1]]), variants[2], ifNA(pop_p[variants[2]]), sum(pop_p)  )
      write(snps, file = file.path(outpath, outfile), append = TRUE)
    }
  }
  
  # FLAG SCRIPT END
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }
  
  invisible(NULL)
}