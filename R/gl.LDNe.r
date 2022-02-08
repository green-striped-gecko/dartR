#' @name gl.LDNe
#' @title Estimate effective population size using the Linkage Disequilibrium 
#' method based on NeEstimator (V2)
#' @description
#' This function is basically a convinience function that runs the LD Ne
#'  estimator using Neestimator2 
#'  (\url{http://www.molecularfisherieslaboratory.com.au/neestimator-software/}) 
#'  within R using the provided genlight object. To be able to do so, the 
#'  software has to be downloaded from their website and the appropriate 
#'  executable Ne2-1 has to be copied into the path as specified in the function 
#'  (see example below).
#' @param x Name of the genlight object containing the SNP data [required].
#' @param outfile File name of the output file [default 'genepopLD.txt'] with 
#' all results from Neestimator 2.
#' @param outpath Path where to save the output file
#' [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath='.'
#'  when calling this function to direct output files to your working directory.
#'  @param neest.path path to the folder of the   NE2-1 file. 
#'  [Please note there are 3 different executables depending on your OS: Ne2-1.exe (=Windows), Ne2-1M (=Mac), Ne2-L (=Linux)]. 
#'  You only need to point to the folder (the function will recognise which OS 
#'  you are running).
#' @param critical (vector of) Critical values that are used to remove alleles 
#' based on their minor allele frequency. This can be done before using the 
#' gl.filter.maf function, therefore the default is set to 0 (no loci are 
#' removed). To run for MAF 0 and MAF 0.05 at the same time specify: critical = 
#' c(0,0.05)
#' @param singleton.rm use this if you want to remove singleton alleles (=TRUE)
#'  [default TRUE].
#' @param mating use formula for Random mating='random' or monogamy= 'monogamy'. 
#' [default 'random'].
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session
#' temporary directory (tempdir) [default FALSE].
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

gl.LDNe <- function(x,
                    outfile = "genepopLD.txt",
                    outpath = tempdir(),
                    neest.path = "./",
                    critical = 0,
                    singleton.rm = TRUE,
                    mating = 'random',
                    plot.out = TRUE, 
                    save2tmp = FALSE,
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
  
  # DO THE JOB
  xx <- gl2genepop(x, outfile = "dummy.gen", outpath = tempdir())
  
  if (singleton.rm == TRUE){
    critical[length(critical) + 1] <- 1
  }
  #copy info file to tempdir
  info <- NA
  info[1] <- "1"
  info[2] <- "./"  #path of input file
  info[3] <- "dummy.gen"  #input file
  info[4] <- 2  #Genepop format
  info[5] <- "./" #path of output file
  info[6] <- outfile #output file
  info[7] <- length(critical)
  info[8] <- paste(critical, collapse = " ")
  mm <- pmatch(mating, c("random", "mono")) - 1
  if (mm == 0 | mm == 1){
    info[9] <- mm
  } else{
    cat(error("  Mating is not either 'random' or 'monogamy'. Please check\n"))
  }
  
  con <- file(file.path(tempdir(), "infodummy"), "w")
  writeLines(info, con)
  close(con)
  
  if (Sys.info()["sysname"] == "Windows"){
    prog <- "Ne2-1.exe"
    cmd <- "Ne2-1.exe i:infodummy"
  }
  
  if (Sys.info()["sysname"] == "Linux"){
    prog <- "Ne2-1L"
    cmd <- "./Ne2-1L i:infodummy"
  }
  
  if (Sys.info()["sysname"] == "Darwin") {
    prog <- "Ne2-1M"
    cmd <- "./Ne2-1M i:infodummy"
  }
  
  #check if file program can be found
  if (file.exists(file.path(neest.path, prog))){
    file.copy(file.path(neest.path, prog),
              to = tempdir(),
              overwrite = TRUE)
  }else{
    cat(error(
        "  Cannot find",
        prog,
        "in the specified folder given by neest.path:",
        neest.path, "\n"
      ))
  }
  
  #change into tempdir (run it there)
  old.path = getwd()
  setwd(tempdir())
  system(cmd)
  res <- read.delim(outfile)
  res <- unlist(lapply(res[,1],function(x){x<- gsub(pattern="Infinite",replacement="Inf",x)}))
  
  pops <- sapply(res[res %like% "Population"], function(x) str_extract(x, "(?<=\\[).*(?=\\])"),USE.NAMES = F)
  pops <- sub("_[^_]+$", "", pops)
  freq <- str_split(res[res %like% "Lowest"],'\\s{3,}')[[1]][-1]
  Estimated_Ne <- lapply(res[res %like% "Estimated"], function(x) strsplit(x, "\\s{3,}")[[1]][-1])
  CI_low_Parametric <- lapply(res[res %like% "* Parametric"], function(x) strsplit(x, "\\s{3,}")[[1]][-1])
  CI_high_Parametric <- lapply(res[grep( "^\\* Parametric",res)+1], function(x) strsplit(x, "\\s{3,}")[[1]][-1])
  CI_low_JackKnife <- lapply(res[res %like% "* JackKnife"], function(x) strsplit(x, "\\s{3,}")[[1]][-1])
  CI_high_JackKnife <- lapply(res[grep( "^\\* JackKnife",res)+1], function(x) strsplit(x, "\\s{3,}")[[1]][-1])
  harmonic_mean <- lapply(res[res %like% "Harmonic"], function(x) strsplit(x, "\\s{3,}")[[1]][-1])
  comparisons <- lapply(res[res %like% "Independent"], function(x) strsplit(x, "\\s{3,}")[[1]][-1])
  overall_r2 <- lapply(res[res %like% "OverAll"], function(x) strsplit(x, "\\s{3,}")[[1]][-1])
  expected_r2 <- lapply(res[res %like% "Expected"], function(x) strsplit(x, "\\s{3,}")[[1]][-1])
  
  pop_list <-  lapply(1:length(pops),function(i){
     df_temp <- as.data.frame(cbind(c("Lowest Allele Frequency Used",
                       "Harmonic Mean Sample Size",
                       "Independent Comparisons",
                       "OverAll r^2",
                       "Expected r^2 Sample",
                       "Estimated Ne^",
                       "CI low Parametric",
                       "CI high Parametric",
                       "CI low JackKnife",
                       "CI high JackKnife"),
                     rbind(freq,
                           as.numeric(harmonic_mean[[i]]),
                           as.numeric(comparisons[[i]]),
                           as.numeric(overall_r2[[i]]),
                           as.numeric(expected_r2[[i]]),
                           as.numeric(Estimated_Ne[[i]]),
                           as.numeric(CI_low_Parametric[[i]]),
                           as.numeric(CI_high_Parametric[[i]]),
                           as.numeric(CI_low_JackKnife[[i]]),
                           as.numeric(CI_high_JackKnife[[i]]))
      ))
     colnames(df_temp) <- c("Statistic",paste("Frequency",1:length(freq)))
     rownames(df_temp) <- 1:nrow(df_temp)
     return(df_temp)
     })
  
  names(pop_list) <- pops
                      
  file.copy(outfile, file.path(outpath, outfile))
  setwd(old.path)
  
  # PLOTS
  pop_list_plot <- lapply(pop_list,function(x){
    setNames(data.frame(t(x[,-1])), x[,1])
  })
  
  pop_list_plot <- lapply(1:length(pops),function(i){
    pop_temp <- pop_list_plot[[i]]
    pop_temp$pop <- pops[i]
    return(pop_temp)
    })
  
  # pop_list_plot[[3]]$pop <- "other"
  
  pop_list_plot <- as.data.frame(rbindlist(pop_list_plot))
  pop_list_plot$`Estimated Ne^` <- as.numeric(  pop_list_plot$`Estimated Ne^` )
  pop_list_plot$`CI low Parametric` <- as.numeric(  pop_list_plot$`CI low Parametric`)
  pop_list_plot$`CI high Parametric` <- as.numeric(  pop_list_plot$`CI high Parametric`)
  pop_list_plot$pop <- as.factor(pop_list_plot$pop)
  pop_list_plot[pop_list_plot==Inf] <- NA
  pop_list_plot <- pop_list_plot[which(pop_list_plot$`Lowest Allele Frequency Used`==freq[1]),]
  pop_list_plot <- unique(pop_list_plot)
  # pop_list_plot <- pop_list_plot[c(1,3),]
  
  # str(pop_list_plot)
  
  
  p3 <-
    ggplot(pop_list_plot ) + 
    geom_point(aes(x =pop,y=`Estimated Ne^`,color= pop,group=pop)) +
   # geom_errorbar(aes(x=pop,ymin=`CI low Parametric`, ymax= `CI high JackKnife`)) +
    # ylim(min, 1) +
    # ylab(" ") + 
     # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    # scale_color_manual(values=c('#999999','#E69F00')) +
    ggtitle("title1")
  

  
  # PRINTING OUTPUTS
  if (plot.out) {
    # using package patchwork
    # p3 <- (p1 / p2) + plot_layout(heights = c(1, 4))
    print(p3)
  }
  print(pop_list, row.names = FALSE)
  
  # SAVE INTERMEDIATES TO TEMPDIR
  
  # creating temp file names
  if (save2tmp) {
    if (plot.out) {
      temp_plot <- tempfile(pattern = "Plot_")
      match_call <-
        paste0(names(match.call()),
               "_",
               as.character(match.call()),
               collapse = "_")
      # saving to tempdir
      saveRDS(list(match_call, p3), file = temp_plot)
      if (verbose >= 2) {
        cat(report("  Saving the ggplot to session tempfile\n"))
      }
    }
    temp_table <- tempfile(pattern = "Table_")
    saveRDS(list(match_call, pop_list), file = temp_table)
    if (verbose >= 2) {
      cat(report("  Saving tabulation to session tempfile\n"))
      cat(
        report(
          "  NOTE: Retrieve output files from tempdir using gl.list.reports() and gl.print.reports()\n"
        )
      )
    }
  }
  
  if(verbose >= 1){
  cat(report("  The results are saved in:", file.path(outpath, outfile), "\n"
  ))
  }
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # RETURN
  return(pop_list)
}
