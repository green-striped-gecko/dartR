gl2stairwayplot2 <- 
function(x, 
         outfile="blueprint.txt", 
         outpath=tempdir(),
         simfolder=NULL,
         minbinsize=1,
         maxbinsize=NULL,
         pct_training=0.67,
         nrand=NULL,
         stairway_plot_dir="stairway_plot_es",
         nreps=200,
         seed=NULL,
         mu=NULL,
         gentime=NULL,
         L = NULL,
         plot_title="Ne against time",
         xmin=0, xmax=0, ymin=0, ymax=0,
         xspacing=2,
         yspacing=2,
         fontsize=12,
         run=FALSE,
         parallel=1,
         verbose=NULL) {
  
  # TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Juliette"
  outfilespec <- file.path(outpath, outfile)
  
  # SET VERBOSITY
  
  if (is.null(verbose)){ 
    if(!is.null(x@other$verbose)){ 
      verbose <- x@other$verbose
    } else { 
      verbose <- 2
    }
  } 
  
  if (verbose < 0 | verbose > 5){
    cat(paste("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n"))
    verbose <- 2
  }
  
  # FLAG SCRIPT START
  
  if (verbose >= 1){
    if(verbose==5){
      cat("Starting",funname,"[ Build =",build,"]\n")
    } else {
      cat("Starting",funname,"\n")
    }
  }
  
  # STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    stop("  Fatal Error: genlight object required!\n")
  }
  
  if (verbose >= 2){
    if (all(x@ploidy == 1)){
      stop("Fatal Error: Detected Presence/Absence (SilicoDArT) data. Please provide a SNP dataset\n")
    } else if (all(x@ploidy == 2)){
      cat("  Processing a SNP dataset\n")
    } else {
      stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
    }
  }
  
  # SCRIPT SPECIFIC ERROR CHECKING
  
  if(is.null(maxbinsize)){
    maxbinsize <- nInd(x)
    if(verbose >= 3){cat("  Max Bin Size not specified, set to",nInd(x),"\n")}
  }
  if(is.null(nrand)){
    nrand <- c(round((nInd(x)-1)/2), round(nInd(x)-1), round((nInd(x)-1)*3/2), round(2*(nInd(x)-1)))
    if(verbose >= 3){cat("  No. of break points not specified, set to:",paste0(nrand, collapse=" "),"\n")}
  }
  if(is.null(stairway_plot_dir)){
    stop("Fatal Error: Directory path for the Stairway Plot 2 executables not specified\n")
  }
  if(is.null(mu)){
    stop("Fatal Error: Mutation rate per site per generation not specified\n")
  }
  if(is.null(gentime)){
    stop("Fatal Error: Generation time (years) not specified\n")
  }
  whether_folded <- "true"
  project_dir <- simfolder
  nseq <- 2*nInd(x)

  
  
  usfs <- colSums(as.matrix(x), na.rm=T)
  
  #fold usfs (in case of uneven individuals, drop the 0s and minbinsize)
  usfs <- ifelse(usfs>nInd(x), nInd(x)*2-usfs, usfs)
  sfs <- table(usfs)[-(1:minbinsize)]

  #fill in empty bins in sfs with zeros
  if (length(sfs) != (nInd(x)-(minbinsize-1) )) {
    
    sfs2 <- rep(0,(nInd(x)-(minbinsize-1) ))
    sfs2[as.numeric(names(sfs))]  <- sfs
    sfs <- sfs2
    
  }
  
  
 #cut at maxbinsize
 #sfs <- sfs[1:maxbinsize] (not necessary)
    
  sfs <- paste(sfs,collapse = " ")
  
  #if total length of sequence is not specificed simply assume nLoc*69 (standard from dart)
  if (is.null(L)) L = nLoc(x)*69 
  
  if (is.null(seed)) {set.seed= as.numeric(Sys.time()) ;seed = round(runif(1)*1e6)}
  
  # DO THE JOB
  
  if (verbose >= 2) {cat(paste("  Extacting SNP data and creating records for each individual\n"))}
  
  #### NEED TO DISCUSS THE GUTS OF IT HERE
  
  
  # Output the results
  
  write.table(paste("#Ne analysis from SNPs for ",as.character(substitute(x))),
              file=outfilespec,row.names=FALSE,col.names=FALSE,
              quote=FALSE)
  write.table(paste("popid:",simfolder,"# id of the population (no white space)"),
              file=outfilespec,
              row.names=FALSE,col.names=FALSE,
              quote=FALSE, sep=" ",append=TRUE)
  write.table(paste("nseq:",nseq,"# number of haploid sequences = 2n"),
              file=outfilespec,
              row.names=FALSE,col.names=FALSE,
              quote=FALSE, sep=" ",append=TRUE)
  write.table(paste("L:",L,"# total number of nucleic sites, including polymorphic and monomorphic"),
              file=outfilespec,
              row.names=FALSE,col.names=FALSE,
              quote=FALSE, sep=" ",append=TRUE)
  write.table(paste("whether_folded:", whether_folded, "# whethr the SFS is folded (true or false)"),
              file=outfilespec,
              row.names=FALSE,col.names=FALSE,
              quote=FALSE, sep=" ",append=TRUE)
  write.table(paste("SFS:", sfs, "# snp frequency spectrum: number of singleton, number of doubleton, etc. (separated by white space)"),
              file=outfilespec,
              row.names=FALSE,col.names=FALSE,
              quote=FALSE, sep=" ",append=TRUE)
  write.table(paste("smallest_size_of_SFS_bin_used_for_estimation:", minbinsize, "# default is 1; to ignore singletons, change this number to 2"),
              file=outfilespec,
              row.names=FALSE,col.names=FALSE,
              quote=FALSE, sep=" ",append=TRUE)
  write.table(paste("largest_size_of_SFS_bin_used_for_estimation:", maxbinsize, "# default is nseq/2 for folded SFS"),
              file=outfilespec,
              row.names=FALSE,col.names=FALSE,
              quote=FALSE, sep=" ",append=TRUE)
  write.table(paste("pct_training:", pct_training, "# proportion of sites for training"),
              file=outfilespec,
              row.names=FALSE,col.names=FALSE,
              quote=FALSE, sep=" ",append=TRUE)
  write.table(paste("nrand:", paste0(nrand, collapse=" "), "# number of random break points for each try (separated by white space)"),
              file=outfilespec,
              row.names=FALSE,col.names=FALSE,
              quote=FALSE, sep=" ",append=TRUE)
  write.table(paste("project_dir:", simfolder, "# project directory"),
              file=outfilespec,
              row.names=FALSE,col.names=FALSE,
              quote=FALSE, sep=" ",append=TRUE)
  write.table(paste("stairway_plot_dir:", stairway_plot_dir, "# directory to the stairway plot files"),
              file=outfilespec,
              row.names=FALSE,col.names=FALSE,
              quote=FALSE, sep=" ",append=TRUE)
  write.table(paste("ninput:", nreps, "# number of input files to be created for each estimation"),
              file=outfilespec,
              row.names=FALSE,col.names=FALSE,
              quote=FALSE, sep=" ",append=TRUE)
  write.table(paste("random_seed:", seed),
              file=outfilespec,
              row.names=FALSE,col.names=FALSE,
              quote=FALSE, sep=" ",append=TRUE)
  write.table(paste("mu:", mu, "# assumed mutation rate per site per generation"),
              file=outfilespec,
              row.names=FALSE,col.names=FALSE,
              quote=FALSE, sep=" ",append=TRUE)
  write.table(paste("year_per_generation:", gentime,"# assumed generation time (in years)"),
              file=outfilespec,
              row.names=FALSE,col.names=FALSE,
              quote=FALSE, sep=" ",append=TRUE)
  write.table(paste("plot_title:", plot_title, "# title of the plot"),
              file=outfilespec,
              row.names=FALSE,col.names=FALSE,
              quote=FALSE, sep=" ",append=TRUE)
  write.table(paste("xrange:", paste0(xmin,',',xmax), "# Time (1k year) range; format: xmin,xmax; 0,0 for default"),
              file=outfilespec,
              row.names=FALSE,col.names=FALSE,
              quote=FALSE, sep=" ",append=TRUE)
  write.table(paste("yrange:", paste0(ymin,",",ymax), "# Ne (1k individual) range; format: ymin,ymax; 0,0 for default"),
              file=outfilespec,
              row.names=FALSE,col.names=FALSE,
              quote=FALSE, sep=" ",append=TRUE)
  write.table(paste("xspacing:", xspacing, "# X axis spacing"),
              file=outfilespec,
              row.names=FALSE,col.names=FALSE,
              quote=FALSE, sep=" ",append=TRUE)
  write.table(paste("yspacing:",yspacing, "# Y axis spacing"),
              file=outfilespec,
              row.names=FALSE,col.names=FALSE,
              quote=FALSE, sep=" ",append=TRUE)
  write.table(paste("fontsize:", fontsize, "# Font size"),
              file=outfilespec,
              row.names=FALSE,col.names=FALSE,
              quote=FALSE, sep=" ",append=TRUE)
  
  if (verbose > 2) {cat(paste("    Stairway Plot 2 blueprint written to",outfilespec,"\n"))}
  
  # FLAG SCRIPT END
  
  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
  
  if (run==TRUE)
  {
    setwd(outpath)
    system(paste0("java -cp stairway_plot_es Stairbuilder ",outfile ))
    
    
  
    
    #run on multiple cores
    if (parallel>1)
    {
      ff <- readLines(paste0(outfile, ".bat"))
      no_cores <- min(parallel, detectCores()-1)
      plan(multisession, workers = no_cores)
      
      index = grep("Stairway_fold_training_testing7", ff)
      
      runs <- ff[index]
      
      runstair <- function(i) {
        system(runs[i])
        return(i)
      }
      
      temp <- future_map(1:length(runs), function(x) runstair(x))
      
      index = grep("MOVE", ff)
      runs <- ff[index]
      runs <- gsub("MOVE /y", "cp " ,runs)
      
      for (i in 1:length(runs)) system(runs[i])
      
      index = grep("Stairpainter", ff)
      system( ff[index])
      er <- system(paste0(outfile, ".plot.bat"))
      
    } else er <- system(paste0(outfile, ".bat"))
    
    
    
    if (er==1)
    {
    cat("Attempt to rerun last step with different settings (lower memory allocation fo the Java Virtual Machine")
    ff <- readLines(paste0(outfile, ".plot.bat"))
    ff <-  gsub("-Xmx4g", "-Xmx1g", ff)
    con <- file(paste0(outfile, ".plot.bat"),"w")
    writeLines(ff, con)
    close(con)
    system(paste0(outfile, ".plot.bat"))
    }
    cat("Check plots (pdf and png files) in folder:", file.path(outpath, simfolder),".\n")
    
  }

  return(NULL)
  
}
