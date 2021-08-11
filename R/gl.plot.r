#' @name gl.plot
#' @title Plotting genlight object as a smear plot (loci by individuals color 
#' coded for scores of 0, 1, 2 and NA)
#' @description 
#' It adds the option to put labels on the individuals and grouping by populations. 
#' If there are too many individuals, it is best to use ind_labels_size = 0.
#'
#' @param x Name of the genlight object containing the SNP or presence/absence
#'  (SilicoDArT) data [required].
#' @param group_pop Group by population [default TRUE].
#' @param ind_labels Labels for individuals [default indNames(x)].
#' @param ind_labels_size Size of the individual labels, if individual labels 
#' are not required set this parameter to 0 [default 10].
#' @param plot_colours Vector with four color names for homozygotes for the 
#' reference allele, heterozygotes, homozygotes for the alternative allele and 
#' for missing values (NA) [default four_colors].
#' @param posi Position of the legend: “left”, “top”, “right”, “bottom” or
#'  "none" [default = "bottom"].
#' @param save2tmp If TRUE, saves plot to the session temporary directory 
#' (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, 
#' progress log ; 3, progress and results summary; 5, full report [default NULL].
#' 
#' @return Returns unaltered genlight object
#' @author Custodian: Luis Mijangos -- Post to \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' gl.plot(bandicoot.gl[1:10,])
#'
#' @seealso \code{\link{gl.filter.callrate}}
#' @family Exploration/visualisation functions
#' @export
#'  

gl.plot <- function (x,
                     group_pop = FALSE,
                     ind_labels = indNames(x), 
                     ind_labels_size = 10,
                     plot_colours = four_colors, 
                     posi = "bottom", 
                     save2tmp = FALSE,
                     verbose = NULL) {
  
  # CHECK IF PACKAGES ARE INSTALLED
  pkg <- "reshape2"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop(error("Package ",pkg," needed for this function to work. Please install it.")) 
  } 
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func=funname,build="Jackson",v=verbose)
  
  # CHECK DATATYPE 
  datatype <- utils.check.datatype(x, verbose=verbose)
  
  # Set a population if none is specified (such as if the genlight object has been generated manually)
  if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
    if (verbose >= 2) {
      cat(warn("  No population assignments detected, 
                             individuals assigned to a single population labelled 'pop1'\n"))
    }
    pop(x) <- array("pop1", dim = nInd(x))
    pop(x) <- as.factor(pop(x))
  }

  # DO THE JOB

  X_temp <- as.data.frame(as.matrix(x))
  colnames(X_temp) <- 1:nLoc(x)
  X_temp$pop <- pop(x)
  X_temp$id <- ind_labels
  X <- reshape2::melt(X_temp,id.vars = c("pop","id"))
  X$value <- as.character(X$value)

  colnames(X) <- c("pop","id","locus","genotype")
  
  loc_labels <- pretty(1:nLoc(x),5)
  
  locus <- id <- genotype <- NA
  
  if(datatype=="SilicoDArT"){
    p3 <- ggplot(X,aes(x=locus,y=id,fill=genotype))+
      geom_raster() +
      scale_fill_discrete(type = four_colors[c(1,3)],na.value=four_colors[4],name="Genotype",labels=c("0","1")) +
      theme_dartR() +
      theme(legend.position=posi,
            axis.text.y = element_text(size = ind_labels_size ))+
      scale_x_discrete(breaks= loc_labels, labels= as.character(loc_labels),name="Loci") +
      ylab("Individuals")
  }
  
  if(datatype=="SNP"){
    p3 <- ggplot(X,aes(x=locus,y=id,fill=genotype))+
      geom_raster() +
      scale_fill_discrete(type = four_colors,na.value=four_colors[4],name="Genotype",labels=c("0","1","2")) +
      theme_dartR() +
      theme(legend.position=posi,
            axis.text.y = element_text(size = ind_labels_size ))+
      scale_x_discrete(breaks= loc_labels, labels= as.character(loc_labels),name="Loci") +
      ylab("Individuals")
  }
    
    if(group_pop == TRUE){
      p3 <- p3 + facet_wrap(~pop, ncol=1,dir="v",scales="free_y")
    }

    
  # PRINTING OUTPUTS
    print(p3)

  # creating temp file names
  if(save2tmp){
    temp_plot <- tempfile(pattern = "Plot_")
    match_call <- paste0(names(match.call()),"_",as.character(match.call()),collapse = "_")
    # saving to tempdir
    saveRDS(list(match_call,p3), file = temp_plot)
    if(verbose>=2){
      cat(report("  Saving the ggplot to session tempfile\n"))
    }
  }

  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # RETURN
  
  invisible(x)
}
