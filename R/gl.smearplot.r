#' @name gl.smearplot
#' @title Smear plot of SNP or presence/absence (SilicoDArT) data
#' @description
#' Each locus is color coded for scores of 0, 1, 2 and NA for SNP data and 0, 1
#' and NA for presence/absence (SilicoDArT) data. Individual labels can be added
#' and individuals can be grouped by population.
#'
#' Plot may become cluttered if ind_labels If there are too many individuals, 
#' it is best to use ind_labels_size = 0.
#'
#' @param x Name of the genlight object containing the SNP or presence/absence
#'  (SilicoDArT) data [required].
#' @param ind_labels If TRUE, individuals are labelled with indNames(x) [default FALSE].
#' @param group_pop If ind_labels is TRUE, group by population [default TRUE].
#' @param ind_labels_size Size of the individual labels [default 10].
#' @param plot_colors Vector with four color names for homozygotes for the
#' reference allele, heterozygotes, homozygotes for the alternative allele and
#' for missing values (NA), e.g. four_colours [default NULL].
#' Can be set to "hetonly", which defines colors to only show heterozygotes in the genlight object
#' @param posi Position of the legend: “left”, “top”, “right”, “bottom” or
#'  'none' [default = 'bottom'].
#' @param save2tmp If TRUE, saves plot to the session temporary directory
#' (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report [default NULL].
#'
#' @return Returns unaltered genlight object
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' gl.smearplot(testset.gl,ind_labels=FALSE)
#' gl.smearplot(testset.gs[1:10,],ind_labels=TRUE)
#' @family Exploration/visualisation functions
#' @export
#'

gl.smearplot <- function(x,
                        ind_labels = FALSE,
                        group_pop = FALSE, 
                        ind_labels_size = 10,
                        plot_colors = NULL,
                        posi = "bottom",
                        save2tmp = FALSE,
                        verbose = NULL) {
    
    # CHECK IF PACKAGES ARE INSTALLED
    pkg <- "reshape2"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      cat(error(
        "Package",
        pkg,
        " needed for this function to work. Please install it.\n"
      ))
      return(-1)
    }
    
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose) 
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbosity = verbose)
    
    # SCRIPT SPECIFIC CHECKS
    if (length(plot_colors)==1) {
      if (plot_colors=="hetonly") plot_colors <- c("#dddddd", "#ff0000", "#dddddd","#dddddd" )
    }
    if(is.null(plot_colors)){
        plot_colors <- c("#a6cee3","#1f78b4","#b2df8a","#dddddd") # = default for plot()
        plot_colors <- c("#1b9e77","#d95f02","#7570b3","#dddddd") # = default for plot()
        
            }
    n10 <- nchar(as.character(nInd(x)))
    lzs <- paste0("%0",as.character(n10),"d")
    
    # Luis approach
    if(ind_labels == TRUE){
      individuals <- indNames(x)
    } else {
      individuals <- seq(1:length(indNames(x)))
    }
    
    # Bernd approach
    # if(ind_labels == TRUE){
    #   individuals <- paste0(sprintf(lzs,1:nInd(x)),"_",indNames(x))
    # } else {
    #     individuals <- paste0(sprintf(lzs,1:nInd(x)))
    # }
    
    # DO THE JOB
    
    X_temp <- as.data.frame(as.matrix(x))
    colnames(X_temp) <- 1:nLoc(x)
    X_temp$id <- individuals
    # converting id to factor using levels parameters
    X_temp$id <- factor(X_temp$id, levels = X_temp$id)
    
    X_temp$pop <- pop(x)
    
    X <- reshape2::melt(X_temp, id.vars = c("pop", "id"))
    X$value <- as.character(X$value)
    X$value <- ifelse(X$value=="NA", NA, X$value)
    colnames(X) <- c("pop", "id", "locus", "genotype")
    loc_labels <- pretty(1:nLoc(x), 5)
    id_labels <- pretty(1:nInd(x), 5)
    
    locus <- id <- genotype <- NA
    
    labels_genotype <- as.character(unique(X$genotype)) 
    labels_genotype[which(is.na(labels_genotype))] <- "Missing data"
    labels_genotype["0"] <- "Homozygote reference\n allele"
    labels_genotype["1"] <- "Heterozygote"
    labels_genotype["2"] <- "Homozygote alternative\n allele"

    
    if (datatype == "SilicoDArT") {
        p3 <-
            ggplot(X, aes(
                x = locus,
                y = id,
                fill = genotype
            )) + geom_raster() + scale_fill_discrete(
                type = plot_colors[c(1, 3)],
                na.value = plot_colors[4],
                name = "Genotype",
                labels = labels_genotype) +
          theme_dartR() + 
          theme(
                legend.position = posi,
                axis.text.y = element_text(size = ind_labels_size)
            ) +
            scale_x_discrete(
                breaks = loc_labels,
                labels = as.character(loc_labels),
                name = "Loci"
            ) + 
            ylab("Individuals")
    }
    
    if (datatype == "SNP") {
        p3 <-
            ggplot(X, aes(
                x = locus,
                y = id,
                fill = genotype
            )) + geom_raster() + 
                scale_fill_discrete(
                type = plot_colors,
                na.value = plot_colors[4],
                name = "Genotype",
                labels = labels_genotype) +
          theme_dartR() + theme(
                legend.position = posi,
                axis.text.y = element_text(size = ind_labels_size)
            ) +
            scale_x_discrete(
                breaks = loc_labels,
                labels = as.character(loc_labels),
                name = "Loci",
                position="bottom"
            ) + 
        ylab("Individuals")
    }
    
    if (ind_labels==TRUE & group_pop == TRUE) {
        p3 <- p3 + facet_wrap(~ pop,
                              ncol = 1,
                              dir = "v",
                              scales = "free")
    }
    
    # PRINTING OUTPUTS
    print(p3)
    
    # creating temp file names
    if (save2tmp) {
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
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    # RETURN
    
    return(p3)
}
