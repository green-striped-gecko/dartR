#' @name gl.report.ld.map
#' @title Calculates pairwise linkage disequilibrium by population
#' @description
#' This function calculates pairwise linkage disequilibrium (LD) by population 
#' using the function \code{\link[snpStats]{ld}} (package snpStats).
#' 
#' If SNPs are not mapped to a reference genome, the parameter
#'  \code{ld_max_pairwise}
#'  should be set as NULL (the default). In this case, the 
#' function will assign the same chromosome ("1") to all the SNPs in the dataset
#'  and assign a sequence from 1 to n loci as the position of each SNP. The 
#'  function will then calculate LD for all possible SNP pair combinations. 
#'  
#' If SNPs are mapped to a reference genome, the parameter 
#' \code{ld_max_pairwise}
#'  should be filled out (i.e. not NULL). In this case, the
#'  information for SNP's position should be stored in the genlight accessor
#'   "@@position" and the SNP's chromosome name in the accessor "@@chromosome"
#'    (see examples). The function will then calculate LD within each chromosome
#'     and for all possible SNP pair combinations within a distance of
#'      \code{ld_max_pairwise}. 
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param ld_max_pairwise Maximum distance in number of base pairs at which LD 
#' should be calculated [default NULL].
#' @param maf Minor allele frequency (by population) threshold to filter out 
#' loci. If a value > 1 is provided it will be interpreted as MAC (i.e. the
#'  minimum number of times an allele needs to be observed) [default 0.05].
#' @param ld_stat The LD measure to be calculated: "LLR", "OR", "Q", "Covar",
#'   "D.prime", "R.squared", and "R". See \code{\link[snpStats]{ld}}
#'    (package snpStats) for details [default "R.squared"].
#' @param ind.limit Minimum number of individuals that a population should
#' contain to take it in account to report loci in LD [default 10].
#' @param stat_keep Name of the column from the slot \code{loc.metrics} to be
#'  used to choose SNP to be kept [default "AvgPIC"].
#' @param ld_threshold_pops LD threshold to report in the plot of "Number of 
#' populations in which the same SNP pair are in LD" [default 0.2].
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot_theme User specified theme [default NULL].
#' @param histogram_colors Vector with two color names for the borders and fill
#' [default NULL].
#' @param boxplot_colors A color palette for box plots by population or a list
#'  with as many colors as there are populations in the dataset
#' [default NULL].
#' @param bins Number of bins to display in histograms [default 50].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session
#' temporary directory (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
#' @details
#' This function reports LD between SNP pairs by population. 
#' The function \code{\link{gl.filter.ld}} filters out the SNPs in LD using as
#' input the results of \code{\link{gl.report.ld.map}}. The actual number of 
#' SNPs to be filtered out depends on the parameters set in the function 
#' \code{\link{gl.filter.ld}}.
#' 
#' Boxplots of LD by population and
#' a histogram showing LD frequency are presented.
#'    
#' @return A dataframe with information for each SNP pair in LD. 
#' @author Custodian: Luis Mijangos -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \dontrun{
#' x <- platypus.gl
#' x$position <- x$other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1
#' x$chromosome <- as.factor(x$other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1)
#' ld_res <- gl.report.ld.map(x,ld_max_pairwise = 1000000)
#' }
#' @seealso \code{\link{gl.filter.ld}}
#' @family report functions
#' @export

gl.report.ld.map <- function(x,
                           ld_max_pairwise = NULL,
                           maf = 0.05,
                           ld_stat = "R.squared",
                           ind.limit = 10,
                           stat_keep = "AvgPIC",
                           ld_threshold_pops = 0.2,
                           plot.out = TRUE,
                           plot_theme = NULL,
                           histogram_colors = NULL,
                           boxplot_colors = NULL,
                           bins = 50,
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
  
  # check if packages are installed
  pkg <- "snpStats"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it."
    ))
  }
  pkg <- "fields"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it."
    ))
  }
  
  # DO THE JOB
  # by default SNPs are mapped to a reference genome
  SNP_map <- TRUE
  
  if(is.null(ld_max_pairwise)){
    x$position <- 1:nLoc(x)
    x$chromosome <- as.factor(rep("1",nLoc(x)))
    # SNPs are not mapped to a reference genome
    SNP_map <- FALSE
  }
  
  x_list <- seppop(x)
  
  df_linkage <- as.data.frame(matrix(nrow = 0, ncol = 11))
  colnames(df_linkage) <- c(
    "pop",
    "chr",
    "pos_loc_a",
    "pos_loc_b",
    "ld_stat",
    "distance",
    "locus_a.snp.name",
    "locus_a.stat_keep",
    "locus_b.snp.name",
    "locus_b.stat_keep",
    "locus_a_b"
  )
  
  for (i in 1:length(x_list)) {

    pop_ld <- x_list[[i]]
    pop_name <- popNames(pop_ld)
    
    if(nInd(pop_ld)<=ind.limit){
        cat(warn(paste("  Skipping population",pop_name,"from analysis because 
                       it has less than",ind.limit,"individuals.\n")))
      next()
    }
    
    if (verbose >= 2) {
      cat(report("  Calculating pairwise LD in population", pop_name, "\n"))
    }
    # ordering SNPs by chromosome and position
    hold <- pop_ld
    pop_ld <- hold[, order(hold$chromosome, hold$position)]
    pop_ld$other$loc.metrics <-
      hold$other$loc.metrics[order(hold$chromosome, hold$position),]
     pop_ld <- gl.recalc.metrics(pop_ld, verbose = 0)
    if (maf > 0) {
      pop_ld <- gl.filter.maf(pop_ld, threshold = maf, verbose = 0)
    }
    gl2plink(
      pop_ld,
      outfile = paste0("gl_plink", "_", pop_name),
      pos_cM = pop_ld$other$loc.metrics[, stat_keep],
      verbose = 0
    )
    
    # Read a pedfile as "SnpMatrix" object using a modified version of the 
    # function read.pedfile from package snpStats
    snp_stats <-
      utils.read.ped(
        file = paste0(tempdir(), "/", "gl_plink", "_", pop_name, ".ped"),
        snps = paste0(tempdir(), "/", "gl_plink", "_", pop_name, ".map") ,
        sep = " ",
        show_warnings = F,
        na.strings = NA
      )
    
    ld_map <- snp_stats$map
    colnames(ld_map) <-
      c("chr",
        "snp.name",
        "stat_keep",
        "loc_bp",
        "allele.1",
        "allele.2")
    ld_map$chr <- pop_ld$chromosome
    genotype <- snp_stats$genotypes
    colnames(genotype@.Data) <- ld_map$loc_bp
    
    chr_list <- as.character(unique(ld_map$chr))
    
    for (chrom in 1:length(chr_list)) {
      ld_loci <- which(ld_map$chr == chr_list[chrom])
      ld_map_loci <- ld_map[ld_loci,]
      genotype_loci <- genotype[, ld_loci]
      # removing loci that have the same location
      dupl_loci <- which(duplicated(ld_map_loci$loc_bp))
      if (length(dupl_loci) > 0) {
        ld_map_loci <- ld_map_loci[-dupl_loci, ]
        genotype_loci <- genotype_loci[, -dupl_loci]
      }
      if (nrow(ld_map_loci) <= 1) {
        next
      }
      
      #if SNPs are mapped to a reference genome
      if(SNP_map==TRUE){
        
        # this is the mean distance between each snp which is used to determine 
        # the depth at which LD analyses are performed
        mean_dis <- mean(diff(ld_map_loci$loc_bp))
        ld_depth_b <- ceiling((ld_max_pairwise / mean_dis)) - 1
        #function to calculate LD
        ld_snps <- snpStats::ld(genotype_loci, depth = ld_depth_b, 
                                stats = ld_stat)
        
        #if SNPs are not mapped to a reference genome 
      }else{
        #function to calculate LD
        ld_snps <- snpStats::ld(genotype_loci,genotype_loci, stats = ld_stat)
        ld_snps[lower.tri(ld_snps,diag = TRUE)] <- 0
        ld_max_pairwise <- nLoc(x)
        
      }
     
      ld_columns <- as.matrix(ld_snps)
      colnames(ld_columns) <- rownames(ld_columns)
      
      ld_columns <- as.data.frame(as.table(as.matrix(ld_columns)))
      # remove cases where LD was not calculated
      ld_columns <- ld_columns[-ld_columns$Freq < 0,] 
      ld_columns$Var1 <- as.numeric(as.character(ld_columns$Var1))
      ld_columns$Var2 <- as.numeric(as.character(ld_columns$Var2))
      # determine the distance at which LD was calculated
      ld_columns$dis <- ld_columns$Var2 - ld_columns$Var1
      # remove pairwise LD results that were calculated at larger distances than 
      # the required in the settings and then filtering and rearranging 
      # dataframes to match each other and then merge them
      df_linkage_temp <- ld_columns[which(ld_columns$dis <= ld_max_pairwise),]
      if (nrow(df_linkage_temp) < 1) {
        next
      }
      ldtb <- data.table(df_linkage_temp , key = "Var1")
      ldtc <- data.table(df_linkage_temp , key = "Var2")
      # this is the location of each snp in cM and in bp
      snp_loc <-
        ld_map_loci[, c("chr", "snp.name", "stat_keep", "loc_bp")]
      dtb <- data.table(snp_loc, key = "loc_bp")
      t_locationb <- ldtb[dtb, c("snp.name", "stat_keep"), nomatch = 0]
      t_locationc <- ldtc[dtb, c("snp.name", "stat_keep"), nomatch = 0]
      
      df_linkage_temp <- df_linkage_temp[order(df_linkage_temp$Var1),]
      df_linkage_temp <- cbind(df_linkage_temp, t_locationb)
      df_linkage_temp <- df_linkage_temp[order(df_linkage_temp$Var2),]
      df_linkage_temp <- cbind(df_linkage_temp, t_locationc)
      df_linkage_temp <- df_linkage_temp[order(df_linkage_temp$Var1),]
      df_linkage_temp <- cbind(chr_list[chrom], df_linkage_temp)
      df_linkage_temp <- cbind(pop_name, df_linkage_temp)
      
      colnames(df_linkage_temp) <- c(
        "pop",
        "chr",
        "pos_loc_a",
        "pos_loc_b",
        "ld_stat",
        "distance",
        "locus_a.snp.name",
        "locus_a.stat_keep",
        "locus_b.snp.name",
        "locus_b.stat_keep"
      )
      
      df_linkage_temp$locus_a_b <-
        paste0(df_linkage_temp$locus_a.snp.name,
               "_",
               df_linkage_temp$locus_b.snp.name)
      
      
      df_linkage <- rbind(df_linkage, df_linkage_temp)
    }
  }
  
  df_ld <- df_linkage
  
  if (plot.out) {
    
    if(is.null(histogram_colors)){
      histogram_colors <- two_colors
    }
    
    if(is.null(plot_theme)){
    plot_theme = theme_dartR()
    }
    
    if(is.null(boxplot_colors)){
      boxplot_colors <- discrete_palette(length(levels(pop(x))))
    }
    
    if (is(boxplot_colors, "function")) {
      boxplot_colors <- boxplot_colors(length(levels(pop(x))))
    }
    
    if (!is(boxplot_colors,"function")) {
      boxplot_colors <- boxplot_colors
    }
    
    # get title for plots
    title1 <- "SNP data - Pairwise LD"
    
    # Boxplot
    p1 <-
      ggplot(df_ld, aes(x=pop,y = ld_stat,color=pop)) +
      geom_boxplot() +
      plot_theme +
      scale_color_manual(values = boxplot_colors) +
      ylab(ld_stat) +
      ggtitle(title1) +
      theme(legend.position = "none")+     
      labs(title = "Pairwise LD by population", color = "") +
      theme(axis.title.x=element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) 
    
    # Histogram
    p2 <-
      ggplot(df_ld, aes(x = ld_stat)) +
      geom_histogram(bins = bins,
                     color = histogram_colors[1],
                     fill = histogram_colors[2]) +
      xlab(ld_stat) +
      ylab("Count") +
      plot_theme
    
    # Number of populations in which the same SNP pairs are in LD
    df_ld_temp <- df_ld[which(df_ld$ld_stat>ld_threshold_pops),]
    ld_pops_tmp <- table(rowSums(as.data.frame.matrix(with(df_ld_temp, 
                                                     table(locus_a_b,pop)))))
    ld_pops <- as.data.frame(cbind(as.numeric(names(ld_pops_tmp)),
                                   unlist(unname(ld_pops_tmp))))
    colnames(ld_pops) <- c("pops","n_loc")
    
    n_loc <- pops <- NULL
    
     p3 <- ggplot(ld_pops,aes(x=pops,y=n_loc))+
      geom_col(color = histogram_colors[1],fill = histogram_colors[2]) +
      ylab("Count")+
      xlab("Number of populations in which the same SNP pair are in LD") +
      plot_theme
  }
  
  # Print out some statistics
  # stats <- summary(ld_stat_res)
  # cat("  Reporting pairwise LD\n")
  # cat("  No. of pairs of loci in LD =", length(ld_stat_res), "\n")
  # cat("  No. of individuals =", nInd(x), "\n")
  # cat("    Minimum      : ", stats[1], "\n")
  # cat("    1st quartile : ", stats[2], "\n")
  # cat("    Median       : ", stats[3], "\n")
  # cat("    Mean         : ", stats[4], "\n")
  # cat("    3r quartile  : ", stats[5], "\n")
  # cat("    Maximum      : ", stats[6], "\n")
  # cat("    Missing Rate Overall: ", round(sum(is.na(as.matrix(
  #   x
  # ))) / (nLoc(x) * nInd(x)), 2), "\n\n")
  
  # PRINTING OUTPUTS
  if (plot.out) {
      # using package patchwork
      p4 <- p1 / p2 / p3
      print(p4)
  }

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
      saveRDS(list(match_call, p4), file = temp_plot)
      if (verbose >= 2) {
        cat(report("  Saving the ggplot to session tempfile\n"))
      }
    }
  }
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # RETURN
  return(invisible(df_ld))
  
}
