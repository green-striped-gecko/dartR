#' @name gl.ld.haplotype
#' @title Visualize patterns of linkage disequilibrium and identification of 
#' haplotypes
#' @description
#' This function plots a Linkage disequilibrium (LD) heatmap, where the colour 
#' shading indicates the strength of LD. Chromosome positions (Mbp) are shown on
#'  the horizontal axis, and haplotypes appear as triangles and delimited by 
#'  dark yellow vertical lines. Numbers identifying each haplotype are shown in 
#'  the upper part of the plot. 
#'  
#'  The heatmap also shows heterozygosity for each SNP. 
#'  
#'  The function identifies haplotypes based on contiguous SNPs that are in 
#'  linkage disequilibrium using as threshold \code{ld_threshold_haplo} and
#' containing more than \code{min_snps} SNPs.
#' 
#' @param x Name of the genlight object containing the SNP data [required].
#' @param pop_name Name of the population to analyse. If NULL all the 
#' populations are analised [default NULL].
#' @param chrom_name Nme of the chromosome to analyse. If NULL all the 
#' chromosomes are analised [default NULL].
#' @param ld_max_pairwise Maximum distance in number of base pairs at which LD
#' should be calculated [default 10000000].
#' @param maf Minor allele frequency (by population) threshold to filter out 
#' loci. If a value > 1 is provided it will be interpreted as MAC (i.e. the
#'  minimum number of times an allele needs to be observed) [default 0.05].
#' @param ld_stat The LD measure to be calculated: "LLR", "OR", "Q", "Covar",
#'   "D.prime", "R.squared", and "R". See \code{\link[snpStats]{ld}}
#'    (package snpStats) for details [default "R.squared"].
#' @param ind.limit Minimum number of individuals that a population should
#' contain to take it in account to report loci in LD [default 10].
#' @param min_snps Minimum number of SNPs that should have a haplotype to call 
#' it [default 10].
#' @param ld_threshold_haplo Minimum LD between adjacent SNPs to call a 
#' haplotype [default 0.5].
#' @param coordinates A vector of two elements with the start and end 
#' coordinates in base pairs to which restrict the 
#' analysis e.g. c(1,1000000) [default NULL]. 
#' @param color_haplo Color palette for haplotype plot. See details
#'  [default "viridis"].
#' @param color_het Color for heterozygosity [default "deeppink"].
#' @param plot.out Specify if heatmap plot is to be produced [default TRUE].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session
#' temporary directory (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
#' @details
#' The information for SNP's position should be stored in the genlight accessor
#'   "@@position" and the SNP's chromosome name in the accessor "@@chromosome"
#'   (see examples). The function will then calculate LD within each chromosome.
#'    
#' The output of the function includes a table with the haplotypes
#'  that were identified and their location.
#'  
#'  Colors of the heatmap (\code{color_haplo}) are based on the function
#'    \code{\link[viridis]{scale_fill_viridis}} from  package \code{viridis}. 
#'    Other color palettes options are "magma", "inferno", "plasma", "viridis",
#'     "cividis", "rocket", "mako" and "turbo".
#' @return A table with the haplotypes that were identified.
#' @family ld functions
#' @examples 
#' require("dartR.data")
#' x <- platypus.gl
#' x <- gl.filter.callrate(x,threshold = 1)
#' x <- gl.keep.pop(x, pop.list = "TENTERFIELD")
#' x$chromosome <- as.factor(x$other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1)
#' x$position <- x$other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1
#' ld_res <- gl.ld.haplotype(x,chrom_name = "NC_041728.1_chromosome_1",
#'                           ld_max_pairwise = 10000000 )
#'      
#' @author Custodian: Luis Mijangos -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @export

gl.ld.haplotype <- function(x,
                            pop_name = NULL,
                            chrom_name = NULL,
                            ld_max_pairwise = 10000000,
                            maf = 0.05,
                            ld_stat = "R.squared",
                            ind.limit = 10,
                            min_snps = 10,
                            ld_threshold_haplo = 0.5,
                            coordinates = NULL,
                            color_haplo = "viridis",
                            color_het = "deeppink",
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
  
  # check if packages are installed
  pkg <- "snpStats"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }
  
  pkg <- "fields"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }
  
  # DO THE JOB
  
  group <- het <- lat <- long <- NULL
  
  if (is.null(chrom_name) == FALSE) {
    chrom_tmp <-  unlist(lapply(chrom_name, function(y) {
      which(x$chromosome == y)
    }))
    
    chrom_tmp <- chrom_tmp[order(chrom_tmp)]
    
    x <- gl.keep.loc(x , loc.list = locNames(x)[chrom_tmp], verbose = 0)
    
  }
  
  if (is.null(pop_name) == FALSE) {
    x <- gl.keep.pop(x, pop.list = pop_name, verbose = 0)
  }
  
  x_list <- seppop(x)
  
  if(is.null(coordinates) == FALSE){
    cat(report("  Restricting the analysis from",coordinates[1],"to",
               coordinates[2],"base pairs\n"))
    x_list <- lapply(x_list,function(y){
      loc_names <- which(y$position >= coordinates[1] & 
                           y$position <= coordinates[2])
      y <- gl.keep.loc(y,loc.list = locNames(y)[loc_names],verbose = 0)
      return(y)
  })
  }
 
haplo_table <- as.data.frame(matrix(nrow = 1,ncol = 10))
colnames(haplo_table) <-  c("population","chromosome","haplotype","start","end",
                            "start_ld_plot","end_ld_plot","midpoint", 
                            "midpoint_ld_plot","labels")

chr_list <- as.character(unique(x$chromosome))

p <- NULL

# p <- rep(list(as.list(rep(NA, length(chr_list)))), length(names(x_list)))
# names(p) <- names(x_list)
# p <- lapply(p, function(x) {
#   names(x) <- paste0("chr_", chr_list)
#   return(x)
# })

  for (pop_n in 1:length(x_list)) {
    pop_ld <- x_list[[pop_n]]
    pop_name <- popNames(pop_ld)
    
    if (nInd(pop_ld) <= ind.limit) {
      cat(warn(
        paste(
          "  Skipping population",
          pop_name,
          "from analysis because it has less than",
          ind.limit,
          "individuals.\n"
        )
      ))
      next()
    }
    
    if (verbose >= 2) {
      cat(report("  Calculating pairwise LD in population", pop_name, "\n"))
    }
    # ordering SNPs by chromosome and position
    hold <- pop_ld
    pop_ld <- hold[, order(hold$chromosome, hold$position)]
    pop_ld$other$loc.metrics <-
      hold$other$loc.metrics[order(hold$chromosome, hold$position), ]
    pop_ld <- gl.recalc.metrics(pop_ld, verbose = 0)
    if (maf > 0) {
      pop_ld <- gl.filter.maf(pop_ld, threshold = maf, verbose = 0)
    }
    gl2plink(pop_ld,
             outfile = paste0("gl_plink", "_", pop_name),
             verbose = 0)
    
    # Read a pedfile as "SnpMatrix" object using a modified version of the
    # function read.pedfile from package snpStats
    snp_stats <-
      utils.read.ped(
        file = paste0(tempdir(), "/", "gl_plink", "_", pop_name, ".ped"),
        snps = paste0(tempdir(), "/", "gl_plink", "_", pop_name, ".map") ,
        sep = " ",
        show_warnings = FALSE,
        na.strings = NA
      )
    
    ld_map <- snp_stats$map
    colnames(ld_map) <-
      c("chr",
        "snp.name",
        "null",
        "loc_bp",
        "allele.1",
        "allele.2")
    ld_map$chr <- pop_ld$chromosome
    genotype <- snp_stats$genotypes
    colnames(genotype@.Data) <- ld_map$loc_bp
    
    for (chrom in 1:length(chr_list)) {
      chr_name <- chr_list[chrom]
      if (verbose >= 2) {
        cat(report("  Analysing chromosome", chr_name, "\n"))
      }
      ld_loci <- which(ld_map$chr == chr_list[chrom])
      ld_map_loci <- ld_map[ld_loci, ]
      genotype_loci <- genotype[, ld_loci]
      # removing loci that have the same location
      dupl_loci <- which(duplicated(ld_map_loci$loc_bp))
      if (length(dupl_loci) > 0) {
        ld_map_loci <- ld_map_loci[-dupl_loci,]
        genotype_loci <- genotype_loci[,-dupl_loci]
      }
      if (nrow(ld_map_loci) <= 1) {
        next
      }
      
      # this is the mean distance between each snp which is used to determine
      # the depth at which LD analyses are performed
      mean_dis <- mean(diff(ld_map_loci$loc_bp))
      ld_depth_b <- ceiling((ld_max_pairwise / mean_dis)) - 1
      
      if(ld_depth_b<5){
        cat(warn("  The maximum distance at which LD should be calculated 
                 (ld_max_pairwise) is too short for chromosome",chr_name,
                 ". Setting this distance to",round(mean_dis*5,0),"bp\n" ))
        ld_depth_b <- 5
      }
      #function to calculate LD
      ld_snps <- snpStats::ld(genotype_loci, depth = ld_depth_b,
                              stats = ld_stat)
      
      ld_matrix_2 <- ld_snps
      colnames(ld_matrix_2) <- rownames(ld_matrix_2)
      rownames(ld_matrix_2) <- 1:nrow(ld_matrix_2)
      colnames(ld_matrix_2) <- 1:ncol(ld_matrix_2)
      ld_columns_2 <- as.data.frame(as.table(as.matrix(ld_matrix_2)))
      #remove cases where LD was not calculated
      ld_columns_2 <- ld_columns_2[-ld_columns_2$Freq < 0, ]
      ld_columns_2$Var1 <- as.numeric(as.character(ld_columns_2$Var1))
      ld_columns_2$Var2 <- as.numeric(as.character(ld_columns_2$Var2))
      ld_columns_2 <- ld_columns_2[complete.cases(ld_columns_2), ]
      raster_haplo <- raster::rasterFromXYZ(ld_columns_2)
      polygon_haplo <-
        raster::rasterToPolygons(
          raster_haplo,
          fun = NULL,
          n = 4,
          na.rm = TRUE,
          digits = 12,
          dissolve = TRUE
        )
      polygon_haplo <- maptools::elide(polygon_haplo, rotate = 45)
      polygon_haplo$id <- rownames(as.data.frame(polygon_haplo))
      #this only has the coordinates
      polygon_haplo.pts <- fortify(polygon_haplo, polygon_haplo = "id")
      # add the attributes back
      polygon_haplo.df <-
        merge(polygon_haplo.pts,
              polygon_haplo,
              by = "id",
              type = 'left')
      width_poly <-  round(max(polygon_haplo.df$long), 0)
      height_poly <-  max(polygon_haplo.df$lat)
      
      reduce_factor <- nrow(ld_map_loci) / width_poly
      enlarge_factor <- width_poly / nrow(ld_map_loci)
      #this is to correct the cell size unit when the polygon figure is more or
      # less than 1000 units
      correction_factor <- (width_poly / 1000)
      
      ld_matrix_3 <- ld_snps
      
      ld_matrix_3 <- as.matrix(ld_matrix_3)
      dimnames(ld_matrix_3) <- NULL
      matrix_rotate <-
        rotate.matrix(x = ld_matrix_3,
                      angle = -45,
                      method = "simple")
      matrix_rotate_2 <- matrix_rotate
      matrix_rotate_2[matrix_rotate_2 == 0] <- NA
      matrix_rotate_3 <-
        t(zoo::na.locf(
          t(matrix_rotate_2),
          fromLast = FALSE,
          na.rm = TRUE
        ))
      
      means_col <- apply(matrix_rotate_3, 2, var, na.rm = TRUE)
      means_col[means_col == 0] <- NA
      means_col <- zoo::na.locf(means_col, fromLast = TRUE)
      means_col <- means_col * 2
      df_col <-
        as.data.frame(matrix(ncol = 2 , nrow = length(means_col)))
      df_col[, 1] <- 1:nrow(df_col)
      df_col[, 2] <- means_col
      df_col[, 2] <-
        scales::rescale(df_col[, 2], to = c(0, height_poly))
      means_col_2 <- means_col
      
      mean_column <- as.numeric(summary(means_col_2, na.rm = TRUE)[1:6])
      mean_column <-
        scales::rescale(mean_column, to = c(0, height_poly))
      
# as the matrix was rotated the position of the snps is not correct anymore.
# the following code reassigns the snp position based on the second row from
# bottom to top of the rotated matrix. the first two and the last snps are
# removed to take in account the snps that are not present in the second row
      second_row_temp <- which(!is.na(matrix_rotate_3[, 1])) - 1
      second_row <-  second_row_temp[length(second_row_temp)]
      second_row_2 <- matrix_rotate_2[second_row, ]
      second_row_3 <- second_row_2
      
      row_snp <- as.numeric(rownames(ld_snps))
      reassign_loc <- row_snp[3:length(row_snp)]
      reassign_loc <- reassign_loc[1:length(reassign_loc) - 1]
      
      element_reassign_loc <- 1
      for (loc in 1:length(second_row_2)) {
        if (is.na(second_row_2[loc])) {
          next
        } else{
          second_row_3[loc] <- reassign_loc[element_reassign_loc]
          element_reassign_loc <-  element_reassign_loc + 1
        }
      }
      
      # putting back the first two snps and the last that were removed
      second_row_3[1:2] <- row_snp[1:2]
      second_row_3[length(second_row_3)] <- row_snp[length(row_snp)]
      #filling the NAs
      second_row_4 <-
        zoo::na.locf(second_row_3, fromLast = TRUE, na.rm = FALSE)
      second_row_4 <-
        zoo::na.locf(second_row_4, fromLast = FALSE, na.rm = FALSE)
      
      second_row_ver_2 <- zoo::na.locf(second_row_2, fromLast = TRUE)
      first_row <- which(!is.na(matrix_rotate_3[, 1]))[1]
      first_row_2 <- matrix_rotate_3[first_row, ]
      first_row_2 <- round(first_row_2, 2)
      
      #this is SNP's heterozygosity to be calculated alone
      het_tmp <- utils.basic.stats(pop_ld)
      snp_het_alone <-
        data.frame(position = pop_ld$position,
                   het = unname(het_tmp$Hs))
      snp_het_alone$position <- 1:nrow(snp_het_alone)
      snp_het_alone[, 1] <-
        scales::rescale(snp_het_alone[, 1], to = c(0, width_poly))
      snp_het_alone[, 2] <-
        scales::rescale(snp_het_alone[, 2], to = c(0, height_poly))
      
      # identifying haplotypes
      if(any(first_row_2>ld_threshold_haplo)){
      haplo_loc_test <- first_row_2 >= ld_threshold_haplo
      haplo_loc_test <- c(haplo_loc_test, FALSE)
      start_haplo <- NULL
      end_haplo <- NULL
      for (i in 1:(length(haplo_loc_test) - 1)) {
        if (haplo_loc_test[i] == haplo_loc_test[i + 1]) {
          next()
        }
        
        if (haplo_loc_test[i] == TRUE) {
          end_haplo_temp <- i
          end_haplo <- c(end_haplo, end_haplo_temp)
        }
        if (haplo_loc_test[i] == FALSE) {
          start_haplo_temp <- i
          start_haplo <- c(start_haplo, start_haplo_temp)
        }
      }
      
      start_haplo <- c(1, start_haplo)
      
      start_haplo_2 <- second_row_4[start_haplo]
      end_haplo_2 <-  second_row_4[end_haplo]
      if(length(start_haplo_2)!=length(end_haplo_2)){
        start_haplo_2 <- start_haplo_2[-length(start_haplo_2)]
      }
      haplo_1_ver_2 <- as.data.frame(cbind(start_haplo_2, end_haplo_2))
      haplo_1_ver_2$size <- (haplo_1_ver_2[, 2] - haplo_1_ver_2[, 1])
      
      n_snps <- as.matrix(haplo_1_ver_2[, 1:2])
      n_snps[, 1] <- n_snps[, 1] + 1
      n_snps <- n_snps[which(n_snps[, 1] != n_snps[, 2]), ]
      
      if(!is.matrix(n_snps)){
      n_snps <- as.matrix(t(n_snps))
      }
      
      n_snps <- n_snps[!duplicated(n_snps[, 1]), ]
      
      if(!is.matrix(n_snps)){
        n_snps <- as.matrix(t(n_snps))
      }
      
      n_snps <- n_snps[!duplicated(n_snps[, 2]), ]
    
      
      df.4.cut <-
        as.data.frame(table(cut(row_snp, breaks = n_snps)), stringsAsFactors =
                        FALSE)
      df.4.cut <- df.4.cut[which(df.4.cut$Freq >= min_snps), ]
      if(nrow(df.4.cut)<1){
        cat(warn(" No haplotypes with more than ",min_snps,"were found. 
                 Try using a lower threshold.\n"))
        next()
      }
      df.4.cut_3 <- gsub("[][()]", "", df.4.cut$Var1, ",")
      df.4.cut_3 <- strsplit(df.4.cut_3, ",")
      df.4.cut_4 <- lapply(df.4.cut_3, as.numeric)
      df.4.cut_4 <- as.data.frame(plyr::laply(df.4.cut_4, rbind))
      df.4.cut_4[, 3] <- (df.4.cut_4[, 2] - df.4.cut_4[, 1])
      
      # this is to calculate the real distance in bp of the polygon figure of LD
      
      real_distance <- c(0, second_row_4)
      real_distance_2 <- diff(real_distance)
      real_distance_3 <- cumsum(real_distance_2)
      real_distance_4 <-
        as.data.frame(cbind(1:length(real_distance_3), real_distance_3))
      
      test_var <- unname(unlist(df.4.cut_4[, 1:2]))
      
      location_test <-
        lapply(test_var, findInterval, vec = as.numeric(paste(
          unlist(real_distance_4$real_distance_3)
        )))
      location_test_2 <-  unlist(location_test)
      test_var_2 <-
        as.data.frame(cbind(1:length(location_test_2), location_test_2))
      
      hap_blocks <-
        as.data.frame(cbind(
          paste0("CHR_", 1:nrow(df.4.cut_4)),
          rep(1, times = nrow(df.4.cut_4)),
          df.4.cut_4[, 1:3],
          df.4.cut$Freq
        ))
      colnames(hap_blocks) <-
        c("BLOCK", "CHR" ,  "BP1" ,  "BP2"   , "SIZE",  "NSNP")
      
      locations_temp <-
        as.data.frame(cbind(hap_blocks$BP1, hap_blocks$BP2))
      locations_temp_2 <- locations_temp
      colnames(locations_temp_2) <- c("start", "end")
      locations_temp_2$start_ld_plot <-
   unlist(lapply(locations_temp_2$start, findInterval, vec = as.numeric(paste(
          unlist(real_distance_4$real_distance_3)
        ))))
      locations_temp_2$end_ld_plot <-
     unlist(lapply(locations_temp_2$end, findInterval, vec = as.numeric(paste(
          unlist(real_distance_4$real_distance_3)
        ))))
      locations_temp_2$midpoint <-
        (locations_temp_2$start + locations_temp_2$end) / 2
      locations_temp_2$midpoint_ld_plot <-
        (locations_temp_2$start_ld_plot + locations_temp_2$end_ld_plot) / 2
      locations_temp_2$labels <-
        paste0(as.character(round(locations_temp_2$start /1000000, 0)), "-", 
               as.character(round(locations_temp_2$end / 1000000, 0)))
      
      haplo_temp_a <-
  locations_temp_2[which(as.numeric(row.names(locations_temp_2)) %% 2 == 1), ]
      haplo_temp_b <-
   locations_temp_2[which(as.numeric(row.names(locations_temp_2)) %% 2 == 0), ]
      
      ticks_breaks <-
        c(locations_temp_2$start_ld_plot,
          locations_temp_2$end_ld_plot)
      ticks_breaks <- ticks_breaks[order(ticks_breaks)]
      ticks_lab <- c(locations_temp_2$start, locations_temp_2$end)
      ticks_lab <- round(ticks_lab / 1000000, 0)
      ticks_lab <- as.character(ticks_lab[order(ticks_lab)])
      
      ticks_joint <- as.data.frame(cbind(ticks_breaks, ticks_lab))
      ticks_joint <- ticks_joint[!duplicated(ticks_joint$ticks_lab), ]
      ticks_joint$ticks_lab <- as.character(ticks_joint$ticks_lab)
      ticks_joint$ticks_breaks <-
        as.numeric(as.character(ticks_joint$ticks_breaks))
      
      colors_plot <-
        c("Heterozygosity" = color_het,
          "Haplotypes limits" = "lightgoldenrod3")
      labels_haplo <- as.character(1:nrow(locations_temp_2))
      
      # for an unknown reason, the name of the fill variable in the geom_polygon
      # changes between Freq and layer. So when there is a error in displaying
      # the graphic, the name of this variable has to be changed for it to work
      
      haplo_table_tmp <- rbind(haplo_temp_a,haplo_temp_b)
      haplo_table_tmp <- cbind(haplotype=as.numeric(rownames(haplo_table_tmp)),
                               haplo_table_tmp)
      haplo_table_tmp <- haplo_table_tmp[order(haplo_table_tmp$haplotype),]
      haplo_table_tmp <- cbind(chromosome=chr_name,haplo_table_tmp)
      haplo_table_tmp <- cbind(population=pop_name,haplo_table_tmp)
      
      haplo_table <- rbind(haplo_table,haplo_table_tmp)
      
      p_temp <- NULL
      
      p_temp <- ggplot() +
        geom_rect(
          aes(
            xmin = haplo_temp_a$start_ld_plot,
            xmax = haplo_temp_a$end_ld_plot,
            ymin = min(polygon_haplo.df$lat) - 30,
            ymax = max(polygon_haplo.df$lat) + 30
          ),
          color = "cornsilk3",
          fill = "cornsilk3"
        ) +
        geom_rect(
          aes(
            xmin = haplo_temp_b$start_ld_plot,
            xmax = haplo_temp_b$end_ld_plot,
            ymin = min(polygon_haplo.df$lat) - 30,
            ymax = max(polygon_haplo.df$lat) + 30
          ),
          color = "cornsilk4",
          fill = "cornsilk4"
        ) +
        geom_polygon(data = polygon_haplo.df,
                     aes(long, lat, group = group, 
                         fill = polygon_haplo.df[, 8])) +
        viridis::scale_fill_viridis(name = ld_stat, option = color_haplo) +
        geom_vline(aes(
          xintercept = c(
            haplo_temp_b$start_ld_plot,
            haplo_temp_b$end_ld_plot,
            haplo_temp_a$start_ld_plot,
            haplo_temp_a$end_ld_plo
          ),
          color = "Haplotypes limits"
        ), size = 1) +
        
        geom_line(
          data = snp_het_alone,
          aes(x = position, y = het, color = "Heterozygosity"),
          inherit.aes = FALSE,
          size = 1 / 2,
          alpha = 1
        ) +
        annotate(
          "text",
          x = locations_temp_2$midpoint_ld_plot,
          y = max(polygon_haplo.df$lat) + 13,
          label = labels_haplo ,
          size = 3,
          color = "black"
        ) +
        annotate(
          "text",
          x = locations_temp_2$midpoint_ld_plot,
          y = min(polygon_haplo.df$lat) - 13,
          label = locations_temp_2$labels ,
          size = 3,
          color = "black"
        ) +
        labs(
          x = "Chromosome location (Mbp)",
          y = "Het",
          title = paste("Population",pop_name,"Chromosome", chr_name, "-", nLoc(pop_ld), "SNPs")
          ) +
        scale_x_continuous(breaks = ticks_joint$ticks_breaks,
                           labels = ticks_joint$ticks_lab) +
        scale_colour_manual(name = "", values = colors_plot) +
        theme_void() +
        theme(
          legend.position = "top",
          legend.text = element_text(size = 10),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_text(hjust = 0.5),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank()
        ) +
        coord_fixed(ratio = 1 / 1)
      
      # p <- c(p,p_temp) 
      
      # p[[pop_n]][[chrom]] <- p_temp
      
      # # PRINTING OUTPUTS
      if (plot.out) {
        print(p_temp)
      }
      
      }else{
    cat(warn("  No haplotypes were identified for chromosome",chr_name,"\n"))
        
        colors_plot <- c("Heterozygosity" = color_het)
        
        p_temp <- NULL
        
        p_temp <- ggplot() +
          geom_polygon(data = polygon_haplo.df,
                       aes(long, lat, group = group, 
                           fill = polygon_haplo.df[, 8])) +
          viridis::scale_fill_viridis(name = ld_stat, option = color_haplo) +
          geom_line(
            data = snp_het_alone,
            aes(x = position, y = het, color = "Heterozygosity"),
            inherit.aes = FALSE,
            size = 1 / 2,
            alpha = 1
          ) +
          labs(
            x = "Chromosome location (Mbp)",
            y = "Het",
            title = paste("Population",pop_name,"Chromosome", chr_name, "-", nLoc(pop_ld), "SNPs")
          ) +
          # scale_x_continuous(breaks = ticks_joint$ticks_breaks,
          #                    labels = ticks_joint$ticks_lab) +
          scale_colour_manual(name = "", values = colors_plot) +
          theme_void() +
          theme(
            legend.position = "top",
            legend.text = element_text(size = 10),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_text(hjust = 0.5),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
          ) +
          coord_fixed(ratio = 1 / 1)
        
        # p <- c(p,p_temp) 
        
        # p[[pop_n]][[chrom]] <- p_temp
        
        # PRINTING OUTPUTS
        if (plot.out) {
          print(p_temp)
        }
        
      }
    }
  }

# # PRINTING OUTPUTS
if (plot.out) {
  print(p)
}

haplo_table <- haplo_table[-1,]
print(haplo_table,row.names = FALSE)

# SAVE INTERMEDIATES TO TEMPDIR creating temp file names
if (save2tmp) {
  if (plot.out) {
    temp_plot <- tempfile(pattern = "Plot_")
    match_call <-
      paste0(names(match.call()),
             "_",
             as.character(match.call()),
             collapse = "_")
    # saving to tempdir
    saveRDS(list(match_call, p), file = temp_plot)
    if (verbose >= 2) {
      cat(report("  Saving the ggplot to session tempfile\n"))
    }
  }
  temp_table <- tempfile(pattern = "Table_")
  saveRDS(list(match_call, haplo_table), file = temp_table)
  if (verbose >= 2) {
    cat(report("  Saving tabulation to session tempfile\n"))
    cat(
      report(
        "  NOTE: Retrieve output files from tempdir using
                    gl.list.reports() and gl.print.reports()\n"
      )
    )
  }
}

# FLAG SCRIPT END

if (verbose >= 1) {
  cat(report("Completed:", funname, "\n"))
}

# RETURN

return(invisible(haplo_table))
  
}
