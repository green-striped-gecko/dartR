#' @name gl.ld.haplotype
#' @title Provides descriptive stats and plots to diagnose linkage
#' disequilibrium patterns
#' @description
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
#' @param save2tmp If TRUE, saves any ggplots and listings to the session
#' temporary directory (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
#' @details
#' This function reports LD between SNP pairs by population. Reports include
#' the number of SNPs that would be filtered out given specific LD thresholds.
#' The function code{\link{gl.filter.ld}} filters out the SNPs in LD using as
#' input the results of gl.report.ld.map. However, the actual number of SNPs to
#' be filtered out depends on the parameters set in the function
#' code{\link{gl.filter.ld}}. Therefore, the number of SNPs to be filtered out
#'  reported by gl.report.ld.map should be used as a guide and not as a definite
#'   number.
#'
#' \enumerate{
#' \item
#'
#' }
#' A bar plot with observed
#' and expected (null expectation) number of significant HWE tests for the same
#' locus in multiple populations (that is, the x-axis shows whether a locus
#' results significant in 1, 2, ..., n populations. The y axis is the count of
#' these occurrences.
#'
#' If SNPs are mapped to a reference genome, the function creates a plot showing
#' the pairwise LD measure against distance in number of base pairs pooled over
#' all the chromosomes and a red line representing the threshold (R.squared =
#' 0.2) that is commonly used to imply that two loci are unlinked (Delourme et
#' al., 2013; Li et al., 2014). Additionally, boxplots of LD by population and
#' a histogram showing LD frequency are presented.
#'
#' If SNPs are not mapped to a reference genome, only boxplots by population and
#'  a histogram showing LD frequency are presented.
#'
#' @references
#' \itemize{
#' \item Delourme, R., Falentin, C., Fomeju, B. F., Boillot, M., Lassalle, G.,
#' André, I., . . . Marty, A. (2013). High-density SNP-based genetic map
#' development and linkage disequilibrium assessment in Brassica napusL. BMC
#' genomics, 14(1), 120.
#' \item Li, X., Han, Y., Wei, Y., Acharya, A., Farmer, A. D., Ho, J., . . .
#' Brummer, E. C. (2014). Development of an alfalfa SNP array and its use to
#' evaluate patterns of population structure and linkage disequilibrium. PLoS
#' One, 9(1), e84329.
#'  }
#' @return A dataframe with information for each SNP pair in LD.
#' @author Custodian: Luis Mijangos -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \dontrun{
#' x <- platypus.gl
#' x$position <- x$other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1
#' x$chromosome <- x$other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1
#' ld_res <- gl.report.ld.map(x,ld_max_pairwise = 1000000)
#' }
#' @seealso \code{\link{gl.filter.ld}}
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
  
  for (i in 1:length(x_list)) {
    pop_ld <- x_list[[i]]
    pop_name <- popNames(pop_ld)
    
    if (nInd(pop_ld) <= ind.limit) {
      cat(warn(
        paste(
          "  Skipping population",
          pop_name,
          "from analysis because
                       it has less than",
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
             # pos_cM = pop_ld$other$loc.metrics[, stat_keep],
             verbose = 0)
    
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
        "null",
        "loc_bp",
        "allele.1",
        "allele.2")
    ld_map$chr <- pop_ld$chromosome
    genotype <- snp_stats$genotypes
    colnames(genotype@.Data) <- ld_map$loc_bp
    
    chr_list <- as.character(unique(ld_map$chr))
    
    for (chrom in 1:length(chr_list)) {
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
          dissolve = T
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
          fromLast = F,
          na.rm = T
        ))
      
      means_col <- apply(matrix_rotate_3, 2, var, na.rm = T)
      means_col[means_col == 0] <- NA
      means_col <- zoo::na.locf(means_col, fromLast = T)
      means_col <- means_col * 2
      df_col <-
        as.data.frame(matrix(ncol = 2 , nrow = length(means_col)))
      df_col[, 1] <- 1:nrow(df_col)
      df_col[, 2] <- means_col
      df_col[, 2] <-
        scales::rescale(df_col[, 2], to = c(0, height_poly))
      means_col_2 <- means_col
      
      mean_column <- as.numeric(summary(means_col_2, na.rm = T)[1:6])
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
      
      ld_threshold <- ld_threshold_haplo
      second_row_ver_2 <- zoo::na.locf(second_row_2, fromLast = T)
      first_row <- which(!is.na(matrix_rotate_3[, 1]))[1]
      first_row_2 <- matrix_rotate_3[first_row, ]
      first_row_2 <- round(first_row_2, 1)
      haplo_loc_test <- first_row_2 >= 1
      haplo_loc_test <- c(haplo_loc_test, F)
      start_haplo <- NULL
      end_haplo <- NULL
      for (i in 1:(length(haplo_loc_test) - 1)) {
        if (haplo_loc_test[i] == haplo_loc_test[i + 1]) {
          next()
        }
        
        if (haplo_loc_test[i] == T) {
          end_haplo_temp <- i
          end_haplo <- c(end_haplo, end_haplo_temp)
        }
        if (haplo_loc_test[i] == F) {
          start_haplo_temp <- i
          start_haplo <- c(start_haplo, start_haplo_temp)
        }
      }
      
      start_haplo <- c(1, start_haplo)
      
      start_haplo_2 <- second_row_4[start_haplo]
      start_haplo_2 <- start_haplo_2[-length(start_haplo_2)]
      end_haplo_2 <-  second_row_4[end_haplo]
      haplo_1_ver_2 <- as.data.frame(cbind(start_haplo_2, end_haplo_2))
      haplo_1_ver_2$size <- (haplo_1_ver_2[, 2] - haplo_1_ver_2[, 1])
      
      n_snps <- as.matrix(haplo_1_ver_2[, 1:2])
      n_snps[, 1] <- n_snps[, 1] + 1
      n_snps <- n_snps[which(n_snps[, 1] != n_snps[, 2]), ]
      n_snps <- n_snps[!duplicated(n_snps[, 1]), ]
      n_snps <- n_snps[!duplicated(n_snps[, 2]), ]
      
      df.4.cut <-
        as.data.frame(table(cut(row_snp, breaks = n_snps)), stringsAsFactors =
                        F)
      df.4.cut <- df.4.cut[which(df.4.cut$Freq >= min_snps), ]
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
      
      colors_haploview <-
        c("Heterozygosity" = "deeppink",
          "Haplotypes limits" = "lightgoldenrod3")
      labels_haplo <- as.character(1:nrow(locations_temp_2))
      
      # for an unknown reason, the name of the fill variable in the geom_polygon
      # changes between Freq and layer. So when there is a error in displaying
      # the graphic, the name of this variable has to be changed for it to work
      
      haploview <- ggplot() +
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
                aes(long, lat, group = group, fill = polygon_haplo.df[, 8])) +
        viridis::scale_fill_viridis(name = "r2 (LD)", option = "viridis") +
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
          inherit.aes = F,
          size = 1 / 5,
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
          title = paste("Chromosome", chrom, "-", nLoc(pop_ld), "SNPs")
        ) +
        scale_x_continuous(breaks = ticks_joint$ticks_breaks,
                           labels = ticks_joint$ticks_lab) +
        scale_colour_manual(name = "", values = colors_haploview) +
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
      
      
      print(haploview)
      
    }
  }
  
}
