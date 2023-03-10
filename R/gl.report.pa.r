#' @name gl.report.pa
#' @title Reports private alleles (and fixed alleles) per pair of populations
#' @description
#' This function reports private alleles in one population compared with a
#' second population, for all populations taken pairwise. It also reports a
#' count of fixed allelic differences and the mean absolute allele frequency
#' differences (AFD) between pairs of populations.
#' @param x Name of the genlight object containing the SNP data [required].
#' @param x2 If two separate genlight objects are to be compared this can be
#' provided here, but they must have the same number of SNPs [default NULL].
#' @param method Method to calculate private alleles: 'pairwise' comparison or
#' compare each population against the rest 'one2rest' [default 'pairwise'].
#' @param loc_names Whether names of loci with private alleles and fixed 
#' differences should reported. If TRUE, loci names are reported using a list
#'  [default FALSE].
#' @param plot.out Specify if Sankey plot is to be produced [default TRUE].
#' @param font_plot Numeric font size in pixels for the node text labels
#' [default 14].
#' @param map.interactive Specify whether an interactive map showing private
#' alleles between populations is to be produced [default FALSE].
#' @param palette_discrete A discrete palette for the color of populations or a
#' list with as many colors as there are populations in the dataset
#'  [default discrete_palette].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session
#' temporary directory (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent, fatal errors only; 1, flag function
#' begin and end; 2, progress log; 3, progress and results summary; 5, full
#' report [default 2 or as specified using gl.set.verbosity].
#' @details
#' Note that the number of paired alleles between two populations is not a
#' symmetric dissimilarity measure.
#'
#' If no x2 is provided, the function uses the pop(gl) hierarchy to determine
#' pairs of populations, otherwise it runs a single comparison between x and
#' x2.
#'
#'\strong{Hint:} in case you want to run comparisons between individuals
#'(assuming individual names are unique), you can simply redefine your
#'population names with your individual names, as below:
#'
#' \code{pop(gl) <- indNames(gl)}
#'
#'\strong{ Definition of fixed and private alleles }
#'
#' The table below shows the possible cases of allele frequencies between
#' two populations (0 = homozygote for Allele 1, x = both Alleles are present,
#' 1 = homozygote for Allele 2).
#'
#'\itemize{
#'\item p: cases where there is a private allele in pop1 compared to pop2 (but
#' not vice versa)
#'\item f: cases where there is a fixed allele in pop1 (and pop2, as those cases
#'are symmetric)
#'}
#'
#'\tabular{ccccc}{
#'\tab \tab \tab \emph{pop1} \tab \cr
#'\tab \tab \strong{0} \tab   \strong{x}  \tab  \strong{1}\cr
#'\tab     \strong{0}\tab -  \tab  p \tab  p,f\cr
#'\emph{pop2} \tab \strong{x}\tab -  \tab- \tab -\cr
#'\tab \strong{1} \tab p,f\tab p \tab   -\cr
#' }
#'
#' The absolute allele frequency difference (AFD) in this function is a simple
#' differentiation metric displaying intuitive properties which provides a
#' valuable alternative to FST. For details about its properties and how it is
#' calculated see Berner (2019).
#' 
#' The function also reports an estimation of the lower bound of the number of
#'  undetected private alleles using the Good-Turing frequency formula,
#'  originally developed for cryptography, which estimates in an ecological 
#'  context the true frequencies of rare species in a single assemblage based on
#'   an incomplete sample of individuals. The approach is described in Chao et 
#'   al. (2017). For this function, the equation 2c is used. This estimate is 
#'   reported in the output table as Chao1 and Chao2. 
#'
#' In this function a Sankey Diagram is used to visualize patterns of private
#' alleles between populations. This diagram allows to display flows (private
#' alleles) between nodes (populations). Their links are represented with arcs
#' that have a width proportional to the importance of the flow (number of
#' private alleles).
#'
#' if save2temp=TRUE, resultant plot(s) and the tabulation(s) are saved to the
#' session's temporary directory.
#'
#' @return A data.frame. Each row shows, for each pair of populations the number
#'  of individuals in each population, the number of loci with fixed differences
#'  (same for both populations) in pop1 (compared to pop2) and vice versa. Same
#'  for private alleles and finally the absolute mean allele frequency
#'  difference between loci (AFD). If loc_names = TRUE, loci names with private
#'   alleles and fixed differences are reported in a list in addition to the 
#'   dataframe. 
#' @author Custodian: Bernd Gruber -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @references \itemize{
#' \item Berner, D. (2019). Allele frequency difference AFD – an intuitive
#' alternative to FST for quantifying genetic population differentiation. Genes,
#'  10(4), 308.
#'  \item Chao, Anne, et al. "Deciphering the enigma of undetected species,
#'  phylogenetic, and functional diversity based on Good-Turing theory." 
#'  Ecology 98.11 (2017): 2914-2929.
#' }
#' @examples
#' out <- gl.report.pa(platypus.gl)
#' @seealso \code{\link{gl.list.reports}},
#'  \code{\link{gl.print.reports}}
#' @family report functions
#' @importFrom tidyr pivot_longer
#' @export

gl.report.pa <- function(x,
                         x2 = NULL,
                         method = "pairwise",
                         loc_names = FALSE,
                         plot.out = TRUE,
                         font_plot = 14,
                         map.interactive = FALSE,
                         palette_discrete = discrete_palette,
                         save2tmp = FALSE,
                         verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jackson",
                   verbosity = verbose)
  
  # CHECK DATATYPE
  datatype1 <- utils.check.datatype(x, verbose = verbose)
  if (!is.null(x2)) {
    datatype2 <- utils.check.datatype(x2, verbose = verbose)
  }
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
  # check if package is installed
  pkg <- "tibble"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }
  
  # check if package is installed
  pkg <- "networkD3"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }
  
  pkg <- "tidyverse"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }
  
  if (!is.null(x2)) {
    pops <- list(pop1 = x, pop2 = x2)
    x <- rbind(x, x2)
  } else {
    if (length(unique(pop(x))) > 1) {
      pops <- seppop(x)
    } else {
      stop(
        error(
          "Only one population provided. Check the @pop slot in your
                    genlight object.\n "
        )
      )
    }
  }
  
  # DO THE JOB For method 'pairwise'
  if (method == "pairwise") {
    pc <- t(combn(length(pops), 2))
    pall <-
      data.frame(
        p1 = pc[, 1],
        p2 = pc[, 2],
        pop1 = names(pops)[pc[, 1]],
        pop2 = names(pops)[pc[, 2]],
        N1 = NA,
        N2 = NA,
        fixed = NA,
        priv1 = NA,
        priv2 = NA,
        Chao1 = NA,
        Chao2= NA,
        totalpriv = NA,
        AFD = NA
      )
    
    pall_loc_names <- rep(list(as.list(rep(NA, 3))), nrow(pc))
    
    names(pall_loc_names) <- paste0(names(pops)[pc[, 1]],
                                    "_",
                                    names(pops)[pc[, 2]])
    
    pall_loc_names <- lapply(pall_loc_names,function(x){
      names(x) <- c("pop1_pop2_pa","pop2_pop1_pa","fd")
      return(x)
    })
    

    for (i in 1:nrow(pc)) {
      i1 <- pall[i, 1]
      i2 <- pall[i, 2]
      
      p1 <- as.matrix(pops[[i1]])
      p2 <- as.matrix(pops[[i2]])
      p1alf <- colMeans(p1, na.rm = TRUE) / 2
      p2alf <- colMeans(p2, na.rm = TRUE) / 2
      
      pall[i, c("N1","N2")] <- c(nrow(p1), nrow(p2))
      # pall[i, "fixed"] <- sum(abs(p1alf - p2alf) == 1, na.rm = T)
      
      # pall[i, "priv1"] <- sum(p2alf == 0 &
      #                     p1alf != 0, na.rm = T) + sum(p2alf == 1 &
      #                                                 p1alf != 1, na.rm = T)
      # pall[i, "priv2"] <- sum(p1alf == 0 &
      #                     p2alf != 0, na.rm = T) + sum(p1alf == 1 &
      #                                                  p2alf != 1, na.rm = T)
      
      # Changing code to get the names of the SNP with fixed differences and 
      # private alleles
      pop1_pop2_pa <- unique(unname(unlist(c(which(p2alf == 0 & p1alf != 0),
                                             which(p2alf == 1 & p1alf != 1)))))
      pop1_pop2_pa_loc_names <- locNames(x)[pop1_pop2_pa]
      pop2_pop1_pa <- unique(unname(unlist(c(which(p1alf == 0 & p2alf != 0),
                                             which(p1alf == 1 & p2alf != 1)))))
      pop2_pop1_pa_loc_names <- locNames(x)[pop2_pop1_pa]
      pop1_pop2_fd <- unique(unname(unlist(which(abs(p1alf - p2alf) == 1))))
      pop1_pop2_fd_loc_names <- locNames(x)[pop1_pop2_fd]
      
      pall[i, "fixed"] <- length(pop1_pop2_fd)
      pall[i, "priv1"] <- length(pop1_pop2_pa)
      pall[i, "priv2"] <- length(pop2_pop1_pa)
      
      pall_loc_names[[i]][["pop1_pop2_pa"]] <- pop1_pop2_pa_loc_names
      pall_loc_names[[i]][["pop2_pop1_pa"]] <- pop2_pop1_pa_loc_names
      pall_loc_names[[i]][["fd"]] <- pop1_pop2_fd_loc_names
      
      pall[i, "totalpriv"] <- pall[i, 8] + pall[i, 9]
      pall[i, "AFD"] <- round(mean(abs(p1alf - p2alf), na.rm = TRUE), 3)
      
      pa_Chao <- utils.pa.Chao(x=x,pop1_m=pops[[i1]],pop2_m=pops[[i2]])
      pall[i,"Chao1"] <- round(pa_Chao[[1]],0)
      pall[i,"Chao2"] <- round(pa_Chao[[2]],0)
      
    }
    
    if (plot.out) {
      mm <- matrix(0, nPop(x), nPop(x))
      
      for (i in 1:nrow(pall)){
        mm[pall[i, 1], pall[i, 2]] <- pall$priv2[i]
      }
      
      for (i in 1:nrow(pall)){
        mm[pall[i, 2], pall[i, 1]] <- pall$priv1[i]
      }
      
      colnames(mm) <- popNames(x)
      rownames(mm) <- popNames(x)
      
      data <- as.data.frame(mm)
      value <- target <- name <- NULL
      data_long <- tibble::rownames_to_column(data, "target")
      data_long <- tibble::as_tibble(data_long)
      data_long <- tidyr::pivot_longer(data_long,
                                       cols=-target,
                                       names_to = "source")
      data_long <- data_long[data_long$value > 0,]
      
      data_long$target <- gsub("\\.", " ", data_long$target)
      data_long$source <- paste0("src_", data_long$source)
      data_long$target <- paste0("trgt_", data_long$target)
      
      nodes <-
        data.frame(name = unique(c(
          data_long$source, data_long$target
        )),
        stringsAsFactors = FALSE)
      nodes <-
        tibble::tibble(name = unique(c(
          data_long$source, data_long$target
        )),
        target = grepl("trgt_", name))
      
      data_long$IDsource <-
        match(data_long$source, nodes$name) - 1
      data_long$IDtarget <-
        match(data_long$target, nodes$name) - 1
      
      nodes$name <- gsub("src_", "", nodes$name)
      nodes$name <- gsub("trgt_", "", nodes$name)
      
      # assigning colors to populations
      if (is(palette_discrete, "function")) {
        colors_pops <- palette_discrete(length(levels(pop(x))))
      }
      
      if (!is(palette_discrete, "function")) {
        colors_pops <- palette_discrete
        if (!any(grepl("#", colors_pops))) {
          colors_pops <- gplots::col2hex(colors_pops)
        }
      }
      
      colors_pops <- paste0("\"", paste0(colors_pops, collapse = "\",\""), "\"")
      
      colorScal <- paste("d3.scaleOrdinal().range([", colors_pops, "])")
      # color links
      data_long$color <-  gsub("src_", "", data_long$source)

      
      p3 <-
        suppressMessages(
          networkD3::sankeyNetwork(
            Links = data_long,
            Nodes = nodes,
            Source = "IDsource",
            Target = "IDtarget",
            LinkGroup = "color",
            Value = "value",
            NodeID = "name",
            sinksRight = FALSE,
            units = "Private alleles",
            colourScale = colorScal,
            nodeWidth = 40,
            fontSize = font_plot,
            nodePadding = 10
          )
        )
    }
  }
  
  # For method 'one2rest'
  
  if (method == "one2rest") {
    pas <- as.data.frame(matrix(nrow = nPop(x), ncol = 11))
    colnames(pas) <-
      c(
        "p1",
        "p2",
        "pop1",
        "pop2",
        "N1",
        "N2",
        "fixed",
        "priv1",
        "priv2",
        "totalpriv",
        "AFD"
      )
    
    pall_loc_names <- rep(list(as.list(rep(NA, 3))), nPop(x))
    
    names(pall_loc_names) <- paste0(popNames(x),
                                    "_",
                                    "Rest")
    
    pall_loc_names <- lapply(pall_loc_names,function(x){
      names(x) <- c("pop1_pop2_pa","pop2_pop1_pa","fd")
      return(x)
    })
    
    for (y in 1:nPop(x)) {
      gl2 <- x
      pop(gl2) <-
        factor(ifelse(pop(gl2) == popNames(x)[y], popNames(x)[y], "zzest"))
      pops <- seppop(gl2)
      pc <- t(combn(length(pops), 2))
      pall <-
        data.frame(
          p1 = pc[, 1],
          p2 = pc[, 2],
          pop1 = names(pops)[pc[, 1]],
          pop2 = names(pops)[pc[, 2]],
          N1 = NA,
          N2 = NA,
          fixed = NA,
          priv1 = NA,
          priv2 = NA,
          totalpriv = NA,
          AFD = NA
        )
      
      for (i in 1:nrow(pc)) {
        i1 <- pall[i, 1]
        i2 <- pall[i, 2]
        
        p1 <- as.matrix(pops[[i1]])
        p2 <- as.matrix(pops[[i2]])
        p1alf <- colMeans(p1, na.rm = T) / 2
        p2alf <- colMeans(p2, na.rm = T) / 2
        
        pop1_pop2_pa <- unique(unname(unlist(c(which(p2alf == 0 & p1alf != 0),
                                               which(p2alf == 1 & p1alf != 1)))))
        pop1_pop2_pa_loc_names <- locNames(x)[pop1_pop2_pa]
        pop2_pop1_pa <- unique(unname(unlist(c(which(p1alf == 0 & p2alf != 0),
                                               which(p1alf == 1 & p2alf != 1)))))
        pop2_pop1_pa_loc_names <- locNames(x)[pop2_pop1_pa]
        pop1_pop2_fd <- unique(unname(unlist(which(abs(p1alf - p2alf) == 1))))
        pop1_pop2_fd_loc_names <- locNames(x)[pop1_pop2_fd]
        
        pall[i, "fixed"] <- length(pop1_pop2_fd)
        pall[i, "priv1"] <- length(pop1_pop2_pa)
        pall[i, "priv2"] <- length(pop2_pop1_pa)
        
        pall_loc_names[[y]][["pop1_pop2_pa"]] <- pop1_pop2_pa_loc_names
        pall_loc_names[[y]][["pop2_pop1_pa"]] <- pop2_pop1_pa_loc_names
        pall_loc_names[[y]][["fd"]] <- pop1_pop2_fd_loc_names
        
        pall[i, 5:6] <- c(nrow(p1), nrow(p2))
        # pall[i, 7] <- sum(abs(p1alf - p2alf) == 1, na.rm = T)
        # 
        # pall[i, 8] <- sum(p2alf == 0 &
        #                     p1alf != 0, na.rm = T) + 
        #   sum(p2alf == 1 & p1alf != 1, na.rm = T)
        # 
        # pall[i, 9] <- sum(p1alf == 0 &
        #                     p2alf != 0, na.rm = T) + 
        #   sum(p1alf == 1 & p2alf != 1, na.rm = T)
        
        pall[i, 10] <- pall[i, 8] + pall[i, 9]
        pall[i, 11] <-
          round(mean(abs(p1alf - p2alf), na.rm = T), 3)
      }
      
      pas[y,] <- pall
    }
    
    pall <- pas
    pall$pop2 <- "Rest"
    
    if (plot.out) {
      # assigning colors to populations
      if (is(palette_discrete, "function")) {
        colors_pops <- palette_discrete(length(levels(pop(x))) + 1)
      }
      
      if (!is(palette_discrete, "function")) {
        colors_pops <- palette_discrete
        # if colors are not in RGB format
        if (grepl("#", colors_pops[1]) == FALSE) {
          colors_pops <- RGB_colors(colors_pops)
        }
      }
      
      data_long_1 <-
        as.data.frame(matrix(nrow = nPop(x), ncol = 6))
      colnames(data_long_1) <-
        c("source",
          "target",
          "value",
          "IDsource",
          "IDtarget",
          "color")
      data_long_1$source <- paste0(pall$pop1, " source")
      data_long_1$target <- "Rest"
      data_long_1$value <- pall$priv1
      data_long_1$IDsource <- (1:nPop(x)) - 1
      data_long_1$IDtarget <- (nPop(x) + 1) - 1
      data_long_1$color <- popNames(x)
      
      data_long_2 <-
        as.data.frame(matrix(nrow = nPop(x), ncol = 6))
      colnames(data_long_2) <-
        c("source",
          "target",
          "value",
          "IDsource",
          "IDtarget",
          "color")
      data_long_2$source <- "Rest"
      data_long_2$target <- paste0(pall$pop1, " target")
      data_long_2$value <- pall$priv2
      data_long_2$IDsource <- (nPop(x) + 1) - 1
      data_long_2$IDtarget <- (nPop(x) + 1):(nPop(x) * 2)
      data_long_2$color <- "Rest"
      
      data_long <- rbind(data_long_1, data_long_2)
      
      nodes <- as.data.frame(matrix(nrow = (nPop(x) * 2) + 1, ncol = 1))
      colnames(nodes) <- c("name")
      nodes$name <- c(data_long_1$source, "Rest", data_long_2$target)
      
      colors_pops <-
        paste0("\"", paste0(colors_pops, collapse = "\",\""), "\"")
      
      colorScal <-
        paste("d3.scaleOrdinal().range([", colors_pops, "])")
      # color links
      data_long$color <-
        gsub("src_", "", data_long$source)
      
      p3 <-
        suppressMessages(
          networkD3::sankeyNetwork(
            Links = data_long,
            Nodes = nodes,
            Source = "IDsource",
            Target = "IDtarget",
            LinkGroup = "color",
            Value = "value",
            NodeID = "name",
            sinksRight = FALSE,
            units = "Private alleles",
            colourScale = colorScal,
            iterations = 0,
            nodeWidth = 40,
            fontSize = 14,
            nodePadding = 20
          )
        )
    }
  }
  
  df <- pall
  
  # PRINTING OUTPUTS
  if (plot.out) {
    if (map.interactive & (method == "pairwise")) {
      labs <- popNames(x)
      gl.map.interactive(x, matrix = mm, symmetric = FALSE)
    }
    # using package patchwork
    print(p3)
  }
  
  if (verbose > 0) {
    print(df)
  }
  
  if (verbose >= 2) {
    cat(report("  Table of private alleles and fixed differences returned\n"))
  }
  
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
      saveRDS(list(match_call, p3), file = temp_plot)
      if (verbose >= 2) {
        cat(report("  Saving the ggplot to session tempfile\n"))
      }
    }
    temp_table <- tempfile(pattern = "Table_")
    saveRDS(list(match_call, df), file = temp_table)
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
  
  if(loc_names==TRUE){
    
    return(invisible(list(table=df,names_loci=pall_loc_names)))
    
  }else{
    
   return(invisible(df))
    
  }
  
}
