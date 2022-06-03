#' @name gl.plot.structure_2
#'
#' @title Plots a STRUCTURE analysis using a genlight object
#'
#' @description
#' This function takes a structure run object (output from
#'  \code{\link{gl.run.structure}}) and plots the typical strcture bar
#'   plot that visiualise the q matrix of a structure run.
#'
#' @param sr Structure run object from \code{\link{gl.run.structure}} [required].
#' @param K The number for K the q matrix should be based on. Needs to
#'  be within you simulated range of K's in your sr structure run object.
#' @param met_clumpp The algorithm to use to infer the correct permutations.
#' One of 'greedy' or 'greedyLargeK' or 'stephens' [default "greedy"].
#' @param iter_clumpp The number of iterations to use if running either 'greedy'
#'  'greedyLargeK' [default 100].
#' @param clumpak [default TRUE].
#' @param plot_theme Theme for the plot. See Details for options
#' [default NULL].
#' @param colors_clusters A color palette for clusters (K) or a list with
#' as many colors as there are clusters (K)
#' [default NULL].
#' @param ind_name Whether to plot individual names [default TRUE].
#' @param border_ind The width of the border line between individuals [default 0.25]. 
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log ; 3, progress and results summary; 5, full report [default
#'   NULL, unless specified using gl.set.verbosity]
#'
#' @details The function outputs a barplot which is the typical output of
#'  structure. This function needs the use of clumpp, which is an external
#'   program that needs to be installed For a evanno plot use gl.evanno.
#'
#' @return a barplot in ggplot format
#'
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' \dontrun{
#' #CLUMPP needs to be installed to be able to run the example
#' #bc <- bandicoot.gl[,1:100]
#' #sr <- gl.run.structure(bc, k.range = 2:5, num.k.rep = 3,
#' exec = './structure.exe')
#' #ev <- gl.evanno(sr)
#' #ev
#' #qmat <- gl.plot.structure(sr, K=3, CLUMPP='d:/structure/')
#' #head(qmat)
#' #gl.map.structure(qmat, bc, scalex=1, scaley=0.5)
#' }
#' @export
#' @seealso \code{gl.run.structure},  \code{clumpp}, \code{gl.plot.structure}
#' @references
#' \itemize{
#' \item Pritchard, J.K., Stephens, M., Donnelly, P. (2000) Inference of
#' population structure using multilocus genotype data. Genetics 155, 945-959.
#' \item Archer, F. I., Adams, P. E. and Schneiders, B. B. (2016) strataG: An R
#' package for manipulating, summarizing and analysing population genetic data.
#' Mol Ecol Resour. doi:10.1111/1755-0998.12559
#' \item Mattias Jakobsson and Noah A. Rosenberg. 2007. CLUMPP: a cluster
#' matching and permutation program for dealing with label switching and
#' multimodality in analysis of population structure. Bioinformatics
#' 23(14):1801-1806. Available at
#' \href{http://web.stanford.edu/group/rosenberglab/clumppDownload.html}{clumpp}
#' }

gl.plot.structure_2 <- function(sr,
                                K = NULL,
                                met_clumpp = "greedyLargeK",
                                iter_clumpp = 100,
                                clumpak = TRUE,
                                plot_theme = NULL,
                                colors_clusters = NULL,
                                ind_name = TRUE,
                                border_ind = 0.15,
                                plot.out = TRUE,
                                verbose = NULL) {
  # DO THE JOB
  
  if (!is(sr, "structure.result")) {
    stop(error(
      "sr is not a structure result object returned from gl.run.structure."
    ))
  }
  
  if (is.null(K)) {
    ks <- range((lapply(sr, function(x) {
      x$summary[1]
    })))
    ks <- ks[1]:ks[2]
  } else{
    ks <- K
  }
  
  res <- list()
  
  for (i in ks) {
    eq.k <- sapply(sr, function(x) {
      x$summary["k"] == i
    })
    
    if (sum(eq.k) == 0) {
      stop(error(paste(
        "No entries for K =", K, "found in 'sr'.\n"
      )))
    }
    
    sr_tmp <- sr[eq.k]
    
    Q_list_tmp <- lapply(sr_tmp, function(x) {
      as.matrix(x[[2]][, 4:ncol(x[[2]])])
    })
    
    # If K = 1
    if (ncol(Q_list_tmp[[1]]) == 1) {
      res <- c(res, as.matrix(Q_list_tmp[1]))
    # If K > 1  
    }else{
      # If just one replicate 
      if(length(Q_list_tmp)==1){
        res_tmp <- Q_list_tmp[[1]]
      # if more than 1 replicate
      }else{
        res_tmp <- clumpp(Q_list_tmp, method = met_clumpp, iter = iter_clumpp)$Q_list
      }
      
      # clumpak method for inferring modes within multiple structure runs as
      # implemented in starmie package
      if (clumpak) {
        # if just one replicate
        if(lenght(res_tmp)==1){
          res_tmp_2 <- res_tmp[[1]]
        # if more than one replicate
        }else{
          simMatrix <- as.matrix(proxy::simil(res_tmp, method = G))
          diag(simMatrix) <- 1
          t <- calcThreshold(simMatrix)
          simMatrix[simMatrix < t] <- 0
          clusters <- mcl(simMatrix, addLoops = TRUE)$Cluster
          res_tmp_2 <- split(res_tmp, clusters)
        }
      }else{
        res_tmp_2 <- res_tmp 
      }
      
      # averaging replicates
      # if there is just one mode
        if(length(res_tmp_2)==1){
          res_tmp_3 <- as.matrix(Reduce("+", res_tmp_2[[1]]) / length(res_tmp_2[[1]]))
      # if there are more than 1 mode  
        }else{
          res_tmp_3 <- lapply(res_tmp_2, function(x) {
            return(Reduce("+", x[[1]]) / length(x[[1]]))
          })
        }
      
      res <- c(res, res_tmp_3)

    }
  }
  
  names(res) <- as.character(1:length(res))
  
  Q_list <- res
  
  #get K labels
  Ks <- unlist(lapply(Q_list, ncol))
  if (length(unique(Ks)) != length(Ks)) {
    #Repeated Ks so label with subnumbering
    Ks <- paste(Ks, ave(Ks, Ks, FUN = seq_along), sep = ".")
  } else{
    Ks <- as.character(Ks)
  }
  
  for (i in 1:length(Q_list)) {
    Q_list_tmp <- data.frame(
      Label = sr[[1]]$q.mat$id,
      Q_list[[i]],
      K = rep(Ks[[i]], nrow(Q_list[[i]])),
      orig.pop = sr[[1]]$q.mat$orig.pop
    )
    n_col <- ncol(Q_list_tmp) - 3
    colnames(Q_list_tmp) <-
      c("Label", paste0(rep("cluster", n_col), 1:n_col), "K", "orig.pop")
    cols_order <- colnames(Q_list_tmp)
    cols_order <- cols_order[grepl("cluster", cols_order)]
    Q_list_tmp$orig.pop <- as.factor(Q_list_tmp$orig.pop)
    Q_list_tmp_list <- split(Q_list_tmp, Q_list_tmp$orig.pop)
    Q_list_tmp_list_2 <- lapply(Q_list_tmp_list, function(x) {
      data.table::setorderv(x, cols = cols_order, order = -1)
    })
    Q_list_tmp <- data.table::rbindlist(Q_list_tmp_list_2)
    Q_list_tmp$ord <- 1:nrow(Q_list_tmp)
    
    # no_cluster <- colSums(Q_list_tmp[, 2:(1 + length(cols_order))]) > 0
    # 
    # no_cluster <- c(TRUE, unname(no_cluster), TRUE, TRUE, TRUE)
    # no_cluster_2 <-  as.vector((1:ncol(Q_list_tmp))[no_cluster])
    # 
    # Q_list_tmp <- Q_list_tmp[, ..no_cluster_2]
    Q_list[[i]] <- Q_list_tmp
  }
  
  if (is.null(plot_theme)) {
    plot_theme <- theme_dartR()
  }
  
  if (is.null(colors_clusters)) {
    colors_clusters <- structure_colors
  }
  
  if (is(colors_clusters, "function")) {
    cols_clusters <- colors_clusters(max(ks))
  }
  
  if (!is(colors_clusters, "function")) {
    cols_clusters <- colors_clusters
  }
  
  # clumpp_plots <- lapply(Q_list, function(z) {
  #   
  #   Q_melt <- data.table::melt(z,id.vars = c("K","Label","orig.pop", "ord"))
  #   Q_melt$orig.pop <- factor(Q_melt$orig.pop, levels = unique(sr[[1]]$q.mat$orig.pop))
  #   Q_melt$variable <- factor(Q_melt$variable, levels = unique(Q_melt$variable))
  #   Q_melt$ord <- as.numeric(Q_melt$ord)
  #   
  #   Q_melt <- as.data.frame(Q_melt)
  #   
  #   p_temp <- ggplot(Q_melt, aes_(x= ~ ord, y = ~value, fill = ~variable)) +
  #     geom_col(color = "black", size = 0.2,width = 1) +
  #      facet_grid(K~orig.pop, scales = "free", space = "free") +
  #     scale_y_continuous(expand = c(0, 0)) +
  #     scale_x_discrete(
  #       breaks = Q_melt$ord,
  #       labels = Q_melt$Label,
  #       expand = c(0, 0)
  #     ) +
  #     scale_fill_manual(values = cols_clusters[1:length(unique(Q_melt$variable))]) +
  #     plot_theme +
  #     theme(
  #       panel.spacing = unit(0, "lines"),
  #       panel.border = element_rect(
  #         color = "black",
  #         fill = NA,
  #         size = 1
  #       ),
  #       strip.background = element_blank(),
  #       strip.text.x = element_text(size = 12, angle = 90),
  #       axis.title.x = element_blank(),
  #       axis.text.x = element_text(
  #         size = 8,
  #         angle = 90,
  #         vjust = 0.5,
  #         hjust = 1
  #       ),
  #       axis.title.y = element_blank(),
  #       axis.text.y = element_blank(),
  #       axis.ticks.y = element_blank() ,
  #       legend.position = "none"
  #     )
  # 
  #    # p_temp
  #   return(p_temp)
  # })
  
  
  # #Melt and append Q matrices
  Q_melt <-
    do.call("rbind",
            lapply(
              Q_list,
              reshape2::melt,
              id.vars = c("Label", "K", "orig.pop", "ord"),
              variable.name = "Cluster"
            ))

  Q_melt$orig.pop <-
    factor(Q_melt$orig.pop, levels = unique(sr[[1]]$q.mat$orig.pop))

  gg <- ggplot(Q_melt, aes_(x= ~ factor(ord), y = ~value, fill = ~Cluster)) +
    geom_col(color = "black", size = border_ind,width = 1) +
    facet_grid(K ~ orig.pop , scales = "free", space = "free") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(
      breaks = unique(Q_melt$ord),
      labels = unique(Q_melt$Label),
      expand = c(0, 0)
    ) +
    scale_fill_manual(values = cols_clusters) +
    plot_theme +
    theme(
      panel.spacing = unit(0, "lines"),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        size = 1
      ),
      strip.background = element_blank(),
      strip.text.x = element_text(size = 12, angle = 90),
      axis.title.x = element_blank(),
      axis.text.x = element_text(
        size = 8,
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank() ,
      legend.position = "none"
    )
  
  if(ind_name==FALSE){
    gg + theme(axis.text.x = element_blank(),
               axis.ticks.x = element_blank())
    
  }
  
  # gg
  
  if(plot.out){
  print(gg)
  }
  
  return(Q_list)
  
}
