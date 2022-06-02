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
#' @param k The number for k the q matrix should be based on. Needs to
#'  be within you simulated range of k's in your sr structure run object.
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
#' #qmat <- gl.plot.structure(sr, k=3, CLUMPP='d:/structure/')
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
                                k = NULL,
                                met_clumpp = "greedyLargeK",
                                iter_clumpp = 100,
                                clumpak = TRUE,
                                plot_theme = NULL,
                                colors_clusters = NULL,
                                verbose = NULL) {
  # DO THE JOB
  
  if (!is(sr, "structure.result")) {
    stop(error(
      "sr is not a structure result object returned from gl.run.structure."
    ))
  }
  
  if (is.null(k)) {
    ks <- range((lapply(sr, function(x) {
      x$summary[1]
    })))
    ks <- ks[1]:ks[2]
  } else{
    ks <- k
  }
  
  res <- list()
  
  for (i in ks) {
    eq.k <- sapply(sr, function(x) {
      x$summary["k"] == i
    })
    
    if (sum(eq.k) == 0) {
      stop(error(paste(
        "No entries for k =", k, "found in 'sr'.\n"
      )))
    }
    
    sr_tmp <- sr[eq.k]
    
    Q_list <- lapply(sr_tmp, function(x) {
      as.matrix(x[[2]][, 4:ncol(x[[2]])])
    })
    
    # clumpak method for inferring modes within multiple structure runs as
    # implemented in starmie package
    
    if (ncol(Q_list[[1]]) == 1) {
      res <- c(res, as.matrix(Q_list[1]))
    } else{
      if (clumpak) {
        simMatrix <- as.matrix(proxy::simil(Q_list, method = G))
        diag(simMatrix) <- 1
        t <- calcThreshold(simMatrix)
        simMatrix[simMatrix < t] <- 0
        clusters <- mcl(simMatrix, addLoops = TRUE)$Cluster
        Q_list_clumpak <- split(Q_list, clusters)
        
        res_tmp <- lapply(Q_list_clumpak,
                          clumpp,
                          method = met_clumpp,
                          iter = iter_clumpp)
        
        # averaging replicates across modes inferred by clumpak method
        res_tmp_2 <- lapply(res_tmp, function(x) {
          if (length(x)[[1]] == 1) {
            return(x[[1]])
          } else{
            Reduce("+", x[[1]]) / length(x[[1]])
          }
        })
      } else{
        res_tmp <- clumpp(Q_list,
                          method = met_clumpp,
                          iter = iter_clumpp)
        
        res_tmp_2 <- res_tmp$Q_list
        
      }
      
      res <- c(res, res_tmp_2)
      
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
      pop = sr[[1]]$q.mat$orig.pop
    )
    n_col <- ncol(Q_list_tmp) - 3
    colnames(Q_list_tmp) <-
      c("Label", paste0(rep("cluster", n_col), 1:n_col), "K", "pop")
    cols_order <- colnames(Q_list_tmp)
    cols_order <- cols_order[grepl("cluster", cols_order)]
    Q_list_tmp$pop <- as.factor(Q_list_tmp$pop)
    Q_list_tmp_list <- split(Q_list_tmp, Q_list_tmp$pop)
    Q_list_tmp_list_2 <- lapply(Q_list_tmp_list, function(x) {
      data.table::setorderv(x, cols = cols_order, order = -1)
    })
    Q_list_tmp <- data.table::rbindlist(Q_list_tmp_list_2)
    Q_list_tmp$ord <- 1:nrow(Q_list_tmp)
    
    no_cluster <- colSums(Q_list_tmp[, 2:(1 + length(cols_order))]) > 0
    
    no_cluster <- c(TRUE, unname(no_cluster), TRUE, TRUE, TRUE)
    no_cluster_2 <-  as.vector((1:ncol(Q_list_tmp))[no_cluster])
    
    Q_list_tmp <- Q_list_tmp[, ..no_cluster_2]
    Q_list[[i]] <- Q_list_tmp
  }
  
  #Melt and append Q matrices
  Q_melt <-
    do.call("rbind",
            lapply(
              Q_list,
              reshape2::melt,
              id.vars = c("Label", "K", "pop", "ord"),
              variable.name = "Cluster"
            ))
  
  Q_melt$pop <-
    factor(Q_melt$pop, levels = unique(sr[[1]]$q.mat$orig.pop))
  
  if (is.null(plot_theme)) {
    plot_theme <- theme_dartR()
  }
  
  if (is.null(colors_clusters)) {
    colors_clusters <- structure_colors
  }
  
  if (is(colors_clusters, "function")) {
    cols_clusters <- colors_clusters(length(levels(Q_melt$Cluster)))
  }
  
  if (!is(colors_clusters, "function")) {
    cols_clusters <- colors_clusters
  }
  
  gg <- ggplot(Q_melt, aes(x = factor(ord), y = value, fill = Cluster)) +
    geom_col(color = "black", size = 0.2,width = 1) +
    facet_grid(K ~ pop , scales = "free", space = "free") +
    scale_y_continuous(expand = c(0, 0)) +
    # coord_cartesian(ylim=c(0,1)) +
    # xlab("Sample ID") +
    # ylab("Proportion of cluster") +
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
  
  return(gg)
  
}
