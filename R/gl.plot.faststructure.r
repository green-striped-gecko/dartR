#' @name gl.plot.faststructure
#'
#' @title Plots fastStructure analysis results (Q-matrix)
#'
#' @description
#' This function takes a fastStructure run object (output from
#'  \code{\link{gl.run.faststructure}}) and plots the typical structure bar
#'   plot that visualize the q matrix of a fastStructure run.
#'
#' @param sr fastStructure run object from \code{\link{gl.run.faststructure}}
#'  [required].
#' @param k.range The number for K of the q matrix that should be plotted. Needs 
#' to be within you simulated range of K's in your sr structure run object. If 
#'  NULL, all the K's are plotted [default NULL].
#' @param met_clumpp The algorithm to use to infer the correct permutations.
#' One of 'greedy' or 'greedyLargeK' or 'stephens' [default "greedyLargeK"].
#' @param iter_clumpp The number of iterations to use if running either 'greedy'
#'  'greedyLargeK' [default 100].
#' @param clumpak Whether use the Clumpak method (see details) [default TRUE].
#' @param plot_theme Theme for the plot. See Details for options
#' [default NULL].
#' @param colors_clusters A color palette for clusters (K) or a list with
#' as many colors as there are clusters (K) [default NULL].
#' @param ind_name Whether to plot individual names [default TRUE].
#' @param border_ind The width of the border line between individuals 
#' [default 0.25]. 
#'
#' @details The function outputs a barplot which is the typical output of
#'  fastStructure.
#'  
#'  This function is based on the methods of CLUMPP and Clumpak as implemented 
#'  in the R package starmie (https://github.com/sa-lee/starmie). 
#'  
#'  The Clumpak method identifies sets of highly similar runs among 
#'  all the replicates of the same K. The method then separates the distinct 
#'  groups of runs representing distinct modes in the space of possible solutions.
#'  
#'  The CLUMPP method permutes the clusters output by independent runs of 
#'  clustering programs such as structure, so that they match up as closely as 
#'  possible.
#'  
#'  This function averages the replicates within each mode identified by the 
#'  Clumpak method.
#'
#' Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'
#' @return List of Q-matrices
#'
#' @author Bernd Gruber & Luis Mijangos (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' \dontrun{
#' t1 <- gl.filter.callrate(platypus.gl,threshold = 1)
#' res <- gl.run.faststructure(t1, exec = "./fastStructure",k.range = 2:3, 
#'                           num.k.rep = 2,output = paste0(getwd(),"/res_str"))
#' qmat <- gl.plot.faststructure(res,k.range=2:3)
#' gl.map.structure(qmat, K=2, t1, scalex=1, scaley=0.5)
#' }
#' @export
#' @seealso \code{gl.run.faststructure}
#' @references
#' \itemize{
#' \item Raj, A., Stephens, M., & Pritchard, J. K. (2014). fastSTRUCTURE: 
#' variational inference of population structure in large SNP data sets. 
#' Genetics, 197(2), 573-589.
#' \item Pritchard, J.K., Stephens, M., Donnelly, P. (2000) Inference of
#' population structure using multilocus genotype data. Genetics 155, 945-959.
#' \item Kopelman, Naama M., et al. "Clumpak: a program for identifying 
#' clustering modes and packaging population structure inferences across K." 
#' Molecular ecology resources 15.5 (2015): 1179-1191.
#' \item Mattias Jakobsson and Noah A. Rosenberg. 2007. CLUMPP: a cluster
#' matching and permutation program for dealing with label switching and
#' multimodality in analysis of population structure. Bioinformatics
#' 23(14):1801-1806. Available at
#' \href{http://web.stanford.edu/group/rosenberglab/clumppDownload.html}{clumpp}
#' }

gl.plot.faststructure <- function(sr,
                                k.range,
                                met_clumpp = "greedyLargeK",
                                iter_clumpp = 100,
                                clumpak = TRUE,
                                plot_theme = NULL,
                                colors_clusters = NULL,
                                ind_name = TRUE,
                                border_ind = 0.15
) {
  
  res <- list()
  
  for (i in k.range) {
    eq.k <- which(names(sr) == as.character(i))
    
    sr_tmp <- sr[[eq.k]]
    
    Q_list_tmp <- lapply(sr_tmp, function(x) {
      as.matrix(x[, 3:ncol(x)])
    })
    
    # If K = 1
    if (ncol(Q_list_tmp[[1]]) == 1) {
      res[[length(res)+1]]  <- c(res, as.matrix(Q_list_tmp[1]))
      # If K > 1  
    }else{
      # If just one replicate 
      if(length(Q_list_tmp)==1){
        res_tmp <- Q_list_tmp[[1]]
        # if more than 1 replicate
      }else{
        res_tmp <- dartR::utils.clumpp(Q_list_tmp, 
                                  method = met_clumpp, 
                                  iter = iter_clumpp)$Q_list
      }
      
      # clumpak method for inferring modes within multiple structure runs as
      # implemented in starmie package
      if (clumpak) {
        # if just one replicate
        if(length(res_tmp)==1){
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
        
        # averaging replicates
        # if there is just one mode
        if(length(res_tmp_2)==1){
          # if there is just one replicate within the mode
          if(length(res_tmp_2[[1]])==1){
            res_tmp_3 <- res_tmp_2[[1]]
            # if there are more than 1 replicate within the mode
          }else{
            res_tmp_3 <- as.matrix(Reduce("+", res_tmp_2[[1]]) / length(res_tmp_2[[1]]))
          }
          # if there are more than 1 mode  
        }else{
          res_tmp_3 <- lapply(res_tmp_2, function(x) {
            # if there is just one replicate within the mode
            if(length(x[1])==1){
              return(x[[1]])
              #if there are more than 1 replicate within the mode
            }else{
              return(Reduce("+", x[1]) / length(x[1]))
            }
          })
        }
        
      }else{
        
        res_tmp_3 <- res_tmp 
      }
      
      # if the object is a list
      if(is.list(res_tmp_3)){
        
        for (y in 1:length(res_tmp_3)) {
          
          res[[length(res)+1]] <- res_tmp_3[y]
          
        }
        
      }else{
        res[[length(res)+1]] <- res_tmp_3
      }
      
    }
  }
  
  #flattening lists
  renquote <- function(l){
    if (is.list(l)){
      lapply(l, renquote) 
    } else {
      enquote(l)
    }
  } 
  
  res <- lapply(unlist(renquote(res)), eval)
  
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
      Label = sr[[1]][[1]]$id,
      Q_list[[i]],
      K = rep(Ks[[i]], nrow(Q_list[[i]])),
      orig.pop = sr[[1]][[1]]$orig.pop
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
    Q_list[[i]] <- Q_list_tmp
    Q_list[[i]] <- as.data.frame(Q_list[[i]])
  }
  
  order_df <- Q_list[[1]][order(Q_list[[1]]$Label),]
  
  Q_list <- lapply(Q_list,function(x){
    tmp <- x[order(x$Label),]
    tmp$ord <- order_df$ord
    return(tmp)
  })
  
  if (is.null(plot_theme)) {
    plot_theme <- dartR::theme_dartR()
  }
  
  if (is.null(colors_clusters)) {
    colors_clusters <- structure_colors
  }
  
  if (is(colors_clusters, "function")) {
    cols_clusters <- colors_clusters(max(k.range))
  }
  
  if (!is(colors_clusters, "function")) {
    cols_clusters <- colors_clusters
  }
  
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
    factor(Q_melt$orig.pop, levels = unique(sr[[1]][[1]]$orig.pop))
  
  p3 <- ggplot(Q_melt, aes_(x= ~ factor(ord), y = ~value, fill = ~Cluster)) +
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
    p3 + theme(axis.text.x = element_blank(),
               axis.ticks.x = element_blank())
    
  }
  
  print(p3)
  
  return(invisible(Q_list))
  
}
