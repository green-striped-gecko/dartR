#' @name gl.plot.structure
#'
#' @title Plots a STRUCTURE analysis using a genlight object
#'
#' @description
#' This function takes a structure run object (output from
#'  \code{\link{gl.run.structure}}) and plots the typical strcture bar
#'   plot that visiualise the q matrix of a structure run.
#'
#' @param sr structure run object from \code{\link{gl.run.structure}} [required].
#' @param k the number for k the q matrix should be based on. Needs to
#'  be within you simulated range of k's in your sr structure run object.
#' @param sort how q matrix is sorted (by population, group proportion etc.)
#'  [default is by the population definition and group proportion from group
#'   1 to k-1, meaning first sorted by group, then by proportion of group 1,
#'    then group 2,.... then group k-1].
#' @param CLUMPP path to the clumpp executable file. Windows: CLUMPP.exe,
#'  macos and linux: CLUMPP (no exe).
#' @param ... additional parameter passed to the clumpp function within
#'  package strataG (\code{clumpp}).
#' @param plot_theme Theme for the plot. See details for options [default
#'  theme_dartR()].
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
#' #sr <- gl.run.structure(bc, k.range = 2:5, num.k.rep = 3, exec = './structure.exe')
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

### @importFrom strataG genind2gtypes structureRun

gl.plot.structure <- function(sr,
                              k,
                              sort = NULL,
                              CLUMPP = "./",
                              ...,
                              plot_theme,
                              verbose) {
  
  ################################################################
  .structureParseQmat <-    function (q.mat.txt, pops) 
  {
    q.mat.txt <- sub("[*]+", "", q.mat.txt)
    q.mat.txt <- sub("[(]", "", q.mat.txt)
    q.mat.txt <- sub("[)]", "", q.mat.txt)
    q.mat.txt <- sub("[|][ ]+$", "", q.mat.txt)
    cols1to4 <- c("row", "id", "pct.miss", 
                  "orig.pop")
    strsplit(q.mat.txt, " ") %>% purrr::map(function(q) {
      q <- q[!q %in% c("", " ", ":")] %>% 
        as.character() %>% rbind() %>% as.data.frame(stringsAsFactors = FALSE)
      stats::setNames(q, c(cols1to4, paste("Group", 1:(ncol(q) - 
                                                         4), sep = ".")))
    }) %>% dplyr::bind_rows() %>% dplyr::mutate_at(dplyr::vars("row", 
                                                               "pct.miss", "orig.pop", dplyr::starts_with("Group.")), 
                                                   as.numeric) %>% dplyr::mutate(orig.pop = if (!is.null(pops)) {
                                                     pops[.data$orig.pop]
                                                   }
                                                   else {
                                                     .data$orig.pop
                                                   })
  }
  
  
  ################################################################
  
    
  clumpp <-   function (sr, k, align.algorithm = "greedy", sim.stat = "g", 
            greedy.option = "ran.order", repeats = 100, order.by.run = 0, 
            label = NULL, delete.files = TRUE) 
  {
    if (k < 2) 
      stop("k must be greater than 1.")
    if (!tolower(align.algorithm) %in% c("full.search", 
                                         "greedy", "large.k")) {
      stop("'align.algorithm' must be either 'full.search', 'greedy', or 'large.k'.")
    }
    if (!tolower(sim.stat) %in% c("g", "g.prime")) {
      stop("'sim.stat' must be either 'g' or 'g.prime'.")
    }
    if (!"structure.result" %in% class(sr)) {
      folder <- sr
      if (!is.null(folder)) {
        if (!file.exists(folder)) 
          dir.create(folder)
        if (!utils::file_test("-d", folder)) {
          stop("'folder' is not a valid folder.")
        }
        label <- file.path(folder, label)
      }
    }
    if (!is.null(label)) 
      label <- gsub(" ", ".", label)
    param.file <- ifelse(is.null(label), "paramfile", paste(label, 
                                                            "paramfile", sep = "_"))
    ind.file <- ifelse(is.null(label), "indfile", paste(label, 
                                                        "indfile", sep = "_"))
    permutation.file <- ifelse(is.null(label), "permutationfile", 
                               paste(label, "permutationfile", sep = "_"))
    out.file <- ifelse(is.null(label), "outfile", paste(label, 
                                                        "outfile", sep = "_"))
    misc.file <- ifelse(is.null(label), "miscfile", paste(label, 
                                                          "miscfile", sep = "_"))
    eq.k <- sapply(sr, function(x) x$summary["k"] == k)
    if (sum(eq.k) == 0) 
      stop(paste("no entries for k =", k, "found in 'sr'."))
    sr <- sr[eq.k]
    align.algorithm <- switch(align.algorithm, full.search = 1, 
                              greedy = 2, large.k = 3)
    sim.stat <- switch(sim.stat, g = 1, g.prime = 2)
    greedy.params <- if (align.algorithm != 1) {
      greedy.option <- switch(greedy.option, all = 1, ran.order = 2)
      gp <- paste("GREEDY_OPTION", greedy.option)
      if (greedy.option != 1) {
        if (is.null(repeats)) {
          stop("'repeats' must be specified if 'greedy.options' is 'greedy' or 'large.k'.")
        }
        gp <- c(gp, paste("REPEATS", repeats), paste("PERMUTATIONFILE", 
                                                     permutation.file))
        perm.mat <- t(sapply(1:repeats, function(i) sample(1:length(sr))))
        utils::write.table(perm.mat, file = permutation.file, 
                           quote = FALSE, sep = " ", row.names = FALSE, 
                           col.names = FALSE)
      }
      gp
    }    else NULL
    param.txt <- c("DATATYPE 0", paste("INDFILE", 
                                       ind.file), paste("OUTFILE", out.file), paste("MISCFILE", 
                                                                                    misc.file), paste("K", k), paste("C", nrow(sr[[1]]$q.mat)), 
                   paste("R", length(sr)), paste("M", align.algorithm), 
                   "W 0", paste("S", sim.stat), greedy.params, 
                   "PRINT_PERMUTED_DATA 0", "PRINT_EVERY_PERM 0", 
                   "PRINT_RANDOM_INPUTORDER 0", "OVERRIDE_WARNINGS 1", 
                   paste("ORDER_BY_RUN", order.by.run))
    param.txt <- paste(param.txt, " ", sep = "")
    write(param.txt, file = param.file)
    id.order <- sr[[1]]$q.mat$id
    q.mat.df <- purrr::map(sr, function(x) {
      q.mat <- x$q.mat
      rownames(q.mat) <- q.mat$id
      q.mat <- q.mat[id.order, ]
      q.mat <- cbind(row = 1:nrow(q.mat), q.mat)
      rownames(q.mat) <- NULL
      cbind(q.mat[, 1:4], sep = rep(":", nrow(q.mat)), 
            q.mat[, 5:ncol(q.mat)])
    }) %>% dplyr::bind_rows()
    q.mat.df$pct.miss <- paste("(", q.mat.df$pct.miss, 
                               ")", sep = "")
    pop.fac <- factor(q.mat.df$orig.pop)
    pops <- levels(pop.fac)
    q.mat.df$orig.pop <- as.numeric(pop.fac)
    id.fac <- factor(q.mat.df$id)
    ids <- levels(id.fac)
    q.mat.df$id <- as.numeric(id.fac)
    utils::write.table(q.mat.df, file = ind.file, quote = FALSE, 
                       sep = " ", row.names = FALSE, col.names = FALSE)
    tryCatch({
      err.code <- system(paste("CLUMPP", param.file))
      if (err.code == 127) {
        stop("You do not have CLUMPP installed.")
      }
      else if (!err.code == 0) {
        stop(paste("Error running CLUMPP. Error code", 
                   err.code, "returned."))
      }
    })
    if (!file.exists(out.file)) {
      stop(paste("CLUMPP exited without creating output file, '", 
                 out.file, "'.", sep = ""))
    }
    out.txt <- scan(out.file, "character", sep = "\n", 
                    quiet = TRUE)
    q.mat <- .structureParseQmat(out.txt, pops)
    q.mat$id <- ids[as.numeric(q.mat$id)]
    q.mat$row <- NULL
    if (delete.files) {
      files <- c(ind.file, param.file, permutation.file, out.file, 
                 misc.file)
      for (f in files) if (file.exists(f)) 
        file.remove(f)
    }
    q.mat
  }  
  
  
        # DO THE JOB run clump
        
        if (!is(sr, "structure.result")) {
            stop(
                error(
                    "sr is not a structure result object returned from gl.run.structure."
                )
            )
        }
        
        # change range of simulated ks in structure object
        ks <- range((lapply(sr, function(x)
            x$summary[1])))
        
        if (is.null(k) | k < ks[1] | k > ks[2]) {
            stop(
                error(
                    "No k provided of k is not in range of the simulated k's in your structure run object"
                )
            )
        }
        
        if (Sys.info()["sysname"] == "Windows") {
            clumppfile <- "CLUMPP.exe"
        } else {
            clumppfile <- "CLUMPP"
            CLUMPP <- "CLUMPP"
        }
        
        if (!file.exists(file.path(CLUMPP, clumppfile))) {
            stop(error(
                "Cannot find clumpp executable. Please provide full path."
            ))
        }
        
        owd <- getwd()
        setwd(CLUMPP)
        q.mat <- clumpp(sr, k = k)
        setwd(owd)
        
        qq <- q.mat[, 4:(k + 3)]
        
        if (is.null(sort)) {
            ll <- data.frame(cbind(as.numeric(factor(
                q.mat$orig.pop
            )), qq[, ]))
            zz <- do.call(order, unname(as.list(ll)))
            
            bb <- t(qq[zz, ])
            colnames(bb) <- q.mat$id
            narg <-
                paste(q.mat$id, q.mat$orig.pop[zz], sep = "_")
            
        } else {
            bb <- t(qq[sort, ])
            colnames(bb) <- q.mat$id[sort]
            narg <- q.mat$id[sort]
        }
        
        # bgg <- reshape2::melt(bb) bgg <- transform(bgg, Var1=factor(Var1, rownames(bb)), Var2=factor(Var2, colnames(bb)))
        
        bbpp <-
            barplot(
                bb,
                col = 1:k,
                las = 2,
                main = paste0("K=", k),
                border = 1:k,
                space = 0,
                names.arg = narg
            )
        # ggplot(bgg, aes(x=Var2, y=value, fill=Var1), )+geom_bar(stat='identity', width = 1)

    return(q.mat)
}
