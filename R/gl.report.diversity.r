#' @name gl.report.diversity
#'
#' @title Calculates diversity indexes for SNPs
#'
#' @description
#'This script takes a genlight object and calculates alpha and beta diversity
#'for q = 0:2. Formulas are taken from Sherwin et al. 2017. The paper describes
#'nicely the relationship between the different q levels and how they relate to
#'population genetic processes such as dispersal and selection.
#'
#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot_theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot_colors A color palette or a list with as many colors as there are 
#' populations in the dataset [default discrete_palette].
#' @param pbar Report on progress. Silent if set to FALSE [default TRUE].
#' @param table Prints a tabular output to the console either 'D'=D values, or
#'  'H'=H values or 'DH','HD'=both or 'N'=no table. [default 'DH'].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session
#' temporary directory (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, unless specified using gl.set.verbosity].
#'
#' @details For all indexes, the entropies (H) and corresponding effective
#'  numbers, i.e. Hill numbers (D), which reflect the number of needed entities
#'  to get the observed values, are calculated. In a nutshell, the alpha indexes
#'  between the different q-values should be similar if there is no deviation
#'  from expected allele frequencies and occurrences (e.g. all loci in HWE &
#'  equilibrium). If there is a deviation of an index, this links to a process
#'  causing it, such as dispersal, selection or strong drift. For a detailed
#'  explanation of all the indexes, we recommend resorting to the literature
#'  provided below. Confidence intervals are +/- 1 standard deviation.
#'
#'\strong{ Function's output }
#'
#'  Plots are saved to the temporal directory (tempdir) and can be accessed with
#'   the function \code{\link{gl.print.reports}} and listed with the function
#'    \code{\link{gl.list.reports}}. Note that they can be accessed only in the
#'     current R session because tempdir is cleared each time that the R session
#'     is closed.
#'
#'  Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'
#' @return A list of entropy indexes for each level of q and equivalent numbers
#'  for alpha and beta diversity.
#'
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr}),
#'  Contributors: William B. Sherwin, Alexander Sentinella
#'
#' @examples
#' div <- gl.report.diversity(bandicoot.gl, table = FALSE, pbar=FALSE)
#' div$zero_H_alpha
#' div$two_H_beta
#' names(div)
#'
#' @family reporting functions
#'
#' @references
#'Sherwin, W.B., Chao, A., Johst, L., Smouse, P.E. (2017). Information Theory
#' Broadens the Spectrum of Molecular Ecology and Evolution. TREE 32(12)
#'  948-963. doi:10.1016/j.tree.2017.09.12
#'
#' @import reshape2
#'
#' @export
#'

### To be done: adjust calculation of betas for population sizes (switch)

gl.report.diversity <- function(x,
                                plot.out = TRUE,
                                pbar = TRUE,
                                table = "DH",
                                plot_theme = theme_dartR(),
                                plot_colors = discrete_palette,
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
    
    # DO THE JOB
    
    if (is.null(pop(x))) {
        pop(x) <- factor(rep("all", nInd(x)))
    }
    
    # split in pops
    pops <- seppop(x)
    
    if (pbar) {
        pb <- txtProgressBar(0, 8, style = 3, width = 20)
        if (verbose >= 2) {
            cat(report(" Counting missing loci...           "))
        }
    }
    # number of missing loci
    nlocpop <-
        lapply(pops, function(x)
            sum(!is.na(colMeans(
                as.matrix(x), na.rm = T
            ))))
    
    if (pbar) {
        setTxtProgressBar(pb, 1)
        if (verbose >= 2) {
            cat(report(" Calculating zero_H/D_alpha ...           "))
        }
        ## 0Halpha (average number of alleles, ignoring missing values and sd)
    }
    zero_H_alpha_es <- lapply(pops, function(x) {
        dummys <- ((colMeans(as.matrix(x), na.rm = T) %% 2) > 0) + 1 - 1
        return(list(
            estH = mean(dummys, na.rm = T),
            sdH = sd(dummys, na.rm = T),
            estD = mean(dummys, na.rm = T) + 1,
            sdD = sd(dummys, na.rm = T)
        ))
    })
    zero_H_alpha <-
        unlist(lapply(zero_H_alpha_es, function(x)
            x[[1]]))
    zero_H_alpha_sd <-
        unlist(lapply(zero_H_alpha_es, function(x)
            x[[2]]))
    zero_D_alpha <-
        unlist(lapply(zero_H_alpha_es, function(x)
            x[[3]]))
    zero_D_alpha_sd <-
        unlist(lapply(zero_H_alpha_es, function(x)
            x[[4]]))
    
    if (pbar) {
        setTxtProgressBar(pb, 2)
        if (verbose >= 2) {
            cat(report(" Calculating one_H/D_alpha ...          "))
        }
    }
    ### one_H_alpha
    shannon <- function(x) {
        x <- x[x > 0]
        p <- x / sum(x)
        - sum(p * log(p))
    }
    
    
    one_H_alpha_es <- lapply(pops, function(x) {
        mat_temp <- as.matrix(x)
        mat <- matrix(nrow = 6, ncol = nLoc(x))
        rownames(mat) <-
            c("AA", "AB", "BB", "A", "B", "shannon")
        mat["AA", ] <- apply(mat_temp, 2, function(y) {
            length(y[which(y == 0)])
        })
        mat["AB", ] <- apply(mat_temp, 2, function(y) {
            length(y[which(y == 1)])
        })
        mat["BB", ] <- apply(mat_temp, 2, function(y) {
            length(y[which(y == 2)])
        })
        mat["A", ] <- 2 * mat["AA", ] + mat["AB", ]
        mat["B", ] <- 2 * mat["BB", ] + mat["AB", ]
        mat_shannon <- mat[c("A", "B"), ]
        
        dummys <- apply(mat_shannon, 2, shannon)
        
        # p <- (2 * p + hets) / 2 q <- (2 * q + hets) / 2 total <- colSums(mat,na.rm = T) p <- colMeans(as.matrix(x), na.rm = T)/2 p <-
        # p[!is.na(p)] #ignore loci with just missing data logp <- ifelse(!is.finite(log(p)), 0, log(p)) log1_p <- ifelse(!is.finite(log(1
        # - p)), 0, log(1 - p)) dummys <- -(p * logp + (1 - p) * log1_p)
        
        return(list(
            estH = mean(dummys),
            sdH = sd(dummys),
            estD = mean(exp(dummys)),
            sdD = sd(exp(dummys)),
            dummys = dummys
        ))
    })
    one_H_alpha <-
        unlist(lapply(one_H_alpha_es, function(x)
            x[[1]]))
    one_H_alpha_sd <-
        unlist(lapply(one_H_alpha_es, function(x)
            x[[2]]))
    one_D_alpha <-
        unlist(lapply(one_H_alpha_es, function(x)
            x[[3]]))
    one_D_alpha_sd <-
        unlist(lapply(one_H_alpha_es, function(x)
            x[[4]]))
    if (pbar) {
        setTxtProgressBar(pb, 3)
        if (verbose >= 2) {
            cat(report(" Calculating two_H/D_alpha ...           "))
        }
    }
    # two_H_alpha
    two_H_alpha_es <- lapply(pops, function(x) {
        p <- colMeans(as.matrix(x), na.rm = T) / 2
        
        p <- p[!is.na(p)]  #ignore loci with just missing data
        dummys <- (1 - (p * p + (1 - p) * (1 - p)))
        
        return(list(
            estH = mean(dummys),
            sdH = sd(dummys),
            estD = mean(1 / (1 - dummys)),
            sdD = sd(1 / (1 - dummys)),
            dummys = dummys
        ))
    })
    two_H_alpha <-
        unlist(lapply(two_H_alpha_es, function(x)
            x[[1]]))
    two_H_alpha_sd <-
        unlist(lapply(two_H_alpha_es, function(x)
            x[[2]]))
    two_D_alpha <-
        unlist(lapply(two_H_alpha_es, function(x)
            x[[3]]))
    two_D_alpha_sd <-
        unlist(lapply(two_H_alpha_es, function(x)
            x[[4]]))
    
    # initiallize betas as NA
    mat_zero_H_beta <- NA
    mat_one_H_beta <- NA
    mat_two_H_beta <- NA
    npops <- length(pops)
    
    if (npops > 1) {
        if (pbar) {
            setTxtProgressBar(pb, 4)
            if (verbose >= 2) {
                cat(report(" Counting pairwise missing loci..."))
            }
        }
        pairs <- t(combn(npops, 2))
        ### pairwise missing loci
        nlocpairpop <- apply(pairs, 1, function(x) {
            pop1 <- pops[[x[1]]]
            pop2 <- pops[[x[2]]]
            pp1 <- colMeans(as.matrix(pop1), na.rm = T) / 2
            pp2 <- colMeans(as.matrix(pop2), na.rm = T) / 2
            index <- !is.na(pp1) & !is.na(pp2)
            return(sum(index))
        })
        mat_nloc_pops <- matrix(NA, nrow = npops, ncol = npops)
        mat_nloc_pops[lower.tri(mat_nloc_pops)] <- nlocpairpop
        colnames(mat_nloc_pops) <-
            rownames(mat_nloc_pops) <- names(pops)
        
        if (pbar) {
            setTxtProgressBar(pb, 5)
            if (verbose >= 2) {
                cat(report(" Calculating zero_H/D_beta ...         "))
            }
        }
        # zero_H_beta
        zero_H_beta_es <- apply(pairs, 1, function(x) {
            pop1 <- pops[[x[1]]]
            pop2 <- pops[[x[2]]]
            pp1 <- colMeans(as.matrix(pop1), na.rm = T) / 2
            
            pp1 <- ifelse(pp1 > 0 & pp1 < 1, 0.5, pp1)
            pp2 <- colMeans(as.matrix(pop2), na.rm = T) / 2
            pp2 <- ifelse(pp2 > 0 & pp2 < 1, 0.5, pp2)
            
            index <- !is.na(pp1) & !is.na(pp2)
            pp1 <- pp1[index]
            pp2 <- pp2[index]
            
            dummys <- abs(pp1 - pp2)
            return(list(
                estH = mean(dummys),
                sdH = sd(dummys),
                estD = mean(dummys) + 1,
                sdD = sd(dummys)
            ))
            
        })
        
        zero_H_beta <-
            unlist(lapply(zero_H_beta_es, function(x)
                x[[1]]))
        zero_H_beta_sd <-
            unlist(lapply(zero_H_beta_es, function(x)
                x[[2]]))
        zero_D_beta <-
            unlist(lapply(zero_H_beta_es, function(x)
                x[[3]]))
        zero_D_beta_sd <-
            unlist(lapply(zero_H_beta_es, function(x)
                x[[4]]))
        
        mat_zero_H_beta <-
            matrix(NA, nrow = npops, ncol = npops)
        mat_zero_H_beta[lower.tri(mat_zero_H_beta)] <-
            zero_H_beta
        mat_zero_H_beta[pairs] <- zero_H_beta_sd
        colnames(mat_zero_H_beta) <-
            rownames(mat_zero_H_beta) <- names(pops)
        
        mat_zero_D_beta <-
            matrix(NA, nrow = npops, ncol = npops)
        mat_zero_D_beta[lower.tri(mat_zero_D_beta)] <-
            zero_D_beta
        mat_zero_D_beta[pairs] <- zero_D_beta_sd
        colnames(mat_zero_D_beta) <-
            rownames(mat_zero_D_beta) <- names(pops)
        
        if (pbar) {
            setTxtProgressBar(pb, 6)
            if (verbose >= 2) {
                cat(report(" Calculating one_H/D_beta ...    "))
            }
        }
        
        # one_H_beta calculate one_H_alpha_all for combined pops
        p <- colMeans(as.matrix(x), na.rm = T) / 2
        i0 <-
            which(!is.na(p))  #ignore loci with just missing data
        logp <- ifelse(!is.finite(log(p)), 0, log(p))
        log1_p <- ifelse(!is.finite(log(1 - p)), 0, log(1 - p))
        one_H_alpha_all <- -(p * logp + (1 - p) * log1_p)
        
        one_H_beta_es <- apply(pairs, 1, function(x) {
            i1 <-
                which(!is.na(colMeans(as.matrix(pops[[x[1]]]), na.rm = T) / 2))
            i2 <-
                which(!is.na(colMeans(as.matrix(pops[[x[2]]]), na.rm = T) / 2))
            tt <- table(c(i0, i1, i2))
            index <- as.numeric(names(tt)[tt == 3])
            dummys <-
                one_H_alpha_all[i0 %in% index] - (one_H_alpha_es[[x[1]]]$dummys[i1 %in% index] + one_H_alpha_es[[x[2]]]$dummys[i2 %in%
                                                                                                                                   index]) / 2
            return(list(
                estH = mean(dummys),
                sdH = sd(dummys),
                estD = mean(exp(dummys)),
                sdD = sd(exp(dummys))
            ))
        })
        
        one_H_beta <-
            unlist(lapply(one_H_beta_es, function(x)
                x[[1]]))
        one_H_beta_sd <-
            unlist(lapply(one_H_beta_es, function(x)
                x[[2]]))
        one_D_beta <-
            unlist(lapply(one_H_beta_es, function(x)
                x[[3]]))
        one_D_beta_sd <-
            unlist(lapply(one_H_beta_es, function(x)
                x[[4]]))
        
        mat_one_H_beta <- matrix(NA, nrow = npops, ncol = npops)
        mat_one_H_beta[lower.tri(mat_one_H_beta)] <- one_H_beta
        mat_one_H_beta[pairs] <- one_H_beta_sd
        colnames(mat_one_H_beta) <-
            rownames(mat_one_H_beta) <- names(pops)
        
        mat_one_D_beta <- matrix(NA, nrow = npops, ncol = npops)
        mat_one_D_beta[lower.tri(mat_one_D_beta)] <- one_D_beta
        mat_one_D_beta[pairs] <- one_D_beta_sd
        colnames(mat_one_D_beta) <-
            rownames(mat_one_D_beta) <- names(pops)
        
        if (pbar) {
            setTxtProgressBar(pb, 7)
            cat(report(" Calculating two_H/D_beta...    "))
        }
        
        p <- colMeans(as.matrix(x), na.rm = T) / 2
        i0 <-
            which(!is.na(p))  #ignore loci with just missing data
        two_H_alpha_all <- (1 - (p * p + (1 - p) * (1 - p)))
        
        two_H_beta_es <- apply(pairs, 1, function(x) {
            i1 <-
                which(!is.na(colMeans(as.matrix(pops[[x[1]]]), na.rm = T) / 2))
            i2 <-
                which(!is.na(colMeans(as.matrix(pops[[x[2]]]), na.rm = T) / 2))
            tt <- table(c(i0, i1, i2))
            index <- as.numeric(names(tt)[tt == 3])
            
            m2Ha <-
                (two_H_alpha_es[[x[1]]]$dummys[i1 %in% index] + two_H_alpha_es[[x[2]]]$dummys[i2 %in% index]) /
                2
            dummys <-
                ((two_H_alpha_all[i0 %in% index] - m2Ha) / (1 - m2Ha)) * (npops / (npops - 1))
            return(list(
                estH = mean(dummys),
                sdH = sd(dummys),
                estD = mean(exp(dummys)),
                sdD = sd(exp(dummys))
            ))
        })
        
        two_H_beta <-
            unlist(lapply(two_H_beta_es, function(x)
                x[[1]]))
        two_H_beta_sd <-
            unlist(lapply(two_H_beta_es, function(x)
                x[[2]]))
        two_D_beta <-
            unlist(lapply(two_H_beta_es, function(x)
                x[[3]]))
        two_D_beta_sd <-
            unlist(lapply(two_H_beta_es, function(x)
                x[[4]]))
        
        mat_two_H_beta <- matrix(NA, nrow = npops, ncol = npops)
        mat_two_H_beta[lower.tri(mat_two_H_beta)] <- two_H_beta
        mat_two_H_beta[pairs] <- two_H_beta_sd
        colnames(mat_two_H_beta) <-
            rownames(mat_two_H_beta) <- names(pops)
        
        mat_two_D_beta <- matrix(NA, nrow = npops, ncol = npops)
        mat_two_D_beta[lower.tri(mat_two_D_beta)] <- two_D_beta
        mat_two_D_beta[pairs] <- two_D_beta_sd
        colnames(mat_two_D_beta) <-
            rownames(mat_two_D_beta) <- names(pops)
    }
    
    if (pbar) {
        setTxtProgressBar(pb, 8)
        cat(report(" Done.                               "))
    }
    
    # PRINTING OUTPUTS
    
    # spectrumplot
    
    fs <- cbind(zero_D_alpha, one_D_alpha, two_D_alpha)
    colnames(fs) <- c("q=0", "q=1", "q=2")
    sds <-
        cbind(zero_D_alpha_sd, one_D_alpha_sd, two_D_alpha_sd)
    up <- fs + sds
    colnames(up) <- c("up_q0", "up_q1", "up_q2")
    low <- fs - sds
    colnames(low) <- c("low_q0", "low_q1", "low_q2")
    
    fs_plot <- reshape2::melt(fs)
    fs_plot_up <- reshape2::melt(up)
    fs_plot_low <- reshape2::melt(low)
    # avoid no visible bindings
    value <- NULL
    fs_final <-
        as.data.frame(cbind(fs_plot, fs_plot_up[, 3], fs_plot_low[, 3]))
    colnames(fs_final) <- c("pop", "q", "value", "up", "low")
    
    if(plot.out==TRUE){
    
    # printing plots and reports assigning colors to populations
    if (is(plot_colors, "function")) {
        colors_pops <- plot_colors(length(levels(pop(x))))
    }
    
    if (!is(plot_colors, "function")) {
        colors_pops <- plot_colors
    }
    
    
    p3 <-
        ggplot(fs_final, aes(x = pop, y = value, fill = pop)) + 
        geom_bar(position = "dodge", stat = "identity",color = "black") + 
        geom_errorbar(aes(ymin = low, ymax = up), width = 0.2) + 
        scale_fill_manual(values = colors_pops) + 
        facet_wrap(~ q, scales = "free_x") + 
        plot_theme + 
        theme(text = element_text(size = 14),
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank()) +
        labs(fill = "Population") + 
        ggtitle("q-profile")
    
    print(p3)
    }
    
    if (!is.na(match(table, c("H", "DH", "HD")))) {
        tt <-
            data.frame(
                nloci = unlist(nlocpop),
                m_0Ha = zero_H_alpha,
                sd_0Ha = zero_H_alpha_sd,
                m_1Ha = one_H_alpha,
                sd_1Ha = one_H_alpha_sd,
                m_2Ha = two_H_alpha,
                sd_2Ha = two_H_alpha_sd
            )
        print(knitr::kable(tt, digits = 3))
        if (npops > 1) {
            cat("\n\npairwise non-missing loci")
            print(knitr::kable(mat_nloc_pops, digits = 3))
            
            cat("\n\n0_H_beta")
            print(knitr::kable(mat_zero_H_beta, digits = 3))
            cat("\n\n1_H_beta")
            print(knitr::kable(mat_one_H_beta, digits = 3))
            cat("\n\n2_H_beta")
            print(knitr::kable(mat_two_H_beta, digits = 3))
        }
    }
    
    if (!is.na(match(table, c("D", "DH", "HD")))) {
        tt <-
            data.frame(
                nloci = unlist(nlocpop),
                m_0Da = zero_D_alpha,
                sd_0Da = zero_D_alpha_sd,
                m_1Da = one_D_alpha,
                sd_1Da = one_D_alpha_sd,
                m_2Da = two_D_alpha,
                sd_2Da = two_D_alpha_sd
            )
        print(knitr::kable(tt, digits = 3))
        if (npops > 1) {
            cat("\n\npairwise non-missing loci")
            print(knitr::kable(mat_nloc_pops, digits = 3))
            cat("\n\n0_D_beta")
            print(knitr::kable(mat_zero_D_beta, digits = 3))
            cat("\n\n1_D_beta")
            print(knitr::kable(mat_one_D_beta, digits = 3))
            cat("\n\n2_D_beta")
            print(knitr::kable(mat_two_D_beta, digits = 3))
        }
    }
    if (npops > 1) {
        out <-
            list(
                nlocpop = unlist(nlocpop),
                nlocpairpop = mat_nloc_pops,
                zero_H_alpha = zero_H_alpha,
                zero_H_alpha_sd = zero_H_alpha_sd,
                one_H_alpha = one_H_alpha,
                one_H_alpha_sd = one_H_alpha_sd,
                two_H_alpha = two_H_alpha,
                two_H_alpha_sd = two_H_alpha_sd,
                zero_D_alpha = zero_D_alpha,
                zero_D_alpha_sd = zero_D_alpha_sd,
                one_D_alpha = one_D_alpha,
                one_D_alpha_sd = one_D_alpha_sd,
                two_D_alpha = two_D_alpha,
                two_D_alpha_sd = two_D_alpha_sd,
                zero_H_beta = mat_zero_H_beta,
                one_H_beta = mat_one_H_beta,
                two_H_beta = mat_two_H_beta,
                zero_D_beta = mat_zero_D_beta,
                one_D_beta = mat_one_D_beta,
                two_D_beta = mat_two_D_beta
            )
    } else {
        out <-
            list(
                nlocpop = unlist(nlocpop),
                zero_H_alpha = zero_H_alpha,
                zero_H_alpha_sd = zero_H_alpha_sd,
                one_H_alpha = one_H_alpha,
                one_H_alpha_sd = one_H_alpha_sd,
                two_H_alpha = two_H_alpha,
                two_H_alpha_sd = two_H_alpha_sd,
                zero_D_alpha = zero_D_alpha,
                zero_D_alpha_sd = zero_D_alpha_sd,
                one_D_alpha = one_D_alpha,
                one_D_alpha_sd = one_D_alpha_sd,
                two_D_alpha = two_D_alpha,
                two_D_alpha_sd = two_D_alpha_sd
            )
    }
    
    # SAVE INTERMEDIATES TO TEMPDIR
    if (save2tmp & plot.out==TRUE) {
        # creating temp file names
        temp_plot <- tempfile(pattern = "Plot_")
        match_call <-
            paste0(names(match.call()),
                   "_",
                   as.character(match.call()),
                   collapse = "_")
        # saving to tempdir
        saveRDS(list(match_call, p3), file = temp_plot)
        if (verbose >= 2) {
            cat(report("  Saving ggplot(s) to the session tempfile\n"))
            cat(
                report(
                    "  NOTE: Retrieve output files from tempdir using gl.list.reports() and gl.print.reports()\n"
                )
            )
        }
    }
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n\n"))
    }
    
    # RETURN
    
    invisible(out)
}
