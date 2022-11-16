#' @name gl.report.heterozygosity
#' @title Reports observed, expected and unbiased heterozygosities and FIS
#' (inbreeding coefficient) by population or by individual from SNP data
#' @description Calculates the observed, expected and unbiased expected (i.e.
#' corrected for sample size) heterozygosities and FIS (inbreeding coefficient)
#' for each population or the observed heterozygosity for each individual in a
#' genlight object.
#'
#' @param x Name of the genlight object containing the SNP [required].
#' @param method Calculate heterozygosity by population (method='pop') or by
#' individual (method='ind') [default 'pop'].
#' @param n.invariant An estimate of the number of invariant sequence tags used
#' to adjust the heterozygosity rate [default 0].
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot_theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot_colors_pop A color palette for population plots or a list with
#' as many colors as there are populations in the dataset
#' [default discrete_palette].
#' @param plot_colors_ind List of two color names for the borders and fill of
#' the plot by individual [default two_colors].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session
#' temporary directory (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, unless specified using gl.set.verbosity].
#'
#' @details
#' Observed heterozygosity for a population takes the proportion of
#' heterozygous loci for each individual then averages over the individuals in
#' that population. The calculations take into account missing values.
#'
#' Expected heterozygosity for a population takes the expected proportion of
#' heterozygotes, that is, expected under Hardy-Weinberg equilibrium, for each
#' locus, then averages this across the loci for an average estimate for the
#' population.
#'
#' Observed heterozygosity for individuals is calculated as the proportion of
#' loci that are heterozygous for that individual.
#'
#' Finally, the loci that are invariant across all individuals in the dataset
#' (that is, across populations), is typically unknown. This can render
#' estimates of heterozygosity analysis specific, and so it is not valid to
#' compare such estimates across species or even across different analyses. This
#' is a similar problem faced by microsatellites. If you have an estimate of the
#' number of invariant sequence tags (loci) in your data, such as provided by
#' \code{\link{gl.report.secondaries}}, you can specify it with the n.invariant
#' parameter to standardize your estimates of heterozygosity.
#'
#' \strong{NOTE}: It is important to realise that estimation of adjusted
#' heterozygosity requires that secondaries not to be removed.
#'
#' Heterozygosities and FIS (inbreeding coefficient) are calculated by locus
#' within each population using the following equations:
#' \itemize{
#' \item Observed heterozygosity (Ho) = number of homozygotes / n_Ind,
#' where n_Ind is the number of individuals without missing data.
#' \item Observed heterozygosity adjusted (Ho.adj) <- Ho * n_Loc /
#'  (n_Loc + n.invariant),
#' where n_Loc is the number of loci that do not have all missing data  and
#' n.invariant is an estimate of the number of invariant loci to adjust
#' heterozygosity.
#' \item Expected heterozygosity (He) = 1 - (p^2 + q^2),
#' where p is the frequency of the reference allele and q is the frequency of
#' the alternative allele.
#' \item Expected heterozygosity adjusted (He.adj) = He * n_Loc /
#' (n_Loc + n.invariant)
#' \item Unbiased expected heterozygosity (uHe) = He * (2 * n_Ind /
#' (2 * n_Ind - 1))
#' \item Inbreeding coefficient (FIS) = 1 - (mean(Ho) / mean(uHe))
#' }
#'
#'\strong{ Function's output }
#'
#' Output for method='pop' is an ordered barchart of observed heterozygosity,
#' expected heterozygosity and FIS (Inbreeding coefficient) across populations
#' together with a table of mean observed and expected heterozygosities and FIS
#' by population and their respective standard deviations (SD).
#' 
#' In the output, it is also reported by population: the number of loci used to
#'  estimate heterozygosity(nLoc), the number of polymorphic loci (polyLoc), 
#'  the number of monomorphic loci (monoLoc) and loci with all missing data
#'   (all_NALoc).
#'
#' Output for method='ind' is a histogram and a boxplot of heterozygosity across
#' individuals.
#'
#'  Plots and table are saved to the session temporary directory (tempdir)
#'
#'  Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'
#' @return A dataframe containing population labels, heterozygosities, FIS,
#' their standard deviations and sample sizes
#'
#' @author Custodian: Luis Mijangos (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' require("dartR.data")
#' df <- gl.report.heterozygosity(platypus.gl)
#' df <- gl.report.heterozygosity(platypus.gl,method='ind')
#' n.inv <- gl.report.secondaries(platypus.gl)
#' gl.report.heterozygosity(platypus.gl, n.invariant = n.inv[7, 2])
#'
#' @seealso \code{\link{gl.filter.heterozygosity}}
#'
#' @family report functions
#' @export

gl.report.heterozygosity <- function(x,
                                     method = "pop",
                                     n.invariant = 0,
                                     plot.out = TRUE,
                                     plot_theme = theme_dartR(),
                                     plot_colors_pop = discrete_palette,
                                     plot_colors_ind = two_colors,
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
    datatype <-
        utils.check.datatype(x, accept = "SNP", verbose = verbose)
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    if (!(method == "pop" | method == "ind")) {
        cat(
            warn(
                "Warning: Method must either be by population or by individual,
                set to method='pop'\n"
            )
        )
        method <- "pop"
    }
    
    if (n.invariant < 0) {
        cat(warn(
            "Warning: Number of invariant loci must be non-negative, set to 
            zero\n"
        ))
        n.invariant <- 0
        if (verbose == 5) {
            cat(
                report(
                    "  No. of invariant loci can be esimated using 
                    gl.report.secondaries\n"
                )
            )
        }
    }
    
    if (any(grepl(x@other$history, pattern = "gl.filter.secondaries") == TRUE) &
        n.invariant > 0) {
        cat(
            warn(
                "  Warning: Estimation of adjusted heterozygosity requires that 
                secondaries not to be removed. A gl.filter.secondaries call was 
                found in the history. This may cause the results to be 
                incorrect\n"
            )
        )
    }
    
    # DO THE JOB
    
    ########### FOR METHOD BASED ON POPULATIONS
    
    if (method == "pop") {
        # Set a population if none is specified (such as if the genlight object
        # has been generated manually)
        if (is.null(pop(x)) |
            is.na(length(pop(x))) | length(pop(x)) <= 0) {
            if (verbose >= 2) {
                cat(
                    warn(
                        "  No population assignments detected,
                             individuals assigned to a single population
                        labelled 'pop1'\n"
                    )
                )
            }
            pop(x) <- array("pop1", dim = nInd(x))
            pop(x) <- as.factor(pop(x))
        }
        
        # Split the genlight object into a list of populations
        sgl <- seppop(x)
        
        # OBSERVED HETEROZYGOSITY
        if (verbose >= 2) {
            cat(
                report(
                    "  Calculating Observed Heterozygosities, averaged across 
                    loci, for each population\n"
                )
            )
        }
        
        # Calculate heterozygosity for each population in the list
        # CP = Carlo Pacioni CP ###
        Ho.loc <-
            lapply(sgl, function(x)
                colMeans(as.matrix(x) == 1, na.rm = TRUE))
        ##########
        
        Ho <-
            unlist(lapply(sgl, function(x)
                mean(
                    colMeans(as.matrix(x) == 1, na.rm = TRUE), na.rm = TRUE
                )))
        
        ### CP ### observed heterozygosity standard deviation
        HoSD <-
            unlist(lapply(sgl, function(x)
                sd(
                    colMeans(as.matrix(x) == 1, na.rm = TRUE), na.rm = TRUE
                )))
        ##########
        
        # Calculate the number of loci that are not all NAs CP ###
        n_loc <-
            unlist(lapply(sgl, function(x)
                sum(!(
                    colSums(is.na(as.matrix(x))) == nrow(as.matrix(x))
                ))))
        ##########
        
        # calculate the number of polymorphic and monomorphic loci by population
        poly_loc <- NULL
        mono_loc <- NULL
        all_na_loc <- NULL
        
        for(y in 1:length(sgl)){

          y_temp <- sgl[[y]]
          hold <- y_temp
          mono_tmp <- gl.alf(y_temp)
          loc.list <- rownames(mono_tmp[which(mono_tmp$alf1==1 |
                                                mono_tmp$alf1 == 0),])
          loc.list_NA <- rownames(mono_tmp[which(is.na(mono_tmp$alf1)),])

          # Remove NAs from list of monomorphic loci and loci with all NAs
          loc.list <- loc.list[!is.na(loc.list)]

          # remove monomorphic loci and loci with all NAs
          if (length(loc.list) > 0) {
            y_temp <- gl.drop.loc(y_temp, loc.list = loc.list, verbose = 0)
          }

          poly_loc <-  c(poly_loc, nLoc(y_temp))  
          mono_loc <- c(mono_loc, (nLoc(hold) - nLoc(y_temp)))
          all_na_loc <- c(all_na_loc, length(loc.list_NA))

        }
     
        # Apply correction CP ###
        Ho.adj <- Ho * n_loc / (n_loc + n.invariant)
        # Manually compute SD for Ho.adj sum of the square of differences from
        # the mean for polymorphic sites plus sum of the square of differences
        #(which is the Ho.adj because Ho=0) from the mean for invariant sites
        Ho.adjSD <-
            sqrt((
                mapply(function(x, Mean)
                    sum((x - Mean) ^ 2, na.rm = TRUE), Ho.loc, Mean = Ho.adj) +
                  n.invariant * Ho.adj ^
                    2
            ) / (n_loc +
                     n.invariant - 1))
        
        ind.count <- function(x) {
            # the loci that are completely missing
            loci.na <-
                which(colSums(is.na(as.matrix(x))) == nrow(as.matrix(x)))
            # the number of samples in the matrix the number of non-genotyped
            # samples remove the loci that are completely missing
            if (length(loci.na) > 0) {
                nind <-
                    mean(nrow(as.matrix(x)) - 
                           colSums(is.na(as.matrix(x)))[-loci.na])
                # the number of samples in the matrix the number of
                # non-genotyped samples
            } else {
                nind <- mean(nrow(as.matrix(x)) - colSums(is.na(as.matrix(x))))
            }
            
            return(nind)
        }
        
        n_ind <- sapply(sgl, ind.count)
        
        ##########
        
        # EXPECTED HETEROZYGOSITY
        if (verbose >= 2) {
            cat(report("  Calculating Expected Heterozygosities\n\n"))
        }
        
        Hexp <- array(NA, length(sgl))
        Hexp.adj <- array(NA, length(sgl))
        
        ### CP ###
        HexpSD <- array(NA, length(sgl))
        uHexp <- array(NA, length(sgl))
        uHexpSD <- array(NA, length(sgl))
        Hexp.adjSD <- array(NA, length(sgl))
        ##########
        
        FIS <- array(NA, length(sgl))
        FISSD <- array(NA, length(sgl))
        
        # For each population
        for (i in 1:length(sgl)) {
            gl <- sgl[[i]]
            gl <- utils.recalc.freqhomref(gl, verbose = 0)
            gl <- utils.recalc.freqhomsnp(gl, verbose = 0)
            gl <- utils.recalc.freqhets(gl, verbose = 0)
            p <- gl@other$loc.metrics$FreqHomRef
            q <- gl@other$loc.metrics$FreqHomSnp
            hets <- gl@other$loc.metrics$FreqHets
            p <- (2 * p + hets) / 2
            q <- (2 * q + hets) / 2
            H <- 1 - (p ^ 2 + q ^ 2)
            
            ### CP ### Unbiased He (i.e. corrected for sample size) hard
            # coded for diploid
            uH <-
                (2 * as.numeric(n_ind[i]) / (2 * as.numeric(n_ind[i]) - 1)) * H
            ##########
            
            Hexp[i] <- mean(H, na.rm = T)
            ### CP ###
            uHexp[i] <- mean(uH, na.rm = T)
            Hexp.adj[i] <-
                Hexp[i] * n_loc[i] / (n_loc[i] + n.invariant)
            HexpSD[i] <- sd(H, na.rm = T)
            uHexpSD[i] <- sd(uH, na.rm = T)
            Hexp.adjSD[i] <-
                sqrt((sum((
                    H - Hexp.adj[i]
                ) ^ 2, na.rm = TRUE) + n.invariant * Hexp.adj[i] ^ 2) / 
                  (n_loc[i] + n.invariant - 1))
            ##########
            FIS_temp <- 1 - (mean(unlist(Ho.loc[i]),na.rm = T) / 
                               mean(uH,na.rm = T))
            FIS[i] <- mean(FIS_temp, na.rm = T)
            #FISSD[i] <- sd(FIS_temp, na.rm = T)
        }
        
        ### CP ###
        df <-
            data.frame(
                pop = popNames(x),
                nInd = n_ind,
                nLoc = n_loc,
                nLoc.adj = n_loc / (n_loc + n.invariant),
                polyLoc = poly_loc ,
                monoLoc = mono_loc , 
                all_NALoc = all_na_loc, 
                Ho = as.numeric(Ho),
                HoSD = HoSD,
                Ho.adj = as.numeric(Ho.adj),
                Ho.adjSD = Ho.adjSD,
                He = round(Hexp, 6),
                HeSD = round(HexpSD, 6),
                uHe = round(uHexp, 6),
                uHeSD = round(uHexpSD, 6),
                He.adj = round(Hexp.adj, 8),
                He.adjSD = round(Hexp.adjSD, 8),
                FIS = FIS
                # FISSD = FISSD
            )
        ##########
        
        if (plot.out) {
            value <- color <- variable <- He.adj <- NULL
            # printing plots and reports assigning colors to populations
            if (is(plot_colors_pop, "function")) {
                colors_pops <- plot_colors_pop(length(levels(pop(x))))
            }
            
            if (!is(plot_colors_pop,"function")) {
                colors_pops <- plot_colors_pop
            }
            
            if (n.invariant == 0) {
                df.ordered <- df
                df.ordered$color <- colors_pops
                df.ordered <- df.ordered[order(df.ordered$Ho),]
                df.ordered$pop <-
                    factor(df.ordered$pop, levels = df.ordered$pop)
                df.ordered <-
                    df.ordered[, c("pop", "nInd", "Ho", "He", "FIS", "color")]
                df.ordered <-
                    reshape2::melt(df.ordered, id = c("pop", "color", "nInd"))
                
                colors_pops_plot <-
                    df.ordered[, c("pop", "color", "variable", "value")]
                colors_pops_plot <-
                    colors_pops_plot[colors_pops_plot$variable == "Ho", ]
                colors_pops_plot <-
               colors_pops_plot[order(as.character(colors_pops_plot$value)), ]
                
                p3 <-
                    ggplot(df.ordered,
                           aes(x = pop,
                               y = value)) + geom_bar(
                                   position = "dodge2",
                                   stat = "identity",
                                   color = "black",
                                   fill = rep(colors_pops_plot$color, each =
                                                  3)
                               ) +
                    scale_x_discrete(labels = paste(
                        df.ordered$pop,
                        round(df.ordered$nInd,
                              0),
                        sep = " | "
                    )) +
                    geom_text(
                        label = rep(c("Ho", "He", "FIS"), length(unique(
                            df.ordered$pop
                        ))),
                        position = position_dodge2(width = 0.9),
                        stat = "identity",
                        vjust = -0.50,
                        size = 4,
                        inherit.aes = TRUE, fontface = "bold"
                    ) +
                    
                    plot_theme +
                    theme(
                        axis.ticks.x = element_blank(),
                        axis.text.x = element_text(
                            angle = 90,
                            hjust = 1,
                            face = "bold",
                            size = 12
                        ),
                        axis.title.x = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.title.y = element_blank(),
                        legend.position = "none"
                    ) + labs(fill = "Population") +
                    ggtitle("Heterozygosities and FIS by Population")
            } else {
                df.ordered <- df
                df.ordered$color <- colors_pops
                df.ordered <-
                    df.ordered[order(df.ordered$Ho.adj),]
                df.ordered$pop <-
                    factor(df.ordered$pop, levels = df.ordered$pop)
                p1 <-
                    ggplot(df.ordered, aes(
                        x = pop,
                        y = Ho.adj,
                        fill = pop
                    )) + geom_bar(position = "dodge",
                                  stat = "identity",
                                  color = "black") +
                    scale_fill_manual(values = df.ordered$color) + 
                  scale_x_discrete(labels = paste(
                        df.ordered$pop,
                        round(df.ordered$nInd, 0),
                        sep = " | "
                    )) + plot_theme + theme(
                        axis.ticks.x = element_blank(),
                        axis.text.x = element_blank(),
                        axis.title.x = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.title.y = element_blank(),
                        legend.position = "none"
                    ) + 
                  labs(fill = "Population") + 
                  ggtitle("Adjusted Observed Heterozygosity by Population")
                
                p2 <-
                    ggplot(df.ordered, aes(
                        x = pop,
                        y = He.adj,
                        fill = pop
                    )) + geom_bar(position = "dodge",
                                  stat = "identity",
                                  color = "black") +
                    scale_fill_manual(values = df.ordered$color) + 
                  scale_x_discrete(labels = paste(
                        df.ordered$pop,
                        round(df.ordered$nInd, 0),
                        sep = " | "
                    )) + plot_theme + theme(
                        axis.ticks.x = element_blank(),
                        axis.text.x = element_text(
                            angle = 90,
                            hjust = 1,
                            face = "bold",
                            size = 12
                        ),
                        axis.title.x = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.title.y = element_blank(),
                        legend.position = "none"
                    ) +
                    labs(fill = "Population") + 
                  ggtitle("Adjusted Expected Heterozygosity by Population")
                
                p3 <- (p1 / p2)
            }
        }
        
        # OUTPUT REPORT
        if (verbose >= 3) {
cat("  Reporting Heterozygosity by Population\n")
cat("\n  No. of loci =", nLoc(x), "\n")
cat("  No. of individuals =", nInd(x), "\n")
cat("  No. of populations =", nPop(x), "\n")
cat("    Minimum Observed Heterozygosity: ", round(min(df$Ho, na.rm = TRUE), 6))
            if (n.invariant > 0) {
          cat("   [Corrected:", round(min(df$Ho.adj, na.rm = TRUE), 6), "]\n")
            } else {
                cat("\n")
            }
cat("    Maximum Observed Heterozygosity: ", round(max(df$Ho, na.rm = TRUE), 6))
            if (n.invariant > 0) {
        cat("   [Corrected:", round(max(df$Ho.adj, na.rm = TRUE), 6), "]\n")
            } else {
                cat("\n")
            }
            cat("    Average Observed Heterozygosity: ",
                round(mean(df$Ho, na.rm = TRUE), 6))
            if (n.invariant > 0) {
      cat("   [Corrected:", round(mean(df$Ho.adj, na.rm = TRUE), 6), "]\n\n")
            } else {
                cat("\n\n")
            }
            cat("    Minimum Unbiased Expected Heterozygosity: ",
                round(min(df$uHe, na.rm = TRUE), 6))
            if (n.invariant > 0) {
           cat("   [Corrected:", round(min(df$He.adj, na.rm = TRUE), 6), "]\n")
            } else {
                cat("\n")
            }
            cat("    Maximum Unbiased Expected Heterozygosity: ",
                round(max(df$uHe, na.rm = TRUE), 6))
            if (n.invariant > 0) {
         cat("   [Corrected:", round(max(df$He.adj, na.rm = TRUE), 6), "]\n")
            } else {
                cat("\n")
            }
            cat("    Average Unbiased Expected Heterozygosity: ",
                round(mean(df$uHe, na.rm = TRUE), 6))
            if (n.invariant > 0) {
      cat("   [Corrected:", round(mean(df$He.adj, na.rm = TRUE), 6), "]\n\n")
            } else {
                cat("\n\n")
            }
            
            if (n.invariant > 0) {
                cat(
                    "  Average correction factor for invariant loci =",
                    mean(n_loc / (n_loc + n.invariant), na.rm = T),
                    "\n"
                )
            } else {
                cat(
     "  Heterozygosity estimates not corrected for uncalled invariant loci\n"
                )
            }
        }
        
        # PRINTING OUTPUTS
        if (plot.out) {
            suppressWarnings(print(p3))
        }
        if (verbose >= 2) {
            if (n.invariant > 0) {
                print(df)
            } else {
                print(df[, c(
                    "pop",
                    "nInd",
                    "nLoc",
                    "polyLoc",
                    "monoLoc", 
                    "all_NALoc", 
                    "Ho",
                    "HoSD",
                    "He",
                    "HeSD",
                    "uHe",
                    "uHeSD",
                    "FIS"
                    # ,
                    # "FISSD"
                )], row.names = FALSE)
            }
        }
    }
    
    ########### FOR METHOD BASED ON INDIVIDUAL
    
    if (method == "ind") {
        if (verbose >= 2) {
            cat(report(
                "  Calculating observed heterozygosity for individuals\n"
            ))
            cat(report(
         "  Note: No adjustment for invariant loci (n.invariant set to 0)\n"
            ))
        }
        # Convert to matrix
        m <- as.matrix(x)
        
        # For each individual determine counts of hets, homs and NAs
        c.na <- array(NA, nInd(x))
        c.hets <- array(NA, nInd(x))
        c.hom0 <- array(NA, nInd(x))
        c.hom2 <- array(NA, nInd(x))
        for (i in 1:nInd(x)) {
            c.na[i] <- sum(is.na(m[i,]))
            c.hets[i] <-
                sum(m[i,] == 1, na.rm = TRUE) / (nLoc(x) - c.na[i])
            c.hom0[i] <-
                sum(m[i,] == 0, na.rm = TRUE) / (nLoc(x) - c.na[i])
            c.hom2[i] <-
                sum(m[i,] == 2, na.rm = TRUE) / (nLoc(x) - c.na[i])
        }
        
        # Join the sample sizes with the heterozygosities
        df <-
            cbind.data.frame(x@ind.names, c.hets, c.hom0, c.hom2)
        names(df) <-
            c("ind.name", "Ho", "f.hom.ref", "f.hom.alt")
        
        # Boxplot
        if (plot.out) {
            upper <- ceiling(max(df$Ho) * 10) / 10
            p1 <-
                ggplot(df, aes(y = Ho)) + 
        geom_boxplot(color = plot_colors_ind[1], fill = plot_colors_ind[2]) + 
              coord_flip() + 
              plot_theme +
                xlim(range = c(-1, 1)) +
              ylim(0, upper) + 
              ylab(" ") + 
         theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()) +
                ggtitle("Observed Heterozygosity by Individual")
            
            # Histogram
            p2 <-
                ggplot(df, aes(x = Ho)) +
geom_histogram(bins =25,color = plot_colors_ind[1],fill =plot_colors_ind[2]) +
              coord_cartesian(xlim = c(0, upper)) +
              xlab("Observed heterozygosity") + 
              ylab("Count") +
              plot_theme
        }
        
        outliers_temp <-
            ggplot_build(p1)$data[[1]]$outliers[[1]]
        outliers <-
            data.frame(ID = as.character(df$ind.name[df$Ho %in% outliers_temp]), 
                       Ho = outliers_temp)
        
        # OUTPUT REPORT
        if (verbose >= 3) {
            cat("Reporting Heterozygosity by Individual\n")
            cat("No. of loci =", nLoc(x), "\n")
            cat("No. of individuals =", nInd(x), "\n")
            cat("  Minimum Observed Heterozygosity: ",
                round(min(df$Ho), 6),
                "\n")
            cat("  Maximum Observed Heterozygosity: ",
                round(max(df$Ho), 6),
                "\n")
            cat("  Average Observed Heterozygosity: ",
                round(mean(df$Ho), 6),
                "\n\n")
            cat("  Results returned as a dataframe\n\n")
            if (nrow(outliers) == 0) {
                cat("  No outliers detected\n\n")
            } else {
                cat("  Outliers detected\n")
                print(outliers)
                cat("\n")
            }
        }
        
        # PRINTING OUTPUTS
        if (plot.out) {
            p3 <- (p1 / p2) + plot_layout(heights = c(1, 4))
            print(p3)
        }
        if (verbose >= 2) {
            print(df, row.names = FALSE)
        }
    }
    
    # SAVE INTERMEDIATES TO TEMPDIR
    if (save2tmp) {
        # creating temp file names
        if (plot.out) {
            temp_plot <- tempfile(pattern = "Plot_")
        }
        temp_table <- tempfile(pattern = "Table_")
        match_call <-
            paste0(names(match.call()),
                   "_",
                   as.character(match.call()),
                   collapse = "_")
        # saving to tempdir
        if (plot.out) {
            saveRDS(list(match_call, p3), file = temp_plot)
            if (verbose >= 2) {
                cat(report("  Saving the ggplot to session tempfile\n"))
            }
        }
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
    
    if (verbose >= 3) {
        cat(report("  Returning a dataframe with heterozygosity values\n"))
    }
    
    # FLAG SCRIPT END
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    # RETURN
    
    invisible(df)
}
