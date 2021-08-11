#' @name gl.report.heterozygosity
#' @title Reports observed and expected heterozygosity by population or by individual from SNP data
#' @description Calculates the observed and expected heterozygosities for each population
#' or the observed heterozygosity for each individual in a genlight object.
#'
#' @param x Name of the genlight object containing the SNP [required].
#' @param method Calculate heterozygosity by population (method='pop') or by individual (method='ind') [default 'pop']
#' @param n.invariant An estimate of the number of invariant sequence tags used to adjust the heterozygosity rate [default 0]
#' @param plot.out Whether produce a plot of the results [default TRUE].
#' @param plot_theme Theme for the plot. See Details for options [default theme_dartR()].
#' @param plot_colours_pop A color palette for population plots [default discrete_palette].
#' @param plot_colours_ind List of two color names for the borders and fill of the plot by individual [default two_colors].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session temporary directory (tempdir) [default FALSE]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default NULL, unless specified using gl.set.verbosity]
#'
#' @details Observed heterozygosity for a population takes the proportion of heterozygous
#' loci for each individual then averages over the individuals in that population. 
#' The calculations take into account missing values.
#' 
#' Expected heterozygosity for a population takes the expected proportion of
#' heterozygotes, that is, expected under Hardy-Weinberg equilibrium, for each locus, then
#' averages this across the loci for an average estimate for the population.
#' The calculations of expected heterozygosity use the unbiassed estimates of Nei, M. (1987) 
#'
#' Output for method='pop' is an ordered barchart of observed heterozygosity across 
#' populations together with a table of observed and expected heterozygosity by population.
#' 
#' Observed heterozygosity for individuals is calculated as the proportion of loci that
#' are heterozygous for that individual.
#' 
#' Finally, the loci that are invariant across all individuals in the dataset (that is,
#' across populations), is typically unknown. This can render estimates of heterozygosity
#' analysis specific, and so it is not valid to compare such estimates across species
#' or even across different analyses. This is a similar problem faced by microsatellites.
#' If you have an estimate of the number of invariant sequence tags (loci) in your data,
#' such as provided by gl.report.secondaries, you can specify it with the n.invariant
#' parameter to standardize your estimates of heterozygosity.
#' 
#'\strong{ Function's output }
#' Output for method='ind' is a histogram and a boxplot of heterozygosity across individuals.
#' 
#'  Plots and table are saved to the session temporary directory (tempdir) 
#'  
#'  Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'
#' @return A dataframe containing population labels, heterozygosities and sample sizes
#'
#' @author Custodian: Luis Mijangos (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#'df <- gl.report.heterozygosity(testset.gl)
#'df <- gl.report.heterozygosity(testset.gl,method='ind')
#'df <- gl.report.heterozygosity(testset.gl,plot.out=FALSE)
#'
#' @seealso \code{\link{gl.filter.heterozygosity}}
#'  
#' @family reporting functions
#' @references Nei, M. and R. K. Chesser (1983). Estimation of fixation indices and gene diversities. Annals of Human Genetics 47:253-259.
#' @importFrom plyr join
#' @export 

gl.report.heterozygosity <- function(x, 
                                     method = "pop", 
                                     n.invariant = 0, 
                                     plot.out = TRUE,
                                     plot_theme = theme_dartR(), 
                                     plot_colours_pop = discrete_palette,
                                     plot_colours_ind = two_colors,
                                     save2tmp = FALSE,
                                     verbose = NULL) {

# SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
# FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func=funname,build="Jackson",v=verbose)
    
# CHECK DATATYPE 
    datatype <- utils.check.datatype(x,accept="SNP", verbose=verbose)

# FUNCTION SPECIFIC ERROR CHECKING

    if (!(method == "pop" | method == "ind")) {
        cat(warn("Warning: Method must either be by population or by individual, set to method='pop'\n"))
        method <- "pop"
    }

    if (n.invariant < 0) {
        cat(warn("Warning: Number of invariant loci must be non-negative, set to zero\n"))
        n.invariant <- 0
        if(verbose==5){cat(report("  No. of invariant loci can be esimated using gl.report.secondaries\n"))}
    }

# DO THE JOB

    ########### FOR METHOD BASED ON POPULATIONS

    if (method == "pop") {

        # Set a population if none is specified (such as if the genlight object has been generated manually)
        if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
            if (verbose >= 2) {
                cat(warn("  No population assignments detected, 
                             individuals assigned to a single population labelled 'pop1'\n"))
            }
            pop(x) <- array("pop1", dim = nInd(x))
            pop(x) <- as.factor(pop(x))
        }

        # Split the genlight object into a list of populations
        sgl <- seppop(x)

        # OBSERVED HETEROZYGOSITY
        if (verbose >= 2) {
            cat(report("  Calculating Observed Heterozygosities, averaged across loci, for each population\n"))
        }

        # Calculate heterozygosity for each population in the list
        Ho <- unlist(lapply(sgl, function(x) mean(colMeans(as.matrix(x) == 1, na.rm = TRUE), na.rm = TRUE)))

        # Calculate the number of loci used
        nl <- unlist(lapply(sgl, function(x) sum(colSums(is.na(as.matrix(x))) == 0)))

        # Apply correction
        Ho.adj <- Ho * nLoc(x)/(nLoc(x) + n.invariant)

        # Calculate sample sizes =Number of individuals
        sums <- data.frame(lapply(sgl, function(x) mean(colSums(as.matrix(x) == 1, na.rm = TRUE), na.rm = TRUE)))
        n <- t(sums/Ho)
        n <- cbind(row.names(n), n, nl, nLoc(x)/(nLoc(x) + n.invariant))

        # Join the sample sizes with the heteozygosities
        df1 <- data.frame(pop = names(Ho), Ho = as.numeric(Ho), Ho.adj = as.numeric(Ho.adj))
        df2 <- data.frame(n)
        names(df2) <- c("pop", "nInd", "nLoc", "nLoc.adj")
        df <- plyr::join(df1, df2, by = "pop")

        # EXPECTED HETEROZYGOSITY
        if (verbose >= 2) {
            cat(report("  Calculating Expected Heterozygosities\n\n"))
        }

        Hexp <- array(NA, length(sgl))
        Hexp.adj <- array(NA, length(sgl))

        # For each population
        for (i in 1:length(sgl)) {
            gl <- sgl[[i]]
            gl <- utils.recalc.freqhomref(gl, verbose = 0)
            gl <- utils.recalc.freqhomsnp(gl, verbose = 0)
            gl <- utils.recalc.freqhets(gl, verbose = 0)
            p <- gl@other$loc.metrics$FreqHomRef
            q <- gl@other$loc.metrics$FreqHomSnp
            hets <- gl@other$loc.metrics$FreqHets
            p <- (2 * p + hets)/2
            q <- (2 * q + hets)/2
            H <- 1 - (p * p + q * q)
            Hexp[i] <- mean(H, na.rm = T)
            Hexp.adj[i] <- Hexp[i] * nLoc(x)/(nLoc(x) + n.invariant)
        }

        df <- data.frame(popNames(x), as.numeric(table(pop(x))), nl, n.invariant, round(df$Ho, 6), round(df$Ho.adj, 6), round(Hexp, 
            6), round(Hexp.adj, 6))
        He <- He.adj <- NULL
        names(df) <- c("pop", "nInd", "nLoc", "nLoc.inv", "Ho", "Ho.adj", "He", "He.adj")

        # printing plots and reports

        if (is.null(n.invariant)) {
            df.ordered <- df[order(df$Ho), ]
            df.ordered$pop <- factor(df.ordered$pop, levels = df.ordered$pop)
            if(plot.out){
            p1 <- ggplot(df.ordered, aes(x = pop, y = Ho, fill = pop)) + 
                geom_bar(position = "dodge", stat = "identity", color = "black") + 
                scale_fill_manual(values = plot_colours_pop(nPop(x))) + 
                scale_x_discrete(labels = paste(df.ordered$pop, df.ordered$nInd,sep = " | ")) + 
                plot_theme + 
                theme(axis.ticks.x = element_blank(), 
                      axis.text.x = element_text(angle = 90,hjust = 1,face = "bold", size = 12), 
                      axis.title.x = element_blank(), 
                      axis.ticks.y = element_blank(), 
                      axis.title.y = element_blank(),
                      legend.position = "none") +
                labs(fill = "Population") +
                ggtitle("Observed Heterozygosity by Population")

            p2 <- ggplot(df.ordered, aes(x = pop, y = He, fill = pop)) + 
                geom_bar(position = "dodge", stat = "identity", color = "black") + 
                scale_fill_manual(values = plot_colours_pop(nPop(x))) + 
                scale_x_discrete(labels = paste(df.ordered$pop, df.ordered$nInd,sep = " | ")) +
                plot_theme + 
                theme(axis.ticks.x = element_blank(), 
                      axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = 12), 
                      axis.title.x = element_blank(),
                      axis.ticks.y = element_blank(), 
                      axis.title.y = element_blank(),
                      legend.position = "none") + 
                labs(fill = "Population") + 
                ggtitle("Expected Heterozygosity by Population")
            }
        } else {
            df.ordered <- df[order(df$Ho.adj), ]
            df.ordered$pop <- factor(df.ordered$pop, levels = df.ordered$pop)
            if(plot.out){
            p1 <- ggplot(df.ordered, aes(x = pop, y = Ho.adj, fill = pop)) + 
                geom_bar(position = "dodge", stat = "identity", 
                color = "black") + scale_fill_manual(values = plot_colours_pop(nPop(x))) + 
                scale_x_discrete(labels = paste(df.ordered$pop, 
                df.ordered$nInd, sep = " | ")) + 
                plot_theme + 
                theme(axis.ticks.x = element_blank(),
                      axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = 12), 
                      axis.title.x = element_blank(), 
                      axis.ticks.y = element_blank(), 
                      axis.title.y = element_blank(), 
                      legend.position = "none") + 
                labs(fill = "Population") + 
                ggtitle("Observed Heterozygosity by Population")

            p2 <- ggplot(df.ordered, aes(x = pop, y = He.adj, fill = pop)) + 
                geom_bar(position = "dodge", stat = "identity", 
                color = "black") + scale_fill_manual(values = plot_colours_pop(nPop(x))) + 
                scale_x_discrete(labels = paste(df.ordered$pop, df.ordered$nInd, sep = " | ")) + 
                plot_theme + 
                theme(axis.ticks.x = element_blank(),
                      axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = 12), 
                      axis.title.x = element_blank(),
                      axis.ticks.y = element_blank(), 
                      axis.title.y = element_blank(),
                      legend.position = "none") + 
                labs(fill = "Population") + 
                ggtitle("Expected Heterozygosity by Population")
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
            cat("    Average Observed Heterozygosity: ", round(mean(df$Ho, na.rm = TRUE), 6))
            if (n.invariant > 0) {
                cat("   [Corrected:", round(mean(df$Ho.adj, na.rm = TRUE), 6), "]\n\n")
            } else {
                cat("\n\n")
            }
            cat("    Miniumum Expected Heterozygosity: ", round(min(df$He, na.rm = TRUE), 6))
            if (n.invariant > 0) {
                cat("   [Corrected:", round(min(df$He.adj, na.rm = TRUE), 6), "]\n")
            } else {
                cat("\n")
            }
            cat("    Maximum Expected Heterozygosity: ", round(max(df$He, na.rm = TRUE), 6))
            if (n.invariant > 0) {
                cat("   [Corrected:", round(max(df$He.adj, na.rm = TRUE), 6), "]\n")
            } else {
                cat("\n")
            }
            cat("    Average Expected Heterozygosity: ", round(mean(df$He, na.rm = TRUE), 6))
            if (n.invariant > 0) {
                cat("   [Corrected:", round(mean(df$He.adj, na.rm = TRUE), 6), "]\n\n")
            } else {
                cat("\n\n")
            }

            if (n.invariant > 0) {
                cat("  Average correction factor for invariant loci =", nLoc(x)/(nLoc(x) + n.invariant), "\n")
            } else {
                cat("  Heterozygosity estimates not corrected for uncalled invariant loci\n")
            }

            if (verbose >= 3) {
                if (n.invariant > 0) {
                  print(df)
                } else {
                  print(df[, c("pop", "nInd", "nLoc", "Ho", "He")])
                }
            }

        }
        
# PRINTING OUTPUTS
        if(plot.out){
            p3 <- (p1/p2)
            print(p3)
        }
        print(df.ordered)
    }

    ########### FOR METHOD BASED ON INDIVIDUAL

    if (method == "ind") {
        if (verbose >= 2) {
            cat(report("  Calculating observed heterozygosity for individuals\n"))
            cat(report("  Note: No adjustment for invariant loci (n.invariant set to 0)\n"))
        }
        # Convert to matrix
        m <- as.matrix(x)

        # For each individual determine counts of hets, homs and NAs
        c.na <- array(NA, nInd(x))
        c.hets <- array(NA, nInd(x))
        c.hom0 <- array(NA, nInd(x))
        c.hom2 <- array(NA, nInd(x))
        for (i in 1:nInd(x)) {
            c.na[i] <- sum(is.na(m[i, ]))
            c.hets[i] <- sum(m[i, ] == 1, na.rm = TRUE)/(nLoc(x) - c.na[i])
            c.hom0[i] <- sum(m[i, ] == 0, na.rm = TRUE)/(nLoc(x) - c.na[i])
            c.hom2[i] <- sum(m[i, ] == 2, na.rm = TRUE)/(nLoc(x) - c.na[i])
        }

        # Join the sample sizes with the heteozygosities
        df <- cbind.data.frame(x@ind.names, c.hets, c.hom0, c.hom2)
        names(df) <- c("ind.name", "Ho", "f.hom.ref", "f.hom.alt")
        
        # Boxplot
        #if(plot.out){
        upper <- ceiling(max(df$Ho)*10)/10
        p1 <- ggplot(df, aes(y = Ho)) + 
          geom_boxplot(color = plot_colours_ind[1], 
                       fill = plot_colours_ind[2]) + 
          coord_flip() + 
          plot_theme + 
          xlim(range = c(-1,1)) + 
          ylim(0,upper) +
          ylab(" ") + 
          theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()) +
          ggtitle("Heterozygosity by Individual")
        
        # Histogram
        p2 <- ggplot(df,aes(x=Ho)) + 
          geom_histogram(bins = 50, color =plot_colours_ind[1], fill = plot_colours_ind[2]) + 
          coord_cartesian(xlim = c(0, upper)) + 
          xlab("Observed heterozygosity") + 
          ylab("Count") + 
          plot_theme
        
        outliers_temp <- ggplot_build(p1)$data[[1]]$outliers[[1]]
        outliers <- data.frame(ID=as.character(df$ind.name[df$Ho %in% outliers_temp]), Ho=outliers_temp)
        #}
        
         # OUTPUT REPORT
        if (verbose >= 3) {
            cat("Reporting Heterozygosity by Individual\n")
            cat("No. of loci =", nLoc(x), "\n")
            cat("No. of individuals =", nInd(x), "\n")
            cat("  Miniumum Observed Heterozygosity: ", round(min(df$Ho), 6), "\n")
            cat("  Maximum Observed Heterozygosity: ", round(max(df$Ho), 6), "\n")
            cat("  Average Observed Heterozygosity: ", round(mean(df$Ho), 6), "\n\n")
            cat("  Results returned as a dataframe\n\n")
            if (nrow(outliers)==0) {
                cat("  No outliers detected\n\n")
            } else {
                cat("  Outliers detected\n")
                print(outliers)
                cat("\n")
            }
        }
        
        # PRINTING OUTPUTS
        if(plot.out){
            p3 <- (p1/p2) + plot_layout(heights = c(1, 4))
            print(p3)
        }
        print(df)
    }
    
    # SAVE INTERMEDIATES TO TEMPDIR  
    if(save2tmp){
    # creating temp file names
    if(plot.out){
        temp_plot <- tempfile(pattern = "Plot_")
        }
    temp_table <- tempfile(pattern = "Table_")
    match_call <- paste0(names(match.call()),"_",as.character(match.call()),collapse = "_")
    # saving to tempdir
    if(plot.out){
        saveRDS(list(match_call,p3), file = temp_plot)
        if(verbose>=2){
          cat(report("  Saving the ggplot to session tempfile\n"))
        }
    }
    saveRDS(list(match_call,df), file = temp_table)
    if(verbose>=2){
        cat(report("  Saving tabulation to session tempfile\n"))
        cat(report("  NOTE: Retrieve output files from tempdir using gl.list.reports() and gl.print.reports()\n"))
    }
    }
    
    if(verbose>=3){
        cat(report("  Returning a dataframe with heterozygosity values\n"))
    }
    
# FLAG SCRIPT END
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
# RETURN
    
    invisible(df)
}
