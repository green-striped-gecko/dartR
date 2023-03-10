#' @name gl.filter.sexlinked
#' @title Filters loci that are sex linked
#' @description
#' Alleles unique to the Y or W chromosome and monomorphic on the X chromosomes
#' will appear in the SNP dataset as genotypes that are heterozygotic in all
#' individuals of the heterogametic sex and homozygous in all individuals of the
#' homogametic sex. This function keeps or drops loci with alleles that behave
#' in this way, as putative sex specific SNP markers.
#'
#' @param x Name of the genlight object containing the SNP or presence/absence
#'  (SilicoDArT) data [required].
#' @param sex Factor that defines the sex of individuals. See explanation in
#' details [default NULL].
#' @param filter Either 'keep' to keep sex linked markers only or 'drop' to drop
#' sex linked markers [required].
#' @param read.depth Additional filter option to keep only loci above a certain
#' read.depth. Default to 0, which means read.depth is not taken into account
#' [default 0].
#' @param t.het Tolerance in the heterogametic sex, that is t.het=0.05 means
#'  that 5\% of the heterogametic sex can be homozygous and still be regarded as
#'   consistent with a sex specific marker [default 0.1].
#' @param t.hom Tolerance in the homogametic sex, that is t.hom=0.05 means that
#' 5\% of the homogametic sex can be heterozygous and still be regarded as
#' consistent with a sex specific marker [default 0.1].
#' @param t.pres Tolerance in presence, that is t.pres=0.05 means that a
#' silicodart marker can be present in either of the sexes and still be regarded
#'  as a sex-linked marker [default 0.1].
#' @param plot.out Creates a plot that shows the heterozygosity of males and
#' females at each loci be regarded as consistent with a sex specific marker 
#' [default TRUE].
#' @param plot_theme Theme for the plot. See Details for options
#'  [default theme_dartR()].
#' @param plot_colors List of three color names for the not sex-linked loci, for
#'  the sex-linked loci and for the area in which sex-linked loci appear 
#'  [default three_colors].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, 
#' progress log; 3, progress and results summary; 5, full report 
#' [default NULL, unless specified using gl.set.verbosity].
#'
#' @details
#' Sex of the individuals for which sex is known with certainty can be provided
#' via a factor (equal to the length of the number of individuals) or to be held
#' in the variable \code{x@other$ind.metrics$sex}.
#' Coding is: M for male, F for female, U or NA for unknown/missing.
#' The script abbreviates the entries here to the first character. So, coding of
#' 'Female' and 'Male' works as well. Character are also converted to upper 
#' cases.
#'
#''\strong{ Function's output }
#'
#' This function creates also a plot that shows the heterozygosity of males and
#' females at each loci for SNP data or percentage of present/absent in the case 
#' of SilicoDArT data.
#'
#'  Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'
#' @return The filtered genlight object (filter = 'keep': sex linked loci,
#' filter='drop', everything except sex linked loci).
#'
#' @author Arthur Georges, Bernd Gruber & Floriaan Devloo-Delva (Post to
#'  \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#'   \donttest{
#' out <- gl.filter.sexlinked(testset.gl, filter='drop')
#' }
#' out <- gl.filter.sexlinked(testset.gs, filter='drop')
#'
#' @family filter functions
#'
#' @export

gl.filter.sexlinked <- function(x,
                                sex = NULL,
                                filter = NULL,
                                read.depth = 0,
                                t.het = 0.1,
                                t.hom = 0.1,
                                t.pres = 0.1,
                                plot.out = TRUE,
                                plot_theme = theme_dartR(),
                                plot_colors = three_colors,
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
    
    # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x, verbose = 0)
    if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {
        cat(warn("  Warning: genlight object contains monomorphic loci\n"))
    }
    
    if (is.null(filter)) {
        stop(
            error(
                "Filter option needs to be set to either 'keep' or 'drop'. 
                Please refer to the help pages if in doubt. 
                [?gl.filter.sexlinked]."
            )
        )
    }
    
    # FLAG SCRIPT START
    
    if (verbose >= 1) {
        if (verbose == 5) {
            cat(report("\n\nStarting", funname, "[ Build =", build, "]\n\n"))
        } else {
            cat(report("\n\nStarting", funname, "\n\n"))
        }
    }
    
    # DO THE JOB
    
    # sex should be provided as it is not a default setting, if not provided it 
    #will be searched here: reproducibility
    if (is.null(sex)) {
        sex <- x@other$ind.metrics$sex
    }
    
    if (is.null(sex)) {
        stop(
            error(
                "No definition for the sex of individuals is provided. If not 
                provided via the function is needs to be at 
                gl@other$ind.metrics$sex."
            )
        )
    }
    
    if (length(sex) != nInd(x)) {
        stop(
            error(
                "The number of individuals and the number of entries defining 
                the sex do not match. Check your genlight object and your sex 
                defining column."
            )
        )
    }
    
    sex <- as.character(sex)
    UP <- toupper(sex)
    UP <- substring(UP, 1, 1)
    sex[UP == "F"] <- "F"
    sex[UP == "M"] <- "M"
    sex <- ifelse(sex == "F" | sex == "M", sex, "U")
    sex[is.na(sex)] <- "U"
    
    ########### FOR SNP data
    
    if (datatype == "SNP")
    {
        # Extract the data for the females
        matf <- as.matrix(x[sex == "F", ])
        # For each individual
        f <- array(data = NA, dim = c(ncol(matf), 3))
        for (i in 1:ncol(matf)) {
            # sum of genotypes 0 for homozygous for reference allele, 1 for
          #heterozygous and 2 for homozygous for alternative allele
            for (j in 1:3) {
                dummy <- sum(matf[, i] == (j - 1), na.rm = T)
                if (is.na(dummy))
                    dummy <- 0
                f[i, j] <- dummy
            }
        }
        dff <- data.frame(f)
        row.names(dff) <- locNames(x)
        # genotypes 0 for homozygous for reference allele is in F0, 1 for
        #heterozygous is in F1 and 2 for homozygous for alternative
        # allele is in F2
        colnames(dff) <- c("F0", "F1", "F2")
        
        # Extract the data for the males
        matm <- as.matrix(x[sex == "M", ])
        # For each individual
        m <- array(data = NA, dim = c(ncol(matm), 3))
        for (i in 1:ncol(matm)) {
            for (j in 1:3) {
                dummy <- sum(matm[, i] == (j - 1), na.rm = T)
                if (is.na(dummy))
                    dummy <- 0
                m[i, j] <- dummy
            }
        }
        dfm <- data.frame(m)
        row.names(dfm) <- locNames(x)
        # genotypes 0 for homozygous for reference allele is in M0, 1 for 
        #heterozygous is in M1 and 2 for homozygous for alternative
        # allele is in M2
        colnames(dfm) <- c("M0", "M1", "M2")
        
        # Combine the two files
        
        df <- cbind(dff, dfm)
        
        df$read.depth <- x@other$loc.metrics$rdepth
        
        # Check for hets in all males, homs in all females (XY); ditto for ZW
        sumf <- df$F0 + df$F1 + df$F2
        summ <- df$M0 + df$M1 + df$M2
        # Pull loci that are 100% homozygous for females and 100% heterozygous 
        #for males
        indexxy <-
            ((df$F0 / (sumf) >= (1 - t.hom) |
                  df$F2 / (sumf) >= (1 - t.hom)) &
                 df$M1 / (summ) >= (1 - t.het))
        # when all loci are homozygous for the reference allele e.g. 
        #df$F0/(sumf), the division is NaN. So, all the NaN's are set as
        # TRUE
        indexxy[is.na(indexxy)] <- TRUE
        
        if (sum(indexxy, na.rm = T) > 0) {
            xy <- cbind(locnr = which(indexxy == TRUE), df[indexxy, ])
            # when F0, F1, F2 or M0, M1, M2 are all 0 due to NAs heterozygosity 
            #is NaN. these cases are removed
            xy$fhet <- xy$F1 / (xy$F0 + xy$F1 + xy$F2)
            xy$mhet <- xy$M1 / (xy$M0 + xy$M1 + xy$M2)
            xy <- xy[complete.cases(xy), ]
        }
        
        if (sum(indexxy, na.rm = T) == 0) {
            xy <- data.frame()
        }
        
        # Pull loci that are 100% homozygous for males and 100% heterozygous for 
        #females
        indexzw <-
            ((df$M0 / (summ) >= (1 - t.hom) |
                  df$M2 / (summ) >= (1 - t.hom)) &
                 df$F1 / (sumf) >= (1 - t.het))
        # when all loci are homozygous for the reference allele e.g.
        #df$M0/(summ), the division is NaN. So all the NaN's are set as
        # TRUE
        indexzw[is.na(indexzw)] <- TRUE
        
        if (sum(indexzw, na.rm = T) > 0) {
            zw <- cbind(locnr = which(indexzw == TRUE), df[indexzw, ])
            # when F0, F1, F2 or M0, M1, M2 are all 0 due to NAs heterozygosity
            #is NaN. these cases are removed
            zw$fhet <- zw$F1 / (zw$F0 + zw$F1 + zw$F2)
            zw$mhet <- zw$M1 / (zw$M0 + zw$M1 + zw$M2)
            zw <- zw[complete.cases(zw), ]
        }
        
        if (sum(indexzw, na.rm = T) == 0) {
            zw <- data.frame()
        }
        
        if (verbose > 0) {
            cat("Number of females:", sum(sex == "F"), "\n")
            cat("Number of males:", sum(sex == "M"), "\n")
            cat("Sex ratio females:(males+females):",
                round(sum(sex == "F") / (
                    sum(sex == "F") + sum(sex == "M")
                ), 2),
                "\n")
        }
        if (nrow(zw) == 0 & verbose > 0) {
            cat(
                important(
 "  No sex linked markers consistent with female heterogamety (ZZ/ZW)\n"
                )
            )
        } else {
            if (verbose > 0) {
cat("  Sex linked loci consistent with female heterogamety (ZZ/ZW)\n\n")
cat(
paste(
"    Threshold proportion for homozygotes in the heterozygotic sex (ZW) is",
                        t.hom,
                        "\n"
                    )
                )
                cat(
                    paste(
"    Threshold proportion for heterozygotes in the homozygotic sex (ZZ) is",
                        t.het,
                        "\n\n"
                    )
                )
cat("- locnr is the location of the locus in the input genlight object\n")
cat("- F0 is the number of homozygous loci for the reference allele in 
    females\n")
cat("- F1 is the number of heterozygous loci in females\n")
cat("- F2 is the number of homozygous loci for the alternative allele in 
    females\n")
cat("- M0 is the number of homozygous loci for the reference allele in males\n")
cat("- M1 is the number of heterozygous loci in males\n")
cat("- M2 is the number of homozygous loci for the alternative allele in 
    males\n")
cat("- fhet is heterozygosity in females\n")
cat("- mhet is heterozygosity in males\n\n")
print(zw)
cat(
important(
"  Note: The most reliable putative markers will have a read depth > 10.\n"
                    )
                )
            }
        }
        
        if (nrow(xy) == 0 & verbose > 0) {
            cat(
                important(
"  No sex linked markers consistent with male heterogamety (XX/XY)\n"
                )
            )
        } else {
            if (verbose > 0) {
cat(" Sex linked loci consistent with male heterogamety (XX/XY)\n\n")
                cat(
                    paste(
"    Threshold proportion for homozygotes in the heterozygotic sex (XY) is",
                        t.hom,
                        "\n"
                    )
                )
                cat(
                    paste(
"    Threshold proportion for heterozygotes in the homozygotic sex (XX) is",
                        t.het,
                        "\n"
                    )
                )
cat("- locnr is the location of the locus in the input genlight object\n")
cat("- F0 is the number of homozygous loci for the reference allele in 
    females\n")
cat("- F1 is the number of heterozygous loci in females\n")
cat("- F2 is the number of homozygous loci for the alternative allele in 
    females\n")
cat("- M0 is the number of homozygous loci for the reference allele in males\n")
cat("- M1 is the number of heterozygous loci in males\n")
cat("- M2 is the number of homozygous loci for the alternative allele in 
    males\n")
cat("- fhet is heterozygosity in females\n")
cat("- mhet is heterozygosity in males\n")
                print(xy)
                cat(
                    important(
"  \nNote: The most reliable putative markers will have a read depth > 10.\n\n"
                    )
                )
            }
        }
        
        if (plot.out) {
            df$fhet <- dff$F1 / (dff$F0 + dff$F1 + dff$F2)
            df$mhet <- dfm$M1 / (dfm$M0 + dfm$M1 + dfm$M2)
            df$xy <- indexxy
            df$zw <- indexzw
            df$test <- df$xy + df$zw
            df_sex_linked <- df[which(df$test == 1), ]
            df_no_sex_linked <- df[which(df$test == 0), ]
            
            gg <-
                ggplot() + geom_rect(
                    aes(
                        xmin = 0,
                        xmax = t.hom,
                        ymin = 1 - t.het,
                        ymax = 1
                    ),
                    fill = three_colors[3],
                    alpha = 1 / 2,
                    color = "black"
                ) +
                geom_text(x = 0, y = 1.03, aes(label = "XX/XY")) + geom_rect(
                    aes(
                        xmin = 1,
                        xmax = 1 - t.het,
                        ymin = 0,
                        ymax = t.hom
                    ),
                    fill = three_colors[3],
                    alpha = 1 / 2,
                    color = "black"
                ) + geom_text(x = 1,
                              y = -0.02,
                              aes(label = "ZZ/ZW")) + geom_point(
                                  data = df_no_sex_linked,
                                  aes(x = fhet,
                                      y = mhet),
                                  alpha = 1 / 3,
                                  size = 2,
                                  color = three_colors[1]
                              ) + geom_point(
                                  data = df_sex_linked,
                                  aes(x = fhet, y = mhet),
                                  alpha = 1 / 3,
                                  size = 3,
                                  color = three_colors[2]
                              ) + xlab("Female Heterozygosity") + 
              ylab("Male Heterozygosity") + xlim(0, 1) + ylim(0, 1) +
                plot_theme
            
            suppressWarnings(print(gg))
        }
        
    }  #end if datatype='SNP'
    
    ########### FOR SilicoDArT data
    
    if (datatype == "SilicoDArT") {
        matf <- as.matrix(x[sex == "F"])
        # For each individual
        f <- array(data = NA, dim = c(ncol(matf), 2))
        for (i in 1:ncol(matf)) {
            for (j in 1:2) {
                dummy <- sum(matf[, i] == (j - 1), na.rm = T)
                if (is.na(dummy))
                    dummy <- 0
                f[i, j] <- dummy
            }
        }
        dff <- data.frame(f)
        row.names(dff) <- locNames(x)
        colnames(dff) <- c("F0", "F1")
        
        # Extract the data for the males
        matm <- as.matrix(x[sex == "M"])
        # For each individual
        m <- array(data = NA, dim = c(ncol(matm), 2))
        for (i in 1:ncol(matm)) {
            for (j in 1:2) {
                dummy <- sum(matm[, i] == (j - 1), na.rm = T)
                if (is.na(dummy))
                    dummy <- 0
                m[i, j] <- dummy
            }
        }
        dfm <- data.frame(m)
        row.names(dfm) <- locNames(x)
        colnames(dfm) <- c("M0", "M1")
        
        # Combine the two files
        
        df <- cbind(dff, dfm)
        
        df$read.depth <- x@other$loc.metrics$AvgReadDepth
        
        # Check for hets in all males, homs in all females (XY); ditto for ZW
        sumf <- df$F0 + df$F1
        summ <- df$M0 + df$M1
        # Pull loci that are 100% present in females and 0% in males
        indexzw <-
            (df$F1 / (sumf) >= (1 - t.pres) &
                 df$M0 / (summ) >= (1 - t.pres))
        # when all loci are homozygous for the reference allele e.g. 
        #df$M0/(summ), the division is NaN. So all the NAN's are set as TRUE
        indexzw[is.na(indexzw)] <- TRUE
        
        if (sum(indexzw, na.rm = T) > 0) {
            zw <- cbind(locnr = which(indexzw == TRUE), df[indexzw, ])
        }
        
        if (sum(indexzw, na.rm = T) == 0) {
            zw <- data.frame()
        }
        
        # Pull loci that are 100% present in males and 0% in females
        indexxy <-
            (df$M1 / (summ) >= (1 - t.pres) &
                 df$F0 / (sumf) >= (1 - t.pres))
        # when all loci are homozygous for the reference allele e.g. 
        #df$F0/(sumf), the division is NaN. So, all the NaN's are set as TRUE
        indexxy[is.na(indexxy)] <- TRUE
        
        if (sum(indexxy, na.rm = T) > 0) {
            xy <- cbind(locnr = which(indexxy == TRUE), df[indexxy, ])
        }
        
        if (sum(indexxy, na.rm = T) == 0) {
            xy <- data.frame()
        }
        
        if (nrow(zw) == 0 & verbose > 0) {
            cat(
                important(
 "  No sex linked markers consistent with female heterogamety (ZZ/ZW)\n"
                )
            )
        } else {
            if (verbose > 0) {
cat("\n  Sex linked loci consistent with female heterogamety (ZZ/ZW)\n")
                cat(paste(
                    "    Threshold proportion for presence/absence is",
                    t.pres,
                    "\n"
                ))
 cat("     - locnr is the location of the locus in the input genlight object\n")
                cat("     - F0 is the number of loci absent in females\n")
                cat("     - F1 is the number of loci present in females\n")
                cat("     - M0 is the number of loci absent in males\n")
                cat("     - M1 is the number of loci present in males\n")
                print(zw)
                cat(
                    important(
 "  Note: The most reliable putative markers will have a read depth > 10.\n"
                    )
                )
            }
        }
        if (nrow(xy) == 0) {
            if (verbose > 0)
                cat(
                    important(
"  No sex linked markers consistent with male heterogamety (XX/XY)\n"
                    )
                )
        } else {
            if (verbose > 0) {
cat("  Sex linked loci consistent with male heterogamety (XX/XY)\n")
                cat(paste(
                    "    Threshold proportion for presence/absence is",
                    t.pres,
                    "\n"
                ))
cat("     - locnr is the location of the locus in the input genlight object\n")
                cat("     - F0 is the number of loci absent in females\n")
                cat("     - F1 is the number of loci present in females\n")
                cat("     - M0 is the number of loci absent in males\n")
                cat("     - M1 is the number of loci present in males\n")
                print(xy)
                cat(
                    important(
"  Note: The most reliable putative markers will have a read depth > 10.\n"
                    )
                )
            }
        }
        
        if (plot.out) {
            fhet <- mhet <- NULL
            df$fhet <- dff$F1 / (dff$F0 + dff$F1)
            df$mhet <- dfm$M1 / (dfm$M0 + dfm$M1)
            df$xy <- indexxy
            df$zw <- indexzw
            df$test <- df$xy + df$zw
            df_sex_linked <- df[which(df$test == 1), ]
            df_no_sex_linked <- df[which(df$test == 0), ]
            
            gg <-
                ggplot() + geom_rect(
                    aes(
                        xmin = 0,
                        xmax = t.pres,
                        ymin = 1 - t.pres,
                        ymax = 1
                    ),
                    fill = three_colors[3],
                    alpha = 1 / 2,
                    color = "black"
                ) +
                geom_text(x = 0, y = 1.03, aes(label = "XX/XY")) + geom_rect(
                    aes(
                        xmin = 1,
                        xmax = 1 - t.pres,
                        ymin = 0,
                        ymax = t.pres
                    ),
                    fill = three_colors[3],
                    alpha = 1 / 2,
                    color = "black"
                ) + geom_text(x = 1,
                              y = -0.02,
                              aes(label = "ZZ/ZW")) + geom_point(
                                  data = df_no_sex_linked,
                                  aes(x = fhet,
                                      y = mhet),
                                  alpha = 1 / 3,
                                  size = 2,
                                  color = three_colors[1]
                              ) + geom_point(
                                  data = df_sex_linked,
                                  aes(x = fhet, y = mhet),
                                  alpha = 1 / 3,
                                  size = 3,
                                  color = three_colors[2]
                              ) + xlab("% present in females") + 
              ylab("% present in males") + xlim(0, 1) + ylim(0, 1) + plot_theme
            
            suppressWarnings(print(gg))
        }
        
    }
    
    # finally filter
    index <- NULL
    
    if (nrow(xy) > 0) {
        index <- xy$locnr[xy$read.depth >= read.depth]
    }
    
    if (nrow(zw) > 0) {
        index2 <- zw$locnr[zw$read.depth >= read.depth]
        if (length(index) > 0) {
            index <- unique(c(index, index2))
        } else {
            index <- index2
        }
    }
    
    
    if (length(index) > 0) {
        if (filter == "drop") {
            index <- -index
        }

        x2 <- x[, index]
        x2@other$loc.metrics <- x@other$loc.metrics[index, ]
      
    }
    
    # case no
    if (length(index) == 0 & filter == "keep") {
        x2 <- NULL
        if (verbose > 0) {
            cat(
                important(
"  No sex-linked loci identified and filter option was 'keep', therefore NULL is
returned.\n"
                )
            )
        }
    }
    
    if (length(index) == 0 & filter == "drop") {
        if (verbose > 0) {
            cat(
                important(
 "  No sex-linked loci identified and filter option was 'drop', therefore the 
 genlight object is returned unchanged.\n"
                )
            )
        }
      x2 <- x
    }
    
    # ADD TO HISTORY
    if (!is.null(x2)) {
        nh <- length(x2@other$history)
        x2@other$history[[nh + 1]] <- match.call()
    }
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    # RETURN
    
    invisible(x2)
}
