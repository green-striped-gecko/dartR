#' @name gl.random.snp
#' @title Randomly changes the allocation of 0's and 2's in a genlight object
#' @description
#' This function samples randomly half of the SNPs and re-codes, in the sampled
#' SNP's, 0's by 2's.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param plot.out Specify if a plot is to be produced [default TRUE].
#' @param save2tmp If TRUE, saves any ggplots to the session temporary directory
#' (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log ; 3, progress and results summary; 5, full report [default NULL,
#' unless specified using gl.set.verbosity].
#'
#' @details
#' DArT calls the most common allele as the reference allele. In a genlight
#' object, homozygous for the reference allele are coded with a '0' and
#' homozygous for the alternative allele are coded with a '2'. This causes some
#' distortions in visuals from time to time.
#'
#' If plot.out = TRUE, two smear plots (pre-randomisation and
#' post-randomisation) are presented using a random subset of individuals (10)
#' and loci (100) to provide an overview of the changes.
#'
#' Resultant ggplots are saved to the session's temporary directory.
#'
#' @return Returns a genlight object with half of the loci re-coded.
#' @author Custodian: Luis Mijangos -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' require("dartR.data")
#' res <- gl.random.snp(platypus.gl[1:5,1:5],verbose = 5)
#'
#' @export

gl.random.snp <- function(x,
                          plot.out = TRUE,
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
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # DO THE JOB
    
    hold <- x
    
    snp_matrix_temp <- as.matrix(x)
    snp_matrix_temp_0 <- snp_matrix_temp == 0
    snp_matrix_temp_2 <- snp_matrix_temp == 2
    
    snp_matrix_temp[snp_matrix_temp_0 == TRUE] <- 2
    snp_matrix_temp[snp_matrix_temp_2 == TRUE] <- 0
    
    random_snps <- sample(1:nLoc(x), nLoc(x) / 2)
    
    snp_matrix <- as.matrix(x)
    snp_matrix[, random_snps] <- snp_matrix_temp[, random_snps]
    
    x@gen <-
        lapply(1:nrow(snp_matrix), function(i)
            new("SNPbin", as.integer(snp_matrix[i, ])))
    
    random_snps <- random_snps[order(random_snps)]
    
    if (verbose == 5) {
        cat(report(paste(
            "The loci that were changed are:",
            paste(random_snps, collapse = ", "),
            "\n"
        )))
    }
    
    if (plot.out) {
        # subsetting objects to provide an overview of the changes
        if (nInd(x) > 10) {
            ind_to_plot <- sample(1:nInd(x), 10)
            x_plot <- x[ind_to_plot, ]
            hold_plot <- hold[ind_to_plot, ]
        } else {
            x_plot <- x
            hold_plot <- hold
        }
        if (nLoc(x_plot) > 100) {
            loc_to_plot <- sample(1:nLoc(x_plot), 100)
            x_plot <- x_plot[, loc_to_plot]
            hold_plot <- hold_plot[, loc_to_plot]
        }
        
        # plot before randomisation
        p1 <-
            gl.smearplot(hold_plot, posi = "none", verbose = 0)
        p1 <-
            p1 + ggtitle("Pre-randomisation") + theme(
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank()
            )
        
        # plot after randomisation
        p2 <- gl.smearplot(x_plot, verbose = 0)
        p2 <- p2 + ggtitle("Post-randomisation")
    }
    
    # PRINTING OUTPUTS
    if (plot.out) {
        # using package patchwork
        p3 <- p1 / p2
        print(p3)
    }
    
    # SAVE INTERMEDIATES TO TEMPDIR
    
    # creating temp file names
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
    }
    
    # ADD TO HISTORY
    nh <- length(x@other$history)
    x@other$history[[nh + 1]] <- match.call()
    
    # FLAG SCRIPT END
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    # RETURN
    invisible(x)
    
}
