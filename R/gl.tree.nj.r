#' Outputs an nj tree to summarize genetic similarity among populations
#'
#' This function is a wrapper for the nj\{ape\} function applied to Euclidian
#' distances calculated from the genlight object.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param outgroup Vector containing the population names that are the outgroups
#'  [deefault NULL].
#' @param type Type of dendrogram "phylogram"|"cladogram"|"fan"|"unrooted"
#'  [default "phylogram"].
#' @param labelsize Size of the labels as a proportion of the graphics default
#'  [default 0.7].
#' @param treefile Name of the file for the tree topology using Newick format 
#' [default NULL].
#' @param verbose Specify the level of verbosity: 0, silent, fatal errors only; 
#' 1, flag function begin and end; 2, progress log; 3, progress and results 
#' summary; 5, full report [default 2].
#' @return A tree file of class phylo.
#' @importFrom stringr str_pad
#' @importFrom ape nj root plot.phylo write.tree
#' @importFrom graphics hist par
#' @export
#' @author Custodian: Arthur Georges (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' # SNP data
#'   gl.tree.nj(testset.gl,type='fan')
#' # Tag P/A data
#'   gl.tree.nj(testset.gs,type='fan')

gl.tree.nj <- function(x,
                       dist.method="euclidean",
                       type = "phylogram",
                       outgroup = NULL,
                       labelsize = 0.7,
                       treefile = NULL,
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
    
<<<<<<< HEAD
    # # Convert gl object to a matrix of allele frequencies, locus by population
    # if (verbose >= 2) {
    #     cat(report(
    #         "  Converting to a matrix of frequencies, locus by populations\n"
    #     ))
    # }
    # t = apply(as.matrix(x), 2, tapply, pop(x), function(e)
    #     mean(e) / 2)
    # # Compute Euclidean distance
    # if (verbose >= 2) {
    #     cat(report("  Computing Euclidean distances\n"))
    # }
    # d <- round(as.matrix(dist(t)), 4)
    # # row.names(d) <- c(paste(row.names(d),' ')) row.names(d) <- substr(row.names(d),1,10)
    # 
    
    d <- gl.dist.pop(x, method=dist.method, output="matrix")
    row.names(d) <- popNames(x)
=======
    # Convert gl object to a matrix of allele frequencies, locus by population
    if (verbose >= 2) {
        cat(report(
            "  Converting to a matrix of frequencies, locus by populations\n"
        ))
    }
    t <-apply(as.matrix(x), 2, tapply, pop(x), function(e)
        mean(e) / 2)
    # Compute Euclidean distance
    if (verbose >= 2) {
        cat(report("  Computing Euclidean distances\n"))
    }
    d <- round(as.matrix(dist(t)), 4)
    # row.names(d) <- c(paste(row.names(d),' ')) row.names(d) <- substr(row.names(d),1,10)
>>>>>>> 06233be3b235c0270e437be3dcab853c38692903
    
    # Plot the distances as an nj tree
    tree <- ape::nj(d)
    if (!is.null(outgroup)) {
        # Function plot.phylo{ape} has the labels all of the same length outgroup <- stringr::str_pad(outgroup, nchar(tree$tip.label[1]),
        # side = c('right'), pad = ' ') # Truncate to 10 characters outgroup <- substr(outgroup,1,10) Root the tree
        tree <- ape::root(tree, outgroup)
        # Plot the tree Save the prior settings for mfrow, oma, mai and pty, and reassign
        op <-
            par(
                mfrow = c(1, 1),
                oma = c(1, 1, 1, 1),
                mai = c(0, 0, 0, 0),
                pty = "m"
            )
        ape::plot.phylo(tree, type = type, cex = labelsize)
    } else {
        # Just plot the tree unrooted
        op <-
            par(
                mfrow = c(1, 1),
                oma = c(1, 1, 1, 1),
                mai = c(0, 0, 0, 0),
                pty = "m"
            )
        ape::plot.phylo(tree, type = type, cex = labelsize)
    }
    
    # Output the tree file
    if (!is.null(treefile)) {
        if (verbose >= 2) {
            cat(report("  Writing the tree topology to", treefile, "\n"))
        }
        ape::write.tree(tree, file = treefile)
    }
    
    # Reset the par options
    par(op)
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(tree)
    
}
