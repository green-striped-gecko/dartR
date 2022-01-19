#' @name gl.select.shapes
#' @title Selects shapes from the base R shape palette and outputs as a vector
#' @description
#' This script draws upon the standard R shape palette to extract a vector of
#' shapes for plotting, where the script that follows has a shape parameter
#' expecting a vector of shapes.
#' @details
#' By default the shape palette will be displayed in full in the graphics window
#' from which shapes can be selected in a subsequent run, and the vector of
#' shapes returned for later use.
#'
#' The select parameter can be used to select shapes from the specified 26
#' shapes available (0-25). For example, select=c(1,1,3) will select shape 1, 1
#' again and 3 to retain in the final vector. This can be useful for fine-tuning
#' shape selection, and matching colors and shapes.
#'
#' @param x Optionally, provide a gl object from which to determine the number
#' of populations [default NULL].
#' @param select Select the shapes to retain in the output vector
#' [default NULL, all shapes shown and returned].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @return A vector with the required number of shapes
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' # SET UP DATASET
#' gl <- testset.gl
#' levels(pop(gl))<-c(rep('Coast',5),rep('Cooper',3),rep('Coast',5),
#' rep('MDB',8),rep('Coast',7),'Em.subglobosa','Em.victoriae')
#' # EXAMPLES
#' shapes <- gl.select.shapes() # Select and display available shapes
#' shapes <- gl.select.shapes(select=c(1,1,1,5,8)) # Select and display a restricted set of shapes
#' shapes <- gl.select.shapes(x=gl,select=c(1,1,1,5,8)) # Select set of shapes and check with no. of pops.
#' @seealso \code{\link{gl.select.colors}}
#' @family Exploration/visualisation functions
#' @export

gl.select.shapes <- function(x = NULL,
                             select = NULL,
                             verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbosity = verbose)
    
    # SCRIPT SPECIFIC ERROR CHECKING
    
    if (is.null(select)) {
        if (verbose >= 2) {
            cat(
                warn(
                    "  Warning: Requires shapes not specified, displaying and returning all available 26 shapes\n"
                )
            )
        }
        select <- c(seq(1:26) - 1)
        nshapes <- 26
    } else {
        if (min(select < 0 | max(select > 25))) {
            stop(error(
                "Fatal Error: specified shapes must be in the range 0-25\n"
            ))
        }
        nshapes <- length(select)
    }
    
    if (!is.null(x)) {
        datatype <- utils.check.datatype(x)
        cat(
            warn(
                "  Specified shapes",
                nshapes,
                "must agree in number with the number of populations",
                nPop(x),
                "in the gl object\n"
            )
        )
        cat(warn("  Setting the number of shapes to number of populations\n"))
        nshapes <- nPop(x)
    }
    
    # DO THE JOB
    
    y <-
        rev(c(rep(1, 6), rep(2, 5), rep(3, 5), rep(4, 5), rep(5, 5)))
    y <- y[1:length(select)]
    x <- c(rep(1:5, 5), 6)
    x <- x[1:length(select)]
    plot(
        x,
        y,
        pch = select,
        cex = 1.5,
        ylim = c(1, 5.5),
        xlim = c(1, 6.5),
        axes = FALSE,
        xlab = "",
        ylab = "",
        bg = "blue"
    )
    text(x, y, labels = select, pos = 3)
    
    if (verbose >= 2) {
        cat(report(
            "  Displaying and returning shapes",
            paste(select, collapse = ", "),
            "\n"
        ))
    }
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(select)
}
