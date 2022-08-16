#' Replays the history and applies it to a genlight object
#'
#' @param x A genlight object (with a history slot) [optional].
#' @param history  If no history is provided the complete history of
#' x is used (recreating the identical object x). If history is a vector it
#' indicates which which part of the history of x is used [\code{c(1,3,4)} uses
#' the first, third and forth entry from \code{x@other$history}]. Or a simple
#' link to a history slot of another genlight object (e.g.
#' code{x2@other$history[c(1,4,5)]}). [optional].
#' @param verbose If set to one then history commands are printed,
#'  which may facilitate reading the output [default 0].
#' @return Returns a genlight object that was created by replaying the provided
#'  applied to the genlight object x. Please note you can 'mix' histories or
#'  part of them and apply them to different genlight objects. If the history
#'  does not contain \code{gl.read.dart}, histories of x and history are
#'  concatenated.
#' @export
#' @importFrom gridExtra grid.table ttheme_default
#' @author Bernd Gruber (bugs? Post to
#' \url{https://groups.google.com/d/forum/dartr}).
#' @details This function basically allows to create a 'template history'
#' (=set of filters) and apply them to any other genlight object. Histories can
#' also be saved and loaded (see. gl.save.history and gl.load.history).
#' @examples
#'\dontrun{
#' dartfile <- system.file('extdata','testset_SNPs_2Row.csv', package='dartR')
#' metadata <- system.file('extdata','testset_metadata.csv', package='dartR')
#' gl <- gl.read.dart(dartfile, ind.metafile = metadata, probar=FALSE)
#' gl2 <- gl.filter.callrate(gl, method='loc', threshold=0.9)
#' gl3 <- gl.filter.callrate(gl2, method='ind', threshold=0.95)
#' #Now 'replay' part of the history 'onto' another genlight object
#' #bc.fil <- gl.play.history(gl.compliance.check(bandicoot.gl),
#' #history=gl3@other$history[c(2,3)], verbose=1)
#' #gl.print.history(bc.fil)
#' }

gl.play.history <- function(x,
                            history = NULL,
                            verbose = 0) {
    if (is.null(history)) {
        hist2 <- x@other$history
    }
    
    if (is.numeric(history)) {
        hist2 <- x@other$history[history]
    }
    
    if (is.list(history)) {
        hist2 <- history
    }
    
    for (i in 1:length(hist2)) {
        glhist <- hist2[[i]]
        narg <-length(glhist)
        ll <- list()
        ll[1:(narg - 1)] <- glhist[2:narg]
        names(ll) <- names(glhist[2:narg])
        
        # check if gl.read.dart or history on different gl
        if (as.character(glhist[[1]]) != "gl.read.dart") {
            if (i == 1) {
                ll[[1]] <- x
            } else {
                ll[[1]] <- gout
            }
        }
        
        # run history one by one
        if (verbose > 0) {
            cat(report("\n################################\n"))
            cat(report("###########Running #############\n"))
            print(glhist)
            cat(report("--------------------------------\n"))
        }
        gout <- do.call(as.character(glhist)[1], ll)
        if (verbose > 0) {
            cat(report("\n###############################\n"))
        }
        flush.console()
    }
    
    # combine or reset histories
    if (as.character(hist2[[1]][[1]]) == "gl.read.dart") {
        gout@other$history <- hist2
    } else {
        gout@other$history <- c(x@other$history, hist2)
    }
    gout
}
