#' @name gl.run.structure
#'
#' @title Runs a STRUCTURE analysis using a genlight object
#'
#' @description
#' This function takes a genlight object and runs a STRUCTURE analysis based on
#' functions from \code{strataG}
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param ... Parameters to specify the STRUCTURE run (check \code{structureRun}
#'  within strataG.
#' for more details). Parameters are passed to the \code{structureRun} function.
#' For example you need to set the k.range and the type of model you would like
#' to run (noadmix, locprior) etc. If those parameter names do not tell you
#' anything, please make sure you familiarize with the STRUCTURE program
#' (Pritchard 2000).
#' @param exec Full path and name+extension where the structure executable is
#' located. E.g. \code{'c:/structure/structure.exe'} under Windows. For Mac and
#' Linux it might be something like \code{'./structure/structure'} if the
#' executable is in a subfolder 'structure' in your home directory
#' [default working directory "."].
#' @param plot.out Create an Evanno plot once finished. Be aware k.range needs
#' to be at least three different k steps [default TRUE].
#' @param plot_theme Theme for the plot. See details for options
#' [default theme_dartR()].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session
#' temporary directory (tempdir) [default FALSE].
#' @param verbose Set verbosity for this function (though structure output
#' cannot be switched off currently) [default NULL]
#' @details The function is basically a convenient wrapper around the beautiful
#' strataG function \code{structureRun} (Archer et al. 2016). For a detailed
#' description please refer to this package (see references below).
#' To make use of this function you need to download STRUCTURE for you system
#' (\bold{non GUI version}) from here
#' \href{https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/html/structure.html}{STRUCTURE}.
#' 
#' \bold{Format note}
#' 
#' For this function to work, make sure that individual and population names 
#' have no spaces. To substitute spaces by underscores you could use the R 
#' function \code{gsub} as below.
#' 
#' \code{popNames(gl) <- gsub(" ","_",popNames(gl))}
#' 
#' \code{indNames(gl) <- gsub(" ","_",indNames(gl)) }
#' 
#' It's also worth noting that Structure truncates individual names at 11 
#' characters. The function will fail if the names of individuals are not unique
#'  after truncation. To avoid this possible problem, a number sequence, as 
#'  shown in the code below, might be used instead of individual names.
#' \code{
#' indNames(gl) <- as.character(1:length(indNames(gl)))
#'}
#' @return An sr object (structure.result list output). Each list entry is a
#' single structurerun output (there are k.range * num.k.rep number of runs).
#' For example the summary output of the first run can be accessed via
#' \code{sr[[1]]$summary} or the q-matrix of the third run via
#' \code{sr[[3]]$q.mat}. To conveniently summarise the outputs across runs
#' (clumpp) you need to run gl.plot.structure on the returned sr object. For
#' Evanno plots run gl.evanno on your sr object.
#'
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' \dontrun{
#' #bc <- bandicoot.gl[,1:100]
#' #sr <- gl.run.structure(bc, k.range = 2:5, num.k.rep = 3, 
#' # exec = './structure.exe')
#' #ev <- gl.evanno(sr)
#' #ev
#' #qmat <- gl.plot.structure(sr, K=3)
#' #head(qmat)
#' #gl.map.structure(qmat, bc, scalex=1, scaley=0.5)
#' }
#' @import patchwork
### @importFrom strataG genind2gtypes structureRun
#' @importFrom dplyr bind_rows mutate_at vars starts_with mutate group_by 
#' ungroup arrange n rename select everything n_distinct bind_rows starts_with
#' @export
### @seealso \code{structureRun}
#' @references
#' \itemize{
#' \item Pritchard, J.K., Stephens, M., Donnelly, P. (2000) Inference of
#' population structure using multilocus genotype data. Genetics 155, 945-959.
#' \item Archer, F. I., Adams, P. E. and Schneiders, B. B. (2016) strataG: An R
#' package for manipulating, summarizing and analysing population genetic data.
#'  Mol Ecol Resour. doi:10.1111/1755-0998.12559
#' }

gl.run.structure <- function(x,
                             ...,
                             exec = ".",
                             plot.out = TRUE,
                             plot_theme = theme_dartR(),
                             save2tmp = FALSE,
                             verbose = NULL) {
    
    pkg <- "purrr"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      cat(error(
        "Package",
        pkg,
        " needed for this function to work. Please install it.\n"
      ))
      return(-1)
    }
        # check that Structure is installed
        structure <- file.exists(exec)
        
        if (!structure) {
            stop(error(
                paste(
"Cannot find Structure executable in the exex path provided:\n",
                    exec,
"\nCheck the help page of ?gl.run.structure on how to download and the exec
parameter to locate it."
                )
            ))
        }
        # SET VERBOSITY
        verbose <- gl.check.verbosity(verbose)
        
        # FLAG SCRIPT START
        funname <- match.call()[[1]]
        utils.flag.start(func = funname,
                         build = "Jody",
                         verbose = verbose)
        
        # CHECK DATATYPE
        datatype <- utils.check.datatype(x, verbose = verbose)
        
        if (datatype != "SNP") {
            stop(error(
                "You need to provide a SNP genlight object (ploidy=2)!"
            ))
        }
        
        # DO THE JOB
        gg <- utils.structure.genind2gtypes(gl2gi(x, verbose = 0))
        
        sr <- utils.structure.run(gg, exec = exec, ...)
        
        ev <- utils.structure.evanno(sr)
        
        pa <- ((ev$plots$mean.ln.k + ev$plots$mean.ln.k) / 
                 (ev$plots$ln.ppk + ev$plots$delta.k)) + plot_theme
        
        # PRINTING OUTPUTS
        if (plot.out) {
            suppressMessages(print(pa))
        }
        
        # SAVE INTERMEDIATES TO TEMPDIR
        if (save2tmp & plot.out) {
            # check for '/' in match.call
            mc <- gsub("/", ":", as.character(funname))
            mc <- gsub(":", "-", mc)
            nmc <- gsub("/", "_over_", names(funname))
            nmc <- gsub(":", "-", nmc)
            
            # creating temp file names
            temp_plot <-
                tempfile(pattern = paste0("Plot", paste0(nmc, "_", mc,
                                                         collapse = "_")))
            
            # saving to tempdir
            saveRDS(pa, file = temp_plot)
            if (verbose >= 2) {
                cat(
                    report(
                        "  Saving the plot in ggplot format to the tempfile as",
                        temp_plot,
                        "using saveRDS\n"
                    )
                )
                cat(
                    report(
                        "  NOTE: Retrieve output files from tempdir using 
                        gl.list.reports() and gl.print.reports()\n"
                    )
                )
            }
        }
        
        # FLAG SCRIPT END
        if (verbose >= 1) {
            cat(report("Completed:", funname, "\n\n"))
        }
        
        # RETURN
        return(sr)
    
}
