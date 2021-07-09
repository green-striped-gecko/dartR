#' @name gl.structure
#'
#' @title Run a STRUCTURE analysis using a genlight object
#'
#' @description 
#' This function takes a genlight object and runs a STRUCTURE analysis based on functions from \code{strataG}
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param ... Parameters to specify the STRUCTURE run (check \code{structureRun} for more details). Parameters are passed to the \code{structureRun} function. For example you need to set the k.range and the type of model you would like to run (noadmix, locprior) etc. If those parameter names do not tell you anything, please make sure you familiarize with the STRUCTURE program (Pritchard 2000).
#' @param exec Full path and name+extension where the structure executable is located. E.g. \code{"c:/structure/structure.exe"} under Windows. For Mac and Linux it might be something like \code{"./structure/structure"} if a the executable is in a subfolder 'structure' in your home directory.
#' @param plot create an Evanno plot once finished. Be aware k.range needs to be at least three different k steps.
#' @param plot_theme Theme for the plot. See details for options [default theme_dartR()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default NULL, unless specified using gl.set.verbosity]
#'
#' @details The function is basically a convenient wrapper around the beautiful
#' strataG function \code{structureRun} (Archer et al. 2016). For a detailed
#' description please refer to this package (see references below).
#' To make use of this function you need to download STRUCTURE for you system
#' (\bold{non GUI version}) from here 
#' \href{https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/html/structure.html}{STRUCTURE}.
#'
#' @return an sr object (structure.result list output). Each list entry is a single structure
#'run output (there are k.range * num.k.rep number of runs). For example the
#' summary output of the first run can be accessed via \code{sr[[1]]$summary} 
#' or the q-matrix of the third run via \code{sr[[3]]$q.mat}. To conveniently
#' summarise the outputs across runs (clumpp) you need to run
#' gl.plotstructure on the returned sr object. For evanno plots run gl.evanno on your sr object.
#'
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' \dontrun{
#' #STRUCTURE needs to be installed to be able to run the example
#' #only the first 100 loci
#' #bc <- bandicoot.gl[,1:100]
#' #sr <- gl.runstructure(bc, k.range = 2:5, num.k.rep = 3, exec = "./structure.exe")
#' }
#' @import patchwork
###@importFrom strataG genind2gtypes structureRun
#'
#' @export
#' @seealso \link[strataG]{structureRun}
#' @references 
#' Pritchard, J.K., Stephens, M., Donnelly, P. (2000) Inference of population structure using multilocus genotype data. Genetics 155, 945-959.
#' 
#' Archer, F. I., Adams, P. E. and Schneiders, B. B. (2016) strataG: An R package for manipulating, summarizing and analysing population genetic data. Mol Ecol Resour. doi:10.1111/1755-0998.12559



gl.runstructure    <- function(x, 
                               ...,
                               exec=".",
                               plot=TRUE,
                               plot_theme=theme_dartR(),
                               verbose = NULL) {
  
  #check for strataG (currently not on CRAN)
  if (!require("strataG")) #not already installed?
  {
  ap <- available.packages()  #check CRAN
  
  oncran <- match("strataG", ap)
  if (is.na(oncran)) {
    warning("package strataG needs to be installed. It is currently not on CRAN, hence I try to install it manually via Github using devtools:\n  devtools::install_github('EricArcher/strataG'")
    devtools::install_github('EricArcher/strataG')
  
  }
  } else {
  
  #check that structure is installed 
    structure <-  file.exists(exec) 
    
    if (!structure) stop(paste("Cannot find structure executable in the exex path provided:\n",exec,"\nCheck the help page of ?gl.runstructure on how to download and the exec parameter to locate it."))
  
  # TRAP COMMAND
  funname <- match.call()[[1]]
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity()
  
  # CHECKS DATATYPE 
  datatype <- utils.check.datatype(x, verbose = verbose)
  if (datatype!="SNP") stop("You need to provide a SNP genlight object (ploidy=2)!")

# FUNCTION SPECIFIC ERROR CHECKING

# FLAG SCRIPT START
  
  if (verbose >= 1) {
    if (verbose == 5) {
      cat(report("\n\nStarting", funname, "[ Build =", 
                 build, "]\n\n"))
    } else {
      cat(report("\n\nStarting", funname, "\n\n"))
    }
  }
  
# DO THE JOB
  gg <- strataG::genind2gtypes(gl2gi(x, verbose=0))
  
  sr <- strataG::structureRun(gg, exec=exec,...)
  
  ev <- strataG::evanno(sr)
  pa <- ((ev$plots$mean.ln.k+ev$plots$mean.ln.k)/(ev$plots$ln.ppk+ev$plots$delta.k))+plot_theme
  # PRINTING OUTPUTS
  if (plot) {
    
    
    suppressMessages(print(pa))
  
    }

  #check for "/" in match.call
  mc <- gsub("/", ":",as.character(funname))
  mc <- gsub(":", "-",mc)
  nmc <- gsub("/", "_over_",names(funname))
  nmc <- gsub(":", "-",nmc)
  
  temp_plot <- tempfile(pattern =paste0("dartR_plot",paste0(nmc,"_",mc,collapse = "_")))
  
  
  # SAVE INTERMEDIATES TO TEMPDIR             
    # creating temp file names
  # saving to tempdir
  saveRDS(pa, file = temp_plot)
  if(verbose>=2){
    cat(report("  Saving the plot in ggplot format to the tempfile as",temp_plot,"using saveRDS\n"))
    cat(report("  NOTE: Retrieve output files from tempdir using gl.list.reports() and gl.print.reports()\n"))
    }

# FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("\n\nCompleted:", funname, "\n\n"))
  }
  # RETURN

  return(sr)
}
}
