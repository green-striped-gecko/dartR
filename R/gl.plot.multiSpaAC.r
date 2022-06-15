#' @name gl.plot.multiSpaAC
#' @title Plot of multiple Spatial autocorrelation analyses 
#' 
#' @description  This function takes multiple  outputs of gl.spatial.autoCorr 
#' and plots them to enable their graphical comparisons. If upper and lower 
#' estimates are present in the outputs, these are plotted as error bars around
#'  the r estimates.
#'  
#' @details The input is a list. If the list has names, these are used 
#' in the legend. If names are absent, these are internally created assuming that each 
#' element of the list represents different populations. 
#' 
#' @param l.spa The list of multiple gl.spatial.autoCorr outputs 
#' @param grp.var.name The name to be used for the grouping variable
#' @inheritParams gl.spatial.autoCorr
#' @return Returns the ggplot object
#' 
#' Custodian: Carlo Pacioni -- Post to
#'   \url{https://groups.google.com/d/forum/dartr}
#' @examples
#'  # Select one pop only
#' plat_Tent <- gl.keep.ind(platypus.gl, 
#' ind.list = platypus.gl@ind.names[pop(platypus.gl) == "TENTERFIELD"], 
#' mono.rm = TRUE)
#' # Compute a simple distance matrix and reverse it so that correlated values 
#' # indicated more similar individuals as we are used to see plots in GenAleEx
#' gd <- 1 - as.matrix(gl.dist.ind(plat_Tent, method = "Simple", plot.out = FALSE))
#' # Replace the diagonal with zeros
#' diag(gd) <- 0
#' # Compute the (approx) geographical distance matrix
#' ggd <- as.matrix(dist(plat_Tent@other$latlon))
#' 
#' # Spatial autocorrelation
#' spa <- gl.spatial.autoCorr(gd, ggd, bins = 3, reps = 100, 
#' permutation = TRUE, bootstrap = TRUE)
#' 
#'  # Does the same for the second pop
#' plat_Sev <- gl.keep.ind(platypus.gl, 
#' ind.list = platypus.gl@ind.names[pop(platypus.gl) == "SEVERN_ABOVE"], 
#' mono.rm = TRUE)
#' # Compute a simple distance matrix and reverse it so that correlated values 
#' # indicated more similar individuals as we are used to see plots in GenAleEx
#' gd_Sev <- 1 - as.matrix(gl.dist.ind(plat_Sev, method = "Simple", plot.out = FALSE))
#' # Replace the diagonal with zeros
#' diag(gd_Sev) <- 0
#' # Compute the (approx) geographical distance matrix
#' ggd_Sev <- as.matrix(dist(plat_Sev@other$latlon))
#' 
#' # Spatial autocorrelation
#' spa_Sev <- gl.spatial.autoCorr(gd, ggd, bins = 3, reps = 100, 
#' permutation = TRUE, bootstrap = TRUE)
#' 
#' # Plot them together 
#' gl.plot.multiSpaAC(list(TENTERFIELD=spa, SEVERN_ABOVE=spa_Sev))
#' 
#' @rawNamespace import(data.table, except = c(melt,dcast))
#' @import ggplot2
#' @export
#' 
gl.plot.multiSpaAC <- function(l.spa, grp.var.name="Pop", 
                                    plot_theme = theme_classic(),
                                    outpath = tempdir(),
                                    out_file_rootname="ac.multiPlot",
                                    verbose=NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jackson",
                   verbosity = verbose)
  # CHECK DATATYPE
  if(!is.list(l.spa))
    stop(error("  The argument 'l.spa' should be a list of gl.spatial.autoCorr outputs\n"))
  
  if(length(names(l.spa)) == 0) names(l.spa) <- paste0("Pop", seq_len(length(l.spa)))
  
  L.r <- U.r <- NULL #avoid global binding error
  
  # DO THE JOB
  
  spa_multi <- rbindlist(l.spa, use.names = TRUE, fill = TRUE, idcol = grp.var.name)
  
p <- ggplot(spa_multi, aes_string("Bin", "r", col=grp.var.name)) + geom_line() + geom_point() + 
  geom_hline(yintercept=0, col="black") +
  xlab("Distance class") + plot_theme
if(length(spa_multi$L.r) > 0) p <- p + geom_errorbar(aes(ymin=L.r, ymax=U.r)) 
print(p)

# SAVE outputs
fn_base <- file.path(outpath, out_file_rootname)
save(p, file = paste(fn_base, "plot.rda", sep = "_"))
ggsave(filename = paste(fn_base, "plot.png", sep = "_"), plot = p, 
       width = 18, height = 12, units = "cm")

if (verbose >= 2) {
  cat(report("  Saving outputs\n"))
  
}
return(p)

}