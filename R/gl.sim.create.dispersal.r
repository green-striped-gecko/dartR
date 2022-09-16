#' @name gl.sim.create_dispersal
#' @title Creates a dispersal file as input for the function gl.sim.WF.run
#' @description
#' This function writes a csv file called "dispersal_table.csv" which contains
#'  the dispersal variables for each pair of populations to be used as input for
#'   the function \code{\link{gl.sim.WF.run}}.
#'
#' The values of the variables can be modified using the columns
#' "transfer_each_gen" and "number_transfers" of this file.
#'
#' See documentation and tutorial for a complete description of the simulations.
#' These documents can be accessed by typing in the R console:
#' browseVignettes(package="dartR”)
#'
#' @param number_pops Number of populations [required].
#' @param dispersal_type One of: "all_connected", "circle" or "line"
#' [default "all_connected"].
#' @param number_transfers Number of dispersing individuals. This value can be .
#' modified by hand after the file has been created [default 1].
#' @param transfer_each_gen Interval of number of generations in which dispersal
#' occur. This value can be modified by hand after the file has been created
#'  [default 1].
#' @param outpath Path where to save the output file. Use outpath=getwd() or
#' outpath='.' when calling this function to direct output files to your working
#'  directory [default tempdir(), mandated by CRAN].
#' @param outfile File name of the output file [default 'dispersal_table.csv'].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @return A csv file containing the dispersal variables for each pair of
#' populations to be used as input for the function \code{\link{gl.sim.WF.run}}.
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' gl.sim.create_dispersal(number_pops=10)
#' @seealso \code{\link{gl.sim.WF.run}}
#' @family simulation functions
#' @export

gl.sim.create_dispersal <- function(number_pops,
                                    dispersal_type = "all_connected",
                                    number_transfers = 1,
                                    transfer_each_gen = 1,
                                    outpath = tempdir(),
                                    outfile = "dispersal_table.csv",
                                    verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jody",
                   verbosity = verbose)
  
  # DO THE JOB
  pops_vector <- 1:number_pops
  
  if (dispersal_type == "all_connected") {
    dispersal_pairs <-
      as.data.frame(expand.grid(pops_vector, pops_vector))
    dispersal_pairs$same_pop <-
      dispersal_pairs$Var1 == dispersal_pairs$Var2
    dispersal_pairs <-
      dispersal_pairs[which(dispersal_pairs$same_pop == FALSE), ]
    colnames(dispersal_pairs) <-
      c("pop1", "pop2", "same_pop")
  }
  
  if (dispersal_type == "line") {
    dispersal_pairs <- as.data.frame(rbind(cbind(head(pops_vector, -1), pops_vector[-1]),
                                           cbind(pops_vector[-1], head(pops_vector, -1))))
    colnames(dispersal_pairs) <- c("pop1", "pop2")
  }
  
  if (dispersal_type == "circle") {
    dispersal_pairs <- as.data.frame(rbind(cbind(
      pops_vector, c(pops_vector[-1], pops_vector[1])
    ),
    cbind(
      c(pops_vector[-1], pops_vector[1]), pops_vector
    )))
    colnames(dispersal_pairs) <- c("pop1", "pop2")
  }
  
  dispersal_pairs <-
    cbind(dispersal_pairs[, 1:2], number_transfers,transfer_each_gen)
  
  dispersal_pairs <- dispersal_pairs[with(dispersal_pairs, order(pop1, pop2)),]
  
  cat(report(
    "  The dispersal table is saved as: ",
    file.path(outpath, outfile, "\n")
  ))
  
  utils::write.table(
    dispersal_pairs,
    file = file.path(outpath, outfile),
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    sep = ","
  )
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
}