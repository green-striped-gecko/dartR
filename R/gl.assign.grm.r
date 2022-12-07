#' @name gl.assign.grm
#' @title Population assignment using grm
#' @description
#' This function takes one individual and estimates
#' their probability of coming from individual populations
#' from multilocus genotype frequencies.
#
#' @param x Name of the genlight object containing the SNP data [required].
#' @param unknown Name of the individual to be assigned to a population [required].
# @param inbreeding_par The inbreeding parameter [default 0].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @details
#' This function is a re-implementation of the function multilocus_assignment
#'  from package gstudio.
#'  Description of the method used in this function can be found at:
#' https://dyerlab.github.io/applied_population_genetics/population-assignment.html
#' @return A \code{data.frame} consisting of assignment probabilities for each
#'  population.
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' require("dartR.data")
#' if ((requireNamespace("rrBLUP", quietly = TRUE)) &(requireNamespace("gplots", quietly = TRUE)) ) {
#' res <- gl.assign.grm(platypus.gl,unknown="T27")
#' }
#' @export
#' @import dartR.data
gl.assign.grm <- function(x,
                          unknown,
                          verbose=NULL){
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

# FLAG SCRIPT START
funname <- match.call()[[1]]
utils.flag.start(func = funname,
                 build = "Jody",
                 verbosity = verbose)

# CHECK DATATYPE
datatype <- utils.check.datatype(x, verbose = verbose)


if (unknown %in% indNames(x) == FALSE) {
  stop(error(
    paste("  Individual", unknown, "is not in the genlight object\n")
  ))
}

# DO THE JOB
# Assign the unknown individual to population 'unknown'
vec <- as.vector(pop(x))
vec[indNames(x) == unknown] <- "unknown"
# Note, population containing the unknown has been reduced in size by 1
pop(x) <- as.factor(vec)  

x_grm <- gl.grm(x,plotheatmap = FALSE,verbose=0)
x_grm[upper.tri(x_grm, diag = TRUE)] <- NA
x_columns <- as.data.frame(as.table(x_grm))
x_columns <- x_columns[which(!is.na(x_columns$Freq)), ]

x_columns_2 <- x_columns[which(x_columns$Var1 == unknown | 
                                 x_columns$Var2 == unknown),]

x_columns_2 <- rbind(x_columns_2,data.frame(Var1=x_columns_2$Var2,Var2=x_columns_2$Var1,Freq=x_columns_2$Freq))

x_merge_1a <- data.frame(Var1=indNames(x),pop_name=as.character(pop(x)))
x_merge_1b <- merge(x_merge_1a,x_columns_2,by="Var1")

x_merge_1b$Var2 <- as.character(x_merge_1b$Var2)

x_merge_1b <- x_merge_1b[which(x_merge_1b$pop_name!="unknown"),]

x_merge_1b$pop_name <- as.factor(x_merge_1b$pop_name)

x_split <- split(x_merge_1b,f=x_merge_1b$pop_name)

rel_list <- unlist(lapply(x_split,function(y){
  return(mean(y$Freq))
}))

rel_list <- rel_list[order(rel_list,decreasing = TRUE)]

res <- names(rel_list)[1]

return(res)

}