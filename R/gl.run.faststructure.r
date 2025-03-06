#' @name gl.run.faststructure
#'
#' @title Runs a faststructure analysis using a genlight object
#'
#' @description
#' This function takes a genlight object and runs a faststructure analysis.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param k.range Range of the number of populations [required].
#' @param num.k.rep Number of replicates [required].
#' @param exec Full path and name+extension where the fastStructure executable
#' is located [default working directory "./fastStructure"].
#' @param output Path to output file [default getwd()].
#' @param tol Convergence criterion [default 10e-6].
#' @param prior Choice of prior: simple or logistic [default "simple"].
#' @param cv Number of test sets for cross-validation, 0 implies no CV step
#'  [default 0].
#' @param seed Seed for random number generator [default NULL].
#' @details
#' Download faststructure binary for your system from here (only runs on Mac or 
#' Linux):
#' 
#' https://github.com/StuntsPT/Structure_threader/tree/master/structure_threader/bins
#' 
#' Move faststructure file to working directory. Make file executable using 
#' terminal app.
#' 
#' \code{system(paste0("chmod u+x ",getwd(), "/faststructure"))}
#' 
#' Download plink binary for your system from here:
#' 
#' https://www.cog-genomics.org/plink/
#' 
#' Move plink file to working directory. Make file executable using 
#' terminal app.
#' 
#' \code{system(paste0("chmod u+x ",getwd(), "/plink"))}
#' 
#' To install fastStructure dependencies follow these directions:
#' https://github.com/rajanil/fastStructure
#'
#' fastStructure performs inference for the simplest, independent-loci,
#' admixture model, with two choices of priors that can be specified using
#' the --prior parameter. Thus, unlike Structure, fastStructure does not require
#' the mainparams and extraparam files. The inference algorithm used by
#'  fastStructure is fundamentally different from that of Structure and
#'  requires the setting of far fewer options.
#'
#'  To identify the number of populations that best approximates the marginal
#'  likelihood of the data, the marginal likelihood is extracted from each run
#'  of K, averaged across replications and plotted.
#'
#' @return A list in which each list entry is a single faststructure run output
#' (there are k.range * num.k.rep number of runs).
#'
#' @author Luis Mijangos (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' \dontrun{
#' t1 <- gl.filter.callrate(platypus.gl,threshold = 1)
#' res <- gl.run.faststructure(t1, exec = "./fastStructure",k.range = 2:3, 
#'                           num.k.rep = 2,output = paste0(getwd(),"/res_str"))
#' qmat <- gl.plot.faststructure(res,k.range=2:3)
#' gl.map.structure(qmat, K=2, t1, scalex=1, scaley=0.5)
#' }
#' @export
#' @references
#' \itemize{
#' \item Raj, A., Stephens, M., & Pritchard, J. K. (2014). fastSTRUCTURE:
#' variational inference of population structure in large SNP data sets.
#' Genetics, 197(2), 573-589.
#' }

gl.run.faststructure <- function(x,
                                 k.range,
                                 num.k.rep,
                                 exec =  "./fastStructure",
                                 output = getwd(),
                                 tol = 10e-6,
                                 prior = "simple",
                                 cv = 0,
                                 seed = NULL) {
  gl2plink(x,
           bed_file = TRUE,
           outpath = output,
           verbose = 0)
  
  for (k_n in k.range) {
    for (rep_n in 1:num.k.rep) {
      if (is.null(seed)) {
        system(
          paste0(
            exec,
            " -K ",
            k_n,
            " --input=",
            output,
            "/gl_plink",
            " --output=",
            output,
            "/genotypes_output",
            "_",
            rep_n,
            " --tol=",
            tol,
            " --prior=",
            prior,
            " --cv=",
            cv
          )
        )
      } else{
        system(
          paste0(
            exec,
            " -K ",
            k_n,
            " --input=",
            output,
            "/gl_plink",
            " --output=",
            output,
            "/genotypes_output",
            "_",
            rep_n,
            " --tol=",
            tol,
            " --prior=",
            prior,
            " --cv=",
            cv,
            " --seed=",
            seed
          )
        )
      }
    }
  }
  
  files_structure <- list.files(path = output,
                                pattern =  "^genotypes_output.+log")
  files_structure_2 <- paste0(output, "/", files_structure)
  n_first_line <- 23
  df_likelihood <-
    as.data.frame(matrix(nrow = dplyr::last(k.range), ncol = num.k.rep))
  for (i in 1:length(files_structure_2)) {
    file_name <- files_structure[i]
    file_name <-
      strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\2', file_name), ' ')
    file_name <- strsplit(file_name[[1]], "\\.")
    k_replicate <- as.numeric(file_name[[1]][1])
    k_run <- as.numeric(file_name[[1]][2])
    likelihood <- as.character(unname(unlist(
      gsubfn::read.pattern(file = files_structure_2[i],
                           pattern = "^Marginal Likelihood = .*")
    )))
    likelihood_2 <-
      as.numeric(substr(likelihood, n_first_line , nchar(likelihood)))
    df_likelihood[k_run, k_replicate] <- likelihood_2
  }
  df_likelihood <-
    df_likelihood[stats::complete.cases(df_likelihood), ]
  df_likelihood_res <- rowMeans(df_likelihood)
  
  p3 <- ggplot() +
    geom_line(aes(x = k.range, y = df_likelihood_res), size = 1) +
    geom_point(aes(x = k.range, y = df_likelihood_res),
               size = 2,
               color = "blue") +
    theme_bw(base_size = 14) +
    theme(legend.title = element_blank()) +
    theme(legend.position = "bottom") +
    xlab("K") +
    ylab("Marginal Likelihood") +
    scale_x_continuous(breaks = round(seq(1, 10, by = 1), 1))
  
  print(p3)
  
  files_structure <- list.files(path = output,
                                pattern =  "^genotypes_output.+log")
  files_structure_2 <- paste0(output, "/", files_structure)
  
  files_q <- list.files(path = output, pattern = "*meanQ")
  files_q_2 <- paste0(output, "/", files_q)
  q_list <-
    rep(list(as.list(rep(NA, num.k.rep))), length(k.range) + 1)
  for (i in 1:length(files_q)) {
    file_name <- files_q[i]
    file_name <-
      strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\2', file_name), ' ')
    file_name <- strsplit(file_name[[1]], "\\.")
    k_replicate <- as.numeric(file_name[[1]][1])
    k_run <- as.numeric(file_name[[1]][2])
    q_df <- read.table(files_q_2[i])
    q_df <- cbind(id = x$ind.names, orig.pop = pop(x), q_df)
    q_list[[k_run]][[k_replicate]] <- q_df
    
  }
  
  names(q_list) <- 1:(length(k.range) + 1)
  
  q_list <- lapply(q_list, function(y) {
    names(y) <- 1:num.k.rep
    return(y)
  })
  
  q_list <- q_list[-(1)]
  
  return(q_list)
  
}
