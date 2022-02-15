#' @name gl.sim.WF.table
#' @title Create reference table for running gl.sim.WF.run
#' @description
#' This function creates a reference table to be used as input for the function
#'  \code{\link{gl.sim.WF.run}}. The created table has seven columns with the 
#'  following information for each locus to be simulated:
#' \itemize{ 
#' \item q - initial frequency.
#' \item h - dominance coefficient.
#' \item s - selection coefficient.
#' \item c - recombination rate.
#' \item loc_bp - chromosome location in base pairs.
#' \item loc_cM - chromosome location in centiMorgans.
#' \item chr_name - chromosome name.
#' } 
#' 
#' The reference table can be further modified as required. 
#' 
#' @param file_var Path of the variables file 'ref_variables.csv' (see details) 
#' [required if interactive_vars = TRUE].
#' @param x Name of the genlight object containing the SNP data to extract
#' values for some simulation variables (see details) [default NULL].
#' @param file_targets_sel Path of the file with the targets for selection (see 
#' details)  [default NULL].
#' @param file_r_map Path of the file with the recombination map (see details)
#' [default NULL].
#' @param interactive_vars Run a shiny app to input interactively the values of
#'  simulations variables [default TRUE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @details
#' Values for the creation of the reference table can be submitted into the function 
#' interactively through a shiny app if interactive_vars = TRUE. Optionally, if 
#' interactive_vars = FALSE, values for variables can be submitted by using the
#' csv file 'ref_variables.csv' which can be found by typing in the R console:
#'  system.file('extdata', 'ref_variables.csv', package ='dartR').
#'  
#' The values of the variables can be modified using the third column (“value”) 
#' of this file. 
#' 
#' If a genlight object is used as input for some of the simulation variables, 
#' this function access the information stored in the slots x$position and 
#' x$chromosome.
#' 
#' Examples of the format required for the recombination map file and the 
#' targets for selection file can be found by typing in the R console:
#' \itemize{ 
#' \item system.file('extdata', 'fly_recom_map.csv', package ='dartR')
#' \item system.file('extdata', 'fly_targets_of_selection.csv', package ='dartR')
#' }
#' 
#' Functions to produce these files are gl.sim.create_rmap() and 
#' gl.sim.create_targets().
#' 
#' To show further information of the variables in interactive mode, it might be
#'  necessary to call first: 'library(shinyBS)' for the information to be 
#'  displayed.
#' @return Returns a list with the reference table used as input for the function
#'  \code{\link{gl.sim.WF.run}} and a table with the values variables used to 
#'  create the reference table.
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' ref_table <- gl.sim.WF.table(file_var=system.file('extdata', 
#' 'ref_variables.csv', package = 'dartR'),interactive_vars = FALSE)
#' res_sim <- gl.sim.WF.run(file_var = system.file('extdata', 
#' 'sim_variables.csv', package ='dartR'),ref_table=ref_table,
#' interactive_vars = FALSE)
#' @seealso \code{\link{gl.sim.WF.run}}
#' @family simulation functions
#' @import shiny
#' @import shinyBS
#' @import shinythemes
#' @import shinyjs
#' @export

gl.sim.WF.table <-
  function(file_var,
           x = NULL,
           file_targets_sel = NULL,
           file_r_map = NULL,
           interactive_vars = TRUE,
           verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbosity = verbose)
    
    if (interactive_vars) {
      pkg <- "shiny"
      if (!(requireNamespace(pkg, quietly = TRUE))) {
        stop(error(
          "Package",
          pkg,
          " needed for this function to work. Please install it."
        ))
      }
      pkg <- "shinyBS"
      if (!(requireNamespace(pkg, quietly = TRUE))) {
        stop(error(
          "Package",
          pkg,
          "needed for this function to work. Please install it."
        ))
      }
      pkg <- "shinythemes"
      if (!(requireNamespace(pkg, quietly = TRUE))) {
        stop(error(
          "Package",
          pkg,
          "needed for this function to work. Please install it."
        ))
      }
      pkg <- "shinyjs"
      if (!(requireNamespace(pkg, quietly = TRUE))) {
        stop(error(
          "Package",
          pkg,
          "needed for this function to work. Please install it."
        ))
      }
    }

    # DO THE JOB
    ##### SIMULATIONS VARIABLES ######
    
    if (interactive_vars) {
      
      ref_vars <- interactive_reference()
      
      ref_vars[ref_vars$variable=="h_distribution" ,"value"] <- 
        paste0("'",ref_vars[ref_vars$variable=="h_distribution" ,"value"],"'")
      
      ref_vars[ref_vars$variable=="chromosome_name" ,"value"] <- 
        paste0("'",ref_vars[ref_vars$variable=="chromosome_name" ,"value"],"'")
      
      ref_vars[ref_vars$variable=="s_distribution" ,"value"] <- 
        paste0("'",ref_vars[ref_vars$variable=="s_distribution" ,"value"],"'")
      
      ref_vars[ref_vars$variable=="q_distribution" ,"value"] <- 
        paste0("'",ref_vars[ref_vars$variable=="q_distribution" ,"value"],"'")
      
      vars_assign <-
        unlist(unname(
          mapply(paste, ref_vars$variable, "<-",
                 ref_vars$value, SIMPLIFY = F)
        ))
      
      eval(parse(text = vars_assign))
      
    } else{
      ref_vars <- suppressWarnings(read.csv(file_var))
      ref_vars <- ref_vars[, 2:3]
      vars_assign <-
        unlist(unname(
          mapply(paste, ref_vars$variable, "<-",
                 ref_vars$value, SIMPLIFY = F)
        ))
      eval(parse(text = vars_assign))
    }
    
    # if the recombination map is provided
    if (!is.null(file_r_map)) {
      map <- read.csv(file_r_map)
      map$Chr <- as.character(map$Chr)
      
      if(!chromosome_name %in% map$Chr){
        cat(error("Chromosome name is not in the recombination map file\n"))
        stop()
      }
      
      map <- map[which(map$Chr == chromosome_name),]
      map <- as.data.frame(map$cM / 1000)
      colnames(map) <- "cM"
      map[is.na(map$cM),] <- 0
      # if the recombination map is not provided
    } else{
      map <- as.data.frame(matrix(nrow = chunk_number))
      map[, 1] <- chunk_recombination / 100
      colnames(map) <- "cM"
    }
    
    chr_length <- (nrow(map) + 1) * map_resolution
    
    # if the targets of selection file is provided
    if (!is.null(file_targets_sel)) {
      targets <- read.csv(file_targets_sel)
      targets$chr_name <- as.character(targets$chr_name)
      
      if(!chromosome_name %in%   targets$chr_name ){
        cat(error("Chromosome name is not in the targets of selection file\n"))
        stop()
      }
      
      targets <-
        targets[which(targets$chr_name == chromosome_name),]
      targets <- targets[!duplicated(targets$start),]
      targets <- targets[!duplicated(targets$end),]
      # if the targets of selection file is not provided
    } else{
      targets <- as.data.frame(matrix(nrow = chunk_number, ncol = 3))
      colnames(targets) <- c("start", "end", "targets")
      targets$start <- head(seq(1, chr_length, map_resolution), -1)
      targets$end <-
        head(seq(map_resolution, chr_length, map_resolution), -1)
      targets$targets <- loci_number_to_simulate / chunk_number
    }
    
    targets$targets <- ceiling(targets$targets * (targets_factor/100))
    
    location_targets <- NULL
    
    for (i in 1:nrow(targets)) {
      location_targets_temp <-
        unlist(mapply(
          FUN = function(a, b) {
            seq(from = a,
                to = b,
                by = 4)
          },
          a = unname(unlist(targets[i, "start"])),
          b = unname(unlist(targets[i, "end"]))
        ))
      location_targets_temp <-
        sample(location_targets_temp, size = targets[i, "targets"])
      location_targets <-
        c(location_targets, location_targets_temp)
    }
    
    # these are the location of the neutral loci in the simulations
    if (real_loc == TRUE & !is.null(x)) {
      location_neutral_loci_analysis_temp <-
        as.data.frame(cbind(as.character(x$chromosome), x$position))
      colnames(location_neutral_loci_analysis_temp) <-
        c("chr", "pos")
      location_neutral_loci_analysis_temp <-
        location_neutral_loci_analysis_temp[location_neutral_loci_analysis_temp$chr ==
                                              chromosome_name, ]
      location_neutral_loci_analysis <-
        as.numeric(location_neutral_loci_analysis_temp[, "pos"])
      
    } else {
      location_neutral_loci_analysis <-
        round(seq(
          map_resolution / (neutral_loci_chunk + 1),
          (nrow(map) * map_resolution),
          map_resolution / neutral_loci_chunk
        ))
      
    }
    
    location_targets <-
      c(location_targets, location_neutral_loci_analysis)
    location_targets <- location_targets[order(location_targets)]
    
    # different transcripts can be located in the same genome location. So,
    # repeated deleterious mutations are deleted
    location_targets <- unique(location_targets)
    
    #this is to fix a bug that crashes the program because the last neutral
    # locus sometimes could be located farther than the last deleterious mutation
    location_targets <- c(location_targets, chr_length)
    loci_number_to_simulate <- length(location_targets)
    
    # the recombination map is produced by cross multiplication. the following
    # lines are the input for doing the cross multiplication.
    recombination_map_temp <- map
    
    recombination_map_temp$midpoint <-
      seq(map_resolution / 2,
          nrow(recombination_map_temp) * map_resolution,
          map_resolution)
    
    recombination_temp <-
      unlist(lapply(location_targets, findInterval, vec = as.numeric(paste(
        unlist(recombination_map_temp$midpoint)
      ))))
    
    # deleterious mutations located below the location in the first row of the
    # recombination map (i.e. 50000) are assigned to row 0, to correct this they
    # are reassigned to row number 1
    recombination_temp[recombination_temp == 0] <- 1
    recombination_2 <-
      recombination_map_temp[recombination_temp, "cM"]
    recombination_map <-
      as.data.frame(cbind(location_targets, recombination_2))
    recombination_map$c <- NA
    # the recombination map is produced by cross multiplication
    #not taking in account the last row for the loop to work
    for (deleterious_row in 1:(nrow(recombination_map) - 1)) {
      recombination_map[deleterious_row, "c"] <-
        ((recombination_map[deleterious_row + 1, "location_targets"] - recombination_map[deleterious_row, "location_targets"]) * recombination_map[deleterious_row, "recombination_2"]) / map_resolution
    }
    
    loci_number <- loci_number_to_simulate
    # In order for the recombination rate to be accurate, we must account for
    # the case when the probability of the total recombination rate is less than 1
    # (i.e. < 100 cM). For this end, the program subtracts from 1 the sum of all
    # the recombination rates and this value inserted in the last row of the
    # recombination_map table. If this row is chosen as the recombination point,
    # recombination does not occur. For example, if a chromosome of 20 cM’s is
    # simulated, the last row of the recombination_map will have a value of 0.8
    # and therefore 80% of the times recombination will not occur.
    # number of recombination events per meiosis
    recom_event <- ceiling(sum(recombination_map[, "c"]))
    recombination_map[loci_number + 1, "c"] <-
      recom_event - sum(recombination_map[, "c"])
    recombination_map[loci_number + 1, "location_targets"] <-
      recombination_map[loci_number, "location_targets"]
    recombination_map[loci_number + 1, "recombination_2"] <-
      recombination_map[loci_number, "recombination_2"]
    recombination_map$accum <- cumsum(recombination_map[, "c"])
    
    location_temp <- location_neutral_loci_analysis
    location_temp <- location_temp[order(location_temp)]
    location <-
      lapply(location_temp, function(x) {
        which(recombination_map$location_targets == x)
      })
    neutral_loci_location <-  as.character(location)
    loci_location <- location
    
    if (s_distribution == "equal") {
      s <- s_gral
    }
    
    if (s_distribution == "gamma") {
      s <- rgamma(loci_number, shape = gamma_shape, scale = gamma_scale)
    }
    
    if (s_distribution == "log_normal") {
      s <- rlnorm(loci_number,
                  meanlog = log(log_mean),
                  sdlog = log(log_sd))
    }
    
    if (h_distribution == "equal") {
      h <- h_gral
    }
    
    if (h_distribution == "normal") {
      h <-
        rnorm(loci_number, mean = dominance_mean, sd = dominance_sd)
    }
    
    # the equation for dominance (h) was taken from Huber 2018 Nature
    if (h_distribution == "equation") {
      h <- 1 / ((1 / intercept) - (-1 * rate * s))
    }
    
    if (q_distribution == "equal") {
      q <- q_gral
    }
    
    if (q_distribution == "equation") {
      a <- s * (1 - (2 * h))
      b <- (h * s) * (1 + mutation_rate)
      c <- rep.int(-(mutation_rate), times = loci_number)
      df_q <- as.data.frame(cbind(a, b, c))
      # q is based on the following equation: (s(1-2h)q^2) + (hs(1+u)q) - u = 0,
      # where u is the mutation rate per generation per site. Taken from Crow &
      # Kimura page 260
      q <-
        mapply(
          q_equilibrium,
          a = df_q$a,
          b = df_q$b,
          c = df_q$c,
          USE.NAMES = F
        )
    }
    
    reference <- as.data.frame(matrix(nrow = loci_number))
    reference$q <- q
    reference$h <- h
    reference$s <- s
    reference$c <- recombination_map[1:loci_number, "c"]
    reference$loc_bp <-
      recombination_map[1:loci_number, "location_targets"]
    reference$loc_cM <- recombination_map[1:loci_number, "accum"]
    reference$chr_name <- chromosome_name
    # setting h and s to 0 in neutral loci
    reference[as.numeric(neutral_loci_location), "s"] <- 0
    reference[as.numeric(neutral_loci_location), "h"] <- 0
    reference[as.numeric(neutral_loci_location), "q"] <- q_neutral
    # NS with very small s have a q > 1. Therefore, we set a maximum q value of
    # 0.5.
    q_more_than_point5 <-
      as.numeric(row.names(reference[reference$q > 0.5, ]))
    reference[q_more_than_point5, "q"] <- 0.5
    # the log normal distribution, with the parameters used in the simulations,
    # generates a few selection coefficients that are > 1. The maximum value of s
    # is set to 0.99
    s_more_than_one <-
      as.numeric(row.names(reference[reference$s > 1, ]))
    reference[s_more_than_one, "s"] <- 0.99
    loci_number <- nrow(reference)
    # setting h and s to 0 in neutral loci
    reference[as.numeric(neutral_loci_location), "s"] <- 0
    reference[as.numeric(neutral_loci_location), "h"] <- 0
    
    ref_res <- list(reference[, -1], ref_vars)
    names(ref_res) <- c("reference", "ref_vars")
    
    return(ref_res)
    
  }
