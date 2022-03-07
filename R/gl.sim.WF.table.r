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
#' \item selection - SNP type.
#' } 
#' 
#' The reference table can be further modified as required. 
#' 
#' See documentation and tutorial for a complete description of the simulations.
#' These documents can be accessed by typing in the R console:
#' browseVignettes(package="dartR”)
#' 
#' @param file_var Path of the variables file 'ref_variables.csv' (see details) 
#' [required if interactive_vars = FALSE].
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
#' @param ... Here you can add separately any variable and its value which will 
#' be changed over the input value supplied by the csv file. See tutorial. 
#' @details
#' Values for the variables to create the reference table can be submitted into the function 
#' interactively through a Shiny app if interactive_vars = TRUE. Optionally, if 
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
#' @rawNamespace import(fields, except = flame)
#' @export

gl.sim.WF.table <- function(file_var, 
                            x = NULL, 
                            file_targets_sel = NULL, 
                            file_r_map = NULL,
                            interactive_vars = TRUE, 
                            verbose = NULL,
                            ...) {
  
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbosity = verbose)
    
    # DO THE JOB
    ##### LOADING VARIABLES ######
    
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
    
    input_list <- list(...)
    
    if(length(input_list>0)){
      
      val_change <- which(ref_vars$variable %in% names(input_list))
      
      ref_vars[val_change,"value"] <- list(input_list)
      vars_assign <-
        unlist(unname(
          mapply(paste, ref_vars$variable, "<-",
                 ref_vars$value, SIMPLIFY = F)
        ))
      eval(parse(text = vars_assign))
    }
    
    ##### LOADING INFORMATION #######
    # recombination map
    if (!is.null(file_r_map)) {
      map <- read.csv(file_r_map)
      map$Chr <- as.character(map$Chr)
      
      if(!chromosome_name %in% map$Chr){
        cat(error("  Chromosome name is not in the recombination map file\n"))
        stop()
      }
      map <- map[which(map$Chr == chromosome_name),]
      #### double check this line #####
      map$cM <- map$cM / 100
      # if there are NA's values converted them to 0
      map[is.na(map$cM),] <- 0
    }else{
      map <- as.data.frame(matrix(nrow = chunk_number))
      map[, 1] <- chunk_cM / 100
      colnames(map) <- "cM"
    }
    # targets of selection 
    targets_temp <- NULL
    if (!is.null(file_targets_sel)) {
      targets_temp <- read.csv(file_targets_sel)
      targets_temp$chr_name <- as.character(targets_temp$chr_name)
      
      if(!chromosome_name %in% targets_temp$chr_name ){
        cat(error("  Chromosome name is not in the targets of selection file\n"))
        stop()
      }
      
      targets_temp <- targets_temp[which(targets_temp$chr_name == chromosome_name),]
    }
    # real dataset 
    location_real_temp <- NULL
    if (real_loc == TRUE & !is.null(x)) {
      location_real_temp <-
        as.data.frame(cbind(as.character(x$chromosome), x$position))
      colnames(location_real_temp) <- c("chr", "pos")
      
      if(!chromosome_name %in% location_real_temp$chr ){
        cat(error("  Chromosome name is not in the genlight object\n"))
        stop()
      }
      
      location_real_temp <-
        location_real_temp[location_real_temp$chr == chromosome_name, ]
      location_real_temp <- 
        as.numeric(location_real_temp[, "pos"])
      
      location_real_temp <- location_real_temp[order(location_real_temp)]
    }
    
    ##### CHROMOSOME LENGTH ####
    if (!is.null(file_r_map)) {
      chr_length <- tail(map$to,1)
    }else{
      chr_length <- chunk_number  * chunk_bp
    }
    
    # if real locations are provided, it is necessary to change the bp per 
    # chromosome chunk. 
    if (real_loc == TRUE & !is.null(x) & is.null(file_r_map)) {
       chr_length <- tail(location_real_temp,1)
       chunk_bp <- chr_length / chunk_number
    }
    
    # if real locations are provided, it is necessary to change the bp per 
    # chromosome chunk. 
    if (!is.null(file_targets_sel) & is.null(file_r_map)) {
      chr_length <- tail(targets_temp$end,1)
      chunk_bp <- chr_length / chunk_number
    }
      
    ##### LOCATIONS ##########
    # real dataset
    location_real_bp <- NULL
    if (real_loc == TRUE & !is.null(x)) {
      location_real_bp <- location_real_temp
      location_real_bp <- location_real_bp[order(location_real_bp)]
      # real locations are rounded to the ten digits and 1 is added to each location
      location_real_bp <- round(location_real_bp,-1)
      location_real_bp <- location_real_bp + 1
    } 
    
    if (real_loc == FALSE & real_freq == TRUE & !is.null(x)) {
      location_real_temp <-
        round(seq(
        chunk_bp / (nLoc(x) + 1),
        (chunk_number * chunk_bp),
        chunk_bp / nLoc(x)
      ))
      
      location_real_bp <- sample(location_real_temp,size =nLoc(x))
      location_real_bp <- location_real_bp[order(location_real_bp)]
      # real locations are rounded to the ten digits and 1 is added to each location
      location_real_bp <- round(location_real_bp,-1)
      location_real_bp <- location_real_bp + 1
    }
    
    # neutral loci simulations
    location_neutral_bp <- NULL
    if(neutral_loci_chunk>0){
      location_neutral_bp <-
        round(seq(
          chunk_bp / (neutral_loci_chunk + 1),
          (chunk_number * chunk_bp),
          chunk_bp / neutral_loci_chunk
        ))
      
      # neutral loci locations are rounded to the ten digits and 2 is added to each location
      location_neutral_bp <- round(location_neutral_bp,-1)
      location_neutral_bp <- location_neutral_bp + 2
    }
    
    # loci under selection
    location_targets_bp <- NULL
    if (!is.null(file_targets_sel) | loci_under_selection>0) {
    if (!is.null(file_targets_sel)) {
      targets <- targets_temp
      targets$targets <- ceiling(targets$targets * (targets_factor/100))
      targets$distance <-  targets$end - targets$start
    } 
    # if the targets of selection file is not provided
    if (is.null(file_targets_sel) ) {
      targets <- as.data.frame(matrix(nrow = chunk_number, ncol = 3))
      colnames(targets) <- c("start", "end", "targets")
      targets$start <- seq(1, chr_length, (chr_length/chunk_number))
      targets$end <- seq((chr_length/chunk_number), chr_length, 
                         (chr_length/chunk_number))
      targets$targets <- round(loci_under_selection / chunk_number)
      targets$distance <-  targets$end - targets$start
    }
    
      sample_resolution <- round(mean(targets$distance)/max(targets$targets)/10)
      
    for (i in 1:nrow(targets)) {
      location_targets_temp <- mapply(
        FUN = function(a, b) {
          seq(from = a,
              to = b,
              by = sample_resolution)
        },
        a = unname(unlist(targets[i, "start"])),
        b = unname(unlist(targets[i, "end"]))
      )
      location_targets_temp <- as.vector(round(location_targets_temp))
      location_targets_temp <- sample(location_targets_temp, 
                                      size = targets[i, "targets"])
      location_targets_bp <- c(location_targets_bp, location_targets_temp)
    }
      location_targets_bp <- location_targets_bp[order(location_targets_bp)]
    # loci under selection locations are rounded to the ten digits and 3 is added to each location
      location_targets_bp <- round(location_targets_bp,-1)
    location_targets_bp <- location_targets_bp + 3
    }
    
    # mutations 
    location_mutations_bp <- NULL
    if (!is.null(file_targets_sel) | (loci_mutation > 0 & mutation == TRUE)) {
    # if the targets of selection file is provided
    if (!is.null(file_targets_sel)) {
      mutations <- targets_temp
      mutations$targets <- ceiling(mutations$targets * (mutations_factor/100))
      mutations$distance <-  mutations$end - mutations$start
      # if the targets of selection file is not provided
    }else{
      mutations <- as.data.frame(matrix(nrow = chunk_number, ncol = 3))
      colnames(mutations) <- c("start", "end", "targets")
      mutations$start <- seq(1, chr_length, (chr_length/chunk_number))
      mutations$end <- seq((chr_length/chunk_number), chr_length,
                           (chr_length/chunk_number))
      mutations$targets <- round(loci_mutation / chunk_number)
      mutations$distance <-  mutations$end - mutations$start
    }
    # the resolution to sample mutations is different to the resolution to sample
    # targets for slection to avoid averlapping locations
    sample_resolution <- round(mean(mutations$distance)/max(mutations$targets)/15)
    
    for (i in 1:nrow(mutations)) {
      location_mutations_temp <- mapply(
        FUN = function(a, b) {
          seq(from = a,
              to = b,
              by = sample_resolution)
        },
        a = unname(unlist(mutations[i, "start"])),
        b = unname(unlist(mutations[i, "end"]))
      )
      location_mutations_temp <- as.vector(round(location_mutations_temp))
      
      location_mutations_temp <- sample(location_mutations_temp, size = mutations[i, "targets"])
      location_mutations_bp <- c(location_mutations_bp, location_mutations_temp)
    }
    location_mutations_bp <- location_mutations_bp[order(location_mutations_bp)]
    # mutation locations are rounded to the ten digits and 4 is added to each location
    location_mutations_bp <- round(location_mutations_bp,-1)
    location_mutations_bp <- location_mutations_bp + 4
    }
    
    location_loci_bp <- c(location_real_bp, location_neutral_bp,
                       location_targets_bp, location_mutations_bp)
    
    location_loci_bp <- location_loci_bp[order(location_loci_bp)]
    
    if(chunk_number > length(location_loci_bp)){
      cat(error("  Number of loci should be more than the number of genome chunks\n"))
      stop()
    }
  
    total_loci <- length(location_loci_bp)
    
    ##### RECOMBINATION MAP #####
    # the recombination map is produced by cross multiplication. the following
    # lines are the input for doing the cross multiplication.
    recombination_map_temp <- map
    recombination_map_temp$midpoint <- seq(chunk_bp / 2, 
                                           chr_length, 
                                           chunk_bp)[1:nrow(recombination_map_temp)]
    recombination_temp <-
      unlist(lapply(location_loci_bp, findInterval, vec = as.numeric(paste(
        unlist(recombination_map_temp$midpoint)
      ))))
    
    # loci located below the location in the first row of the
    # recombination map are assigned to row 0, to correct this, they
    # are reassigned to row number 1
    recombination_temp[recombination_temp == 0] <- 1
    recombination_2 <- recombination_map_temp[recombination_temp, "cM"]
    recombination_map <- as.data.frame(cbind(location_loci_bp, recombination_2))
    recombination_map$c <- NA
    # the recombination map is produced by cross multiplication
    #not taking in account the last row for the loop to work
    for (target_row in 1:(nrow(recombination_map) - 1)) {
      recombination_map[target_row, "c"] <-
        ((recombination_map[target_row + 1, "location_loci_bp"] - recombination_map[target_row, "location_loci_bp"]) * recombination_map[target_row, "recombination_2"]) / chunk_bp
    }
    # The last element of the recombination rate column must be zero,
    # otherwise the recombination function crashes.
    recombination_map[nrow(recombination_map), "c"] <- 0
    recombination_map$accum <- cumsum(recombination_map[, "c"])
    
    # getting the row of the recombination map for neutral loci
    location_neutral_row <- NULL
    if(neutral_loci_chunk>0){
    location_neutral_row <- lapply(location_neutral_bp, function(x) {
      which(recombination_map$location_loci == x)
      })
    location_neutral_row <- unname(unlist(location_neutral_row))
    }
    
    # getting the row of the recombination map for real loci
    location_real_row <- NULL
    if(real_loc==TRUE | real_freq == TRUE){
      location_real_row <- lapply(location_real_bp, function(x) {
        which(recombination_map$location_loci == x)
        })
      location_real_row <- unname(unlist(location_real_row))
    }
    
    # getting the row of the recombination map for loci under selection
    location_targets_row <- NULL
    if(neutral_loci_chunk>0){
      location_targets_row <- lapply(location_targets_bp, function(x) {
        which(recombination_map$location_loci == x)
      })
      location_targets_row <- unname(unlist(location_targets_row))
    }
    
    # getting the row of the recombination map for mutations
    location_mutations_row <- NULL
    if(mutation==T){
      location_mutations_row  <- lapply(location_mutations_bp, function(x) {
        which(recombination_map$location_loci == x)
      })
      location_mutations_row <- unname(unlist(location_mutations_row))
    }
    
    ##### REFERENCE TABLE ########
    # selection coefficient
    if(s_distribution== "exponential_gamma" ){
      s_advantageous <- rexp(total_loci,rate = exp_rate) * -1
      s_deleterious <- rgamma(total_loci, shape = gamma_shape, scale = gamma_scale)
      s_temp <- c(sample(s_advantageous,
                         size=round(total_loci*(percent_adv/100))),
             sample(s_deleterious
                    ,size=round(total_loci*(1-percent_adv/100))))
      s <- sample(s_temp,size =total_loci,replace = T )
    }
    
    if (s_distribution == "equal") {
      s <- s_gral
    }
    
    if (s_distribution == "gamma") {
      s <- rgamma(total_loci, shape = gamma_shape, scale = gamma_scale)
    }
    
    if (s_distribution == "log_normal") {
      s <- rlnorm(total_loci,
                  meanlog = log(log_mean),
                  sdlog = log(log_sd))
    }
    
    # dominance
    if (h_distribution == "equal") {
      h <- h_gral
    }
    
    if (h_distribution == "normal") {
      h <- rnorm(total_loci, mean = dominance_mean, sd = dominance_sd)
    }
    
    # the equation for dominance (h) was taken from Huber 2018 Nature
    if (h_distribution == "equation") {
      h <- 1 / ((1 / intercept) - (-1 * rate * abs(s)))
    }
    
    # initial frequency
    if (q_distribution == "equal") {
      q <- q_gral
    }
    
    if (q_distribution == "equation") {
      a <- abs(s) * (1 - (2 * h))
      b <- (h * abs(s)) * (1 + mutation_rate)
      c <- rep.int(-(mutation_rate), times = total_loci)
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
    
    reference <- as.data.frame(matrix(nrow = total_loci))
    reference$q <- q
    reference$h <- h
    reference$s <- s
    reference$c <- recombination_map[1:total_loci, "c"]
    reference$loc_bp <- recombination_map[1:total_loci, "location_loci_bp"]
    reference$loc_cM <- recombination_map[1:total_loci, "accum"]
    reference$chr_name <- chromosome_name
    
    # setting h and s to 0 in neutral loci and loci from real dataset
    reference[as.numeric(location_real_row), "s"] <- 0
    reference[as.numeric(location_real_row), "h"] <- 0
    reference[as.numeric(location_real_row), "q"] <- q_neutral
    reference[as.numeric(location_neutral_row), "s"] <- 0
    reference[as.numeric(location_neutral_row), "h"] <- 0
    reference[as.numeric(location_neutral_row), "q"] <- q_neutral
    # NS with very small s have a q > 1. Therefore, we set a maximum q value of
    # 0.5.
    q_more_than_point5 <- as.numeric(row.names(reference[reference$q > 0.5, ]))
    reference[q_more_than_point5, "q"] <- 0.5
    # the log normal distribution, with the parameters used in the simulations,
    # generates a few selection coefficients that are > 1. The maximum value of s
    # is set to 0.99
    s_more_than_one <- as.numeric(row.names(reference[reference$s > 1, ]))
    reference[s_more_than_one, "s"] <- 0.99
    # the exponential distribution, with the parameters used in the simulations,
    # generates a few selection coefficients that are < -1. The minimum value of s
    # is set to -0.5
    s_less_minus_one <- as.numeric(row.names(reference[reference$s < - 0.5, ]))
    reference[s_less_minus_one, "s"] <- -0.5

    # setting q to 0 in loci available to mutation
    reference[location_mutations_row, "q"] <- 0
  
    reference <- reference[, -1]
    
    mutation_loci_adv <- which(reference$q == 0 & reference$s <0)
    mutation_loci_del <- which(reference$q == 0 & reference$s >0)
    deleterious <- which(reference$q > 0 & reference$s > 0 )
    advantageous <- which(reference$q > 0 & reference$s < 0 )

    reference[location_neutral_row, "selection"] <- "neutral"
    reference[location_real_row, "selection"] <- "real"
    reference[mutation_loci_adv, "selection"] <- "mutation_adv"
    reference[mutation_loci_del, "selection"] <- "mutation_del"
    reference[deleterious, "selection"] <- "deleterious"
    reference[advantageous, "selection"] <- "advantageous"
  
    ref_res <- list(reference, ref_vars)
    names(ref_res) <- c("reference", "ref_vars")
    ##### END ######
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
      cat(report("Completed:", funname, "\n"))
    }
    
    # RETURN
    
    return(invisible(ref_res))
    
  }
