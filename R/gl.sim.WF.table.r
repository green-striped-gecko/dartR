#' @name gl.sim.WF.table
#' @title Creates the reference table for running gl.sim.WF.run
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
#' \item type - SNP type.
#' } 
#' 
#' The reference table can be further modified as required. 
#' 
#' See documentation and tutorial for a complete description of the simulations.
#' These documents can be accessed at http://georges.biomatix.org/dartR 
#' 
#' @param file_var Path of the variables file 'ref_variables.csv' (see details) 
#' [required if interactive_vars = FALSE].
#' @param x Name of the genlight object containing the SNP data to extract
#' values for some simulation variables (see details) [default NULL].
#' @param file_targets_sel Path of the file with the targets for selection (see 
#' details) [default NULL].
#' @param file_r_map Path of the file with the recombination map (see details)
#' [default NULL].
#' @param interactive_vars Run a shiny app to input interactively the values of
#'  simulation variables [default TRUE].
#' @param seed Set the seed for the simulations [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @param ... Any variable and its value can be added separately within the 
#' function, will be changed over the input value supplied by the csv file. See 
#' tutorial. 
#' @details
#' Values for the variables to create the reference table can be submitted into 
#' the function interactively through a Shiny app if interactive_vars = TRUE. 
#' Optionally, if interactive_vars = FALSE, values for variables can be 
#' submitted by using the csv file 'ref_variables.csv' which can be found by 
#' typing in the R console:
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
#' \dontrun{
#' #uncomment to run 
#' res_sim <- gl.sim.WF.run(file_var = system.file('extdata', 
#' 'sim_variables.csv', package ='dartR'),ref_table=ref_table,
#' interactive_vars = FALSE)
#' }
#' @seealso \code{\link{gl.sim.WF.run}}
#' @family simulation functions
#' @rawNamespace import(fields, except = flame)
#' @export

gl.sim.WF.table <- function(file_var, 
                            x = NULL, 
                            file_targets_sel = NULL, 
                            file_r_map = NULL,
                            interactive_vars = TRUE, 
                            seed = NULL,
                            verbose = NULL,
                            ...) {
  
  # setting the seed
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbosity = verbose)
    
    # DO THE JOB
    ##### LOADING VARIABLES ######
    
    if (interactive_vars==TRUE) {
      
      ref_vars <- interactive_reference()
      
      ref_vars[ref_vars$variable=="chromosome_name" ,"value"] <-
        paste0("'",ref_vars[ref_vars$variable=="chromosome_name" ,"value"],"'")
      
      ref_vars[ref_vars$variable=="h_distribution_del" ,"value"] <-
        paste0("'",ref_vars[ref_vars$variable=="h_distribution_del" ,"value"],"'")

      ref_vars[ref_vars$variable=="s_distribution_del" ,"value"] <-
        paste0("'",ref_vars[ref_vars$variable=="s_distribution_del" ,"value"],"'")

      ref_vars[ref_vars$variable=="q_distribution_del" ,"value"] <-
        paste0("'",ref_vars[ref_vars$variable=="q_distribution_del" ,"value"],"'")
      
      ref_vars[ref_vars$variable=="h_distribution_adv" ,"value"] <-
        paste0("'",ref_vars[ref_vars$variable=="h_distribution_adv" ,"value"],"'")
      
      ref_vars[ref_vars$variable=="s_distribution_adv" ,"value"] <-
        paste0("'",ref_vars[ref_vars$variable=="s_distribution_adv" ,"value"],"'")
      
      ref_vars[ref_vars$variable=="q_distribution_adv" ,"value"] <-
        paste0("'",ref_vars[ref_vars$variable=="q_distribution_adv" ,"value"],"'")

      ref_vars <- ref_vars[order(ref_vars$variable),]
      
      vars_assign <-
        unlist(unname(
          mapply(paste, ref_vars$variable, "<-",
                 ref_vars$value, SIMPLIFY = F)
        ))
      
      eval(parse(text = vars_assign))
      
    } else {
      ref_vars <- suppressWarnings(read.csv(file_var))
      ref_vars <- ref_vars[, 2:3]
      
      ref_vars <- ref_vars[order(ref_vars$variable),]
      
      vars_assign <-
        unlist(unname(
          mapply(paste, ref_vars$variable, "<-",
                 ref_vars$value, SIMPLIFY = F)
        ))
      eval(parse(text = vars_assign))
    }
    
    input_list <- list(...)
    
    if(length(input_list)>0){
      
      ref_vars <- ref_vars[order(ref_vars$variable),]
      
      input_list <- input_list[order(names(input_list))]
      
      val_change <- which(ref_vars$variable %in% names(input_list))
      
      ref_vars[val_change,"value"] <- unlist(input_list)
      
      ref_vars[ref_vars$variable=="chromosome_name" ,"value"] <-
        paste0("'",ref_vars[ref_vars$variable=="chromosome_name" ,"value"],"'")
      
      ref_vars[ref_vars$variable=="h_distribution_del" ,"value"] <-
        paste0("'",ref_vars[ref_vars$variable=="h_distribution_del" ,"value"],"'")
      
      ref_vars[ref_vars$variable=="s_distribution_del" ,"value"] <-
        paste0("'",ref_vars[ref_vars$variable=="s_distribution_del" ,"value"],"'")
      
      ref_vars[ref_vars$variable=="q_distribution_del" ,"value"] <-
        paste0("'",ref_vars[ref_vars$variable=="q_distribution_del" ,"value"],"'")
      
      ref_vars[ref_vars$variable=="h_distribution_adv" ,"value"] <-
        paste0("'",ref_vars[ref_vars$variable=="h_distribution_adv" ,"value"],"'")
      
      ref_vars[ref_vars$variable=="s_distribution_adv" ,"value"] <-
        paste0("'",ref_vars[ref_vars$variable=="s_distribution_adv" ,"value"],"'")
      
      ref_vars[ref_vars$variable=="q_distribution_adv" ,"value"] <-
        paste0("'",ref_vars[ref_vars$variable=="q_distribution_adv" ,"value"],"'")
      
      vars_assign <-
        unlist(unname(
          mapply(paste, ref_vars$variable, "<-",
                 ref_vars$value, SIMPLIFY = F)
        ))
      eval(parse(text = vars_assign))
    }
    
    s_distribution_del <- gsub('\"', "", s_distribution_del, fixed = TRUE)
    s_distribution_adv <- gsub('\"', "", s_distribution_adv, fixed = TRUE)
    h_distribution_del <- gsub('\"', "", h_distribution_del, fixed = TRUE)
    h_distribution_adv <- gsub('\"', "", h_distribution_adv, fixed = TRUE)
    q_distribution_del <- gsub('\"', "", q_distribution_del, fixed = TRUE)
    q_distribution_adv <- gsub('\"', "", q_distribution_adv, fixed = TRUE)
    
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
    if( (real_loc==TRUE | real_freq==TRUE) && is.null(x)){
      cat(error(" The real dataset to extract information is missing\n"))
      stop()
    }
    location_real_temp <- NULL
    if (real_loc == TRUE & !is.null(x)) {
      location_real_temp <-
        as.data.frame(cbind(as.character(x$chromosome), x$position))
      colnames(location_real_temp) <- c("chr", "pos")
      
      if(!chromosome_name %in% location_real_temp$chr ){
        cat(error("  Chromosome name is not in the genlight object\n"))
        stop()
      }
      
      location_real_temp <- location_real_temp[location_real_temp$chr == chromosome_name, ]
      location_real_temp <- as.numeric(location_real_temp[, "pos"])
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
    if(chunk_neutral_loci>0){
      location_neutral_bp <-
        round(seq(
          chunk_bp / (chunk_neutral_loci + 1),
          (chunk_number * chunk_bp),
          chunk_bp / chunk_neutral_loci
        ))
      
      # neutral loci locations are rounded to the ten digits and 2 is added to each location
      location_neutral_bp <- round(location_neutral_bp,-1)
      location_neutral_bp <- location_neutral_bp + 2
    }
    
    # deleterious loci
    location_deleterious_bp <- NULL
    if (!is.null(file_targets_sel) | loci_deleterious>0) {
    if (!is.null(file_targets_sel)) {
      del <- targets_temp
      del$targets <- ceiling(del$targets * (deleterious_factor/100))
      del$distance <-  del$end - del$start
    } 
    # if the targets of selection file is not provided
    if (is.null(file_targets_sel) ) {
      del <- as.data.frame(matrix(nrow = chunk_number, ncol = 3))
      colnames(del) <- c("start", "end", "targets")
      del$start <- seq(1, chr_length, (chr_length/chunk_number))
      del$end <- seq((chr_length/chunk_number), chr_length, 
                         (chr_length/chunk_number))
      if(loci_deleterious<chunk_number){
        row_targets <- sample(1:chunk_number,size = loci_deleterious)
        del[row_targets,"targets"] <- 1
        del[is.na(del$targets),"targets"] <- 0
      }else{
        del$targets <- round(loci_deleterious / chunk_number)
      }
      del$distance <-  del$end - del$start
    }
        sample_resolution <- round(mean(del$distance)/max(del$targets)/10)
    for (i in 1:nrow(del)) {
      location_deleterious_temp <- mapply(
        FUN = function(a, b) {
          seq(from = a,
              to = b,
              by = sample_resolution)
        },
        a = unname(unlist(del[i, "start"])),
        b = unname(unlist(del[i, "end"]))
      )
      location_deleterious_temp <- as.vector(round(location_deleterious_temp))
      location_deleterious_temp <- sample(location_deleterious_temp, 
                                      size = del[i, "targets"])
      location_deleterious_bp <- c(location_deleterious_bp, location_deleterious_temp)
    }
      location_deleterious_bp <- location_deleterious_bp[order(location_deleterious_bp)]
    # deleterious loci locations are rounded to the ten digits and 3 is added to each location
      location_deleterious_bp <- round(location_deleterious_bp,-1)
    location_deleterious_bp <- location_deleterious_bp + 3
    loci_deleterious <- length(location_deleterious_bp)
    }
      
    # advantageous loci
    location_advantageous_bp <- NULL
    if (loci_advantageous>0) {
      location_advantageous_bp <- sample(1:chr_length,loci_advantageous)
      location_advantageous_bp <- location_advantageous_bp[order(location_advantageous_bp)]
       # advantageous loci locations are rounded to the ten digits and 4 is added to each location
        location_advantageous_bp <- round(location_advantageous_bp,-1)
        location_advantageous_bp <- location_advantageous_bp + 4
    }
    
    # mutations 
    loci_mutation <- loci_mut_neu + loci_mut_del + loci_mut_adv
    
    location_mutations_bp <- NULL
    if (!is.null(file_targets_sel) | (loci_mutation > 0)) {
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
      mutations$targets <- ceiling(loci_mutation / chunk_number)
      mutations$distance <-  mutations$end - mutations$start
    }
    # the resolution to sample mutations is different to the resolution to sample
    # deleterious and advantageous loci to avoid overlapping locations
    sample_resolution <- round(mean(mutations$distance)/max(mutations$targets)/20)
    
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
    # mutation locations are rounded to the ten digits and 5 is added to each location
    location_mutations_bp <- round(location_mutations_bp,-1)
    location_mutations_bp <- location_mutations_bp + 5
    location_mutations_bp <- location_mutations_bp[sample(1:length(location_mutations_bp),size=loci_mutation)]
    }
    
    location_deleterious_bp <- unique(location_deleterious_bp)
    loci_deleterious <- length(location_deleterious_bp)
    
    location_advantageous_bp <- unique(location_advantageous_bp)
    loci_advantageous <- length(location_advantageous_bp)
    
    location_mutations_bp <- unique(location_mutations_bp)
    loci_mutation <- length(location_mutations_bp)

    location_loci_bp <- c(location_real_bp, location_neutral_bp,
                       location_deleterious_bp, location_advantageous_bp,
                       location_mutations_bp)
    
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
    if(chunk_neutral_loci>0){
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
    
    # getting the row of the recombination map for deleterious loci
    location_deleterious_row <- NULL
    if(loci_deleterious>0){
      location_deleterious_row <- lapply(location_deleterious_bp, function(x) {
        which(recombination_map$location_loci == x)
      })
      location_deleterious_row <- unname(unlist(location_deleterious_row))
    }
    
    # getting the row of the recombination map for deleterious loci
    location_advantageous_row <- NULL
    if(loci_advantageous>0){
      location_advantageous_row <- lapply(location_advantageous_bp, function(x) {
        which(recombination_map$location_loci == x)
      })
      location_advantageous_row <- unname(unlist(location_advantageous_row))
    }
    
    # getting the row of the recombination map for mutations
    location_mutations_row <- NULL
    if(loci_mutation>0){
      location_mutations_row  <- lapply(location_mutations_bp, function(x) {
        which(recombination_map$location_loci == x)
      })
      location_mutations_row <- unname(unlist(location_mutations_row))
    }
    
    ##### REFERENCE TABLE ########
    mut_tmp <- location_mutations_row
    if(loci_mut_neu>0){
    neu_tmp <- sample(1:length(mut_tmp),size=loci_mut_neu)
    location_mut_neu_row <- mut_tmp[neu_tmp]
    location_mut_neu_row <- location_mut_neu_row[order(location_mut_neu_row)]
    mut_tmp <- mut_tmp[-neu_tmp]
    }
    
    if(loci_mut_del>0){
    del_tmp <- sample(1:length(mut_tmp),size=loci_mut_del)
    location_mut_del_row <- mut_tmp[del_tmp]
    location_mut_del_row <- location_mut_del_row[order(location_mut_del_row)]
    mut_tmp <- mut_tmp[-del_tmp]
    }
    
    if(loci_mut_adv>0){
    adv_tmp <- sample(1:length(mut_tmp),size=loci_mut_adv)
    location_mut_adv_row <- mut_tmp[adv_tmp]
    location_mut_adv_row <- location_mut_adv_row[order(location_mut_adv_row)]
    mut_tmp <- mut_tmp[-adv_tmp]
    }
    
    s <- rep(NA,total_loci)
    h <- rep(NA,total_loci)
    q <- rep(NA,total_loci)
    
    ##### SELECTION COEFFICIENT
    # neutral loci simulations
    if(chunk_neutral_loci>0){
      s[location_neutral_row] <- 0
    }
    #real locations
    if(real_loc == TRUE){
    s[location_real_row] <- 0
    }
    # deleterious
    if(loci_deleterious>0){
    if (s_distribution_del == "equal") {
      s[location_deleterious_row] <- s_del
    }
    
    if (s_distribution_del == "gamma") {
      s[location_deleterious_row] <- rgamma(loci_deleterious, 
                                            shape = gamma_shape, 
                                            scale = gamma_scale)
    }
    
    if (s_distribution_del == "log_normal") {
      s[location_deleterious_row] <- rlnorm(loci_deleterious,
                  meanlog = log(log_mean),
                  sdlog = log(log_sd))
    }
    }
    #advantageous
    if(loci_advantageous>0){
      if (s_distribution_adv == "equal") {
        s[location_advantageous_row] <- s_adv * -1
      }
      if(s_distribution_adv=="exponential"){
        s[location_advantageous_row] <- rexp(loci_advantageous,
                                             rate = exp_rate) * -1
      }
    }
    # neutral mutations 
    if(loci_mut_neu>0){
      s[location_mut_neu_row] <- 0
    }
    # deleterious mutations
    if(loci_mut_del>0){
      if (s_distribution_del == "equal") {
        s[location_mut_del_row] <- s_del
      }
      
      if (s_distribution_del == "gamma") {
        s[location_mut_del_row] <- rgamma(loci_mut_del, 
                                              shape = gamma_shape, 
                                              scale = gamma_scale)
      }
      
      if (s_distribution_del == "log_normal") {
        s[location_mut_del_row] <- rlnorm(loci_mut_del,
                                              meanlog = log(log_mean),
                                              sdlog = log(log_sd))
      }
    }
    # advantageous mutations
    if(loci_mut_adv>0){
      if (s_distribution_adv == "equal") {
        s[location_mut_adv_row] <- s_adv * -1
      }
      if(s_distribution_adv=="exponential"){
        s[location_mut_adv_row] <- rexp(loci_mut_adv,
                                             rate = exp_rate) * -1
      }
      
    }

    ##### DOMINANCE
    # neutral loci simulations
    if(chunk_neutral_loci>0){
      h[location_neutral_row] <- 0
    }
    #real locations
    if(real_loc == TRUE){
      h[location_real_row] <- 0
    }
    # deleterious
    if(loci_deleterious>0){
      if (h_distribution_del == "equal") {
        h[location_deleterious_row] <- h_del
      }
      if (h_distribution_del == "normal") {
        h[location_deleterious_row] <- rnorm(loci_deleterious, 
                                             mean = h_mean_del,
                                             sd = h_sd_del)
      }
      # the equation for dominance (h) was taken from Huber 2018 Nature
      if (h_distribution_del == "equation") {
        h[location_deleterious_row] <- 1 / ((1 / h_intercept_del) - (-1 * h_rate_del * abs(s[location_deleterious_row])))
      }
    }
    #advantageous
    if(loci_advantageous>0){
      if (h_distribution_adv == "equal") {
        h[location_advantageous_row] <- h_adv
      }
      if (h_distribution_adv == "normal") {
        h[location_advantageous_row] <- rnorm(loci_advantageous, 
                                             mean = h_mean_adv,
                                             sd = h_sd_adv)
      }
      # the equation for dominance (h) was taken from Huber 2018 Nature
      if (h_distribution_adv == "equation") {
        h[location_advantageous_row] <- 1 / ((1 / h_intercept_adv) - (-1 * h_rate_adv * abs(s[location_advantageous_row])))
      }
    }
    # neutral mutations 
    if(loci_mut_neu>0){
      h[location_mut_neu_row] <- 0
    }
    # deleterious mutations
    if(loci_mut_del>0){
      if (h_distribution_del == "equal") {
        h[location_mut_del_row] <- h_del
      }
      if (h_distribution_del == "normal") {
        h[location_mut_del_row] <- rnorm(loci_mut_del, 
                                             mean = h_mean_del,
                                             sd = h_sd_del)
      }
      # the equation for dominance (h) was taken from Huber 2018 Nature
      if (h_distribution_del == "equation") {
        h[location_mut_del_row] <- 1 / ((1 / h_intercept_del) - (-1 * h_rate_del * abs(s[location_mut_del_row])))
      }
    }
    # advantageous mutations
    if(loci_mut_adv>0){
      if (h_distribution_adv == "equal") {
        h[location_mut_adv_row] <- h_adv
      }
      if (h_distribution_adv == "normal") {
        h[location_mut_adv_row] <- rnorm(loci_mut_adv, 
                                              mean = h_mean_adv,
                                              sd = h_sd_adv)
      }
      # the equation for dominance (h) was taken from Huber 2018 Nature
      if (h_distribution_adv == "equation") {
        h[location_mut_adv_row] <- 1 / ((1 / h_intercept_adv) - (-1 * h_rate_adv * abs(s[location_mut_adv_row])))
      }
    }

    ###### INITIAL FREQUENCY
    # neutral loci simulations
    if(chunk_neutral_loci>0){
      q[location_neutral_row] <- q_neutral
    }
    # real locations 
    # the initial frequency of real locations is set in the function
    # gl.sim.WF.run
    if(real_loc==TRUE){
      q[location_real_row] <- q_neutral
    }
    # deleterious
    if(loci_deleterious>0){
      if (q_distribution_del == "equal") {
        q[location_deleterious_row] <- q_del
      }
      if (q_distribution_del == "equation") {
        a <- abs(s[location_deleterious_row]) * 
          (1 - (2 * h[location_deleterious_row]))
        b <- (h[location_deleterious_row] * 
                abs(s[location_deleterious_row])) * 
          (1 + q_equation_del)
        c <- rep.int(-(q_equation_del), times = loci_deleterious)
        df_q <- as.data.frame(cbind(a, b, c))
        # q is based on the following equation: (s(1-2h)q^2) + (hs(1+u)q) - u = 0,
        # where u is the mutation rate per generation per site. Taken from Crow &
        # Kimura page 260
        q[location_deleterious_row] <-
          mapply(
            q_equilibrium,
            a = df_q$a,
            b = df_q$b,
            c = df_q$c,
            USE.NAMES = F
          )
      }
    }
    # advantageous
    if(loci_advantageous>0){
      if (q_distribution_adv == "equal") {
        q[location_advantageous_row] <- q_adv
      }
      if (q_distribution_adv == "equation") {
        a <- abs(s[location_advantageous_row]) * 
          (1 - (2 * h[location_advantageous_row]))
        b <- (h[location_advantageous_row] * 
                abs(s[location_advantageous_row])) * 
          (1 + q_equation_adv)
        c <- rep.int(-(q_equation_adv), times = loci_advantageous)
        df_q <- as.data.frame(cbind(a, b, c))
        # q is based on the following equation: (s(1-2h)q^2) + (hs(1+u)q) - u = 0,
        # where u is the mutation rate per generation per site. Taken from Crow &
        # Kimura page 260
        q[location_advantageous_row] <-
          mapply(
            q_equilibrium,
            a = df_q$a,
            b = df_q$b,
            c = df_q$c,
            USE.NAMES = F
          )
      }
    }
    #neutral mutations
    if(loci_mut_neu>0){
      q[location_mut_neu_row] <- 0
    }
    if(loci_mut_del>0){
      q[location_mut_del_row] <- 0
    }
    if(loci_mut_adv>0){
      q[location_mut_adv_row] <- 0 
    }
    
    # Producing reference table 
    reference <- as.data.frame(matrix(nrow = total_loci))
    reference$q <- q
    reference$h <- h
    reference$s <- s
    reference$c <- recombination_map[1:total_loci, "c"]
    reference$loc_bp <- recombination_map[1:total_loci, "location_loci_bp"]
    reference$loc_cM <- recombination_map[1:total_loci, "accum"]
    reference$chr_name <- chromosome_name
    reference$type <- NA
    
    reference <- reference[, -1]
    
    if(real_loc==TRUE){
      reference[location_real_row, "type"] <- "real"
    }
    if(chunk_neutral_loci>0){
      reference[location_neutral_row, "type"] <- "neutral"
    }
    if(loci_deleterious>0){
      reference[location_deleterious_row, "type"] <- "deleterious"
    }
    if(loci_advantageous>0){
      reference[location_advantageous_row, "type"] <- "advantageous"
    }
    if(loci_mut_neu>0){
      reference[location_mut_neu_row, "type"] <- "mutation_neu"
    }
    if(loci_mut_del>0){
      reference[location_mut_del_row, "type"] <- "mutation_del"
    }
    if(loci_mut_adv>0){
      reference[location_mut_adv_row, "type"] <- "mutation_adv"
    }
    
    # NS with very small s have a q > 1. Therefore, we set a maximum q value of
    # 0.5.
    if(q_distribution_del != "equal"){
    q_more_than_point5 <- as.numeric(row.names(reference[reference$q > 0.5, ]))
    reference[q_more_than_point5, "q"] <- 0.5
    }
    
    # the log normal distribution, with the parameters used in the simulations,
    # generates a few selection coefficients that are > 1. The maximum value of s
    # is set to 0.99
    if(s_distribution_del != "equal"){
    s_more_than_one <- as.numeric(row.names(reference[reference$s > 1, ]))
    reference[s_more_than_one, "s"] <- 0.99
    # the exponential distribution, with the parameters used in the simulations,
    # generates a few selection coefficients that are < -1. The minimum value of s
    # is set to -0.5
    s_less_minus_one <- as.numeric(row.names(reference[reference$s < - 0.5, ]))
    reference[s_less_minus_one, "s"] <- -0.5
    }
    
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
