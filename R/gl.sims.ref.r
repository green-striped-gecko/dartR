
gl.sims.ref <- function(file_ref=system.file('extdata','ref_variables.csv',package='dartR'),
                        file_targets_sel = NULL,
                        file_r_map = NULL,
                        verbose = NULL){
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jody",
                   verbosity = verbose)
  
  # DO THE JOB
  
  ##### SIMULATIONS VARIABLES ######
  ref_vars <- suppressWarnings(read.csv(file_ref))
  cat(important("  Use the second column ('value') to change the values of the 
                simulation variables and then clic on the quit icon\n"))
  cat(important("  The description of the variables and their range of values 
                are shown in the fourth and fifth colums respectively\n"))
  ref_vars <- edit(ref_vars)
  vars_assign <- unlist(unname(mapply(paste, ref_vars$variable, "<-", 
                                      ref_vars$value, SIMPLIFY = F)))
  eval(parse(text = vars_assign))

###### recombination map
  
  ###### input targets of selection
  if(!is.null(file_targets_sel)){
    targets <- read.csv(file_targets_sel)
    targets$chr_name <- as.character(targets$chr_name)
    targets <- targets[which(targets$chr_name == chromosome_name), ]
    targets <- targets[!duplicated(targets$start), ]
    targets <- targets[!duplicated(targets$end), ]
  }
  
if(!is.null(file_r_map)){
  map <- read_csv(file_r_map)
  map$Chr <- as.character(map$Chr)
  map <- map[which(map$Chr == chromosome_name), ]
  map <- as.data.frame(map$cM / 1000)
  colnames(map) <- "cM"
  map[is.na(map$cM), ] <- 0
  chr_length <- (nrow(map) + 1) * map_resolution
  recombination_map_temp <- map
  location_targets <- NULL
  targets$targets <- targets$ns * targets_factor
    
    transcripts <- as.data.frame(targets)
    for (i in 1:nrow(transcripts)) {
      location_targets_temp <-
        unlist(mapply(
          FUN = function(a, b) {
            seq(from = a,
                to = b,
                by = 4)
          },
          a = unname(unlist(transcripts[i, "start"])),
          b = unname(unlist(transcripts[i, "end"]))
        ))
      location_targets_temp <-
        sample(location_targets_temp, size = transcripts[i, "targets"])
      location_targets <-
        c(location_targets, location_targets_temp)
    }
    # here are added the location of the loci genotyped in the fly experiment
    if (experiment_loci == T) {
      # these are the location of the neutral loci in the simulations
      location_neutral_loci_analysis <-
        c(seq(
          map_resolution / 2,
          (nrow(recombination_map_temp) * map_resolution),
          map_resolution
        ),
        location_msats_experiment)
      location_neutral_loci_analysis <-
        location_neutral_loci_analysis[order(location_neutral_loci_analysis)]
      location_targets <-
        c(location_targets, location_neutral_loci_analysis)
    }
    if (experiment_loci == F) {
      location_neutral_loci_analysis <-
        seq(map_resolution / 2,
            (nrow(recombination_map_temp) * map_resolution),
            map_resolution)
      location_targets <-
        c(location_targets, location_neutral_loci_analysis)
    }
    location_targets <-
      location_targets[order(location_targets)]
    # different transcripts can be located in the same genome location So,
    # repeated deleterious mutations are deleted, however transcripts that are in
    # the same place are taken in account to calculate fitness in the fitness 
    # function
    location_targets <- unique(location_targets)
    
    if (experiment_freq == T) {
      loc_exp_loci <- location_neutral_loci_analysis
      loc_exp_loci <- loc_exp_loci[order(loc_exp_loci)]
      loc_exp_loci <- which(loc_exp_loci %% 10000 != 0)
      loc_exp_loci_2 <-
        unlist(lapply(location_msats_experiment, function(x) {
          which(x == location_targets)
        }))
    }
    
    chromosome_length <-
      (nrow(recombination_map_temp) + 1) * map_resolution
    #this is to fix a bug that crashes the program because the last neutral 
    # locus sometimes could be located farther than the last deleterious mutation
    location_targets <- c(location_targets, chromosome_length)
    loci_number_to_simulate <- length(location_targets)
    # the recombination map is produced by cross multiplication. the following 
    # lines are the input for doing the cross multiplication.
    recombination_map_temp$midpoint <-
      seq(
        map_resolution / 2,
        nrow(recombination_map_temp) * map_resolution,
        map_resolution
      )
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
    recombination_map[loci_number + 1, "c"] <- recom_event - sum(recombination_map[, "c"])
    recombination_map[loci_number + 1, "location_targets"] <- recombination_map[loci_number, "location_targets"]
    recombination_map[loci_number + 1, "recombination_2"] <- recombination_map[loci_number, "recombination_2"]
    recombination_map$accum <- cumsum(recombination_map[, "c"])
    
    location_temp <- location_neutral_loci_analysis
    location_temp <- location_temp[order(location_temp)]
    location <-
      lapply(location_temp, function(x) {
        which(recombination_map$location_targets == x)
      })
    neutral_loci_location <-  as.character(location)
    loci_location <- location
    
    s <- rlnorm(loci_number,
                meanlog = log(log_mean),
                sdlog = log(log_sd))
    # the equation for dominance (h) was taken from Huber 2018 Nature
    h <- rnorm(loci_number, mean = dominance_mean, sd = sqrt(0.001))
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
    reference <- as.data.frame(matrix(nrow = loci_number))
    reference$q <- q
    reference$h <- h
    reference$s <- s
    # NS with very small s have a q > 1. Therefore, we set a maximum q value of 
    # 0.5.
    q_more_than_point5 <-
      as.numeric(row.names(reference[reference$q > 0.5,]))
    reference[q_more_than_point5, "q"] <- 0.5
    # the log normal distribution, with the parameters used in the simulations,
    # generates a few selection coefficients that are > 1. the maximum value of s 
    # is set to 0.99
    s_more_than_one <-
      as.numeric(row.names(reference[reference$s > 1,]))
    reference[s_more_than_one, "s"] <- 0.99
    loci_number <- nrow(reference)
    # setting h and s to 0 in neutral loci
    reference[as.numeric(neutral_loci_location),"s"] <- 0
    reference[as.numeric(neutral_loci_location),"h"] <- 0
    reference$location <- recombination_map[1:loci_number, "location_targets"]
    
    # one is subtracted from the recombination map to account for the last row that
    # was added in the recombination map to avoid that the recombination function crashes
    plink_map <- as.data.frame(matrix(nrow = loci_number, ncol = 4))
    plink_map[, 1] <- chromosome_name
    plink_map[, 2] <- rownames(recombination_map[-(loci_number + 1), ])
    plink_map[, 3] <- recombination_map[-(loci_number + 1), "accum"]
    plink_map[, 4] <- recombination_map[-(loci_number + 1), "location_targets"]
    
}else{
    chromosome_length <- chunk_number * map_resolution
    # The recombination map is generated by creating a table with as many rows 
    # as loci_number_to_simulate and one column, which is filled with the 
    # recombination rate between loci (chunk_recombination).
    map <- as.data.frame(matrix(nrow = chunk_number))
    map[, 1] <- chunk_recombination / 100
    map[, 2] <- 1
    colnames(map) <- c("cM", "loc_under_sel")
    # Adjusting the total number of loci to simulate by adding one more locus to
    # each window because each window (row) contains one SNP. This is done 
    # with the aim that the number of adaptive loci and the number of neutral 
    # loci are completely independent from each other. The neutral locus is 
    # located at the middle of the window
    location_neutral_loci_analysis <-
      seq(map_resolution / 2, (nrow(map) * map_resolution), map_resolution)
    # the number of loci under selection and the number of neutral loci are 
    # added to obtain the total number of loci per each chunk_number. 
    loci_number_to_simulate <- loci_number_to_simulate + nrow(map)
    map$loci_density <-
      ceiling(map$loc_under_sel / (sum(map$loc_under_sel) / loci_number_to_simulate))
    # To determine the location of each locus in the chromosome, first the cM’s 
    # and map_resolution are divided between the total number of loci. Then, each 
    # locus is placed along the chromosome sequentially using the accumulative 
    # sum of cM’s and number of base pairs for each locus.
    recombination_map <- NULL
    for (value in 1:nrow(map)) {
      temp <- NULL
      temp <-
        as.data.frame(rep((map[value, "cM"] / map[value, "loci_density"]),
                          map[value, "loci_density"]))
      distance <- map_resolution / nrow(temp)
      temp$loc_bp <- distance * (1:nrow(temp))
      if (value == 1) {
        temp$loc_bp <- ceiling(temp$loc_bp)
      } else{
        temp$loc_bp <-
          ceiling(temp$loc_bp + ((value * map_resolution) - map_resolution))
      }
      recombination_map <- rbind(recombination_map, temp)
    }
    loci_number <- nrow(recombination_map)
    recombination_map$accum <- cumsum(recombination_map[, 1])
    colnames(recombination_map) <-
      c("c", "location_loci", "accum")
    
    location <-
      lapply(location_neutral_loci_analysis,
             findInterval,
             vec = as.numeric(paste(unlist(
               recombination_map$location_loci
             ))))
    location[[1]] <- 1
    
    neutral_loci_location <-  as.character(location)
    
    loci_number <- nrow(recombination_map)
    reference <- as.data.frame(matrix(nrow = loci_number))
    reference$q <- q_gral
    reference$h <- h_gral
    reference$s <- s_gral
    reference$c <- recombination_map[1:loci_number, "c"]
    reference$loc_bp <- recombination_map[1:loci_number, "location_loci"]
    reference$loc_cM <- recombination_map[1:loci_number, "accum"]
    reference$chr_name <- chromosome_name
    # setting h and s to 0 in neutral loci
    reference[as.numeric(neutral_loci_location),"s"] <- 0
    reference[as.numeric(neutral_loci_location),"h"] <- 0
    reference[as.numeric(neutral_loci_location),"q"] <- q_neutral
    
}

return(reference[,-1])

}





