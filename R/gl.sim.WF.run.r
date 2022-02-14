#' @name gl.sim.WF.run
#' @title Wright-Fisher simulations
#' @description
#' Forward simulattions
#' @param file_var [required if interactive_vars = TRUE].
#' @param ref_table [required].
#' @param x Name of the genlight object containing the SNP data [default NULL].
#' @param number_iterations [default 1].
#' @param every_gen [default 5].
#' @param store_phase1 [default FALSE].
#' @param interactive_vars [default TRUE].
#' @param seed [default NULL].
#' @param parallel [default FALSE].
#' @param n.cores [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' 
#' @details
#' 
#' system.file('extdata', 'sim_variables.csv', package ='dartR')
#' 
#' @return Returns unaltered genlight object
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' 
#' @seealso \code{\link{gl.sim.WF.table}}
#' @family simulation functions
#' @import shiny
#' @import shinyBS
#' @import shinythemes
#' @import shinyjs
#' @export

gl.sim.WF.run <-
  function(file_var,
           ref_table,
           x = NULL,
           number_iterations = 1,
           every_gen = 5,
           store_phase1 = FALSE,
           interactive_vars = TRUE,
           seed = NULL,
           parallel = FALSE,
           n.cores = NULL,
           verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbosity = verbose)
    
    # check if package is installed
    pkg <- "data.table"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      stop(error(
        "Package",
        pkg,
        "needed for this function to work. Please install it."
      ))
    }
    
    if (interactive_vars) {
      pkg <- "shiny"
      if (!(requireNamespace(pkg, quietly = TRUE))) {
        stop(error(
          "Package",
          pkg,
          "needed for this function to work. Please install it."
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
      interactive_sim_run( env_fun = environment())
      
    } else {
      sim_vars <- suppressWarnings(read.csv(file_var))
      sim_vars <- sim_vars[, 2:3]
      vars_assign <-
        unlist(unname(
          mapply(paste, sim_vars$variable, "<-",
                 sim_vars$value, SIMPLIFY = F)
        ))
      eval(parse(text = vars_assign))
    }
    
    reference <- ref_table$reference
    ref_vars <- ref_table$ref_vars
    
    q_neutral <- as.numeric(ref_vars[ref_vars$variable=="q_neutral","value"])
    
    # setting the seed
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    if (phase1 == FALSE) {
      gen_number_phase1 <- 0
    }
    
    # This is the total number of generations
    number_generations <-
      gen_number_phase1 + gen_number_phase2
    
    # This is the list to store the final genlight objects
    gen_store <-
      c(seq(1, number_generations, every_gen), number_generations)
    final_res <-
      rep(list(as.list(rep(
        NA, length(gen_store)
      ))), number_iterations)
    
    loci_number <- nrow(reference)
    recombination_map <- reference[, c("c", "loc_bp", "loc_cM")]
    # The last element of the first column must be zero, otherwise the
    # recombination function crashes.
    recombination_map[nrow(recombination_map), 1] <- 0
    # In order for the recombination rate to be accurate, we must account for
    # the case when the probability of the total recombination rate is less than
    # 1 (i.e. < 100 cM) or more than 1 (> 100 cM). For the first case, the program
    # subtracts from 1 the sum of all the recombination rates and this value
    # inserted in the last row of the recombination_map table. If this row is
    # chosen as the recombination point, recombination does not occur. For
    # example, if a chromosome of 20 cM’s is simulated, the last row of the
    # recombination_map will have a value of 0.8 and therefore 80% of the times
    # recombination will not occur. For the second case, having more than 100 cM,
    # means that more than 1 recombination event occurs. So, one recombination
    # event is perform for each 100 cM. Then, the program subtracts the number of
    # recombination events from the sum of all the recombination rates and this
    # value inserted in the last row of the recombination_map table, in the same
    # way as in the first case.
    # number of recombination events per meiosis
    recom_event <- ceiling(sum(recombination_map[, "c"]))
    # filling the probability of recombination when the total recombination rate
    # is less than an integer (recom_event) and placing it at the end of the
    # recombination map
    
    recombination_map$loc_cM <- cumsum(recombination_map$c)
    
    recombination_map[loci_number + 1, 1] <-
      recom_event - sum(recombination_map[, 1])
    recombination_map[loci_number + 1, 2] <-
      recombination_map[loci_number, 2]
    recombination_map[loci_number + 1, 3] <-
      recombination_map[loci_number, 3]
    
    neutral_loci_location <- which(reference$s == 0)
    
    reference$selection <- "under_selection"
    reference[as.numeric(neutral_loci_location), "selection"] <-
      "neutral"
    
    # one is subtracted from the recombination map to account for the last row that
    # was added in the recombination map to avoid that the recombination function crashes
    plink_map <-
      as.data.frame(matrix(nrow = nrow(reference), ncol = 4))
    plink_map[, 1] <- reference$chr_name
    plink_map[, 2] <- rownames(reference)
    plink_map[, 3] <- reference$loc_cM
    plink_map[, 4] <- reference$loc_bp
    
    population_size_phase2 <-
      as.numeric(unlist(strsplit(population_size_phase2, " ")))
    
    population_size_phase1 <-
      as.numeric(unlist(strsplit(population_size_phase1, " ")))
    
    if(phase1 == TRUE & number_pops_phase1!=number_pops_phase2 ){
      cat(error("Number of populations in phase 1 and phase 2 must be the same\n"))
      stop()
    }
    
    if(length(population_size_phase2)!=number_pops_phase2){
      cat(error("Number of entries for population sizes do not agree with the number of populations for phase 2\n"))
      stop()
    }
    
    if(length(population_size_phase1)!=number_pops_phase1 & phase1==TRUE){
      cat(error("Number of entries for population sizes do not agree with the number of populations for phase 1\n"))
      stop()
    }
    
    
    
    if (phase1 == TRUE & real_pops == FALSE) {
      number_pops <- number_pops_phase1
    }
    if (phase1 == TRUE & real_pops == TRUE & !is.null(x)) {
      number_pops <- nPop(x)
    }
    if (phase1 == FALSE & real_pops == FALSE) {
      number_pops <- number_pops_phase2
    }
    if (phase1 == FALSE & real_pops == TRUE & !is.null(x)) {
      number_pops <- nPop(x)
    }
    
    if (real_freq == TRUE & !is.null(x)) {
      pop_list_freq_temp <- seppop(x)
      loc_to_keep <-
        locNames(pop_list_freq_temp[[1]])[which(pop_list_freq_temp[[1]]$chromosome == chromosome_name)]
      pop_list_freq_temp <-
        lapply(pop_list_freq_temp,
               gl.keep.loc,
               loc.list = loc_to_keep,
               verbose = 0)
      pop_list_freq <- lapply(pop_list_freq_temp, gl.alf)
      
    } else{
      pop_list_freq <- rep(NA, number_pops)
    }
    

    
    ##### ANALYSIS VARIABLES #####
    # This is to calculate the density of mutations per centimorgan. The density
    # is based on the number of heterozygous loci in each individual. Based on HW
    # equation (p^2+2pq+q^2), the proportion of heterozygotes (2pq) for each locus
    # is calculated and then averaged. This proportion is then multiplied by the
    # number of loci and divided by the length of the chromosome in centiMorgans.
    # According to Haddrill 2010, the mean number of heterozygous deleterious
    # mutations per fly is 5,000 deleterious mutations per individual, with an
    # estimated mean selection coefficient (sh) of 1.1X10^-5. The chromosome
    # arm 2L has 17% of the total number of non-synonymous mutations and is 55/2
    # cM long (cM are divided by two because there is no recombination in males),
    # with these parameters the density per cM is (5000*0.17)/(55/2) = 30.9
    freq_deleterious <-
      reference[-as.numeric(neutral_loci_location),]
    freq_deleterious_b <-
      mean(2 * (freq_deleterious$q) * (1 - freq_deleterious$q))
    density_mutations_per_cm <-
      (freq_deleterious_b * nrow(freq_deleterious)) /
      (recombination_map[loci_number, "loc_cM"] * 100)
    
    dispersal_rate_phase2 <-
      (number_transfers_phase2 / transfer_each_gen_phase2) / (population_size_phase2)
    Fst_expected_phase2 <-
      1 / ((4 * Ne_fst_phase2 * dispersal_rate_phase2) * ((2 / (2 - 1)) ^ 2) + 1)
    mi_expected_phase2 <-
      (0.22 / (sqrt(
        2 * Ne_fst_phase2 * dispersal_rate_phase2
      ))) - (0.69 / ((2 * Ne_fst_phase2) * sqrt(dispersal_rate_phase2)))
    rate_of_loss_phase2 <- 1 - (1 / (2 * Ne_phase2))
    
    dispersal_rate_phase1 <-
      (number_transfers_phase1 / transfer_each_gen_phase1) / (population_size_phase1)
    Fst_expected_phase1 <-
      1 / ((4 * Ne_fst_phase1 * dispersal_rate_phase1) * ((2 / (2 - 1)) ^ 2) + 1)
    mi_expected_phase1 <-
      (0.22 / (sqrt(
        2 * Ne_fst_phase1 * dispersal_rate_phase1
      ))) - (0.69 / ((2 * Ne_fst_phase1) * sqrt(dispersal_rate_phase1)))
    rate_of_loss_phase1 <- 1 - (1 / (2 * Ne_phase1))
    
    ##### START ITERATION LOOP #####
    for (iteration in 1:number_iterations) {
      if (iteration %% 1 == 0) {
        cat(report(" Starting iteration =", iteration, "\n"))
      }
      
      ##### VARIABLES PHASE 1 #######
      if (phase1 == TRUE) {
        if (real_pop_size == TRUE & !is.null(x)) {
          population_size_phase1 <- unname(unlist(table(pop(x))))
          # converting odd population sizes to even
          population_size_phase1 <-
            (population_size_phase1 %% 2 != 0) + population_size_phase1
        } else{
          population_size <- population_size_phase1
        }
        
        selection <- selection_phase1
        dispersal <- dispersal_phase1
        dispersal_type <- dispersal_type_phase1
        number_transfers <- number_transfers_phase1
        transfer_each_gen <- transfer_each_gen_phase1
        variance_offspring <- variance_offspring_phase1
        number_offspring <- number_offspring_phase1
        
        store_values <- store_phase1
        
        if(store_phase1==TRUE){
        # counter to store values every generation
        gen <- 0
        }
        
        # pick which sex is going to be transferred first
        if (number_transfers >= 2) {
          maletran <- TRUE
          femaletran <- TRUE
        } else if (number_transfers == 1) {
          maletran <- TRUE
          femaletran <- FALSE
        }
        
      } else{
        if (real_pop_size == TRUE) {
          population_size_phase2 <- unname(unlist(table(pop(x))))
          population_size_phase2 <-
            (population_size_phase2 %% 2 != 0) + population_size_phase2
        } else{
          population_size <- population_size_phase2
        }
        
      }
      
      # tic("initialisation")
      ##### INITIALISE POPS #####
      if (verbose >= 2) {
        cat(report("  Initialising populations\n"))
      }
      
      pops_vector <- 1:number_pops
      pop_list <- lapply(pops_vector, function(y) {
        initialise(
          pop_number = y,
          pop_size = population_size[y],
          refer = reference,
          q_neu = q_neutral,
          n_l_loc = neutral_loci_location,
          r_freq = pop_list_freq[[y]]
        )
      })
      # toc()
      
      #if there is just one population set dispersal to FALSE
      if (length(pop_list) == 1) {
        dispersal <- FALSE
      }
      
      ##### START GENERATION LOOP #####
      for (generation in 1:number_generations) {
        if (phase1 == TRUE & generation == 1) {
          cat(report(" Starting phase 1\n"))
        }
        
        if (generation %% 5 == 0) {
          cat(report("  Starting generation =", generation, "\n"))
        }
        
        ##### VARIABLES PHASE 2 #######
        if (generation == (gen_number_phase1 + 1)) {
          cat(report(" Starting phase 2\n"))
          
          selection <- selection_phase2
          dispersal <- dispersal_phase2
          dispersal_type <- dispersal_type_phase2
          number_transfers <- number_transfers_phase2
          transfer_each_gen <- transfer_each_gen_phase2
          variance_offspring <- variance_offspring_phase2
          number_offspring <- number_offspring_phase2
          dispersal <- dispersal_phase2
          
          store_values <- TRUE
          
          # pick which sex is going to be transferred first
          if (number_transfers >= 2) {
            maletran <- TRUE
            femaletran <- TRUE
          } else if (number_transfers == 1) {
            maletran <- TRUE
            femaletran <- FALSE
          }
          
          # counter to store genlight objects
          count_store <- 0
          # counter to store values every generation
          gen <- 0
          
          if (phase1 == TRUE) {
            if (same_line == TRUE) {
              # pop_list_temp is used because pop_list is used to sample populations
              pop_sample <- sample(pops_vector, 1)
              
              pop_list_temp <- lapply(pops_vector, function(x) {
                pop_temp <-
                  rbind(pop_list[[pop_sample]][sample(
                    which(pop_list[[pop_sample]]$V1 == "Male"),
                    size =  population_size[x] / 2,
                    replace = T
                  ),],
                  pop_list[[pop_sample]][sample(
                    which(pop_list[[pop_sample]]$V1 == "Female"),
                    size = population_size[x] / 2,
                    replace = T
                  ),])
                pop_temp$V2 <- x
                return(pop_temp)
              })
              pop_list <- pop_list_temp
            }
            
            if (same_line == FALSE) {
              pop_list <- lapply(pops_vector, function(x) {
                pop_temp <-
                  rbind(pop_list[[x]][sample(
                    which(pop_list[[x]]$V1 == "Male"),
                    size =  population_size[x] / 2,
                    replace = T
                  ),],
                  pop_list[[x]][sample(
                    which(pop_list[[x]]$V1 == "Female"),
                    size = population_size[x] / 2,
                    replace = T
                  ),])
                pop_temp$V2 <- x
                return(pop_temp)
              })
              
            }
          }
        }
        
        # generation counter
        if (store_values == TRUE) {
          gen <- gen + 1
        }
        
        # tic("dispersal")
        ##### DISPERSAL ######
        # dispersal is symmetrical
        if (dispersal == TRUE) {
          if (dispersal_type == "all_connected") {
            dispersal_pairs <-
              as.data.frame(expand.grid(pops_vector, pops_vector))
            dispersal_pairs$same_pop <-
              dispersal_pairs$Var1 == dispersal_pairs$Var2
            dispersal_pairs <-
              dispersal_pairs[which(dispersal_pairs$same_pop == FALSE),]
            colnames(dispersal_pairs) <-
              c("pop1", "pop2", "same_pop")
          }
          
          if (dispersal_type == "line") {
            dispersal_pairs <- as.data.frame(rbind(
              cbind(head(pops_vector,-1), pops_vector[-1]),
              cbind(pops_vector[-1], head(pops_vector,-1))
            ))
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
          
          # defining the population size of each population
          if (real_pop_size == TRUE) {
            dispersal_pairs$size_pop1 <-
              unname(unlist(table(pop(x))))[dispersal_pairs$pop1]
            dispersal_pairs$size_pop2 <-
              unname(unlist(table(pop(x))))[dispersal_pairs$pop2]
          } else{
            dispersal_pairs$size_pop1 <- population_size[dispersal_pairs$pop1]
            dispersal_pairs$size_pop2 <-
              population_size[dispersal_pairs$pop2]
          }
          
          for (dis_pair in 1:nrow(dispersal_pairs)) {
            res <- migration(
              population1 = pop_list[[dispersal_pairs[dis_pair, "pop1"]]],
              population2 = pop_list[[dispersal_pairs[dis_pair, "pop2"]]],
              gen = generation,
              size_pop1 = dispersal_pairs$size_pop1[dis_pair],
              size_pop2 = dispersal_pairs$size_pop2[dis_pair],
              trans_gen = transfer_each_gen,
              male_tran = maletran,
              female_tran = femaletran,
              n_transfer = number_transfers
            )
            
            pop_list[[dispersal_pairs[dis_pair, "pop1"]]] <-
              res[[1]]
            pop_list[[dispersal_pairs[dis_pair, "pop1"]]]$V2 <-
              dispersal_pairs[dis_pair, "pop1"]
            pop_list[[dispersal_pairs[dis_pair, "pop2"]]] <-
              res[[2]]
            pop_list[[dispersal_pairs[dis_pair, "pop2"]]]$V2 <-
              dispersal_pairs[dis_pair, "pop2"]
            maletran <- res[[3]]
            femaletran <- res[[4]]
          }
        }
        # toc()
        
        # tic(reproduction)
        ##### REPRODUCTION #########
        offspring_list <- lapply(pops_vector, function(x) {
          reproduction(
            pop = pop_list[[x]],
            pop_number = x,
            pop_size = population_size[x],
            var_off = variance_offspring,
            num_off = number_offspring,
            r_event = recom_event,
            recom = recombination,
            r_males = recombination_males,
            r_map_1 = recombination_map,
            n_loc = loci_number
          )
        })
        # toc()
        
        # tic("selection")
        ##### SELECTION #####
        if (selection == TRUE) {
          offspring_list <- lapply(pops_vector, function(x) {
            selection_fun(
              offspring = offspring_list[[x]],
              reference_pop = reference,
              sel_model = natural_selection_model,
              g_load = genetic_load
            )
          })
        }
        # toc()
        
        # tic("sampling_next_gen")
        ##### SAMPLING NEXT GENERATION ########
        # testing whether any population became extinct, if so break the
        test_extinction <- unlist(lapply(pops_vector, function(x) {
          length(which(offspring_list[[x]]$V1 == "Male")) < population_size / 2 |
            length(which(offspring_list[[x]]$V1 == "Female")) < population_size / 2
        }))
        
        if (any(test_extinction == TRUE)) {
          pops_extinct <- which(test_extinction == TRUE)
          cat(
            important(
              "  Population",
              pops_extinct,
              "became EXTINCT at generation",
              generation,
              "\n"
            )
          )
          cat(important(
            "  Breaking this iteration and passing to the next iteration",
            "\n"
          ))
          
          s_vars_temp <- sim_vars
          s_vars_temp <-
            setNames(data.frame(t(s_vars_temp[,-1])), s_vars_temp[, 1])
          s_vars_temp$generation <- generation
          s_vars_temp$iteration <- iteration
          s_vars_temp$seed <- seed
          s_vars_temp$del_ind_cM <- density_mutations_per_cm
          
          s_vars_temp$dispersal_rate_phase1 <-
            paste(dispersal_rate_phase1, collapse = " ")
          s_vars_temp$Fst_expected_phase1 <-
            paste(Fst_expected_phase1, collapse = " ")
          s_vars_temp$mi_expected_phase1 <-
            paste(mi_expected_phase1, collapse = " ")
          s_vars_temp$rate_of_loss_phase1 <- rate_of_loss_phase1
          
          s_vars_temp$dispersal_rate_phase2 <-
            paste(dispersal_rate_phase2, collapse = " ")
          s_vars_temp$Fst_expected_phase2 <-
            paste(Fst_expected_phase2, collapse = " ")
          s_vars_temp$mi_expected_phase2 <-
            paste(mi_expected_phase2, collapse = " ")
          s_vars_temp$rate_of_loss_phase2 <- rate_of_loss_phase2
          
          final_res[[iteration]][[count_store]] <-
            store(
              p_vector = pops_vector,
              p_size = population_size,
              p_list = pop_list,
              n_loc_1 = loci_number,
              paral = parallel,
              n_cores = n.cores,
              ref = reference,
              p_map = plink_map,
              s_vars = s_vars_temp
            )
          
          break()
        }
        
        if (selection == FALSE) {
          pop_list <- lapply(pops_vector, function(x) {
            rbind(offspring_list[[x]][sample(which(offspring_list[[x]]$V1 == "Male"),
                                             size =  population_size[x] / 2),],
                  offspring_list[[x]][sample(which(offspring_list[[x]]$V1 == "Female"),
                                             size = population_size[x] / 2),])
          })
        }
        
        if (selection == TRUE &
            natural_selection_model == "absolute") {
          pop_list <- lapply(pops_vector, function(x) {
            rbind(offspring_list[[x]][sample(which(offspring_list[[x]]$V1 == "Male"),
                                             size =  population_size[x] / 2),],
                  offspring_list[[x]][sample(which(offspring_list[[x]]$V1 == "Female"),
                                             size = population_size[x] / 2),])
          })
        }
        
        if (selection == TRUE &
            natural_selection_model == "relative") {
          # We modeled selection as Lesecque et al. 2012: offspring are randomly
          # selected to become parents of the next generation in proportion to
          # their relative fitness, for example, if we had four individuals
          # with fitness (W) of 0.1, 0.2, 0.3, and 0.2 the first individual
          # would be selected on average 0.1/(0.1+0.2+0.3+0.2)=0.125 of the time
          # to become parent of the next generation. The vector of probabilities
          # used in sample is multiplied by two because in the selection function
          # (selection_fun), the proportional relative fitness was calculated for
          # all offspring together, and below the males and females are separated
          # in groups, with the objective that exactly the parents of the next
          # generation are half males and half females
          
          pop_list <- lapply(pops_vector, function(x) {
            males_pop <-
              offspring_list[[x]][which(offspring_list[[x]]$V1 == "Male"),]
            females_pop <-
              offspring_list[[x]][which(offspring_list[[x]]$V1 == "Female"),]
            
            rbind(males_pop[sample(
              row.names(males_pop),
              size = (population_size[x] / 2),
              prob = (males_pop$relative_fitness * 2)
            ), ],
            females_pop[sample(
              row.names(females_pop),
              size = (population_size[x] / 2),
              prob = (females_pop$relative_fitness * 2)
            ), ])
          })
        }
        # toc()
        
        # tic("store")
        ##### STORE VALUES ########
        if (generation %in% gen_store & exists("count_store")) {
          # counter to store genlight objects
          count_store <- count_store + 1
          s_vars_temp <- rbind(ref_vars, sim_vars)
          s_vars_temp <-
            setNames(data.frame(t(s_vars_temp[,-1])), s_vars_temp[, 1])
          s_vars_temp$generation <- generation
          s_vars_temp$iteration <- iteration
          s_vars_temp$seed <- seed
          s_vars_temp$del_ind_cM <- density_mutations_per_cm
          
          s_vars_temp$dispersal_rate_phase1 <-
            paste(dispersal_rate_phase1, collapse = " ")
          s_vars_temp$Fst_expected_phase1 <-
            paste(Fst_expected_phase1, collapse = " ")
          s_vars_temp$mi_expected_phase1 <-
            paste(mi_expected_phase1, collapse = " ")
          s_vars_temp$rate_of_loss_phase1 <- rate_of_loss_phase1
          
          s_vars_temp$dispersal_rate_phase2 <-
            paste(dispersal_rate_phase2, collapse = " ")
          s_vars_temp$Fst_expected_phase2 <-
            paste(Fst_expected_phase2, collapse = " ")
          s_vars_temp$mi_expected_phase2 <-
            paste(mi_expected_phase2, collapse = " ")
          s_vars_temp$rate_of_loss_phase2 <- rate_of_loss_phase2
          
          final_res[[iteration]][[count_store]] <-
            store(
              p_vector = pops_vector,
              p_size = population_size,
              p_list = pop_list,
              n_loc_1 = loci_number,
              paral = parallel,
              n_cores = n.cores,
              ref = reference,
              p_map = plink_map,
              s_vars = s_vars_temp
            )
        }
        # toc()
      }
    }
    
    # naming list elements
    names(final_res) <- paste0("iteration_", 1:number_iterations)
    
    final_res <- lapply(final_res, function(x) {
      names(x) <- paste0("generation_", gen_store)
      return(x)
    })
    
    # removing NA's from results
    final_res <- lapply(final_res, function(x)
      x[!is.na(x)])
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
      cat(report("Completed:", funname, "\n"))
    }
    
    # RETURN
    
    return(invisible(final_res))
    
  }
