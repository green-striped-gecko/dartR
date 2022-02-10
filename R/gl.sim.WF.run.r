#' @export

gl.sim.WF.run <-
  function(file_var = system.file('extdata', 'sim_variables.csv', package =
                                    'dartR'),
           ref_table,
           number_iterations = 1,
           every_gen = 5,
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
    }
    
    # DO THE JOB
    
    ##### SIMULATIONS VARIABLES ######
    if (interactive_vars) {
      
      ui <- fluidPage(
        
        h5(em("Hover over input boxes to display more information about the variable")),

        h3(strong("General variables")),
        
        fluidRow(
          
          column(
            3,
            numericInput(
              "number_pops",
              strong("Number of populations to simulate"),
              value = 2,
              min = 0
            ),
            bsTooltip(id = "number_pops", 
                      title = "Here is some text with your instructions"),
          )
        ),
        
        h4("Post-adaptation phase variables"),
        
        fluidRow(
          
          column(
            3,
            numericInput(
              "population_size_dispersal",
              "Census Population size of the post-adaptation phase (must be even)",
              value = 100,
              min = 0
            )
          ),
          
          column(
            3,
            numericInput(
              "gen_number_dispersal",
              strong("Number of generations of the post-adaptation phase"),
              value = 100,
              min = 0
            )
          )
        ),
        
        h4("Pre-adaptation phase variables"),
        
        fluidRow(
          
          column(
            3,
            radioButtons(
              "pre_adaptation",
              strong("Whether pre-adaptation phase occur"),
              choices = list("TRUE" = TRUE,
                             "FALSE" = FALSE),
              selected = FALSE
            )
          ),
          
          column(
            3,
            numericInput(
              "population_size_pre_adaptation",
              "Census Population size of the pre-adaptation phase (must be even)",
              value = 100,
              min = 0
            )
          ),
          
          column(
            3,
            numericInput(
              "gen_number_pre_adaptation",
              "Number of generations of the pre-adaptation phase",
              value = 100,
              min = 0
            )
          )),
        
        fluidRow(
          
          column(
          3,
          radioButtons(
            "same_line",
            "Whether the post-adaptation populations are sampled from the same pre-adaptation population or from different pre-adaptation populations",
            choices = list("Same population" = TRUE,
                           "Different populations" = FALSE),
            selected = FALSE
          )
        ),
        
        column(
          3,
          radioButtons(
            "dispersal_pre_adaptation",
            "Whether dispersal occurs in the pre-adaptation phase",
            choices = list("TRUE" = TRUE, "FALSE" = FALSE),
            selected = FALSE
          )
        )
        ),
        
        hr(),
        
        h3(strong("Initialisation variables")),
        
        fluidRow(
          
          column(
          3,
          textInput(
            "chromosome_name",
            "Name of the chromosome to be simulated",
            value = "NA"
          )
        ),
        
        column(
          3,
          numericInput(
            "real_freq",
            "Frequency of each one of the alleles of the loci from a real dataset",
            value = NA,
            min = 0
          )
        )
        ),
        
        hr(),

        h3(strong("Dispersal variables")),
        
        fluidRow(
          
          column(
            3,
            radioButtons(
              "dispersal_dispersal",
              strong("Whether dispersal occurs in the post-adaptation phase"),
              choices = list("TRUE" = TRUE, "FALSE" = FALSE),
              selected = TRUE
            )
          ),
        
        column(
          3,
          numericInput(
            "number_transfers",
            "Number of dispersing individuals in each dispersal event",
            value = 1,
            min = 0
          )
        ),
        
        column(
          3,
          numericInput(
            "transfer_each_gen",
            "Interval of number of generations in which a dispersal event occurs",
            value = 1,
            min = 0
          )
        )
        ),
        
        fluidRow(
          
          column(
            3,
            radioButtons(
              "dispersal_type",
              "Type of dispersal",
              choices = list(
                "All connected" = "all_connected",
                "Circle" = "circle",
                "Line" = "line"
              ),
              selected = "all_connected"
            )
          )
        ),
        
        hr(),
        
        h3(strong("Reproduction variables")),
        
        fluidRow(
          
          column(
            3,
            numericInput(
              "number_offspring",
              "Mean number offspring per mating",
              value = 10,
              min = 0
            )
          ),
          
          column(
            3,
            numericInput(
              "variance_offspring",
              "Coefficient that determines the variance in the number of offspring per mating",
              value = 1000000,
              min = 0
            )
          )
        ),
        
        hr(),
        
        h3(strong("Recombination variables")),
        
        fluidRow(
          
          column(
            3,
            radioButtons(
              "recombination",
              "Whether recombination occurs",
              choices = list("TRUE" = TRUE,
                             "FALSE" = FALSE),
              selected = TRUE
            )
          ),
          
          column(
            3,
            radioButtons(
              "recombination_males",
              "Whether recombination occurs in males and females (or only females)",
              choices = list("TRUE" = TRUE,
                             "FALSE" = FALSE),
              selected = TRUE
            )
          )
        ),
        
        hr(),
        
        h3(strong("Selection variables")),
        
        fluidRow(
          
          column(
            3,
            radioButtons(
              "selection",
              "Whether selection occurs",
              choices = list("TRUE" = TRUE,
                             "FALSE" = FALSE),
              selected = FALSE
            )
          ),
          
          column(
            3,
            radioButtons(
              "natural_selection_model",
              "Selection model to use",
              choices = list("Relative" = "relative",
                             "Absolute" = "absolute"),
              selected = "relative"
            )
          ),
          
          column(
            3,
            numericInput(
              "genetic_load",
              "Approximation of the genetic load of the proportion of the genome that is simulated. This variable is used in the absolute fitness model",
              value = 0.8,
              min = 0
            )
          )
          ),
        
        hr(),
        
        h3(strong("Analysis variables")),
        
        fluidRow(
          
          column(
            3,
            numericInput(
              "Ne",
              "Ne value to be used in the equation of the expected rate of loss of heterozygosity",
              value = 50,
              min = 0
            )
          ),
          
          column(
            3,
            numericInput(
              "Ne_fst",
              "Ne value to be used in the equation of the expected FST",
              value = 50,
              min = 0
            )
          )
          ),
        
        hr(),
        
        fluidRow(
          column(
          9,
          actionButton(
            "close",
            label = h4(strong("RUN")),
            icon = icon("play"),
            width = "100%",
            class = "btn-success"
          )
        )
        ),
        
        br()
        
      )
      
      server <- function(input, output) {
        observeEvent(input$Ne, {
          Ne <<- input$Ne
        })
        observeEvent(input$Ne_fst, {
          Ne_fst <<- input$Ne_fst
        })
        observeEvent(input$dispersal_dispersal, {
          dispersal_dispersal <<- input$dispersal_dispersal
        })
        observeEvent(input$dispersal_pre_adaptation, {
          dispersal_pre_adaptation <<- input$dispersal_pre_adaptation
        })
        observeEvent(input$dispersal_type, {
          dispersal_type <<- input$dispersal_type
        })
        observeEvent(input$number_transfers, {
          number_transfers <<- input$number_transfers
        })
        observeEvent(input$transfer_each_gen, {
          transfer_each_gen <<- input$transfer_each_gen
        })
        observeEvent(input$gen_number_dispersal, {
          gen_number_dispersal <<- input$gen_number_dispersal
        })
        observeEvent(input$gen_number_pre_adaptation, {
          gen_number_pre_adaptation <<- input$gen_number_pre_adaptation
        })
        observeEvent(input$number_pops, {
          number_pops <<- input$number_pops
        })
        observeEvent(input$population_size_dispersal, {
          population_size_dispersal <<- input$population_size_dispersal
        })
        observeEvent(input$population_size_pre_adaptation, {
          population_size_pre_adaptation <<- input$population_size_pre_adaptation
        })
        observeEvent(input$pre_adaptation, {
          pre_adaptation <<- input$pre_adaptation
        })
        observeEvent(input$same_line, {
          same_line <<- input$same_line
        })
        observeEvent(input$chromosome_name, {
          chromosome_name <<- input$chromosome_name
        })
        observeEvent(input$real_freq, {
          real_freq <<- input$real_freq
        })
        observeEvent(input$recombination, {
          recombination <<- input$recombination
        })
        observeEvent(input$recombination_males, {
          recombination_males <<- input$recombination_males
        })
        observeEvent(input$number_offspring, {
          number_offspring <<- input$number_offspring
        })
        observeEvent(input$variance_offspring, {
          variance_offspring <<- input$variance_offspring
        })
        observeEvent(input$genetic_load, {
          genetic_load <<- input$genetic_load
        })
        
        observeEvent(input$natural_selection_model, {
          natural_selection_model <<- input$natural_selection_model
        })
        observeEvent(input$selection, {
          selection <<- input$selection
        })
        observeEvent(input$close, {
          sim_vars_temp <<- as.data.frame(cbind(
            c("Ne",
              "Ne_fst",
              "dispersal_dispersal",
              "dispersal_pre_adaptation",
              "dispersal_type",
              "number_transfers",
              "transfer_each_gen",
              "gen_number_dispersal",
              "gen_number_pre_adaptation",
              "number_pops",
              "population_size_dispersal",
              "population_size_pre_adaptation",
              "pre_adaptation",
              "same_line",
              "chromosome_name",
              "real_freq",
              "recombination",
              "recombination_males",
              "number_offspring",
              "variance_offspring",
              "genetic_load",
              "natural_selection_model",
              "selection"),
            c(Ne,
              Ne_fst,
              dispersal_dispersal,
              dispersal_pre_adaptation,
              dispersal_type,
              number_transfers,
              transfer_each_gen,
              gen_number_dispersal,
              gen_number_pre_adaptation,
              number_pops,
              population_size_dispersal,
              population_size_pre_adaptation,
              pre_adaptation,
              same_line,
              chromosome_name,
              real_freq,
              recombination,
              recombination_males,
              number_offspring,
              variance_offspring,
              genetic_load,
              natural_selection_model,
              selection)
          ))
          
          colnames(sim_vars_temp) <- c("variable","value")
          
          sim_vars <<- sim_vars_temp
          
          stopApp()
          
        })
        
      }
      
      runApp(shinyApp(ui = ui, server = server))
      
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
    
    # setting the seed
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    if (pre_adaptation == F) {
      gen_number_pre_adaptation <- 0
    }
    
    # This is the total number of generations
    number_generations <-
      gen_number_pre_adaptation + gen_number_dispersal
    
    #this is the list to store the final genlight objects
    gen_store <-
      c(seq(1, number_generations, every_gen), number_generations)
    final_res <-
      rep(list(as.list(rep(
        NA, length(gen_store)
      ))), number_iterations)
    
    # pick which sex is going to be transferred first
    if (number_transfers >= 2) {
      maletran <- TRUE
      femaletran <- TRUE
    } else if (number_transfers == 1) {
      maletran <- TRUE
      femaletran <- FALSE
    }
    
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
    freq_deleterious <- reference[-as.numeric(neutral_loci_location), ]
    freq_deleterious_b <-
      mean(2 * (freq_deleterious$q) * (1 - freq_deleterious$q))
    density_mutations_per_cm <-
      (freq_deleterious_b * nrow(freq_deleterious)) /
      (recombination_map[loci_number, "loc_cM"] * 100)
    
    dispersal_rate <-
      (number_transfers / transfer_each_gen) / (population_size_dispersal)
    Fst_expected <-
      1 / ((4 * Ne_fst * dispersal_rate) * ((2 / (2 - 1)) ^ 2) + 1)
    mi_expected <-
      (0.22 / (sqrt(2 * Ne_fst * dispersal_rate))) - (0.69 / ((2 * Ne_fst) * sqrt(dispersal_rate)))
    rate_of_loss <- 1 - (1 / (2 * Ne))
    
    ##### START ITERATION LOOP #####
    for (iteration in 1:number_iterations) {
      if (iteration %% 1 == 0) {
        cat(report("  Iteration =", iteration, "\n"))
      }
      ##### VARIABLES PRE_ADAPTATION PHASE #######
      if (pre_adaptation == TRUE) {
        population_size <- population_size_pre_adaptation
        dispersal <- dispersal_pre_adaptation
        store_values <- FALSE
      } else{
        population_size <- population_size_dispersal
        dispersal <- dispersal_dispersal
        store_values <- TRUE
      }
      
      ##### INITIALISE POPS #####
      pops_vector <- 1:number_pops
      pop_list <- lapply(pops_vector, function(x) {
        initialise(
          pop_number = x,
          pop_size = population_size,
          refer = reference,
          n_l_loc = neutral_loci_location,
          r_freq = real_freq
        )
      })
      
      #if there is just one population set dispersal to FALSE
      if (length(pop_list) == 1) {
        dispersal <- FALSE
        dispersal_dispersal <- FALSE
      }
      
      ##### START GENERATION LOOP #####
      for (generation in 1:number_generations) {
        if (generation %% 5 == 0) {
          cat(report("   Generation =", generation, "\n"))
        }
        ##### VARIABLES DISPERSAL PHASE #######
        if (generation == (gen_number_pre_adaptation + 1)) {
          cat(report("  Starting post-adaptation phase\n"))
          population_size <- population_size_dispersal
          dispersal <- dispersal_dispersal
          store_values <- TRUE
          # counter to store genlight objects
          count_store <- 0
          # counter to store values every generation
          gen_dispersal <- 0
          
          if (pre_adaptation == TRUE) {
            if (same_line == TRUE) {
              # pop_list_temp is used because pop_list is used to sample populations
              pop_sample <- sample(pops_vector, 1)
              
              pop_list_temp <- lapply(pops_vector, function(x) {
                pop_temp <-
                  rbind(pop_list[[pop_sample]][sample(which(pop_list[[pop_sample]]$V1 == "Male"),
                                                      size =  population_size / 2), ],
                        pop_list[[pop_sample]][sample(which(pop_list[[pop_sample]]$V1 == "Female"),
                                                      size = population_size / 2), ])
                pop_temp$V2 <- x
                return(pop_temp)
              })
              pop_list <- pop_list_temp
            }
            
            if (same_line == FALSE) {
              pop_list <- lapply(pops_vector, function(x) {
                pop_temp <-
                  rbind(pop_list[[x]][sample(which(pop_list[[x]]$V1 == "Male"),
                                             size =  population_size / 2), ],
                        pop_list[[x]][sample(which(pop_list[[x]]$V1 == "Female"),
                                             size = population_size / 2), ])
                pop_temp$V2 <- x
                return(pop_temp)
              })
              
            }
          }
        }
        # generation counter
        if (store_values == TRUE) {
          gen_dispersal <- gen_dispersal + 1
        }
        ##### DISPERSAL ######
        # dispersal is symmetrical
        if (dispersal == TRUE) {
          if (dispersal_type == "all_connected") {
            dispersal_pairs <-
              as.data.frame(expand.grid(pops_vector, pops_vector))
            dispersal_pairs$same_pop <-
              dispersal_pairs$Var1 == dispersal_pairs$Var2
            dispersal_pairs <-
              dispersal_pairs[which(dispersal_pairs$same_pop == FALSE), ]
            colnames(dispersal_pairs) <- c("pop1", "pop2", "same_pop")
            
            for (dis_pair in 1:nrow(dispersal_pairs)) {
              res <- migration(
                population1 = pop_list[[dispersal_pairs[dis_pair, "pop1"]]],
                population2 = pop_list[[dispersal_pairs[dis_pair, "pop2"]]],
                generation = generation,
                pop_size = population_size,
                trans_gen = transfer_each_gen,
                male_tran = maletran,
                female_tran = femaletran,
                n_transfer = number_transfers
              )
              
              pop_list[[dispersal_pairs[dis_pair, "pop1"]]] <- res[[1]]
              pop_list[[dispersal_pairs[dis_pair, "pop1"]]]$V2 <-
                dispersal_pairs[dis_pair, "pop1"]
              pop_list[[dispersal_pairs[dis_pair, "pop2"]]] <- res[[2]]
              pop_list[[dispersal_pairs[dis_pair, "pop2"]]]$V2 <-
                dispersal_pairs[dis_pair, "pop2"]
              maletran <- res[[3]]
              femaletran <- res[[4]]
            }
          }
          
          if (dispersal_type == "line") {
            dispersal_pairs <-
              as.data.frame(rbind(
                cbind(head(pops_vector, -1), pops_vector[-1]),
                cbind(pops_vector[-1], head(pops_vector, -1))
              ))
            colnames(dispersal_pairs) <- c("pop1", "pop2")
            
            for (dis_pair in 1:nrow(dispersal_pairs)) {
              res <- migration(
                population1 = pop_list[[dispersal_pairs[dis_pair, "pop1"]]],
                population2 = pop_list[[dispersal_pairs[dis_pair, "pop2"]]],
                generation = generation,
                pop_size = population_size,
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
          
          if (dispersal_type == "circle") {
            dispersal_pairs <-
              as.data.frame(rbind(cbind(
                pops_vector, c(pops_vector[-1], pops_vector[1])
              ),
              cbind(
                c(pops_vector[-1], pops_vector[1]), pops_vector
              )))
            colnames(dispersal_pairs) <- c("pop1", "pop2")
            
            for (dis_pair in 1:nrow(dispersal_pairs)) {
              res <- migration(
                population1 = pop_list[[dispersal_pairs[dis_pair, "pop1"]]],
                population2 = pop_list[[dispersal_pairs[dis_pair, "pop2"]]],
                generation = generation,
                pop_size = population_size,
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
        }
        
        ##### REPRODUCTION #########
        offspring_list <- lapply(pops_vector, function(x) {
          reproduction(
            pop = pop_list[[x]],
            pop_number = x,
            pop_size = population_size,
            var_off = variance_offspring,
            num_off = number_offspring,
            r_event = recom_event,
            recom = recombination,
            r_males = recombination_males,
            r_map_1 = recombination_map,
            n_loc = loci_number
          )
        })
        
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
              "Population",
              pops_extinct,
              "became EXTINCT at generation",
              generation,
              "\n"
            )
          )
          cat(important(
            "Breaking this iteration and passing to the next iteration",
            "\n"
          ))
          
          s_vars_temp <- sim_vars[, 2:3]
          s_vars_temp <-
            setNames(data.frame(t(s_vars_temp[, -1])), s_vars_temp[, 1])
          s_vars_temp$generation <- generation
          s_vars_temp$iteration <- iteration
          s_vars_temp$seed <- seed
          s_vars_temp$del_ind_cM <- density_mutations_per_cm
          s_vars_temp$disp_rate <- dispersal_rate
          s_vars_temp$fst_exp <- Fst_expected
          s_vars_temp$mi_exp <- mi_expected
          s_vars_temp$rate_loss <- rate_of_loss
          
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
                                             size =  population_size / 2), ],
                  offspring_list[[x]][sample(which(offspring_list[[x]]$V1 == "Female"),
                                             size = population_size / 2), ])
          })
        }
        
        if (selection == TRUE &
            natural_selection_model == "absolute") {
          pop_list <- lapply(pops_vector, function(x) {
            rbind(offspring_list[[x]][sample(which(offspring_list[[x]]$V1 == "Male"),
                                             size =  population_size / 2), ],
                  offspring_list[[x]][sample(which(offspring_list[[x]]$V1 == "Female"),
                                             size = population_size / 2), ])
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
              offspring_list[[x]][which(offspring_list[[x]]$V1 == "Male"), ]
            females_pop <-
              offspring_list[[x]][which(offspring_list[[x]]$V1 == "Female"), ]
            
            rbind(males_pop[sample(
              row.names(males_pop),
              size = (population_size / 2),
              prob = (males_pop$relative_fitness * 2)
            ),],
            females_pop[sample(
              row.names(females_pop),
              size = (population_size / 2),
              prob = (females_pop$relative_fitness * 2)
            ),])
          })
        }
        
        ##### STORE VALUES ########
        if (generation %in% gen_store) {
          # counter to store genlight objects
          count_store <- count_store + 1
          s_vars_temp <- rbind(ref_vars,sim_vars)
          s_vars_temp <-
            setNames(data.frame(t(s_vars_temp[, -1])), s_vars_temp[, 1])
          s_vars_temp$generation <- generation
          s_vars_temp$iteration <- iteration
          s_vars_temp$seed <- seed
          s_vars_temp$del_ind_cM <- density_mutations_per_cm
          s_vars_temp$disp_rate <- dispersal_rate
          s_vars_temp$fst_exp <- Fst_expected
          s_vars_temp$mi_exp <- mi_expected
          s_vars_temp$rate_loss <- rate_of_loss
          
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
      }
    }
    
    # naming list elements
    names(final_res) <- paste0("iteration_",1:number_iterations)
    
    final_res <- lapply(final_res,function(x){
      names(x) <- paste0("generation_",gen_store)
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
