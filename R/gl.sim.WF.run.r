#' @name gl.sim.WF.run
#' @title Runs Wright-Fisher simulations
#' @description
#' This function simulates populations made up of diploid organisms that 
#' reproduce in non-overlapping generations. Each individual has a pair of 
#' homologous chromosomes that contains interspersed selected and neutral loci. 
#' For the initial generation, the genotype for each individual’s chromosomes is
#'  randomly drawn from distributions at linkage equilibrium and in 
#'  Hardy-Weinberg equilibrium. 
#' 
#' See documentation and tutorial for a complete description of the simulations.
#' These documents can be accessed at http://georges.biomatix.org/dartR 
#' 
#' Take into account that the simulations will take a little bit longer the
#'  first time you use the function gl.sim.WF.run() because C++ functions must
#'   be compiled.
#' @param file_var Path of the variables file 'sim_variables.csv' (see details) 
#' [required if interactive_vars = FALSE].
#' @param ref_table Reference table created by the function 
#' \code{\link{gl.sim.WF.table}} [required].
#' @param x Name of the genlight object containing the SNP data to extract
#' values for some simulation variables (see details) [default NULL].
#' @param file_dispersal Path of the file with the dispersal table created with
#'  the function \code{\link{gl.sim.create_dispersal}} [default NULL]. 
#' @param number_iterations Number of iterations of the simulations [default 1].
#' @param every_gen Generation interval at which simulations should be stored in
#'  a genlight object [default 10].
#' @param sample_percent Percentage of individuals, from the total population, 
#' to sample and save in the genlight object every generation [default 50].
#' @param store_phase1 Whether to store simulations of phase 1 in genlight
#'  objects [default FALSE].
#' @param interactive_vars Run a shiny app to input interactively the values of
#'  simulations variables [default TRUE].
#' @param seed Set the seed for the simulations [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @param ... Any variable and its value can be added separately within the 
#' function, will be changed over the input value supplied by the csv file. See 
#' tutorial. 
#' @details
#' Values for simulation variables can be submitted into the function 
#' interactively through a shiny app if interactive_vars = TRUE. Optionally, if 
#' interactive_vars = FALSE, values for variables can be submitted by using the
#' csv file 'sim_variables.csv' which can be found by typing in the R console:
#'  system.file('extdata', 'sim_variables.csv', package ='dartR').
#'  
#' The values of the variables can be modified using the third column (“value”) 
#' of this file. 
#' 
#' The output of the simulations can be analysed seemingly with other dartR 
#' functions.
#' 
#' If a genlight object is used as input for some of the simulation variables, 
#' this function access the information stored in the slots x$position and 
#' x$chromosome.
#' 
#' To show further information of the variables in interactive mode, it might be
#'  necessary to call first: 'library(shinyBS)' for the information to be 
#'  displayed.
#' 
#' The main characteristics of the simulations are:
#' \itemize{ 
#' \item Simulations can be parameterised with real-life genetic 
#' characteristics such as the number, location, allele frequency and the 
#' distribution of fitness effects (selection coefficients and dominance) of 
#' loci under selection. 
#' \item Simulations can recreate specific life histories and demographics, such
#'  as source populations, dispersal rate, number of generations, founder 
#'  individuals, effective population size and census population size.
#' \item Each allele in each individual is an agent (i.e., each allele is 
#' explicitly simulated).
#' \item Each locus can be customisable regarding its allele frequencies, 
#' selection coefficients, and dominance.
#' \item The number of loci, individuals, and populations to be simulated is 
#' only limited by computing resources.
#' \item Recombination is accurately modeled, and it is possible to use real 
#' recombination maps as input.
#' \item The ratio between effective population size and census population size 
#' can be easily controlled.
#' \item The output of the simulations are genlight objects for each generation 
#' or a subset of generations.
#' \item Genlight objects can be used as input for some simulation variables.
#' }
#' @return Returns genlight objects with simulated data.
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \dontrun{
#' ref_table <- gl.sim.WF.table(file_var=system.file('extdata', 
#' 'ref_variables.csv', package = 'dartR'),interactive_vars = FALSE)
#' res_sim <- gl.sim.WF.run(file_var = system.file('extdata', 
#' 'sim_variables.csv', package ='dartR'),ref_table=ref_table,
#' interactive_vars = FALSE)
#' }
#' @seealso \code{\link{gl.sim.WF.table}}
#' @family simulation functions
#' @import stats
#' @import shiny
#' @export

gl.sim.WF.run <-
  function(file_var,
           ref_table,
           x = NULL,
           file_dispersal = NULL,
           number_iterations = 1,
           every_gen = 10,
           sample_percent = 50,
           store_phase1 = FALSE,
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
    
    # CHECK IF PACKAGES ARE INSTALLED
    pkg <- "stringi"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      cat(error(
        "Package",
        pkg,
        " needed for this function to work. Please install it.\n"
      ))
      return(-1)
    }
    
    # DO THE JOB
    
    ##### SIMULATIONS VARIABLES ######

    if (interactive_vars==TRUE) {
      
      sim_vars <- interactive_sim_run()
      
     sim_vars[sim_vars$variable=="population_size_phase2" ,"value"] <-
paste0("'",sim_vars[sim_vars$variable=="population_size_phase2" ,"value"],"'")

     sim_vars[sim_vars$variable=="population_size_phase1" ,"value"] <-
paste0("'",sim_vars[sim_vars$variable=="population_size_phase1" ,"value"],"'")

     sim_vars[sim_vars$variable=="dispersal_type_phase2" ,"value"] <- 
paste0("'",sim_vars[sim_vars$variable=="dispersal_type_phase2" ,"value"],"'")
     
     sim_vars[sim_vars$variable=="dispersal_type_phase1" ,"value"] <- 
paste0("'",sim_vars[sim_vars$variable=="dispersal_type_phase1" ,"value"],"'")
     
     sim_vars[sim_vars$variable=="natural_selection_model" ,"value"] <- 
paste0("'",sim_vars[sim_vars$variable=="natural_selection_model" ,"value"],"'")
     
     sim_vars <- sim_vars[order(sim_vars$variable),]
     
     vars_assign <-
       unlist(unname(
         mapply(paste, sim_vars$variable, "<-",
                sim_vars$value, SIMPLIFY = F)
       ))
     
     eval(parse(text = vars_assign))
      
    } else {
      sim_vars <- suppressWarnings(read.csv(file_var))
      sim_vars <- sim_vars[, 2:3]
      
      sim_vars <- sim_vars[order(sim_vars$variable),]
      
      vars_assign <-
        unlist(unname(
          mapply(paste, sim_vars$variable, "<-",
                 sim_vars$value, SIMPLIFY = F)
        ))
      eval(parse(text = vars_assign))
    }
    
    input_list <- list(...)
    
    if(length(input_list)>0){
      
      sim_vars <- sim_vars[order(sim_vars$variable),]
      
      input_list <- input_list[order(names(input_list))]
      
      val_change <- which(sim_vars$variable %in% names(input_list))
      
      sim_vars[val_change,"value"] <- unlist(input_list)
      
      sim_vars[sim_vars$variable=="population_size_phase2" ,"value"] <-
paste0("'",sim_vars[sim_vars$variable=="population_size_phase2" ,"value"],"'")
      
      sim_vars[sim_vars$variable=="population_size_phase1" ,"value"] <-
paste0("'",sim_vars[sim_vars$variable=="population_size_phase1" ,"value"],"'")
      
      sim_vars[sim_vars$variable=="dispersal_type_phase2" ,"value"] <- 
paste0("'",sim_vars[sim_vars$variable=="dispersal_type_phase2" ,"value"],"'")
      
      sim_vars[sim_vars$variable=="dispersal_type_phase1" ,"value"] <- 
paste0("'",sim_vars[sim_vars$variable=="dispersal_type_phase1" ,"value"],"'")
      
      sim_vars[sim_vars$variable=="natural_selection_model" ,"value"] <- 
paste0("'",sim_vars[sim_vars$variable=="natural_selection_model" ,"value"],"'")
      
      vars_assign <-
        unlist(unname(
          mapply(paste, sim_vars$variable, "<-",
                 sim_vars$value, SIMPLIFY = F)
        ))
      eval(parse(text = vars_assign))
    }
    
    reference <- ref_table$reference
    ref_vars <- ref_table$ref_vars
    
    neutral_loci_location <- which(reference$type == "neutral" |
                                     reference$type == "real")
    adv_loci <- 
      which(reference$type=="mutation_adv" | reference$type=="advantageous")
    mutation_loci_adv <- which(reference$type == "mutation_adv" )
    mutation_loci_del <- which(reference$type == "mutation_del")
    mutation_loci_neu <- which(reference$type == "mutation_neu")
    
    mutation_loci_location <- 
      c(mutation_loci_adv,mutation_loci_del,mutation_loci_neu)
    mutation_loci_location <- 
      mutation_loci_location[order(mutation_loci_location)]
    
    real <- which(reference$type == "real")
    
    q_neutral <- as.numeric(ref_vars[ref_vars$variable=="q_neutral","value"])
    
    # this is the option of real_freq from the reference table 
    real_freq_table <- ref_vars[ref_vars$variable=="real_freq","value"]
    if(real_freq_table != real_freq){
      cat(error("  The value for the real_freq parameter was set differently in 
                the simulations and in the creation of the reference table. 
                They should be the same. Please check it\n"))
      stop()
    }
    
    # this is the option of real_loc from the reference table 
    real_loc_table <- ref_vars[ref_vars$variable=="real_loc","value"]
    if(real_loc_table != real_loc){
      cat(error("  The value for the real_loc parameter was set differently in 
                the simulations and in the creation of the reference table. 
                They should be the same. Please check it\n"))
      stop()
    }
    
    if( (real_pops ==TRUE | real_pop_size ==TRUE | real_loc ==TRUE | 
         real_freq==TRUE) && is.null(x)){
      cat(error(" The real dataset to extract information is missing\n"))
      stop()
    }
    
    if (phase1 == FALSE) {
      gen_number_phase1 <- 0
    }
    
    # This is the total number of generations
    number_generations <- gen_number_phase1 + gen_number_phase2
    
    # This is the list to store the final genlight objects
    gen_store <- c(seq(1, number_generations, every_gen), number_generations)
    final_res <- rep(list(as.list(rep(NA, length(gen_store)))), 
                     number_iterations)
    
    loci_number <- nrow(reference)
    recombination_map <- reference[, c("c", "loc_bp", "loc_cM")]
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
    recom_event <- ceiling(sum(recombination_map[, "c"],na.rm = TRUE))
    # filling the probability of recombination when the total recombination rate
    # is less than an integer (recom_event) and placing it at the end of the
    # recombination map

    recombination_map[loci_number + 1, 1] <- recom_event - sum(recombination_map[, 1])
    recombination_map[loci_number + 1, 2] <- recombination_map[loci_number, 2]
    recombination_map[loci_number + 1, 3] <- recombination_map[loci_number, 3]
    
# one is subtracted from the recombination map to account for the last row that
# was added in the recombination map to avoid that the recombination function 
    # crashes
    plink_map <- as.data.frame(matrix(nrow = nrow(reference), ncol = 4))
    plink_map[, 1] <- reference$chr_name
    plink_map[, 2] <- rownames(reference)
    plink_map[, 3] <- reference$loc_cM
    plink_map[, 4] <- reference$loc_bp
    
    dispersal_type_phase2 <- gsub('\"', "", dispersal_type_phase2, fixed = TRUE)
    dispersal_type_phase1 <- gsub('\"', "", dispersal_type_phase1, fixed = TRUE)
    natural_selection_model <- gsub('\"', "", natural_selection_model, 
                                    fixed = TRUE)
    chromosome_name <- gsub('\"', "", chromosome_name, fixed = TRUE)

    population_size_phase2 <- gsub('\"', "", population_size_phase2,
                                   fixed = TRUE)
    
    population_size_phase2 <- 
      as.numeric(unlist(strsplit(population_size_phase2, " ")))
      
    population_size_phase1 <- 
      gsub('\"', "", population_size_phase1, fixed = TRUE)
    
    population_size_phase1 <-
      as.numeric(unlist(strsplit(population_size_phase1, " ")))
    
    local_adap <- gsub('\"', "", local_adap, fixed = TRUE)
    
    local_adap <- as.numeric(unlist(strsplit(local_adap, " ")))
    
    clinal_adap <- gsub('\"', "", clinal_adap, fixed = TRUE)
    
    clinal_adap <- as.numeric(unlist(strsplit(clinal_adap, " ")))

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
      number_pops <- number_pops_phase2 <- nPop(x)
    }
  
    if (real_freq == TRUE & !is.null(x) & real_loc == TRUE) {
      pop_list_freq_temp <- seppop(x)
      loc_to_keep <-
locNames(pop_list_freq_temp[[1]])[which(pop_list_freq_temp[[1]]$chromosome == 
                                          chromosome_name)]
      pop_list_freq_temp <-
        lapply(pop_list_freq_temp,
               gl.keep.loc,
               loc.list = loc_to_keep,
               verbose = 0)
      pop_list_freq <- lapply(pop_list_freq_temp, gl.alf)
      
    } 
    
    if (real_freq == TRUE & !is.null(x) & real_loc == FALSE) {
      pop_list_freq_temp <- seppop(x)
      pop_list_freq <- lapply(pop_list_freq_temp, gl.alf)
      }
    
    if (real_freq == FALSE) {
      pop_list_freq <- rep(NA, number_pops)
    }

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
    freq_deleterious <- reference[-as.numeric(neutral_loci_location),]
    freq_deleterious_b <- mean(2 * (freq_deleterious$q) *
                                 (1 - freq_deleterious$q))
    density_mutations_per_cm <- (freq_deleterious_b * nrow(freq_deleterious)) /
      (recombination_map[loci_number, "loc_cM"] * 100)
   
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
        
      } else {
        
        if (real_pop_size == TRUE & !is.null(x)) {
          population_size_phase2 <- unname(unlist(table(pop(x))))
          # converting odd population sizes to even
          population_size_phase2 <-
            (population_size_phase2 %% 2 != 0) + population_size_phase2
          population_size <- population_size_phase2

        } else{
          
          population_size <- population_size_phase2
          
        }
        
      }
      
      if(phase1 == TRUE & number_pops_phase1!=number_pops_phase2 ){
cat(error("  Number of populations in phase 1 and phase 2 must be the same\n"))
        stop()
      }
      
      if(length(population_size_phase2)!=number_pops_phase2){
        cat(error("  Number of entries for population sizes do not agree with 
                  the number of populations for phase 2\n"))
        stop()
      }
      
      if(length(population_size_phase1)!=number_pops_phase1 & phase1==TRUE){
        cat(error("  Number of entries for population sizes do not agree with 
                  the number of populations for phase 1\n"))
        stop()
      }
      
      ##### INITIALISE POPS #####
     #tic("initialisation")
      if (verbose >= 2) {
        cat(report("  Initialising populations\n"))
      }
      
      # make chromosomes
      #to hack package checking...
      make_chr <- function(){}  
      
      Rcpp::cppFunction(plugins="cpp11",
                        
        'StringVector make_chr(int j, NumericVector q) {
    StringVector out(j);
    int size = 1;
    IntegerVector x = IntegerVector::create(1,0);
    bool rep = false;
for (int i = 0; i < j; i++) {
std::ostringstream temp;
for (int z = 0; z < q.length(); z++) {
NumericVector p = NumericVector::create(q[z],1-q[z]);
   temp << sample(x, size, rep, p);
  }
      out[i] = temp.str();
    }
    return out;
  }'
      )
      
      chr_temp <- make_chr(j=sum(population_size)*2,q=reference$q)
      chr_pops_temps <- split(chr_temp, rep(1:number_pops, 
                                            (c(population_size)*2)))
      chr_pops <- lapply(chr_pops_temps,split,c(1:2))
       
      pops_vector <- 1:number_pops
      
      pop_list <- as.list(pops_vector)
      
      for(pop_n in pops_vector){
        pop <- as.data.frame(matrix(ncol = 4, nrow = population_size[pop_n]))
        pop[, 1] <- rep(c("Male", "Female"), each = population_size[pop_n] / 2)
        pop[, 2] <- pop_n # second column stores population number
        pop[, 3] <- chr_pops[[pop_n]][1]
        pop[, 4] <- chr_pops[[pop_n]][2]
     
        if(real_freq == TRUE & real_loc == TRUE){
        for (individual_pop in 1:population_size[pop_n]) {
          q_prob_t <- pop_list_freq[[pop_n]]$alf1
          q_prob_t2 <- cbind(q_prob_t,1-q_prob_t)
          q_prob <- split(q_prob_t2, row(q_prob_t2))

          stringi::stri_sub_all(pop[individual_pop, 3], from=real,length = 1) <- 
              mapply(function(y){sample(x=c(1,0),size=1,prob=y,replace=FALSE)},
                     q_prob,
                     USE.NAMES = FALSE)
            
          stringi::stri_sub_all(pop[individual_pop, 4], from=real,length = 1) <- 
            mapply(function(y){sample(x=c(1,0),size=1,prob=y,replace=FALSE)},
                   q_prob,
                   USE.NAMES = FALSE)
            
          }
        }
        
        if(real_freq == TRUE & real_loc == FALSE){
          for (individual_pop in 1:population_size[pop_n]) {
            q_prob_t <- pop_list_freq[[pop_n]]$alf1
            q_prob_t2 <- cbind(q_prob_t,1-q_prob_t)
            q_prob <- split(q_prob_t2, row(q_prob_t2))
            
stringi::stri_sub_all(pop[individual_pop, 3], from=real,length = 1) <- 
              mapply(function(y){sample(x=c(1,0),size=1,prob=y,replace=FALSE)},
                     q_prob,
                     USE.NAMES = FALSE)
            
stringi::stri_sub_all(pop[individual_pop, 4], from=real,length = 1) <- 
              mapply(function(y){sample(x=c(1,0),size=1,prob=y,replace=FALSE)},
                     q_prob,
                     USE.NAMES = FALSE)
            
          }
        }
        
        pop_list[[pop_n]] <- pop
      }
      
      #if there is just one population set dispersal to FALSE
      if (length(pop_list) == 1) {
        dispersal <- FALSE
      }
     #toc()
      
      ##### START GENERATION LOOP #####
      for (generation in 1:number_generations) {
        if (phase1 == TRUE & generation == 1) {
          cat(report(" Starting phase 1\n"))
        }
        if (generation %% 10 == 0) {
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
                    replace = TRUE
                  ),],
                  pop_list[[pop_sample]][sample(
                    which(pop_list[[pop_sample]]$V1 == "Female"),
                    size = population_size[x] / 2,
                    replace = TRUE
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
                    replace = TRUE
                  ),],
                  pop_list[[x]][sample(
                    which(pop_list[[x]]$V1 == "Female"),
                    size = population_size[x] / 2,
                    replace = TRUE
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
        
        ##### DISPERSAL ######
       #tic("dispersal")
        if(number_pops==1){
          dispersal <- FALSE
        }
        # dispersal is symmetrical
        if (dispersal == TRUE) {
          
          if(is.null(file_dispersal)){
            
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
            
            dispersal_pairs$number_transfers <- number_transfers
            dispersal_pairs$transfer_each_gen <- transfer_each_gen
            
          }else{
            # if dispersal file is provided
            dispersal_pairs <- suppressWarnings(read.csv(file_dispersal))
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
              trans_gen = dispersal_pairs$transfer_each_gen[dis_pair],
              n_transfer = dispersal_pairs$number_transfers[dis_pair],
              male_tran = maletran,
              female_tran = femaletran
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
       #toc()

        ##### REPRODUCTION #########
       #tic("reproduction")
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
       #toc()

        ##### MUTATION #####
       #tic("mutation")
        if(mutation==TRUE){
          
          for(off_pop in 1:length(offspring_list)){
            
            offspring_pop <- offspring_list[[off_pop]]
            
            offspring_pop$runif <- runif(nrow(offspring_pop))
            
            for (offspring_ind in 1:nrow(offspring_pop)) {
              
              if (length(mutation_loci_location) == 0) {
                cat(important("  No more locus to mutate\n"))
                break()
              }
              
              if (offspring_pop[offspring_ind, "runif"] < mut_rate) {
                locus_to_mutate <- sample(mutation_loci_location, 1)
                mutation_loci_location <-
                  mutation_loci_location[-which(mutation_loci_location == 
                                                  locus_to_mutate)]
                chromosomes <- c(offspring_pop[offspring_ind, 3], 
                                 offspring_pop[offspring_ind, 4])
                chr_to_mutate <- sample(1:2, 1)
                chr_to_mutate_b <- chromosomes[chr_to_mutate]
                substr(chr_to_mutate_b, as.numeric(locus_to_mutate),
                       as.numeric(locus_to_mutate)) <- "1"
                chromosomes[chr_to_mutate] <- chr_to_mutate_b
                offspring_pop[offspring_ind, 3] <- chromosomes[1]
                offspring_pop[offspring_ind, 4] <- chromosomes[2]
                
              } else{
                next()
              }
            }
            offspring_list[[off_pop]] <- offspring_pop
          }
        }
       #toc()

        ##### SELECTION #####
       #tic("selection")
        if (selection == TRUE) {
          if(!is.null(local_adap)){
            pops_local <- setdiff(pops_vector, local_adap)
            reference_local <- replicate(length(pops_vector), 
                                         reference[,c("s","h")], 
                                         simplify = FALSE)
            
            reference_local[pops_local] <- lapply(reference_local[pops_local], 
                                                  function(y){
              y[adv_loci,"s"] <- 0
              return(y)
            })

            offspring_list <- lapply(pops_vector, function(x) {
              selection_fun(
                offspring = offspring_list[[x]],
                h = reference_local[[x]][,"h"],
                s = reference_local[[x]][,"s"],
                sel_model = natural_selection_model,
                g_load = genetic_load
              )
            })
            
          } else if(!is.null(clinal_adap)){
            pops_clinal <- seq(clinal_adap[1],clinal_adap[2],1)
            reference_clinal_temp <- replicate(length(pops_vector), 
                                               reference[,c("s","h")], 
                                               simplify = FALSE)
            
            clinal_s <- 1 - c(0,(1:(length(pops_clinal)-1)*
                                   (clinal_strength/100)))
            
            reference_clinal <- lapply(pops_clinal,function(y){
              reference_clinal_temp[[y]][adv_loci,"s"] <-  
                reference_clinal_temp[[y]][adv_loci,"s"] * clinal_s[y]
              return(reference_clinal_temp[[y]])
            })
            
            offspring_list <- lapply(pops_vector, function(x) {
              selection_fun(
                offspring = offspring_list[[x]],
                h = reference_clinal[[x]][,"h"],
                s = reference_clinal[[x]][,"s"],
                sel_model = natural_selection_model,
                g_load = genetic_load
              )
            })
            
          } else {
          offspring_list <- lapply(pops_vector, function(x) {
            selection_fun(
              offspring = offspring_list[[x]],
              h = reference[,"h"],
              s = reference[,"s"],
              sel_model = natural_selection_model,
              g_load = genetic_load
            )
          })
          }
        }
       #toc()
        
        ##### SAMPLING NEXT GENERATION ########
        #tic("sampling_next_gen")
        # testing whether any population became extinct, if so break the
        # iteration and pass to the next 
        test_extinction <- unlist(lapply(pops_vector, function(x) {
        length(which(offspring_list[[x]]$V1 == "Male")) < population_size / 2 |
        length(which(offspring_list[[x]]$V1 == "Female")) < population_size / 2
        }))

        if (any(test_extinction == TRUE)) {
          
          cat(
            important(
              " One Population became EXTINCT at generation",
              generation,
              "\n"
            )
          )
          cat(important(
            "  Breaking this iteration and passing to the next iteration",
            "\n"
          ))
          
          if (sample_percent != 100) {
            
          population_size_temp <- round(population_size * (sample_percent/100))
            # Converting odd even population sizes to even 
            population_size_temp <- (population_size_temp %% 2 != 0) +
              population_size_temp
            pop_list_temp <- lapply(pops_vector, function(x) {
              rbind(pop_list[[x]][sample(which(pop_list[[x]]$V1 == "Male"),
                                         size =  population_size_temp[x] / 2),],
                    pop_list[[x]][sample(which(pop_list[[x]]$V1 == "Female"),
                                         size = population_size_temp[x] / 2),])
            })
            
          }else{
            
            population_size_temp <- population_size
            pop_list_temp <-  pop_list
            
          }
          
          # formatting the values of the variables to be saved in the genlight
          # object
          s_vars_temp <- rbind(ref_vars, sim_vars)
          s_vars_temp <- setNames(data.frame(t(s_vars_temp[,-1])),
                                  s_vars_temp[, 1])
          s_vars_temp$generation <- generation
          s_vars_temp$iteration <- iteration
          s_vars_temp$seed <- seed
          s_vars_temp$del_ind_cM <- density_mutations_per_cm
          s_vars_temp$sample_percent <- sample_percent
          s_vars_temp$file_dispersal <- file_dispersal
          
          s_vars_temp$dispersal_rate_phase1 <-
            paste(dispersal_rate_phase1, collapse = " ")
          
          s_vars_temp$dispersal_rate_phase2 <-
            paste(dispersal_rate_phase2, collapse = " ")
         
          
          final_res[[iteration]][[count_store]] <-
            store(
              p_vector = pops_vector,
              p_size = population_size_temp,
              p_list = pop_list_temp,
              n_loc_1 = loci_number,
              ref = reference,
              p_map = plink_map,
              s_vars = s_vars_temp
            )
          
          if(real_pops==TRUE){
            
            popNames(final_res[[iteration]][[count_store]]) <- popNames(x)
            
          }
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
       #toc()

        ##### RECYCLE MUTATIONS ###########
        #tic("mutation_2")
        # making available to mutation those loci in which deleterious alleles 
        # have been eliminated from all populations
        if(mutation==T){
          
          pops_merge <- rbindlist(pop_list)
          pops_seqs <- c(pops_merge$V3,pops_merge$V4)
          
          # make frequencies
          #to hack package checking...
          make_freqs <- function(){}  

Rcpp::cppFunction(plugins="cpp11",
                  
"NumericVector make_freqs(StringVector seqs) {
  int seqN = seqs.length();
  int locN = strlen(seqs(0));
  NumericMatrix freq_mat = NumericMatrix(seqN,locN);
  NumericVector out(locN);
  for (int i = 0; i < seqN; i++) {
    for (int j = 0; j < locN; j++) {
      freq_mat(i,j) = seqs(i)[j] - '0';
    }
  }
     for (int y = 0; y < locN; y++){
          out[y] = sum(freq_mat(_,y));
          }
  return out;
}"
                  
)

          freqs <- make_freqs(pops_seqs)
          deleterious_eliminated <- which(freqs==0)
          mutation_loci_location <- 
            union(mutation_loci_location,deleterious_eliminated)
        
        }
       #toc()

        ##### STORE VALUES ########
        #tic("store")
        if (generation %in% gen_store & exists("count_store")) {
          # counter to store genlight objects
          count_store <- count_store + 1
          
          # subsampling individuals 
          
          if (sample_percent < 100) {
            
            population_size_temp <- round(population_size * 
                                            (sample_percent/100))
            # Converting odd population sizes to even 
            population_size_temp <- (population_size_temp %% 2 != 0) + 
              population_size_temp
            pop_list_temp <- lapply(pops_vector, function(x) {
              rbind(pop_list[[x]][sample(which(pop_list[[x]]$V1 == "Male"),
                                      size =  population_size_temp[x] / 2),],
                    pop_list[[x]][sample(which(pop_list[[x]]$V1 == "Female"),
                                      size = population_size_temp[x] / 2),])
            })
            
          }else{
            
          population_size_temp <- population_size
          pop_list_temp <-  pop_list
            
          }
          
          # formatting the values of the variables to be saved in the genlight
          # object
          s_vars_temp <- rbind(ref_vars, sim_vars)
          s_vars_temp <- setNames(data.frame(t(s_vars_temp[,-1])), 
                                  s_vars_temp[, 1])
          s_vars_temp$generation <- generation
          s_vars_temp$iteration <- iteration
          s_vars_temp$seed <- seed
          s_vars_temp$del_ind_cM <- density_mutations_per_cm
          s_vars_temp$sample_percent <- sample_percent
          s_vars_temp$file_dispersal <- file_dispersal
          
          if(dispersal==TRUE){
          s_vars_temp$number_transfers_phase2 <- 
            paste(dispersal_pairs$number_transfers, collapse = " ")  
          s_vars_temp$transfer_each_gen_phase2 <-
            paste(dispersal_pairs$transfer_each_gen, collapse = " ") 
          }
          
          final_res[[iteration]][[count_store]] <-
            store(
              p_vector = pops_vector,
              p_size = population_size_temp,
              p_list = pop_list_temp,
              n_loc_1 = loci_number,
              ref = reference,
              p_map = plink_map,
              s_vars = s_vars_temp,
              g = generation
            )
          
          if(real_pops==TRUE){
            
           popNames(final_res[[iteration]][[count_store]]) <- popNames(x)
            
          }else{
            
  popNames(final_res[[iteration]][[count_store]]) <- as.character(pops_vector)
            
          }
          
        }
       #toc()
      }
    }
    
    # naming list elements
    names(final_res) <- paste0("iteration_", 1:number_iterations)
    
    final_res <- lapply(final_res, function(x) {
      names(x) <- paste0("generation_", gen_store)
      return(x)
    })
    
    # removing NA's from results
    final_res <- lapply(final_res, function(x){
      x[!is.na(x)]
      })
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
      cat(report("Completed:", funname, "\n"))
    }
    
    # RETURN
    
    return(invisible(final_res))
    
  }
