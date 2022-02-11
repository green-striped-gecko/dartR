ped <- function(df_ped,
                n_loc = loci_number) {
  chromosome1 <- df_ped[3]
  chromosome2 <- df_ped[4]
  split_seqs <- strsplit(c(chromosome1, chromosome2), split = "")
  genotypes <- as.data.frame(matrix(nrow = n_loc, ncol = 2))
  genotypes$V1 <- split_seqs[[1]]
  genotypes$V2 <- split_seqs[[2]]
  genotypes_final <-
    paste0(paste(genotypes[, 1], genotypes[, 2]), collapse = " ")
  return(genotypes_final)
}

delta <- function(a, b, c) {
  # # Constructing delta
  b ^ 2 - 4 * a * c
}

q_equilibrium <- function(a, b, c) {
  # Constructing Quadratic Formula
  x_1 = (-b + sqrt(delta(a, b, c))) / (2 * a)
  return(x_1)
}

migration <- function(population1,
                      population2,
                      generation,
                      pop_size = population_size,
                      trans_gen = transfer_each_gen,
                      male_tran = maletran,
                      female_tran = femaletran,
                      n_transfer = number_transfers) {
  if (generation != 1 & generation %% trans_gen == 0) {
    if (n_transfer == 1) {
      if (male_tran) {
        malepoptran <- sample(c(1:(pop_size / 2)), size = 1)
        temppop1 <- population1[malepoptran, ]
        temppop2 <- population2[malepoptran, ]
        population2[malepoptran, ] <- temppop1
        population1[malepoptran, ] <- temppop2
      }
      if (female_tran) {
        fempoptran <-
          sample(c(((pop_size / 2) + 1):pop_size), size = 1)
        temppop1 <- population1[fempoptran, ]
        temppop2 <- population2[fempoptran, ]
        population2[fempoptran, ] <- temppop1
        population1[fempoptran, ] <- temppop2
      }
    }
    if (n_transfer >= 2) {
      size_malepoptran <- ceiling(n_transfer / 2)
      size_femalepoptran <- floor(n_transfer / 2)
      if (male_tran) {
        malepoptran <- sample(c(1:(pop_size / 2)), size = size_malepoptran)
        temppop1 <- population1[malepoptran, ]
        temppop2 <- population2[malepoptran, ]
        population2[malepoptran, ] <- temppop1
        population1[malepoptran, ] <- temppop2
      }
      if (female_tran) {
        fempoptran <-
          sample(c(((pop_size / 2) + 1):pop_size), size = size_femalepoptran)
        temppop1 <- population1[fempoptran, ]
        temppop2 <- population2[fempoptran, ]
        population2[fempoptran, ] <- temppop1
        population1[fempoptran, ] <- temppop2
      }
    }
    # if necessary flip transfer
    if (n_transfer == 1) {
      male_tran <- !male_tran
      female_tran <- !female_tran
    }
  }
  return (list(population1, population2, male_tran, female_tran))
}

selection_fun <- function(offspring,
                          reference_pop,
                          sel_model = natural_selection_model,
                          g_load = genetic_load) {
  offspring$fitness <-
    apply(offspring, 1, fitness, ref = reference_pop)
  if (sel_model == "absolute") {
    offspring$random_deviate <-
      runif(nrow(offspring), min = 0, max = g_load)
    offspring$alive <- offspring$fitness > offspring$random_deviate
    offspring <- offspring[which(offspring$alive == TRUE),]
  }
  if (sel_model == "relative") {
    fitnes_proportion <- sum(offspring$fitness)
    offspring$relative_fitness <-
      offspring$fitness / fitnes_proportion
    offspring$relative_fitness[offspring$relative_fitness < 0] <-
      0
  }
  return(offspring)
}

# this is the function to calculate fitness
fitness <- function(df_fitness,
                    ref) {
  #remove neutral loci from chromosomes
  chromosome1 <- gsub("[-^1-9]", "0", df_fitness[3])
  chromosome2 <- gsub("[-^1-9]", "0", df_fitness[4])
  split_seqs <- strsplit(c(chromosome1, chromosome2), split = "")
  ref$hom <-
    (split_seqs[[1]] == split_seqs[[2]] & split_seqs[[1]] == "a")
  ref$het_chr1 <- (split_seqs[[1]] == "a" & split_seqs[[2]] == "A")
  ref$het_chr2 <- (split_seqs[[1]] == "A" & split_seqs[[2]] == "a")
  ref$hom <- as.numeric(ref$hom)
  sum_hom <- sum(as.numeric(ref$hom))
  ref$het_chr1 <- as.numeric(ref$het_chr1)
  ref$het_chr2 <- as.numeric(ref$het_chr2)
  ref$het_sel_chr1 <- ref$h * ref$s * ref$het_chr1
  ref$het_sel_chr2 <- ref$h * ref$s * ref$het_chr2
  ref$hom_sel <- ref$s * ref$hom
  ref$tot_sel <- ref$het_sel_chr1 + ref$het_sel_chr2 + ref$hom_sel
  ref$fitness <- 1 - (ref$tot_sel)
  # keeping only loci with NS
  deleterious <- ref[ref$fitness < 1, ]
  net_fitness <- prod(deleterious$fitness)
  
  return(net_fitness)
}

# this is the function for recombination
recomb <- function(r_chromosome1,
                   r_chromosome2, 
                   r_map,
                   loci) {
  chiasma <-
    as.numeric(sample(row.names(r_map), size = 1, prob = r_map[, "c"]))
  if (chiasma < (loci + 1)) {
    split_seqs <- strsplit(c(r_chromosome1, r_chromosome2), split = "")
    r_chr1 <-
      paste0(c(split_seqs[[1]][1:chiasma], split_seqs[[2]][(chiasma + 1):loci]), collapse = "")
    r_chr2 <-
      paste0(c(split_seqs[[2]][1:chiasma], split_seqs[[1]][(chiasma + 1):loci]), collapse = "")
    return(list(r_chr1, r_chr2))
  } else{
    return(list(r_chromosome1, r_chromosome2))
  }
}

initialise <- function(pop_number, 
                       pop_size = population_size, 
                       refer = reference,
                       n_l_loc = neutral_loci_location, 
                       r_freq = real_freq) {
  
    pop <- as.data.frame(matrix(ncol = 4, nrow = pop_size))
    pop[, 1] <- rep(c("Male", "Female"), each = pop_size / 2)
    pop[, 2] <- pop_number # second column stores population
    for (individual_pop in 1:pop_size) {
      chromosome1 <-
        paste0(mapply(sample_alleles, q = refer$q,  USE.NAMES = F), collapse = "")
      chromosome2 <-
        paste0(mapply(sample_alleles, q = refer$q,  USE.NAMES = F), collapse = "")
      for (element in n_l_loc) {
        substr(chromosome1, start = element, stop = element) <-
          sample(as.character(c(1:2)),
                 size = 1,
                 prob = c(rep(1 / 2, 2)))
        substr(chromosome2, start = element, stop = element) <-
          sample(as.character(c(1:2)),
                 size = 1,
                 prob = c(rep(1 / 2, 2)))
      }
      if (!is.na(r_freq)) {
        counter_msats <- 1
        for (element in 1:length(loc_exp_loci_2)) {
          substr(
            chromosome1,
            start = as.numeric(loc_exp_loci_2[element]),
            stop = as.numeric(loc_exp_loci_2[element])
          ) <-
            sample(as.character(c(1:length(
              msats_freq[[counter_msats]]
            ))), size = 1, prob = msats_freq[[counter_msats]])
          substr(
            chromosome2,
            start = as.numeric(loc_exp_loci_2[element]),
            stop = as.numeric(loc_exp_loci_2[element])
          ) <-
            sample(as.character(c(1:length(
              msats_freq[[counter_msats]]
            ))), size = 1, prob = msats_freq[[counter_msats]])
          counter_msats <- counter_msats + 1
        }
      }
      pop[individual_pop, 3] <- chromosome1
      pop[individual_pop, 4] <- chromosome2
    }
    return(pop)
  }

sample_alleles <- function(alleles = c("a", "A"), q) {
  sample(alleles, size = 1, prob = c(q, 1 - q))
}

reproduction <- function(pop,
                           pop_number,
                           pop_size = population_size,
                           var_off = variance_offspring,
                           num_off = number_offspring,
                           r_event = recom_event,
                           recom = recombination,
                           r_males = recombination_males,
                           r_map_1 = recombination_map,
                           n_loc = loci_number) {
  parents_matrix <-
    as.data.frame(matrix(nrow = pop_size / 2, ncol = 2))
  male_pool <-
    sample(rownames(pop[1:(pop_size / 2),]), size = pop_size / 2)
  parents_matrix[, 1] <- sample(male_pool, size = pop_size / 2)
  parents_matrix[, 2] <-
    sample(rownames(pop[((pop_size / 2) + 1):pop_size,]), size = pop_size / 2)
  offspring <- NULL
  for (parent in 1:dim(parents_matrix)[1]) {
    pairing_offspring <-  rnbinom(1, size = var_off, mu = num_off)
    offspring_temp <-
      as.data.frame(matrix(nrow = pairing_offspring, ncol = 6))
    if (pairing_offspring < 1) {
      next
    }
    offspring_temp[, 1] <-
      sample(c("Male", "Female"), size = pairing_offspring, replace = TRUE) # sex
    offspring_temp[, 2] <- pop_number # source population
    male_chromosomes <-
      list(pop[parents_matrix[parent, 1], 3], pop[parents_matrix[parent, 1], 4])
    female_chromosomes <-
      list(pop[parents_matrix[parent, 2], 3], pop[parents_matrix[parent, 2], 4])
    for (offs in 1:pairing_offspring) {
      males_recom_events <- rpois(1, r_event)
      females_recom_events <- rpois(1, r_event)
      #recombination in males
      if (recom == TRUE &
          r_males == TRUE & males_recom_events > 1) {
        for (event in males_recom_events) {
          male_chromosomes <-
            recomb(male_chromosomes[[1]],
                   male_chromosomes[[2]],
                   r_map = r_map_1,
                   loci = n_loc)
        }
        offspring_temp[offs, 3] <-
          male_chromosomes[[sample(c(1, 2), 1)]]
      } else{
        offspring_temp[offs, 3] <- male_chromosomes[[sample(c(1, 2), 1)]]
      }
      #recombination in females
      if (recom == TRUE & females_recom_events > 1) {
        for (event in females_recom_events) {
          female_chromosomes <-
            recomb(
              female_chromosomes[[1]],
              female_chromosomes[[2]],
              r_map = r_map_1,
              loci = n_loc
            )
        }
        offspring_temp[offs, 4] <-
          female_chromosomes[[sample(c(1, 2), 1)]]
      } else{
        offspring_temp[offs, 4] <- female_chromosomes[[sample(c(1, 2), 1)]]
      }
    }
    offspring_temp[, 5] <- parents_matrix[parent, 1] #id father
    offspring_temp[, 6] <- parents_matrix[parent, 2] #id mother
    offspring <- rbind(offspring, offspring_temp)
  }
  return(offspring)
}

store <- function(p_vector = pops_vector,
                  p_size = population_size,
                  p_list = pop_list,
                  n_loc_1 = loci_number,
                  paral = parallel,
                  n_cores = n.cores,
                  ref = reference,
                  p_map = plink_map,
                  s_vars = s_vars_temp){
  pop_names <- rep(paste0("pop",p_vector),p_size)
  pop_names <- pop_names[order(pop_names)]
  df_genotypes <- rbindlist(p_list)
  df_genotypes$V1[df_genotypes$V1 == "Male"]   <- 1
  df_genotypes$V1[df_genotypes$V1 == "Female"] <- 2
  df_genotypes[, 2] <- pop_names
  df_genotypes$id <- paste0(unlist(unname(df_genotypes[, 2])), "_", rep(1:p_size,length(p_vector)))
  plink_ped <- apply(df_genotypes, 1, ped, n_loc = n_loc_1)
  # converting allele names to numbers
  plink_ped <- gsub("a", "1", plink_ped) 
  plink_ped <- gsub("A", "2", plink_ped)
  plink_ped <-
    lapply(plink_ped, function(x) {
      gsub(" ", "", strsplit(x, '(?<=([^ ]\\s){2})', perl = TRUE)[[1]])
    })
  plink_ped_2 <- lapply(plink_ped, function(x) {
    x[x == "22"] <- 2
    x[x == "11"] <- 0
    x[x == "21"] <- 1
    x[x == "12"] <- 1
    return(x)
  })
  
  if (paral && is.null(n_cores)) {
    n_cores <- parallel::detectCores()
  }
  
  loc.names <- 1:nrow(ref)
  n.loc <- length(loc.names)
  misc.info <- lapply(1:6, function(i)
    NULL)
  names(misc.info) <-
    c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
  res <- list()
  temp <-
    as.data.frame(cbind(
      df_genotypes[,2],
      df_genotypes[,7],
      df_genotypes[,5],
      df_genotypes[,6],
      df_genotypes[,1],
      1
    ))
  
  for (i in 1:6) {
    misc.info[[i]] <- temp[, i]
  }
  txt <-
    lapply(plink_ped_2, function(e)
      suppressWarnings(as.integer(e)))
  if (paral) {
    res <-
      c(
        res,
        parallel::mclapply(txt, function(e)
          new("SNPbin", snp = e, ploidy = 2L), mc.cores = n_cores, mc.silent = TRUE, mc.cleanup = TRUE, mc.preschedule = FALSE)
      )
  } else {
    res <-
      c(res, lapply(txt, function(e)
        new(
          "SNPbin", snp = e, ploidy = 2L
        )))
  }
  
  res <- new("genlight", res, ploidy = 2L,parallel = paral)
  
  indNames(res) <- misc.info$IID
  pop(res) <- misc.info$FID
  locNames(res) <- loc.names
  misc.info <- misc.info[c("SEX", "PHENOTYPE", "PAT", "MAT")]
  names(misc.info) <- tolower(names(misc.info))
  misc.info$sex[misc.info$sex == 1] <- "m"
  misc.info$sex[misc.info$sex == 2] <- "f"
  misc.info$sex <- factor(misc.info$sex)
  misc.info$phenotype[misc.info$phenotype == 1] <- "control"
  misc.info$phenotype[misc.info$phenotype == 2] <- "case"
  misc.info$phenotype <- factor(misc.info$phenotype)
  res$other$ind.metrics <- as.data.frame(misc.info)
  loc_metrics_temp <- as.data.frame(cbind(p_map, ref[, 2:5]))
  colnames(loc_metrics_temp) <-
    c("chr", "loc_id", "loc_cM", "loc_bp", "q", "h", "s","selection")
  res$other$loc.metrics <- loc_metrics_temp
  res$other$sim.vars <- s_vars
  return(res)
  
}

interactive_reference <- function(){

ui <- fluidPage(
  
  h5(
    em(
      "Hover over input boxes to display more information about the variable"
    )
  ),
  
  h3(strong("General variables")),
  
  fluidRow(
    column(
      3,
      numericInput(
        "chunk_number",
        "Number of rows of the recombination map",
        value = 100,
        min = 0
      )
    ),
    
    column(
      3,
      numericInput(
        "loci_number_to_simulate",
        "Number of loci under selection",
        value = 1000,
        min = 0
      )
    ),
    
    column(
      3,
      numericInput(
        "real_loc",
        "Locations of the loci from a real dataset",
        value = NULL,
        min = 0
      )
    )
  ),
  
  fluidRow(column(
    3,
    textInput(
      "chromosome_name",
      "Name of the chromosome to be simulated",
      value = "NA"
    )
  )),
  
  hr(),
  
  h3(strong("Alelle frequency variables")),
  
  fluidRow(
    column(
      3,
      sliderInput(
        "q_neutral",
        "Initial frequencies of neutral alleles",
        value = 0.5,
        min = 0,
        max = 1
      )
    ),
    
    column(
      3,
      radioButtons(
        "q_distribution",
        strong(
          "How the initial allele frequency of the deleterious allele (q) should be determined"
        ),
        choices = list("All equal" = "equal", "From equation" = "equation"),
        selected = "equal"
      )
    ),
    
    column(
      3,
      numericInput(
        "mutation_rate",
        "Mutation rate per generation per site. Value only used in the equation to determine q",
        value = 5 * 10 ^ -5
      )
    )
  ),
  
  fluidRow(column(
    3,
    sliderInput(
      "q_gral",
      "Initial frequencies of deleterious alleles",
      value = 0.15,
      min = 0,
      max = 1
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
  )),
  
  hr(),
  
  h3(strong("Recombination variables")),
  
  fluidRow(
    column(
      3,
      numericInput(
        "map_resolution",
        "Resolution of the recombination map (bp)",
        value = 100000,
        min = 0
      )
    ),
    
    column(
      3,
      numericInput(
        "chunk_recombination",
        "Recombination rate (cM) per region of size chunk_number",
        value = 1,
        min = 0
      )
    )
  ),
  
  hr(),
  
  h3(strong("Selection variables")),
  
  h4("Selection coefficient variables"),
  
  fluidRow(
    
    column(
      3,
      radioButtons(
        "s_distribution",
        strong(
          "Name of the distribution to use to sample the values of the selection coefficient (s) for each locus under selection"
        ),
        choices = list("All equal" = "equal", 
                       "Gamma distribution" = "gamma",
                       "Log normal distribution" = "log_normal"),
        selected = "equal"
      )
    )
  ),
  
  h5(strong("All equal variables")),
  
  fluidRow(
    
    column(
      3,
      numericInput(
        "s_gral",
        "Selection coefficient of deleterious alleles",
        value = 0.001,
        min = 0
      )
    )
  ),
  
  h5(strong("Gamma distribution variables")),
  
  fluidRow(
    
    column(
      3,
      numericInput(
        "gamma_scale",
        "Scale of the gamma distribution",
        value = 0.03,
        min = 0
      )
    ),
    
    column(
      3,
      numericInput(
        "gamma_shape",
        "Shape of the gamma distribution",
        value = 0.25,
        min = 0
      )
    )
    
  ),
  
  h5(strong("Log normal distribution variables")),
  
  fluidRow(
    
    column(
      3,
      numericInput(
        "log_mean",
        "Mean of the log normal distribution",
        value = 0.002,
        min = 0
      )
    ),
    
    column(
      3,
      numericInput(
        "log_sd",
        "Standard deviation of the log normal distribution",
        value = 4,
        min = 0
      )
    )
  ),
  
  hr(),
  
  h4("Dominance coefficient variables"),
  
  fluidRow(
    
    column(
      3,
      radioButtons(
        "h_distribution",
        strong(
          "Name of the distribution to use to sample the values of the dominance coefficient (h) for each locus under selection"
        ),
        choices = list("All equal" = "equal", 
                       "Normal distribution" = "normal",
                       "From equation" = "equation"),
        selected = "equal"
      )
    )
  ),
  
  h5(strong("All equal variables")),
  
  fluidRow(
    
    column(
      3,
      sliderInput(
        "h_gral",
        "Dominance coefficient of deleterious alleles",
        value = 0.25,
        min = 0,
        max = 1
      )
    )
  ),
  
  h5(strong("Normal distribution variables")),
  
  fluidRow(
    
    column(
      3,
      sliderInput(
        "dominance_mean",
        "Mean of the normal distribution from where h values are sampled",
        value = 0.25,
        min = 0,
        max = 1
      )
    )
  ),
  
  h5(strong("Equation variables")),
  
  fluidRow(
    
    column(
      3,
      sliderInput(
        "intercept",
        "Value for the intercept of the equation",
        value = 0.5,
        min = 0,
        max = 1
      )
    ),
    
    column(
      3,
      numericInput(
        "rate",
        "Value for the variable rate of the equation",
        value = 500,
        min = 0
      )
    )
  ),
  
  hr(),
  
  h4("Targets of selection variables"),
  
  fluidRow(
    
    column(
      3,
      numericInput(
        "targets_factor",
        "Factor to sample the number of loci under selection from the input file 'targets_of_selection.csv'",
        value = 0.05,
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
  
  observeEvent(input$loci_number_to_simulate, {
    loci_number_to_simulate <<- input$loci_number_to_simulate
  })
  observeEvent(input$mutation_rate, {
    mutation_rate <<- input$mutation_rate
  })
  observeEvent(input$map_resolution, {
    map_resolution <<- input$map_resolution
  })
  observeEvent(input$dominance_mean, {
    dominance_mean <<- input$dominance_mean
  })
  observeEvent(input$gamma_scale, {
    gamma_scale <<- input$gamma_scale
  })
  observeEvent(input$gamma_shape, {
    gamma_shape <<- input$gamma_shape
  })
  observeEvent(input$h_gral, {
    h_gral <<- input$h_gral
  })
  observeEvent(input$intercept, {
    intercept <<- input$intercept
  })
  observeEvent(input$log_mean, {
    log_mean <<- input$log_mean
  })
  observeEvent(input$log_sd, {
    log_sd <<- input$log_sd
  })
  observeEvent(input$q_gral, {
    q_gral <<- input$q_gral
  })
  observeEvent(input$rate, {
    rate <<- input$rate
  })
  observeEvent(input$s_gral, {
    s_gral <<- input$s_gral
  })
  observeEvent(input$targets_factor, {
    targets_factor <<- input$targets_factor
  })
  observeEvent(input$chromosome_name, {
    chromosome_name <<- input$chromosome_name
  })
  observeEvent(input$chunk_number, {
    chunk_number <<- input$chunk_number
  })
  observeEvent(input$chunk_recombination, {
    chunk_recombination <<- input$chunk_recombination
  })
  observeEvent(input$real_loc, {
    real_loc <<- input$real_loc
  })
  observeEvent(input$real_freq, {
    real_freq <<- input$real_freq
  })
  observeEvent(input$q_neutral, {
    q_neutral <<- input$q_neutral
  })
  observeEvent(input$q_distribution, {
    q_distribution <<- input$q_distribution
  })
  observeEvent(input$h_distribution, {
    h_distribution <<- input$h_distribution
  })
  observeEvent(input$s_distribution, {
    s_distribution <<- input$s_distribution
  })
  
  observeEvent(input$close, {
    
    ref_vars_temp <<- as.data.frame(cbind(
      c(
        "real_freq",
        "real_loc",
        "loci_number_to_simulate",
        "mutation_rate",
        "map_resolution",
        "dominance_mean",
        "gamma_scale",
        "gamma_shape",
        "h_gral",
        "intercept",
        "log_mean",
        "log_sd",
        "q_gral",
        "rate",
        "s_gral",
        "targets_factor",
        "chromosome_name",
        "chunk_number",
        "chunk_recombination",
        "q_neutral",
        "q_distribution",
        "h_distribution",
        "s_distribution"
      ),
      c(
        real_freq,
        real_loc,
        loci_number_to_simulate,
        mutation_rate,
        map_resolution,
        dominance_mean,
        gamma_scale,
        gamma_shape,
        h_gral,
        intercept,
        log_mean,
        log_sd,
        q_gral,
        rate,
        s_gral,
        targets_factor,
        chromosome_name,
        chunk_number,
        chunk_recombination,
        q_neutral,
        q_distribution,
        h_distribution,
        s_distribution
      )
    ))
    
    colnames(ref_vars_temp) <- c("variable", "value")
    
    ref_vars <<- ref_vars_temp
    
    stopApp()
    
  })
  
}

runApp(shinyApp(ui = ui, server = server))
}

interactive_sim_run <- function(){
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
}