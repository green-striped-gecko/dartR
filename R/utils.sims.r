###############################################################################
######################## CONVERT TO PED FORMAT ################################
###############################################################################

ped <- function(df_ped, n_loc) {
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

###############################################################################
######################## CALCULATION OF Q #####################################
###############################################################################

delta <- function(a, b, c) {
  # # Constructing delta
  d <- b ^ 2 - 4 * a * c
  
  return(d)
  
}

q_equilibrium <- function(a, b, c) {
  # Constructing Quadratic Formula
  x_1 <- (-b + sqrt(delta(a, b, c))) / (2 * a)
  
  return(x_1)
  
}

###############################################################################
################################ MIGRATION ####################################
###############################################################################

migration <-
  function(population1,
           population2,
           gen,
           size_pop1,
           size_pop2,
           trans_gen,
           male_tran,
           female_tran,
           n_transfer) {
    if (gen != 1 & gen %% trans_gen == 0) {
      if (n_transfer == 1) {
        if (male_tran) {
          malepoptran_pop1 <- sample(c(1:(size_pop1 / 2)), size = 1)
          malepoptran_pop2 <- sample(c(1:(size_pop2 / 2)), size = 1)
          temppop1 <- population1[malepoptran_pop1, ]
          temppop2 <- population2[malepoptran_pop2, ]
          population2[malepoptran_pop2, ] <- temppop1
          population1[malepoptran_pop1, ] <- temppop2
        }
        if (female_tran) {
          fempoptran_pop1 <-
            sample(c(((
              size_pop1 / 2
            ) + 1):size_pop1), size = 1)
          fempoptran_pop2 <-
            sample(c(((
              size_pop2 / 2
            ) + 1):size_pop2), size = 1)
          temppop1 <- population1[fempoptran_pop1, ]
          temppop2 <- population2[fempoptran_pop2, ]
          population2[fempoptran_pop2, ] <- temppop1
          population1[fempoptran_pop1, ] <- temppop2
        }
      }
      if (n_transfer >= 2) {
        size_malepoptran <- ceiling(n_transfer / 2)
        size_femalepoptran <- floor(n_transfer / 2)
        if (male_tran) {
          malepoptran_pop1 <-
            sample(c(1:(size_pop1 / 2)), size = size_malepoptran)
          malepoptran_pop2 <-
            sample(c(1:(size_pop2 / 2)), size = size_malepoptran)
          temppop1 <- population1[malepoptran_pop1, ]
          temppop2 <- population2[malepoptran_pop2, ]
          population2[malepoptran_pop2, ] <- temppop1
          population1[malepoptran_pop1, ] <- temppop2
        }
        if (female_tran) {
          fempoptran_pop1 <-
            sample(c(((
              size_pop1 / 2
            ) + 1):size_pop1), size = size_femalepoptran)
          fempoptran_pop2 <-
            sample(c(((
              size_pop2 / 2
            ) + 1):size_pop2), size = size_femalepoptran)
          temppop1 <- population1[fempoptran_pop1, ]
          temppop2 <- population2[fempoptran_pop2, ]
          population2[fempoptran_pop2, ] <- temppop1
          population1[fempoptran_pop1, ] <- temppop2
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

###############################################################################
########################## SELECTION ##########################################
###############################################################################

selection_fun <-
  function(offspring,
           reference_pop,
           sel_model,
           g_load) {
    offspring$fitness <-
      apply(offspring, 1, fitness, ref = reference_pop)
    if (sel_model == "absolute") {
      offspring$random_deviate <-
        runif(nrow(offspring), min = 0, max = g_load)
      offspring$alive <-
        offspring$fitness > offspring$random_deviate
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

###############################################################################
########################## FITNESS ############################################
###############################################################################

fitness <- function(df_fitness, ref) {
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

###############################################################################
########################## RECOMBINATION ######################################
###############################################################################

recomb <- function(r_chromosome1, r_chromosome2, r_map, loci) {
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

###############################################################################
########################## POP INITIALISATION #################################
###############################################################################

initialise <-
  function(pop_number,
           pop_size,
           refer,
           n_l_loc,
           r_freq,
           q_neu) {
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
                 prob = c(q_neu, 1 - q_neu))
        substr(chromosome2, start = element, stop = element) <-
          sample(as.character(c(1:2)),
                 size = 1,
                 prob = c(q_neu, 1 - q_neu))
      }
      if (!anyNA(r_freq)) {
        counter_snp <- 1
        for (element in 1:nrow(r_freq)) {
          substr(chromosome1,
                 start = as.numeric(n_l_loc[element]),
                 stop = as.numeric(n_l_loc[element])) <-
            sample(as.character(c(1:2)),
                   size = 1,
                   prob = unname(unlist(r_freq[counter_snp,])))
          
          substr(chromosome2,
                 start = as.numeric(n_l_loc[element]),
                 stop = as.numeric(n_l_loc[element])) <-
            sample(as.character(c(1:2)),
                   size = 1,
                   prob = unname(unlist(r_freq[counter_snp,])))
          
          counter_snp <- counter_snp + 1
          
        }
      }
      pop[individual_pop, 3] <- chromosome1
      pop[individual_pop, 4] <- chromosome2
    }
    
    return(pop)
    
  }

###############################################################################
########################## SAMPLE ALLELES #####################################
###############################################################################

sample_alleles <- function(alleles = c("a", "A"), q) {
  s_alleles <- sample(alleles, size = 1, prob = c(q, 1 - q))
  
  return(s_alleles)
  
}

###############################################################################
########################## REPRODUCTION #######################################
###############################################################################

reproduction <-
  function(pop,
           pop_number,
           pop_size,
           var_off,
           num_off,
           r_event,
           recom,
           r_males,
           r_map_1,
           n_loc) {
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

###############################################################################
########################## STORE GENLIGHT #####################################
###############################################################################

store <-
  function(p_vector,
           p_size,
           p_list,
           n_loc_1,
           paral,
           n_cores,
           ref,
           p_map,
           s_vars) {
    pop_names <- rep(paste0("pop", p_vector), p_size)
    pop_names <- pop_names[order(pop_names)]
    df_genotypes <- rbindlist(p_list)
    df_genotypes$V1[df_genotypes$V1 == "Male"]   <- 1
    df_genotypes$V1[df_genotypes$V1 == "Female"] <- 2
    df_genotypes[, 2] <- pop_names
    df_genotypes$id <-
      paste0(unlist(unname(df_genotypes[, 2])), "_", unlist(lapply(p_size, function(x) {
        1:x
      })))
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
      as.data.frame(
        cbind(
          df_genotypes[, 2],
          df_genotypes[, 7],
          df_genotypes[, 5],
          df_genotypes[, 6],
          df_genotypes[, 1],
          1
        )
      )
    
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
            new("SNPbin", snp = e, ploidy = 2L), mc.cores = n_cores,
            mc.silent = TRUE, mc.cleanup = TRUE, mc.preschedule = FALSE)
        )
    } else {
      res <-
        c(res, lapply(txt, function(e)
          new(
            "SNPbin", snp = e, ploidy = 2L
          )))
    }
    
    res <- new("genlight", res, ploidy = 2L, parallel = paral)
    
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
      c("chr", "loc_id", "loc_cM", "loc_bp", "q", "h", "s", "selection")
    res$other$loc.metrics <- loc_metrics_temp
    res$other$sim.vars <- s_vars
    
    return(res)
    
  }

###############################################################################
######### SHINY APP FOR THE VARIABLES OF THE REFRENCE TABLE ###################
###############################################################################
#' @name interactive_reference
#' @title Shiny app for the input of the reference table for the simulations
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @param env_fun Environment of the calling function
#' @import shiny
#' @import shinyBS
#' @import shinythemes
#' @import shinyjs

interactive_reference <- function(env_fun) {
  
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
  
  ui <- fluidPage(
    
    useShinyjs(),
    
    theme = shinytheme("darkly"),
    
    h5(
      em(
        "Enlarge window for a better visualisation of the variables"
      )
    ),
    
    h5(
      em(
        "Click on the window and hover over an input box to display more information about the variable"
      )
    ),
    
    h5(
      em(
        "Title of input box is the variable's name as used in documentation, tutorials and code of simulations"
      )
    ),
    
    h3(strong("General variables")),
    
    fluidRow(
      
      column(
        4,
        numericInput(
          "chunk_number",  
          tags$div(tags$i(HTML("chunk_number<br/>")),
            "Number of rows of recombination map"
          ),
          value = 10,
          min = 0
        ), 
        bsTooltip(id = "chunk_number",
                        title = "Information pending")
      ),
      
      column(
        4,
        numericInput(
          "loci_number_to_simulate", 
          tags$div(tags$i(HTML("loci_number_to_simulate<br/>")),
            "Number of loci under selection"
          ),
          value = 100,
          min = 0
        ),     
        bsTooltip(id = "loci_number_to_simulate",
                         title = "Information pending")
      ),
      
      column(
        4,
        radioButtons(
          "real_loc",
          tags$div(tags$i(HTML("real_loc<br/>")),
            "Should the number and location of neutral loci be based on the genlight object?"
          ) ,
          choices = list("TRUE" = TRUE,
                         "FALSE" = FALSE),
          selected = FALSE
        ),    
        bsTooltip(id = "real_loc",
                         title = "Information pending")
      )
    ),
    
    fluidRow(
      
      column(
      4,
      numericInput(
        "neutral_loci_chunk",  
        tags$div(tags$i(HTML("neutral_loci_chunk<br/>")),
        "Number of neutral loci per chromosome chunk"
      ),
        value = 1,
        min = 0
      ),    
      bsTooltip(id = "neutral_loci_chunk",
                title = "Information pending")
    )
    
    ),
    
    fluidRow(
      
      column(
      4,
      textInput(
        "chromosome_name", 
        tags$div(tags$i(HTML("chromosome_name<br/>")),
        "Name of the chromosome to be simulated"),
        value = "2L"
      ),    
      bsTooltip(id = "chromosome_name",
                title = "Information pending")
    )
    
    ),
    
    hr(),
    
    h3(strong("Alelle frequency variables")),
    
    fluidRow(
      
      column(
        4,
        sliderInput(
          "q_neutral", 
          tags$div(tags$i(HTML("q_neutral<br/>")),
          "Initial frequencies of neutral alleles"),
          value = 0.5,
          min = 0,
          max = 1
        ),    
        bsTooltip(id = "q_neutral",
                  title = "Information pending")
      ),
      
      column(
        4,
        radioButtons(
          "q_distribution",   
          tags$div(tags$i(HTML("q_distribution<br/>")),
            "How the initial allele frequency of the deleterious allele (q) should be determined"),
          choices = list("All equal" = "equal", "From equation" = "equation"),
          selected = "equal"
        ),    
        bsTooltip(id = "q_distribution",
                  title = "Information pending")
      ),
      
      column(
        4,
        numericInput(
          "mutation_rate", 
          tags$div(tags$i(HTML("mutation_rate<br/>")),
          "Mutation rate per generation per site. Value only used in the equation to determine q"),
          value = 5 * 10 ^ -5
        ),    
        bsTooltip(id = "mutation_rate",
                  title = "Information pending")
      )
      
    ),
    
    fluidRow(
      
      column(
      4,
      sliderInput(
        "q_gral", 
        tags$div(tags$i(HTML("q_gral<br/>")),
        "Initial frequencies of deleterious alleles"),
        value = 0.15,
        min = 0,
        max = 1
      ),    
      bsTooltip(id = "q_gral",
                title = "Information pending")
    ),
    
    column(
      4,
      radioButtons(
        "real_freq", 
        tags$div(tags$i(HTML("real_freq<br/>")),
          "Should the allele frequency of neutral loci be based on the genlight object?"),
        choices = list("TRUE" = TRUE,
                       "FALSE" = FALSE),
        selected = FALSE
      ),    
      bsTooltip(id = "real_freq",
                title = "Information pending")
    )
    
    ),
    
    hr(),
    
    h3(strong("Recombination variables")),
    
    fluidRow(
      
      column(
      4,
      numericInput(
        "map_resolution",
        tags$div(tags$i(HTML("map_resolution<br/>")),
        "Resolution of the recombination map (bp)"),
        value = 100000,
        min = 0
      ),    
      bsTooltip(id = "map_resolution",
                title = "Information pending")
    ),
    
    column(
      4,
      numericInput(
        "chunk_recombination",
        tags$div(tags$i(HTML("chunk_recombination<br/>")),
        "Recombination rate (cM) per region of size chunk_number"),
        value = 1,
        min = 0
      ),    
      bsTooltip(id = "chunk_recombination",
                title = "Information pending")
    )
    
    ),
    
    hr(),
    
    h3(strong("Selection variables")),
    
    h4("Selection coefficient variables"),
    
    fluidRow(
      
      column(
      4,
      radioButtons(
        "s_distribution",
        tags$div(tags$i(HTML("s_distribution<br/>")),
          "Distribution to sample selection coefficients"),
        choices = list(
          "All equal" = "equal",
          "Gamma distribution" = "gamma",
          "Log normal distribution" = "log_normal"
        ),
        selected = "equal"
      ),    
      bsTooltip(id = "s_distribution",
                title = "Information pending")
    )
    
    ),
    
    h5(strong("All equal variables")),
    
    fluidRow(
      
      column(
      4,
      numericInput(
        "s_gral",
        tags$div(tags$i(HTML("s_gral<br/>")),
        "Selection coefficient"),
        value = 0.001,
        min = 0
      ),    
      bsTooltip(id = "s_gral",
                title = "Information pending")
    )
    
    ),
    
    h5(strong("Gamma distribution variables")),
    
    fluidRow(
      
      column(
      4,
      numericInput(
        "gamma_scale",
        tags$div(tags$i(HTML("gamma_scale<br/>")),
        "Scale of gamma distribution"),
        value = 0.03,
        min = 0
      ),    
      bsTooltip(id = "gamma_scale",
                title = "Information pending")
    ),
    
    column(
      4,
      numericInput(
        "gamma_shape",
        tags$div(tags$i(HTML("gamma_shape<br/>")),
        "Shape of gamma distribution"),
        value = 0.25,
        min = 0
      ),    
      bsTooltip(id = "gamma_shape",
                title = "Information pending")
    )
    
    ),
    
    h5(strong("Log normal distribution variables")),
    
    fluidRow(
      
      column(
      4,
      numericInput(
        "log_mean",
        tags$div(tags$i(HTML("log_mean<br/>")),
        "Mean of log normal distribution"),
        value = 0.002,
        min = 0
      ),    
      bsTooltip(id = "log_mean",
                title = "Information pending")
    ),
    
    column(
      4,
      numericInput(
        "log_sd",
        tags$div(tags$i(HTML("log_sd<br/>")),
        "Standard deviation of log normal distribution"),
        value = 4,
        min = 0
      ),    
      bsTooltip(id = "log_sd",
                title = "Information pending")
    )
    
    ),
    
    hr(),
    
    h4("Dominance coefficient variables"),
    
    fluidRow(
      
      column(
      4,
      radioButtons(
        "h_distribution",
        tags$div(tags$i(HTML("h_distribution<br/>")),
          "Distribution to sample dominance coefficient"),
        choices = list(
          "All equal" = "equal",
          "Normal distribution" = "normal",
          "From equation" = "equation"
        ),
        selected = "equal"
      ),    
      bsTooltip(id = "h_distribution",
                title = "Information pending")
    )
    
    ),
    
    h5(strong("All equal variables")),
    
    fluidRow(
      
      column(
      4,
      sliderInput(
        "h_gral",
        tags$div(tags$i(HTML("h_gral<br/>")),
        "Dominance coefficient"),
        value = 0.25,
        min = 0,
        max = 1
      ),    
      bsTooltip(id = "h_gral",
                title = "Information pending")
    )
    
    ),
    
    h5(strong("Normal distribution variables")),
    
    fluidRow(
      
      column(
      4,
      sliderInput(
        "dominance_mean",
        tags$div(tags$i(HTML("dominance_mean<br/>")),
        "Mean of normal distribution"),
        value = 0.25,
        min = 0,
        max = 1
      ),    
      bsTooltip(id = "dominance_mean",
                title = "Information pending")
    ),
    
    column(
      4,
      numericInput(
        "dominance_sd",
        tags$div(tags$i(HTML("dominance_sd<br/>")),
        "Standard deviation of normal distribution"),
        value = sqrt(0.001),
        min = 0
      ),    
      bsTooltip(id = "dominance_sd",
                title = "Information pending")
    )
    
    ),
    
    h5(strong("Equation variables")),
    
    fluidRow(
      
      column(
      4,
      sliderInput(
        "intercept",
        tags$div(tags$i(HTML("intercept<br/>")),
        "Value for the intercept of the equation"),
        value = 0.5,
        min = 0,
        max = 1
      ),    
      bsTooltip(id = "intercept",
                title = "Information pending")
    ),
    
    column(
      4,
      numericInput(
        "rate",
        tags$div(tags$i(HTML("rate<br/>")),
        "Value for the variable rate of the equation"),
        value = 500,
        min = 0
      ),    
      bsTooltip(id = "rate",
                title = "Information pending")
    )
    
    ),
    
    hr(),
    
    h4("Targets of selection variables"),
    
    fluidRow(
      
      column(
      4,
      numericInput(
        "targets_factor",
        tags$div(tags$i(HTML("targets_factor<br/>")),
        "Factor to sample the number of loci under selection from the input file 'targets_of_selection.csv'"),
        value = 0.01,
        min = 0
      ),    
      bsTooltip(id = "targets_factor",
                title = "Information pending")
    )
    
    ),
    
    hr(),
    
    fluidRow(column(
      9,
      actionButton(
        "close",
        label = h4(strong("RUN")),
        icon = icon("play"),
        width = "100%",
        class = "btn-success"
      )
    )),
    
    br()
    
  )
  
  server <- function(input, output, session) {
    
    observeEvent(input$loci_number_to_simulate, {
      assign("loci_number_to_simulate",input$loci_number_to_simulate,envir = env_fun)
    })
    observeEvent(input$mutation_rate, {
      assign("mutation_rate",input$mutation_rate,envir = env_fun)
    })
    observeEvent(input$map_resolution, {
      assign("map_resolution",input$map_resolution,envir = env_fun)
    })
    observeEvent(input$dominance_mean, {
      assign("dominance_mean",input$dominance_mean,envir = env_fun)
    })
    observeEvent(input$gamma_scale, {
      assign("gamma_scale",input$gamma_scale,envir = env_fun)
    })
    observeEvent(input$gamma_shape, {
      assign("gamma_shape",input$gamma_shape,envir = env_fun)
    })
    observeEvent(input$h_gral, {
      assign("h_gral",input$h_gral,envir = env_fun)
    })
    observeEvent(input$intercept, {
      assign("intercept",input$intercept,envir = env_fun)
    })
    observeEvent(input$log_mean, {
      assign("log_mean",input$log_mean,envir = env_fun)
    })
    observeEvent(input$log_sd, {
      assign("log_sd",input$log_sd,envir = env_fun)
    })
    observeEvent(input$q_gral, {
      assign("q_gral",input$q_gral,envir = env_fun)
    })
    observeEvent(input$rate, {
      assign("rate",input$rate,envir = env_fun)
    })
    observeEvent(input$s_gral, {
      assign("s_gral",input$s_gral,envir = env_fun)
    })
    observeEvent(input$targets_factor, {
      assign("targets_factor",input$targets_factor,envir = env_fun)
    })
    observeEvent(input$chromosome_name, {
      assign("chromosome_name",input$chromosome_name,envir = env_fun)
    })
    observeEvent(input$chunk_number, {
      assign("chunk_number",input$chunk_number,envir = env_fun)
    })
    observeEvent(input$chunk_recombination, {
      assign("chunk_recombination",input$chunk_recombination,envir = env_fun)
    })
    observeEvent(input$real_loc, {
      assign("real_loc",input$real_loc,envir = env_fun)
    })
    observeEvent(input$real_freq, {
      assign("real_freq",input$real_freq,envir = env_fun)
    })
    observeEvent(input$q_neutral, {
      assign("q_neutral",input$q_neutral,envir = env_fun)
    })
    observeEvent(input$q_distribution, {
      assign("q_distribution",input$q_distribution,envir = env_fun)
    })
    observeEvent(input$h_distribution, {
      assign("h_distribution",input$h_distribution,envir = env_fun)
    })
    observeEvent(input$s_distribution, {
      assign("s_distribution",input$s_distribution,envir = env_fun)
    })
    observeEvent(input$neutral_loci_chunk, {
      assign("neutral_loci_chunk",input$neutral_loci_chunk,envir = env_fun)
    })
    observeEvent(input$dominance_sd, {
      assign("dominance_sd",input$dominance_sd,envir = env_fun)
    })
    
    observeEvent(input$close, {
      assign("ref_vars_temp", as.data.frame(cbind(
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
          "s_distribution",
          "neutral_loci_chunk",
          "dominance_sd"
        ),
        c(
          input$real_freq,
          input$real_loc,
          input$loci_number_to_simulate,
          input$mutation_rate,
          input$map_resolution,
          input$dominance_mean,
          input$gamma_scale,
          input$gamma_shape,
          input$h_gral,
          input$intercept,
          input$log_mean,
          input$log_sd,
          input$q_gral,
          input$rate,
          input$s_gral,
          input$targets_factor,
          input$chromosome_name,
          input$chunk_number,
          input$chunk_recombination,
          input$q_neutral,
          input$q_distribution,
          input$h_distribution,
          input$s_distribution,
          input$neutral_loci_chunk,
          input$dominance_sd
        )
      )))
      
      colnames(ref_vars_temp) <- c("variable", "value")
      
      assign("ref_vars",ref_vars_temp,envir = env_fun)
      
      stopApp()
      
    })
    
  }
  
  runApp(shinyApp(ui, server))
  # runGadget(shinyApp(ui, server), viewer =  paneViewer())
  
}

###############################################################################
######### SHINY APP FOR THE VARIABLES OF THE SIMULATIONS ######################
###############################################################################
#' @name interactive_sim_run
#' @title Shiny app for the input of the simulations variables
#' @param env_fun Environment of the calling function
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @import shiny
#' @import shinyBS
#' @import shinythemes
#' @import shinyjs

interactive_sim_run <- function(env_fun) {
  
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
  
  ui <- fluidPage(
    
    h5(
      em(
        "Enlarge window for a better visualisation of the variables"
      )
    ),
    
    h5(
      em(
        "Click on the window and hover over an input box to display more information about the variable"
      )
    ),
    
    h5(
      em(
        "Title of input box is the variable's name as used in documentation, tutorials and code of simulations"
      )
    ),
    
    
    h3(strong("Real dataset variables")),
    
    fluidRow(
      
      column(
        4,
        radioButtons(
          "real_pops",
          tags$div(tags$i(HTML("real_pops<br/>")),
            "Should the number of populations in the simulations be based on the genlight object?"),
          choices = list("TRUE" = TRUE,
                         "FALSE" = FALSE),
          selected = FALSE
        ),    
        bsTooltip(id = "real_pops",
                  title = "Information pending")
      ),
      
      column(
        4,
        radioButtons(
          "real_pop_size",
          tags$div(tags$i(HTML("real_pop_size<br/>")),
            "Should the census population size in the simulations be based on the genlight object?"),
          choices = list("TRUE" = TRUE,
                         "FALSE" = FALSE),
          selected = FALSE
        ),    
        bsTooltip(id = "real_pop_size",
                  title = "Information pending")
      ),
      
      column(
        4,
        radioButtons(
          "real_freq",
          tags$div(tags$i(HTML("real_freq<br/>")),
            "Should the initial frequency of neutral alleles in the simulations be based on the genlight object?"),
          choices = list("TRUE" = TRUE,
                         "FALSE" = FALSE),
          selected = FALSE
        ),    
        bsTooltip(id = "real_freq",
                  title = "Information pending")
      )
      
    ),
    
    hr(),
    
    h3(strong("Phase 2 variables")),
    
    fluidRow(
      
      column(
        4,
        numericInput(
          "number_pops_phase2",
          tags$div(tags$i(HTML(
            "number_pops_phase2<br/>"
          )),
          "Number of populations in phase 2"),
          value = 2,
          min = 0
        ),
        bsTooltip(id = "number_pops_phase2",
                  title = "Information pending")
      ),
      
      column(
        4,
        textInput(
          "population_size_phase2",
          tags$div(
            tags$i(HTML("population_size_phase2<br/>")),
            "Census population size of phase 2 (must be even and space delimited)"),
          value = c("10", "10")
        ),
        bsTooltip(id = "population_size_phase2",
                  title = "Information pending")
      ),
      
      column(
        4,
        numericInput(
          "gen_number_phase2",
          tags$div(tags$i(HTML(
            "gen_number_phase2<br/>"
          )),
          "Number of generations of phase 2"),
          value = 10,
          min = 0
        ),
        bsTooltip(id = "gen_number_phase2",
                  title = "Information pending")
      )
    ),
    
    fluidRow(
      
      column(
        4,
        radioButtons(
          "dispersal_phase2",
          tags$div(tags$i(HTML(
            "dispersal_phase2<br/>"
          )),
          "Whether dispersal occurs in phase 2"),
          choices = list("TRUE" = TRUE, "FALSE" = FALSE),
          selected = TRUE
        ),
        bsTooltip(id = "dispersal_phase2",
                  title = "Dispersal between populations is symmetric and constant across generations. Dispersal rate (m) is the fraction of individuals in a population that is composed of dispersers or the probability that a randomly chosen individual in this generation came from a population different from the one in which it was found in the preceding generation (Holsinger, 2020, p. 93). Dispersal rate is calculated as (number_transfers / transfer_each_gen) / pop_size.")
      ),
      
      column(
        4,
        radioButtons(
          "dispersal_type_phase2",
          tags$div(tags$i(HTML(
            "dispersal_type_phase2<br/>"
          )),
          "Type of dispersal for phase 2"),
          choices = list(
            "All connected" = "all_connected",
            "Circle" = "circle",
            "Line" = "line"
          ),
          selected = "all_connected"
        ),
        bsTooltip(id = "dispersal_type_phase2",
                  title = "Information pending")
      ),
      
      column(
        4,
        numericInput(
          "number_transfers_phase2",
          tags$div(
            tags$i(HTML("number_transfers_phase2<br/>")),
            "Number of dispersing individuals in each dispersal event in phase 2"),
          value = 1,
          min = 0
        ),
        bsTooltip(id = "number_transfers_phase2",
                  title = "Information pending")
      )
    ),
    
    fluidRow(
      column(
        4,
        numericInput(
          "transfer_each_gen_phase2",
          tags$div(
            tags$i(HTML("transfer_each_gen_phase2<br/>")),
            "Interval of number of generations in which a dispersal event occurs in phase 2"),
          value = 1,
          min = 0
        ),
        bsTooltip(id = "transfer_each_gen_phase2",
                  title = "Information pending")
      ),
      
      column(
        4,
        numericInput(
          "variance_offspring_phase2",
          tags$div(
            tags$i(HTML("variance_offspring_phase2<br/>")),
            "Coefficient that determines the variance in the number of offspring per mating for phase 2"
          ),
          value = 1000000,
          min = 0
        ),
        bsTooltip(id = "variance_offspring_phase2",
                  title = "This variable controls the variance of the negative binomial distribution that is used to determine the number of offspring that each mating pair produces")
      ),
      
      column(
        4,
        numericInput(
          "number_offspring_phase2",
          tags$div(tags$i(HTML(
            "number_offspring_phase2<br/>"
          )),
          "Mean number offspring per mating in phase 2"),
          value = 10,
          min = 0
        ),
        bsTooltip(id = "number_offspring_phase2",
                  title = "This variable controls the mean of the negative binomial distribution. This variable allows to control the number of offspring per mating which is convenient when there is a need that each pair of parents produce enough offspring in each generation for the population not to become extinct. However, in the simulations the mean number of offspring per mating each generation is effectively equal to two for two reasons: a) the population size remains constant from generation to generation, this means that  on average, across generations and replicates, two offspring per pair of parents are selected to become the parents of the next generation; and b) there is no variance in reproductive success (whether or not an individual gets to reproduce at all) because all individuals reproduce once")
      )
    ),
    
    fluidRow(
      column(
        4,
        radioButtons(
          "selection_phase2",
          tags$div(tags$i(HTML(
            "selection_phase2<br/>"
          )),
          "Whether selection occurs in phase 2"),
          choices = list("TRUE" = TRUE,
                         "FALSE" = FALSE),
          selected = FALSE
        ),
        bsTooltip(id = "selection_phase2",
                  title = "Selection is directional (replaces one allele by another) and multiplicative")
      ),
      
      column(
        4,
        numericInput(
          "Ne_phase2",
          tags$div(
            tags$i(HTML("Ne_phase2<br/>")),
            "Ne value to be used in the equation of the expected rate of loss of heterozygosity for phase 2"
          ),
          value = 50,
          min = 0
        ),
        bsTooltip(id = "Ne_phase2",
                  title = "Equation: He_t = He_0 (1-1 / 2 * Ne)^t, where He_0 is heterozygosity at generation 0 and t is the number of generations")
      ),
      
      column(
        4,
        numericInput(
          "Ne_fst_phase2",
          tags$div(
            tags$i(HTML("Ne_fst_phase2<br/>")),
            "Ne value to be used in the equation of the expected FST for phase 2"
          ),
          value = 50,
          min = 0
        ),
        bsTooltip(id = "Ne_fst_phase2",
                  title = "FST = 1/(4*Ne*m(n/(n-1))^2+1), where Ne is effective populations size of each individual subpopulation, m is dispersal rate and n the number of subpopulations (Takahata, 1983).")
      )
    ),
    
    hr(),
    
    h3(strong("Phase 1 variables")),
    
    fluidRow(
      
      column(
      4,
      radioButtons(
        "phase1",
        tags$div(tags$i(HTML("phase1<br/>")),
   "Whether phase 1 occur"),
        choices = list("TRUE" = TRUE,
                       "FALSE" = FALSE),
        selected = FALSE
      ),    
   bsTooltip(id = "phase1",
             title = "Information pending")
    )
   
    ),
    
    fluidRow(    
      
    column(
      4,
      radioButtons(
        "same_line",
        tags$div(tags$i(HTML("same_line<br/>")),
        "Whether phase 2 populations are sampled from the same phase 1 population or from different phase 1 populations"),
        choices = list(
          "Same population" = TRUE,
          "Different populations" = FALSE
        ),
        selected = FALSE
      ),    
      bsTooltip(id = "same_line",
                title = "Information pending")
    )
    
    ),
    
    fluidRow(
      
      column(
        4,
        numericInput(
          "number_pops_phase1",
          tags$div(tags$i(HTML(
            "number_pops_phase1<br/>"
          )),
          "Number of populations in phase 1"),
          value = 2,
          min = 0
        ),
        bsTooltip(id = "number_pops_phase1",
                  title = "Information pending")
      ),
      
      column(
        4,
        textInput(
          "population_size_phase1",
          tags$div(
            tags$i(HTML("population_size_phase1<br/>")),
            "Census population size of phase 1 (must be even and space delimited)"
          ),
          value = c("10", "10")
        ),
        bsTooltip(id = "population_size_phase1",
                  title = "Information pending")
      ),
      
      column(
        4,
        numericInput(
          "gen_number_phase1",
          tags$div(tags$i(HTML(
            "gen_number_phase1<br/>"
          )),
          "Number of generations of phase 1"),
          value = 10,
          min = 0
        ),
        bsTooltip(id = "gen_number_phase1",
                  title = "Information pending")
      )
    ),
    
    fluidRow(
      
      column(
        4,
        radioButtons(
          "dispersal_phase1",
          tags$div(tags$i(HTML(
            "dispersal_phase1<br/>"
          )),
          "Whether dispersal occurs in phase 1"),
          choices = list("TRUE" = TRUE, "FALSE" = FALSE),
          selected = TRUE
        ),
        bsTooltip(id = "dispersal_phase1",
                  title = "Dispersal between populations is symmetric and constant across generations. Dispersal rate (m) is the fraction of individuals in a population that is composed of dispersers or the probability that a randomly chosen individual in this generation came from a population different from the one in which it was found in the preceding generation (Holsinger, 2020, p. 93). Dispersal rate is calculated as (number_transfers / transfer_each_gen) / pop_size.")
      ),
      
      column(
        4,
        radioButtons(
          "dispersal_type_phase1",
          tags$div(tags$i(HTML(
            "dispersal_type_phase1<br/>"
          )),
          "Type of dispersal for phase 1"),
          choices = list(
            "All connected" = "all_connected",
            "Circle" = "circle",
            "Line" = "line"
          ),
          selected = "all_connected"
        ),
        bsTooltip(id = "dispersal_type_phase1",
                  title = "Information pending")
      ),
      
      column(
        4,
        numericInput(
          "number_transfers_phase1",
          tags$div(
            tags$i(HTML("number_transfers_phase1<br/>")),
            "Number of dispersing individuals in each dispersal event in phase 1"
          ),
          value = 1,
          min = 0
        ),
        bsTooltip(id = "number_transfers_phase1",
                  title = "Information pending")
      )
    ),
    
    fluidRow(
      column(
        4,
        numericInput(
          "transfer_each_gen_phase1",
          tags$div(
            tags$i(HTML("transfer_each_gen_phase1<br/>")),
            "Interval of number of generations in which a dispersal event occurs in phase 1"
          ),
          value = 1,
          min = 0
        ),
        bsTooltip(id = "transfer_each_gen_phase1",
                  title = "Information pending")
      ),
      
      column(
        4,
        numericInput(
          "variance_offspring_phase1",
          tags$div(
            tags$i(HTML("variance_offspring_phase1<br/>")),
            "Coefficient that determines the variance in the number of offspring per mating for phase 1"
          ),
          value = 1000000,
          min = 0
        ),
        bsTooltip(id = "variance_offspring_phase1",
                  title = "This variable controls the variance of the negative binomial distribution that is used to determine the number of offspring that each mating pair produces")
      ),
      
      column(
        4,
        numericInput(
          "number_offspring_phase1",
          tags$div(tags$i(HTML(
            "number_offspring_phase1<br/>"
          )),
          "Mean number offspring per mating in phase 1"),
          value = 10,
          min = 0
        ),
        bsTooltip(id = "number_offspring_phase1",
                  title = "This variable controls the mean of the negative binomial distribution. This variable allows to control the number of offspring per mating which is convenient when there is a need that each pair of parents produce enough offspring in each generation for the population not to become extinct. However, in the simulations the mean number of offspring per mating each generation is effectively equal to two for two reasons: a) the population size remains constant from generation to generation, this means that  on average, across generations and replicates, two offspring per pair of parents are selected to become the parents of the next generation; and b) there is no variance in reproductive success (whether or not an individual gets to reproduce at all) because all individuals reproduce once")
      )
    ),
    
    fluidRow(
      column(
        4,
        radioButtons(
          "selection_phase1",
          tags$div(tags$i(HTML(
            "selection_phase1<br/>"
          )),
          "Whether selection occurs in phase 1"),
          choices = list("TRUE" = TRUE,
                         "FALSE" = FALSE),
          selected = FALSE
        ),
        bsTooltip(id = "selection_phase1",
                  title = "Selection is directional (replaces one allele by another) and multiplicative")
      ),
      
      column(
        4,
        numericInput(
          "Ne_phase1",
          tags$div(
            tags$i(HTML("Ne_phase1<br/>")),
            "Ne value to be used in the equation of the expected rate of loss of heterozygosity for phase 1"
          ),
          value = 50,
          min = 0
        ),
        bsTooltip(id = "Ne_phase1",
                  title = "Equation: He_t = He_0 (1-1 / 2 * Ne)^t, where He_0 is heterozygosity at generation 0 and t is the number of generations")
      ),
      
      column(
        4,
        numericInput(
          "Ne_fst_phase1",
          tags$div(
            tags$i(HTML("Ne_fst_phase1<br/>")),
            "Ne value to be used in the equation of the expected FST for phase 1"
          ),
          value = 50,
          min = 0
        ),
        bsTooltip(id = "Ne_fst_phase1",
                  title = "FST = 1/(4*Ne*m(n/(n-1))^2+1), where Ne is effective populations size of each individual subpopulation, m is dispersal rate and n the number of subpopulations (Takahata, 1983).")
      )
    ),
    
    hr(),
    
    h3(strong("Initialisation variables")),
    
    fluidRow(
      
      column(
      4,
      textInput(
        "chromosome_name",
        tags$div(tags$i(HTML("chromosome_name<br/>")),
        "Name of the chromosome to be simulated"),
        value = "2L"
      ),    
      bsTooltip(id = "chromosome_name",
                title = "Information pending")
    )
    
    ),
    
    hr(),
    
    h3(strong("Recombination variables")),
    
    fluidRow(
      
      column(
      4,
      radioButtons(
        "recombination",
        tags$div(tags$i(HTML("recombination<br/>")),
        "Whether recombination occurs"),
        choices = list("TRUE" = TRUE,
                       "FALSE" = FALSE),
        selected = TRUE
      ),    
      bsTooltip(id = "recombination",
                title = "Information pending")
    ),
    
    column(
      4,
      radioButtons(
        "recombination_males",
        tags$div(tags$i(HTML("recombination_males<br/>")),
        "Whether recombination occurs in males and females (or only females)"),
        choices = list("TRUE" = TRUE,
                       "FALSE" = FALSE),
        selected = TRUE
      ),    
      bsTooltip(id = "recombination_males",
                title = "Information pending")
    )
    
    ),
    
    hr(),
    
    h3(strong("Selection variables")),
    
    fluidRow(
      
      column(
      4,
      radioButtons(
        "natural_selection_model",
        tags$div(tags$i(HTML("natural_selection_model<br/>")),
        "Selection model to use"),
        choices = list("Relative" = "relative",
                       "Absolute" = "absolute"),
        selected = "relative"
      ),    
      bsTooltip(id = "natural_selection_model",
                title = "Information pending")
    ),
    
    column(
      4,
      numericInput(
        "genetic_load",
        tags$div(tags$i(HTML("genetic_load<br/>")),
        "Approximation of the genetic load of the proportion of the genome that is simulated. This variable is used in the absolute fitness model"),
        value = 0.8,
        min = 0
      ),    
      bsTooltip(id = "genetic_load",
                title = "Information pending")
    )
    
    ),
    
    hr(),
    
    fluidRow(column(
      9,
      actionButton(
        "close",
        label = h4(strong("RUN")),
        icon = icon("play"),
        width = "100%",
        class = "btn-success"
      )
    )),
    
    hr()
    
  )
  
  server <- function(input, output, session) {
    # phase 2
    
    observeEvent(input$number_pops_phase2, {
      assign("number_pops_phase2",input$number_pops_phase2,envir = env_fun)
    })
    observeEvent(input$population_size_phase2, {
      assign("population_size_phase2",input$population_size_phase2,envir = env_fun)
    })
    observeEvent(input$gen_number_phase2, {
      assign("gen_number_phase2",input$gen_number_phase2,envir = env_fun)
    })
    observeEvent(input$dispersal_phase2, {
      assign("dispersal_phase2",input$dispersal_phase2,envir = env_fun)
    })
    observeEvent(input$dispersal_type_phase2, {
      assign("dispersal_type_phase2",input$dispersal_type_phase2,envir = env_fun)
    })
    observeEvent(input$number_transfers_phase2, {
      assign("number_transfers_phase2",input$number_transfers_phase2,envir = env_fun)
    })
    observeEvent(input$transfer_each_gen_phase2, {
      assign("transfer_each_gen_phase2",input$transfer_each_gen_phase2,envir = env_fun)
    })
    observeEvent(input$variance_offspring_phase2, {
      assign("variance_offspring_phase2",input$variance_offspring_phase2,envir = env_fun)
    })
    observeEvent(input$number_offspring_phase2, {
      assign("number_offspring_phase2",input$number_offspring_phase2,envir = env_fun)
    })
    observeEvent(input$selection_phase2, {
      assign("selection_phase2",input$selection_phase2,envir = env_fun)
    })
    observeEvent(input$Ne_phase2, {
      assign("Ne_phase2",input$Ne_phase2,envir = env_fun)
    })
    observeEvent(input$Ne_fst_phase2, {
      assign("Ne_fst_phase2",input$Ne_fst_phase2,envir = env_fun)
    })
    
    # phase 1
    
    observeEvent(input$phase1, {
      assign("phase1",input$phase1,envir = env_fun)
      
      toggleElement(
        id = "same_line",
        condition = input$phase1 == TRUE
        )
      
      toggleElement(
        id = "number_pops_phase1",
        condition = input$phase1 == TRUE
      )
      
      toggleElement(
        id = "population_size_phase1",
        condition = input$phase1 == TRUE
      )

      toggleElement(
        id = "gen_number_phase1",
        condition = input$phase1 == TRUE
      )

      toggleElement(
        id = "dispersal_phase1",
        condition = input$phase1 == TRUE
      )

      toggleElement(
        id = "dispersal_type_phase1",
        condition = input$phase1 == TRUE
      )

      toggleElement(
        id = "number_transfers_phase1",
        condition = input$phase1 == TRUE
      )

      toggleElement(
        id = "transfer_each_gen_phase1",
        condition = input$phase1 == TRUE
      )

      toggleElement(
        id = "variance_offspring_phase1",
        condition = input$phase1 == TRUE
      )

      toggleElement(
        id = "number_offspring_phase1",
        condition = input$phase1 == TRUE
      )

      toggleElement(
        id = "selection_phase1",
        condition = input$phase1 == TRUE
      )

      toggleElement(
        id = "Ne_phase1",
        condition = input$phase1 == TRUE
      )

      toggleElement(
        id = "Ne_fst_phase1",
        condition = input$phase1 == TRUE
      )
      
    })
    observeEvent(input$same_line, {
      assign("same_line",input$same_line,envir = env_fun)
    })
    observeEvent(input$number_pops_phase1, {
      assign("number_pops_phase1",input$number_pops_phase1,envir = env_fun)
    })
    observeEvent(input$population_size_phase1, {
      assign("population_size_phase1",input$population_size_phase1,envir = env_fun)
    })
    observeEvent(input$gen_number_phase1, {
      assign("gen_number_phase1",input$gen_number_phase1,envir = env_fun)
    })
    observeEvent(input$dispersal_phase1, {
      assign("dispersal_phase1",input$dispersal_phase1,envir = env_fun)
    })
    observeEvent(input$dispersal_type_phase1, {
      assign("dispersal_type_phase1",input$dispersal_type_phase1,envir = env_fun)
    })
    observeEvent(input$number_transfers_phase1, {
      assign("number_transfers_phase1",input$number_transfers_phase1,envir = env_fun)
    })
    observeEvent(input$transfer_each_gen_phase1, {
      assign("transfer_each_gen_phase1",input$transfer_each_gen_phase1,envir = env_fun)
    })
    observeEvent(input$variance_offspring_phase1, {
      assign("variance_offspring_phase1",input$variance_offspring_phase1,envir = env_fun)
    })
    observeEvent(input$number_offspring_phase1, {
      assign("number_offspring_phase1",input$number_offspring_phase1,envir = env_fun)
    })
    observeEvent(input$selection_phase1, {
      assign("selection_phase1",input$selection_phase1,envir = env_fun)
    })
    observeEvent(input$Ne_phase1, {
      assign("Ne_phase1",input$Ne_phase1,envir = env_fun)
    })
    observeEvent(input$Ne_fst_phase1, {
      assign("Ne_fst_phase1",input$Ne_fst_phase1,envir = env_fun)
    })
    
    # real dataset
    
    observeEvent(input$real_pops, {
      assign("real_pops",input$real_pops,envir = env_fun)
    })
    observeEvent(input$real_pop_size, {
      assign("real_pop_size",input$real_pop_size,envir = env_fun)
    })
    observeEvent(input$real_freq, {
      assign("real_freq",input$real_freq,envir = env_fun)
    })
    
    # recombination
    
    observeEvent(input$recombination, {
      assign("recombination",input$recombination,envir = env_fun)
    })
    observeEvent(input$recombination_males, {
      assign("recombination_males",input$recombination_males,envir = env_fun)
    })
    
    # selection
    
    observeEvent(input$genetic_load, {
      assign("genetic_load",input$genetic_load,envir = env_fun)
    })
    observeEvent(input$natural_selection_model, {
      assign("natural_selection_model",input$natural_selection_model,envir = env_fun)
    })
    
    # intialisation
    
    observeEvent(input$chromosome_name, {
      assign("chromosome_name",input$chromosome_name,envir = env_fun)
    })
    
    observeEvent(input$close, {
      assign("sim_vars_temp", as.data.frame(cbind(
        c(
          "number_pops_phase2",
          "population_size_phase2",
          "gen_number_phase2",
          "dispersal_phase2",
          "dispersal_type_phase2",
          "number_transfers_phase2",
          "transfer_each_gen_phase2",
          "variance_offspring_phase2",
          "number_offspring_phase2",
          "selection_phase2",
          "Ne_phase2",
          "Ne_fst_phase2",
          "number_pops_phase1",
          "population_size_phase1",
          "gen_number_phase1",
          "dispersal_phase1",
          "dispersal_type_phase1",
          "number_transfers_phase1",
          "transfer_each_gen_phase1",
          "variance_offspring_phase1",
          "number_offspring_phase1",
          "selection_phase1",
          "Ne_phase1",
          "Ne_fst_phase1",
          "phase1",
          "same_line",
          "real_freq",
          "real_pops",
          "real_pop_size",
          "recombination",
          "recombination_males",
          "genetic_load",
          "natural_selection_model",
          "chromosome_name"
        ),
        c(
          input$number_pops_phase2,
          input$population_size_phase2,
          input$gen_number_phase2,
          input$dispersal_phase2,
          input$dispersal_type_phase2,
          input$number_transfers_phase2,
          input$transfer_each_gen_phase2,
          input$variance_offspring_phase2,
          input$number_offspring_phase2,
          input$selection_phase2,
          input$Ne_phase2,
          input$Ne_fst_phase2,
          input$number_pops_phase1,
          input$population_size_phase1,
          input$gen_number_phase1,
          input$dispersal_phase1,
          input$dispersal_type_phase1,
          input$number_transfers_phase1,
          input$transfer_each_gen_phase1,
          input$variance_offspring_phase1,
          input$number_offspring_phase1,
          input$selection_phase1,
          input$Ne_phase1,
          input$Ne_fst_phase1,
          input$phase1,
          input$same_line,
          input$real_freq,
          input$real_pops,
          input$real_pop_size,
          input$recombination,
          input$recombination_males,
          input$genetic_load,
          input$natural_selection_model,
          input$chromosome_name
        )
      )))
      
      colnames(sim_vars_temp) <- c("variable", "value")
      
      assign("sim_vars",sim_vars_temp,envir = env_fun)
      
      stopApp()
      
    })
    
  }

   # runGadget(shinyApp(ui, server), viewer =  paneViewer())
  # runGadget(shinyApp(ui, server), viewer = dialogViewer("", width = 600, height = 600))
   runApp(shinyApp(ui, server))
  

}
