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
    parents_matrix <- as.data.frame(matrix(nrow = pop_size / 2, ncol = 2))
    parents_matrix[, 1] <- sample(rownames(pop[1:(pop_size / 2),]), 
                                  size = pop_size / 2)
    parents_matrix[, 2] <- sample(rownames(pop[((pop_size / 2) + 1):pop_size,]),
                                  size = pop_size / 2)
    offspring <- NULL
    for (parent in 1:dim(parents_matrix)[1]) {
      pairing_offspring <- rnbinom(1, size = var_off, mu = num_off)
      offspring_temp <- as.data.frame(matrix(nrow = pairing_offspring, ncol = 6))
      if (pairing_offspring < 1) {
        next
      }
      offspring_temp[, 1] <- sample(c("Male", "Female"), 
                                    size = pairing_offspring, 
                                    replace = TRUE) # sex
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
              recomb(
                chr1 = male_chromosomes[[1]],
                chr2 = male_chromosomes[[2]],
                r_map = r_map_1,
                loci = n_loc)
          }
          offspring_temp[offs, 3] <- male_chromosomes[[sample(c(1, 2), 1)]]
        } else{
          offspring_temp[offs, 3] <- male_chromosomes[[sample(c(1, 2), 1)]]
        }
        #recombination in females
        if (recom == TRUE & females_recom_events > 1) {
          for (event in females_recom_events) {
            female_chromosomes <-
              recomb(
                chr1 = female_chromosomes[[1]],
                chr2 = female_chromosomes[[2]],
                r_map = r_map_1,
                loci = n_loc
              )
          }
          offspring_temp[offs, 4] <- female_chromosomes[[sample(c(1, 2), 1)]]
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
########################## RECOMBINATION ######################################
###############################################################################

recomb <- function(chr1, 
                   chr2,
                   r_map, 
                   loci) {
  chiasma <- as.numeric(sample(row.names(r_map), size = 1, prob = r_map[, "c"]))
  if (chiasma < (loci + 1)) {
    r_chr1 <- paste0(substr(chr1,1,chiasma), substr(chr2, chiasma + 1,loci))
    r_chr2 <- paste0(substr(chr2,1,chiasma), substr(chr1, chiasma + 1,loci))
    return(list(r_chr1, r_chr2))
  } else{
    return(list(chr1, chr2))
  }
}

###############################################################################
########################## SELECTION ##########################################
###############################################################################

selection_fun <-
  function(offspring,
           h,
           s,
           sel_model,
           g_load) {
    
    
    # make genotypes
    #to hack package checking...
    make_fit <- function(){}  
    
    Rcpp::cppFunction(plugins="cpp11",
                      
"NumericVector make_fit(StringMatrix seqs, NumericVector h, NumericVector s){
int loc_number = strlen(seqs(0,0));
int indN = seqs.nrow();
NumericVector out(indN);
for (int i = 0; i < indN; i++) {
NumericVector fit_ind(loc_number);
fit_ind.fill(1);
  for (int loc = 0; loc < loc_number; loc++){
    char chr1 = seqs(i,0)[loc], chr2 = seqs(i,1)[loc];
    if (chr1 == chr2 && chr1=='1')
    fit_ind[loc] = 1 - s[loc];
    if (chr1 != chr2)
    fit_ind[loc] = 1 - (h[loc] * s[loc]);
  }
  out[i] = algorithm::prod(fit_ind.begin(), fit_ind.end());
}
  return(out);
}"
    )

    offspring$fitness <- make_fit(as.matrix(offspring[,3:4]),h,s)
    if (sel_model == "absolute") {
      offspring$random_deviate <- runif(nrow(offspring), min = 0, max = g_load)
      offspring$alive <- offspring$fitness > offspring$random_deviate
      offspring <- offspring[which(offspring$alive == TRUE),]
    }
    if (sel_model == "relative") {
      fitnes_proportion <- sum(offspring$fitness)
      offspring$relative_fitness <- offspring$fitness / fitnes_proportion
      offspring$relative_fitness[offspring$relative_fitness < 0] <- 0
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
           ref,
           p_map,
           s_vars,
           g) {
    pop_names <- rep(as.character(p_vector), p_size)
    df_genotypes <- rbindlist(p_list)
    df_genotypes$V1[df_genotypes$V1 == "Male"]   <- 1
    df_genotypes$V1[df_genotypes$V1 == "Female"] <- 2
    df_genotypes[, 2] <- pop_names
    df_genotypes$id <-
      paste0(unlist(unname(df_genotypes[, 2])), "_", unlist(lapply(p_size, function(x) {
        1:x
      })))
    
    # make genotypes
    #to hack package checking...
    make_geno <- function(){}  
    
    Rcpp::cppFunction(plugins="cpp11",
                      
'List make_geno(StringMatrix mat) {
    int ind = mat.nrow();
    int loc = strlen(mat(0,0));
    List out(ind);
for (int i = 0; i < ind; i++) {
 std::string chr1 (mat(i,0));
 std::string chr2 (mat(i,1));
 StringVector temp(loc);
for (int j = 0; j < loc; j++) {
 StringVector geno = StringVector::create(chr1[j],chr2[j]);
    temp[j] = collapse(geno);
  }
      out[i] = temp;
    }
    return out;
  }'
    )
    
    plink_temp <- as.matrix(df_genotypes[,3:4])
    plink_ped <- make_geno(plink_temp)
    plink_ped_2 <- lapply(plink_ped, function(x) {
      x[x == "11"] <- 2
      x[x == "00"] <- 0
      x[x == "01"] <- 1
      x[x == "10"] <- 1
      
      return(x)
      
    })
    
    loc.names <- 1:nrow(ref)
    n.loc <- length(loc.names)
    misc.info <- lapply(1:6, function(i){NULL})
    names(misc.info) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
    res <- list()
    temp <-
      as.data.frame(
        cbind(
          df_genotypes[, "V2"],
          df_genotypes[, "id"],
          df_genotypes[, "V5"],
          df_genotypes[, "V6"],
          df_genotypes[, "V1"],
          1
        )
      )
    
    for (i in 1:6) {
      misc.info[[i]] <- temp[, i]
    }
    txt <-
      lapply(plink_ped_2, function(e)
        suppressWarnings(as.integer(e)))
      res <-
        c(res, lapply(txt, function(e)
          new(
            "SNPbin", snp = e, ploidy = 2L
          )))
    
    res <- new("genlight", res, ploidy = 2L)
    
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
    res$other$loc.metrics <- ref
    chromosome(res) <- res$other$loc.metrics$chr_name
    position(res) <- res$other$loc.metrics$loc_bp
    res$loc.all <- rep("G/C",nLoc(res))
    res$other$sim.vars <- s_vars
    res <- utils.reset.flags(res,verbose=0)
    
    
    return(res)
    
  }

###############################################################################
######### SHINY APP FOR THE VARIABLES OF THE REFERENCE TABLE ###################
###############################################################################
#' @name interactive_reference
#' @title Shiny app for the input of the reference table for the simulations
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

interactive_reference <- function() {
  
  pkg <- "shinyBS"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }
  
  pkg <- "shinythemes"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }
  
  pkg <- "shinyjs"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }
  
  pkg <- "shinyWidgets"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }
  
  ui <- fluidPage(
    
    shinyjs::useShinyjs(),
    
    theme = shinythemes::shinytheme("superhero"),

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
    
    h3(strong("Basic variables")),
    
    fluidRow(
      
      column(
        4,
        numericInput(
          "chunk_number",  
          tags$div(tags$i(HTML("chunk_number<br/>")),
            "Number of chromosome chunks"
          ),
          value = 100,
          min = 0
        ), 
        shinyBS::bsTooltip(id = "chunk_number",
                        title = "The lenght of the chromosome is chunk_number * chunk_bp and chunk_number * chunk_cM")
      ),
      
      column(
        4,
        numericInput(
          "chunk_bp",
          tags$div(tags$i(HTML("chunk_bp<br/>")),
                   "Number of basepairs (bp) per chromosome chunk"),
          value = 100000,
          min = 0
        ),    
        shinyBS::bsTooltip(id = "chunk_bp",
                           title = "This variable can be also seen as the resolution of the recombination map")
      ),
      
      column(
        4,
        numericInput(
          "chunk_cM",
          tags$div(tags$i(HTML("chunk_cM<br/>")),
                   "Number of centiMorgans (cM) per chromosome chunk"),
          value = 10,
          min = 0
        ),    
        shinyBS::bsTooltip(id = "chunk_cM",
                           title = "1 cM corresponds to one percent of probability that two loci will be separated by a recombination event in each meiosis")
      )
      
      ),
    
    fluidRow(
      
      column(
        4,
        numericInput(
          "chunk_neutral_loci",  
          tags$div(tags$i(HTML("chunk_neutral_loci<br/>")),
                   "Number of loci with neutral alleles per chromosome chunk"
          ),
          value = 1,
          min = 0
        ),    
        shinyBS::bsTooltip(id = "chunk_neutral_loci",
                  title = "Neutral alleles have no effect on fitness and are evenly distributed in each genome chunk")
      ),
      
      column(
        4,
        sliderInput(
          "q_neutral", 
          tags$div(tags$i(HTML("q_neutral<br/>")),
                   "Initial frequency for all neutral alleles"),
          value = 0.5,
          min = 0,
          max = 1
        ),    
        shinyBS::bsTooltip(id = "q_neutral",
                           title = "Loci with neutral alleles are bi-allelic")
      )
      
    ),
    
    fluidRow(
      
      column(
        4,
        numericInput(
          "loci_deleterious", 
          tags$div(tags$i(HTML("loci_deleterious<br/>")),
                   "Number of loci with deleterious alleles"
          ),
          value = 0,
          min = 0
        ),     
        shinyBS::bsTooltip(id = "loci_deleterious",
                           title = "Deleterious alleles have detrimental effects on fitness")
      ),
      
      column(
        4,
        numericInput(
          "loci_advantageous", 
          tags$div(tags$i(HTML("loci_advantageous<br/>")),
                   "Number of loci with advantageous alleles"
          ),
          value = 0,
          min = 0
        ),     
        shinyBS::bsTooltip(id = "loci_advantageous",
                           title = "Advantageous alleles have beneficial effects on fitness")
      )
      
    ),
    
    hr(),
    
    fluidRow(
      
      column(
        4,
        shinyWidgets::awesomeCheckbox(inputId = "s_coeff",
                        label = h5("Selection coefficient (s) variables"),
                        value = FALSE,
                        status = "success")
      )
      
    ),
    
    fluidRow(
    
    column(
      4,
      radioButtons(
        "s_distribution_del",
        tags$div(tags$i(HTML("s_distribution_del<br/>")),
                 "Distribution to sample selection coefficients for deleterious alleles"),
        choices = list(
          "All equal" = "equal",
          "Gamma distribution" = "gamma",
          "Log normal distribution" = "log_normal"
          ),
        selected = "equal"
      ),    
      shinyBS::bsTooltip(id = "s_distribution_del",
                         title = "s ranges from 0 to 1, where s = 0 means that allele has no effect on fitness and s = 1 means allele is lethal")
    ),
    
    column(
      4,
      radioButtons(
        "s_distribution_adv",
        tags$div(tags$i(HTML("s_distribution_adv<br/>")),
                 "Distribution to sample selection coefficients for advantageous alleles"),
        choices = list(
          "All equal" = "equal",
          "Exponential distribution" = "exponential"
        ),
        selected = "equal"
      ),    
      shinyBS::bsTooltip(id = "s_distribution_adv",
                         title = "s ranges from 0 to 1, where s = 0 means that allele has no effect on fitness")
    )
    
    ),
    
    fluidRow(
      
      column(
        4,
        numericInput(
          "s_del",
          tags$div(tags$i(HTML("s_del<br/>")),
                   "Selection coefficient for all deleterious alleles"),
          value = 0.001,
          min = 0
        ),    
        shinyBS::bsTooltip(id = "s_del",
                           title = "")
      ),
      
      column(
        4,
        numericInput(
          "s_adv",
          tags$div(tags$i(HTML("s_adv<br/>")),
                   "Selection coefficient for all advantageous alleles"),
          value = 0.001,
          min = 0
        ),    
        shinyBS::bsTooltip(id = "s_adv",
                           title = "")
      )
      
    ),
    
    fluidRow(
      
      column(
        4,
        numericInput(
          "gamma_scale",
          tags$div(tags$i(HTML("gamma_scale<br/>")),
                   "Scale of gamma distribution for deleterious alleles"),
          value = 0.03,
          min = 0
        ),    
        shinyBS::bsTooltip(id = "gamma_scale",
                           title = "See Huber et al., 2017")
      ),
      
      column(
        4,
        numericInput(
          "exp_rate",
          tags$div(tags$i(HTML("exp_rate<br/>")),
                   "Mean of the exponential distribution for advantageous alleles"),
          value = 16,
          min = 0
        ),    
        shinyBS::bsTooltip(id = "exp_rate",
                           title = "See Tataru et al., 2017")
      )
      
    ),
    
    fluidRow(

      column(
        4,
        numericInput(
          "gamma_shape",
          tags$div(tags$i(HTML("gamma_shape<br/>")),
                   "Shape of gamma distribution for deleterious alleles"),
          value = 0.25,
          min = 0
        ),    
        shinyBS::bsTooltip(id = "gamma_shape",
                           title = "See Huber et al., 2017")
      )
      
    ),
    
    fluidRow(
      
      column(
        4,
        numericInput(
          "log_mean",
          tags$div(tags$i(HTML("log_mean<br/>")),
                   "Mean of log normal distribution for deleterious alleles"),
          value = 0.002,
          min = 0
        ),    
        shinyBS::bsTooltip(id = "log_mean",
                           title = "See Charlesworth, 2015")
      )
      
      ),
    
    fluidRow(
      
      column(
        4,
        numericInput(
          "log_sd",
          tags$div(tags$i(HTML("log_sd<br/>")),
                   "Standard deviation of log normal distribution for deleterious alleles"),
          value = 4,
          min = 0
        ),    
        shinyBS::bsTooltip(id = "log_sd",
                           title = "See Charlesworth, 2015")
      )
      
    ),
    
    hr(),
    
    fluidRow(
      
      column(
        4,
        shinyWidgets::awesomeCheckbox(inputId = "dominance",
                        label = h5("Dominance (h) variables"),
                        value = FALSE,
                        status = "success")
      )
      
    ),
    
    fluidRow(
      
      column(
        4,
        radioButtons(
          "h_distribution_del",
          tags$div(tags$i(HTML("h_distribution_del<br/>")),
                   "Distribution to sample dominance for deleterious alleles"),
          choices = list(
            "All equal" = "equal",
            "Normal distribution" = "normal",
            "From equation" = "equation"
          ),
          selected = "equal"
        ),    
        shinyBS::bsTooltip(id = "h_distribution_del",
                           title = "h ranges from 0 to 1, where h = 0 is completely recessive and h = 1 is completely dominant")
      ),
      
      column(
        4,
        radioButtons(
          "h_distribution_adv",
          tags$div(tags$i(HTML("h_distribution_adv<br/>")),
                   "Distribution to sample dominance for advantageous alleles"),
          choices = list(
            "All equal" = "equal",
            "Normal distribution" = "normal",
            "From equation" = "equation"
          ),
          selected = "equal"
        ),    
        shinyBS::bsTooltip(id = "h_distribution_adv",
                           title = "")
      )
      
    ),
    
    fluidRow(
      
      column(
        4,
        sliderInput(
          "h_del",
          tags$div(tags$i(HTML("h_del<br/>")),
                   "Dominance for all deleterious alleles"),
          value = 0.25,
          min = 0,
          max = 1
        ),    
        shinyBS::bsTooltip(id = "h_del",
                           title = "")
      ),
      
      column(
        4,
        sliderInput(
          "h_adv",
          tags$div(tags$i(HTML("h_adv<br/>")),
                   "Dominance for all advantageous alleles"),
          value = 0.25,
          min = 0,
          max = 1
        ),    
        shinyBS::bsTooltip(id = "h_adv",
                           title = "")
      )
      
    ),
    
    fluidRow(
      
      column(
        4,
        sliderInput(
          "h_mean_del",
          tags$div(tags$i(HTML("h_mean_del<br/>")),
                   "Mean of normal distribution for deleterious alleles"),
          value = 0.25,
          min = 0,
          max = 1
        ),    
        shinyBS::bsTooltip(id = "h_mean_del",
                           title = "See Charlesworth, 2015")
      ),
      
      column(
        4,
        sliderInput(
          "h_mean_adv",
          tags$div(tags$i(HTML("h_mean_adv<br/>")),
                   "Mean of normal distribution for advantageous alleles"),
          value = 0.25,
          min = 0,
          max = 1
        ),    
        shinyBS::bsTooltip(id = "h_mean_adv",
                           title = "See Charlesworth, 2015")
      )
      
      ),
    
    fluidRow(
      
      column(
        4,
        numericInput(
          "h_sd_del",
          tags$div(tags$i(HTML("h_sd_del<br/>")),
                   "Standard deviation of normal distribution for deleterious alleles"),
          value = sqrt(0.001),
          min = 0
        ),    
        shinyBS::bsTooltip(id = "h_sd_del",
                           title = "See Charlesworth, 2015")
      ),
      
      column(
        4,
        numericInput(
          "h_sd_adv",
          tags$div(tags$i(HTML("h_sd_adv<br/>")),
                   "Standard deviation of normal distribution for advantageous alleles"),
          value = sqrt(0.001),
          min = 0
        ),    
        shinyBS::bsTooltip(id = "h_sd_adv",
                           title = "See Charlesworth, 2015")
      )
      
    ),
    
    fluidRow(
      
      column(
        4,
        sliderInput(
          "h_intercept_del",
          tags$div(tags$i(HTML("h_intercept_del<br/>")),
                   "Value for the intercept of the equation for deleterious alleles"),
          value = 0.5,
          min = 0,
          max = 1
        ),    
        shinyBS::bsTooltip(id = "h_intercept_del",
                           title = "See Huber et al., 2018")
      ),
      
      column(
        4,
        sliderInput(
          "h_intercept_adv",
          tags$div(tags$i(HTML("h_intercept_adv<br/>")),
                   "Value for the intercept of the equation for advantageous alleles"),
          value = 0.5,
          min = 0,
          max = 1
        ),    
        shinyBS::bsTooltip(id = "h_intercept_adv",
                           title = "See Huber et al., 2018")
      )
      
      ),
    
    fluidRow(
      
      column(
        4,
        numericInput(
          "h_rate_del",
          tags$div(tags$i(HTML("h_rate_del<br/>")),
                   "Value for the variable rate of the equation for deleterious alleles"),
          value = 500,
          min = 0
        ),    
        shinyBS::bsTooltip(id = "h_rate_del",
                           title = "See Huber et al., 2018")
      ),
      
      column(
        4,
        numericInput(
          "h_rate_adv",
          tags$div(tags$i(HTML("h_rate_adv<br/>")),
                   "Value for the variable rate of the equation for advantageous alleles"),
          value = 500,
          min = 0
        ),    
        shinyBS::bsTooltip(id = "h_rate_adv",
                           title = "See Huber et al., 2018")
      )
      
    ),
    
    hr(),
    
    fluidRow(
      
      column(
        4,
        shinyWidgets::awesomeCheckbox(inputId = "frequency",
                        label = h5("Initial frequency (q) variables"),
                        value = FALSE,
                        status = "success")
      )
      
    ),
    
    fluidRow(
      
      column(
        4,
        shinyWidgets::prettyRadioButtons(
          inputId = "q_distribution_del",
          label =  tags$div(tags$i(HTML("q_distribution_del<br/>")),
                            "Method to determine initial allele frequency for deleterious alleles"),
          choices = list("All equal" = "equal", "From equation" = "equation"),
          selected = "equal"
        ),    
        shinyBS::bsTooltip(id = "q_distribution_del",
                           title = "")
        
      ),
      
      column(
        4,
        shinyWidgets::prettyRadioButtons(
          inputId = "q_distribution_adv",
          label =  tags$div(tags$i(HTML("q_distribution_adv<br/>")),
                            "Method to determine initial allele frequency for advantageous alleles"),
          choices = list("All equal" = "equal", "From equation" = "equation"),
          selected = "equal"
        ),    
        shinyBS::bsTooltip(id = "q_distribution_adv",
                           title = "")
        
      )
      
    ),
    
    fluidRow(
      
      column(
        4,
        sliderInput(
          "q_del", 
          tags$div(tags$i(HTML("q_del<br/>")),
                   "Initial frequencies for all deleterious alleles"),
          value = 0.05,
          min = 0,
          max = 1
        ),    
        shinyBS::bsTooltip(id = "q_del",
                           title = "")
      ),
      
      column(
        4,
        sliderInput(
          "q_adv", 
          tags$div(tags$i(HTML("q_adv<br/>")),
                   "Initial frequencies for all advantageous alleles"),
          value = 0.05,
          min = 0,
          max = 1
        ),    
        shinyBS::bsTooltip(id = "q_adv",
                           title = "")
      )
      
      ),
    
    fluidRow(
      
      column(
        4,
        numericInput(
          "q_equation_del", 
          tags$div(tags$i(HTML("q_equation_del<br/>")),
                   "Mutation rate per generation per site (u) to be used in the equation that approximates the mean frequency of a recessive deleterious variant sampled from a large population in mutation-selection equilibrium"),
          value = 5 * 10 ^ -5
        ),    
        shinyBS::bsTooltip(id = "q_equation_del",
                           title = "Equation: s(1-2h)q^2 + hs(1+u)q - u = 0, where u is the mutation rate per generation per site; see Crow & Kimura, 1970, p. 260")
      ),
      
      column(
        4,
        numericInput(
          "q_equation_adv", 
          tags$div(tags$i(HTML("q_equation_adv<br/>")),
                   "Mutation rate per generation per site (u) to be used in the equation that approximates the mean frequency of a recessive deleterious variant sampled from a large population in mutation-selection equilibrium"),
          value = 5 * 10 ^ -5
        ),    
        shinyBS::bsTooltip(id = "q_equation_adv",
                           title = "Equation: s(1-2h)q^2 + hs(1+u)q - u = 0, where u is the mutation rate per generation per site; see Crow & Kimura, 1970, p. 260")
      )
      
    ),
    
    hr(),
    
    h5(
      em(
        "The number of loci available to mutation must be defined before the start of simulations. To ensure that there are enough loci available to mutation, it is necessary to consider the values defined in the variables: mutation rate per genome per generation, population size and number of generations."
      )
    ),
    
    fluidRow(
      
      column(
        4,
        shinyWidgets::awesomeCheckbox(inputId = "mutation",
                        label = h5("Mutation variables"),
                        value = FALSE,
                        status = "success")
      )
      
    ),
    
    fluidRow(
      
      column(
        4,
        numericInput(
          "loci_mut_neu", 
          tags$div(tags$i(HTML("loci_mut_neu<br/>")),
                   "Number of loci with neutral alleles available to mutation"
          ),
          value = 0,
          min = 0
        ),     
        shinyBS::bsTooltip(id = "loci_mut_neu",
                           title = "")
      ),
      
      column(
        4,
        numericInput(
          "loci_mut_del", 
          tags$div(tags$i(HTML("loci_mut_del<br/>")),
                   "Number of loci with deleterious alleles available to mutation "
          ),
          value = 0,
          min = 0
        ),     
        shinyBS::bsTooltip(id = "loci_mut_del",
                           title = "s and h values are taken from the values set for deleterious alleles")
      ),
      
      column(
        4,
        numericInput(
          "loci_mut_adv", 
          tags$div(tags$i(HTML("loci_mut_adv<br/>")),
                   "Number of loci with advantageous alleles available to mutation "
          ),
          value = 0,
          min = 0
        ),     
        shinyBS::bsTooltip(id = "loci_mut_adv",
                           title = "s and h values are taken from the values set for advantageous alleles")
      )
      
    ),
    
    hr(),
    
    fluidRow(
      
      column(
        4,
        shinyWidgets::awesomeCheckbox(inputId = "real_v",
                        label = h5("Real information variables"),
                        value = FALSE,
                        status = "success")
      )
      
    ),
    
    fluidRow(
    
    column(
      4,
      radioButtons(
        "real_loc",
        tags$div(tags$i(HTML("real_loc<br/>")),
                 "Extract location of neutral loci from genlight object"
        ) ,
        choices = list("TRUE" = TRUE,
                       "FALSE" = FALSE),
        selected = FALSE
      ),    
      shinyBS::bsTooltip(id = "real_loc",
                         title = "If real_loc = TRUE, the last SNP is used as the length of the chromosome. If real_freq = FALSE and real_loc = TRUE, the initial frequency is taken from q_neutral")
    ),
    
    column(
      4,
      radioButtons(
        "real_freq",
        tags$div(tags$i(HTML("real_freq<br/>")),
                 "Extract allele frequencies for neutral loci from genlight object"),
        choices = list("TRUE" = TRUE,
                       "FALSE" = FALSE),
        selected = FALSE
      ),    
      shinyBS::bsTooltip(id = "real_freq",
                         title = "If real_freq = TRUE and real_loc = FALSE, the locations of the genlight object are randomly choose but following their order in the genlight object")
    ),
    
    column(
      4,
      textInput(
        "chromosome_name", 
        tags$div(tags$i(HTML("chromosome_name<br/>")),
                 "Chromosome name from where to extract location, alllele frequency, recombination map and targets of selection (if provided)"),
        value = "1"
      ),    
      shinyBS::bsTooltip(id = "chromosome_name",
                         title = "")
    )
    
    ),
    
    fluidRow(
    
      column(
        4,
        numericInput(
          "deleterious_factor",
          tags$div(tags$i(HTML("deleterious_factor<br/>")),
                   "Percentage of the number targets from the input file 'targets_of_selection.csv' to use for loci with deleterious alleles"),
          value = 1,
          min = 0
        ),    
        shinyBS::bsTooltip(id = "deleterious_factor",
                  title = "")
      ),
      
      column(
        4,
        numericInput(
          "mutations_factor",
          tags$div(tags$i(HTML("mutations_factor<br/>")),
                   "Percentage of the number targets from the input file 'targets_of_selection.csv' to use for mutations "),
          value = 1,
          min = 0
        ),    
        shinyBS::bsTooltip(id = "mutations_factor",
                           title = "")
      )
    
      ),
    
    hr(),
    
    fluidRow(column(
      12,
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
    
    #SELECTION COEFFICIENT DELETERIOUS
    
    toListen_s_del <- reactive({
      list(input$s_coeff,input$s_distribution_del,isolate(input$loci_deleterious))
    })
    
    observeEvent(toListen_s_del(), {
      
      if(input$s_coeff == TRUE & input$loci_deleterious>0){
        shinyjs::show('s_distribution_del')
      }else{
        shinyjs::hide('s_distribution_del')
      }
      
      if(input$s_distribution_del == "equal" & input$s_coeff == TRUE &
         input$loci_deleterious>0){
        shinyjs::show('s_del')
      }else{
        shinyjs::hide('s_del')
      }
      
      if(input$s_distribution_del == "gamma" & input$s_coeff == TRUE &
         input$loci_deleterious>0){
        shinyjs::show('gamma_scale')
        shinyjs::show('gamma_shape')
      }else{
        shinyjs::hide('gamma_scale')
        shinyjs::hide('gamma_shape')
      }
      
      if(input$s_distribution_del == "log_normal" & input$s_coeff == TRUE &
         input$loci_deleterious>0){
        shinyjs::show('log_mean')
        shinyjs::show('log_sd')
      }else{
        shinyjs::hide('log_mean')
        shinyjs::hide('log_sd')
      }
    
    })
    
    # DOMINANCE DELETERIOUS
    
    toListen_h_del <- reactive({
      list(input$dominance,input$h_distribution_del,isolate(input$loci_deleterious))
    })
    
    observeEvent(toListen_h_del(), {
      
      if(input$dominance == TRUE & input$loci_deleterious>0){
        shinyjs::show('h_distribution_del')
      }else{
        shinyjs::hide('h_distribution_del')
      }
      
      if(input$h_distribution_del == "equal" & input$dominance == TRUE &
         input$loci_deleterious>0){
        shinyjs::show('h_del')
      }else{
        shinyjs::hide('h_del')
      }
      
      if(input$h_distribution_del == "normal" & input$dominance == TRUE &
         input$loci_deleterious>0){
        shinyjs::show('h_mean_del')
        shinyjs::show('h_sd_del')
      }else{
        shinyjs::hide('h_mean_del')
        shinyjs::hide('h_sd_del')
      }
      
      if(input$h_distribution_del == "equation" & input$dominance == TRUE &
         input$loci_deleterious>0){
        shinyjs::show('h_intercept_del')
        shinyjs::show('h_rate_del')
      }else{
        shinyjs::hide('h_intercept_del')
        shinyjs::hide('h_rate_del')
      }
      
    })
    
    #### FREQUENCY DELETERIOUS
    
    toListen_freq_del <- reactive({
      list(input$frequency,input$q_distribution_del,isolate(input$loci_deleterious))
           })

    observeEvent(toListen_freq_del(), {

      if(input$frequency == TRUE & input$loci_deleterious>0){
        shinyjs::show('q_distribution_del')
      }else{
        shinyjs::hide('q_distribution_del')
      }

      if(input$frequency == TRUE & input$q_distribution_del == "equal" &
         input$loci_deleterious>0){
        shinyjs::show('q_del')
      }else{
        shinyjs::hide('q_del')
      }

      if(input$frequency == TRUE && input$q_distribution_del == "equation" &
         input$loci_deleterious>0){
        shinyjs::show('q_equation_del')
      }else{
        shinyjs::hide('q_equation_del')
      }

    })
    
    # SELECTION COEFFICIENT ADVANTAGEOUS
    
    toListen_s_adv<- reactive({
      list(input$s_coeff,input$s_distribution_adv,isolate(input$loci_advantageous))
    })
    
    observeEvent(toListen_s_adv(), {
      
      if(input$s_coeff == TRUE & input$loci_advantageous>0){
        shinyjs::show('s_distribution_adv')
      }else{
        shinyjs::hide('s_distribution_adv')
      }
      
      if(input$s_distribution_adv == "equal" & input$s_coeff == TRUE &
         input$loci_advantageous>0){
        shinyjs::show('s_adv')
      }else{
        shinyjs::hide('s_adv')
      }
      
      if(input$s_distribution_adv == "exponential" & input$s_coeff == TRUE &
         input$loci_advantageous>0){
        shinyjs::show('exp_rate')
      }else{
        shinyjs::hide('exp_rate')
      }
      
    })
    
    # DOMINANCE ADVANTAGEOUS
    
    toListen_h_adv <- reactive({
      list(input$dominance,input$h_distribution_adv,isolate(input$loci_advantageous))
    })
    
    observeEvent(toListen_h_adv(), {
      
      if(input$dominance == TRUE & input$loci_advantageous>0){
        shinyjs::show('h_distribution_adv')
      }else{
        shinyjs::hide('h_distribution_adv')
      }
      
      if(input$h_distribution_adv == "equal" & input$dominance == TRUE &
         input$loci_advantageous>0){
        shinyjs::show('h_adv')
      }else{
        shinyjs::hide('h_adv')
      }
      
      if(input$h_distribution_adv == "normal" & input$dominance == TRUE &
         input$loci_advantageous>0){
        shinyjs::show('h_mean_adv')
        shinyjs::show('h_sd_adv')
      }else{
        shinyjs::hide('h_mean_adv')
        shinyjs::hide('h_sd_adv')
      }
      
      if(input$h_distribution_adv == "equation" & input$dominance == TRUE &
         input$loci_advantageous>0){
        shinyjs::show('h_intercept_adv')
        shinyjs::show('h_rate_adv')
      }else{
        shinyjs::hide('h_intercept_adv')
        shinyjs::hide('h_rate_adv')
      }
      
    })
    
    #### FREQUENCY ADVANTAGEOUS
    
    toListen_freq_adv <- reactive({
      list(input$frequency,input$q_distribution_adv,isolate(input$loci_advantageous))
    })
    
    observeEvent(toListen_freq_adv(), {
      
      if(input$frequency == TRUE & input$loci_advantageous>0){
        shinyjs::show('q_distribution_adv')
      }else{
        shinyjs::hide('q_distribution_adv')
      }
      
      if(input$frequency == TRUE & input$q_distribution_adv == "equal" &
         input$loci_advantageous>0){
        shinyjs::show('q_adv')
      }else{
        shinyjs::hide('q_adv')
      }
      
      if(input$frequency == TRUE && input$q_distribution_adv == "equation" &
         input$loci_advantageous>0){
        shinyjs::show('q_equation_adv')
      }else{
        shinyjs::hide('q_equation_adv')
      }
      
    })
    
    observeEvent(input$mutation, {
      
      if(input$mutation == TRUE){
        shinyjs::show('loci_mut_neu')
        shinyjs::show('loci_mut_del')
        shinyjs::show('loci_mut_adv')
      }else{
        shinyjs::hide('loci_mut_neu')
        shinyjs::hide('loci_mut_del')
        shinyjs::hide('loci_mut_adv')   
        }
      
    })
    
    observeEvent(input$real_v, {
      
      if(input$real_v == TRUE){
        shinyjs::show('real_loc')
        shinyjs::show('real_freq')
        shinyjs::show('chromosome_name')
        shinyjs::show('deleterious_factor')
        shinyjs::show('mutations_factor')
      }else{
        shinyjs::hide('real_loc')
        shinyjs::hide('real_freq')
        shinyjs::hide('chromosome_name')
        shinyjs::hide('deleterious_factor')
        shinyjs::hide('mutations_factor')
      }
      
    })

    observeEvent(input$close, {
      
      ref_vars_temp <- as.data.frame(cbind(
        c("chunk_number",
          "chunk_bp",
          "chunk_cM",
          "chunk_neutral_loci",
          "loci_deleterious",
          "loci_mut_del",
          "q_neutral",
          "q_distribution_del",
          "q_del",
          "q_equation_del",
          "s_distribution_del",
          "s_del",
          "exp_rate",
          "gamma_scale",
          "gamma_shape",
          "log_mean",
          "log_sd",
          "h_distribution_del",
          "h_del",
          "h_mean_del",
          "h_sd_del",
          "h_intercept_del",
          "h_rate_del",
          "real_loc",
          "real_freq",
          "chromosome_name",
          "deleterious_factor",
          "mutations_factor",
          "q_adv",
          "loci_advantageous",
          "h_adv",
          "s_adv",
          "s_distribution_adv",
          "h_distribution_adv",
          "h_mean_adv",
          "h_sd_adv",
          "h_intercept_adv",
          "h_rate_adv",
          "q_distribution_adv",
          "q_equation_adv",
          "loci_mut_adv",
          "loci_mut_neu"
        ),
        c(input$chunk_number,
          input$chunk_bp,
          input$chunk_cM,
          input$chunk_neutral_loci,
          input$loci_deleterious,
          input$loci_mut_del,
          input$q_neutral,
          input$q_distribution_del,
          input$q_del,
          input$q_equation_del,
          input$s_distribution_del,
          input$s_del,
          input$exp_rate,
          input$gamma_scale,
          input$gamma_shape,
          input$log_mean,
          input$log_sd,
          input$h_distribution_del,
          input$h_del,
          input$h_mean_del,
          input$h_sd_del,
          input$h_intercept_del,
          input$h_rate_del,
          input$real_loc,
          input$real_freq,
          input$chromosome_name,
          input$deleterious_factor,
          input$mutations_factor,
          input$q_adv,
          input$loci_advantageous,
          input$h_adv,
          input$s_adv,
          input$s_distribution_adv,
          input$h_distribution_adv,
          input$h_mean_adv,
          input$h_sd_adv,
          input$h_intercept_adv,
          input$h_rate_adv,
          input$q_distribution_adv,
          input$q_equation_adv,
          input$loci_mut_adv,
          input$loci_mut_neu
        )
      ))
      
      colnames(ref_vars_temp) <- c("variable", "value")
      
      stopApp(ref_vars_temp)
      
    })
    
  }
  
  runApp(shinyApp(ui, server))

}

###############################################################################
######### SHINY APP FOR THE VARIABLES OF THE SIMULATIONS ######################
###############################################################################
#' @name interactive_sim_run
#' @title Shiny app for the input of the simulations variables
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

interactive_sim_run <- function() {
  
  pkg <- "shinyBS"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }
  
  pkg <- "shinythemes"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }
  
  pkg <- "shinyjs"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }
  
  pkg <- "shinyWidgets"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }
  
  ui <- fluidPage(
    
    shinyjs::useShinyjs(),
    
    theme = shinythemes::shinytheme("superhero"),
    
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

    h3(strong("Basic variables")),
    
    fluidRow(
      
      column(
        4,
        numericInput(
          "number_pops_phase2",
          tags$div(tags$i(HTML(
            "number_pops_phase2<br/>"
          )),
          "Number of populations"),
          value = 2,
          min = 0
        ),
        shinyBS::bsTooltip(id = "number_pops_phase2",
                  title = "")
      ),
      
      column(
        4,
        numericInput(
          "gen_number_phase2",
          tags$div(tags$i(HTML(
            "gen_number_phase2<br/>"
          )),
          "Number of generations"),
          value = 10,
          min = 0
        ),
        shinyBS::bsTooltip(id = "gen_number_phase2",
                           title = "Generations are not overlapping (i.e., parents and offspring do not coexist)")
      )
      
    ),
    
    hr(),
    
    h4("Genetic drift"),
    
    fluidRow(
      
      column(
        4,
        textInput(
          "population_size_phase2",
          tags$div(
            tags$i(HTML("population_size_phase2<br/>")),
            "Census population size of each population (must be even and space delimited)"),
          value = c("10", "10")
        ),
        textOutput("equal"),
        tags$head(tags$style("#equal{color: red}")),
        shinyBS::bsTooltip(id = "population_size_phase2",
                  title = "Population size remain constant across generations, i.e., in each generation the entire population is replaced by sampling the same number of offspring as there were parents in the previous generation")
      ),
      
      column(
        4,
        numericInput(
          "variance_offspring_phase2",
          tags$div(
            tags$i(HTML("variance_offspring_phase2<br/>")),
            "Coefficient that determines the variance in the number of offspring per mating"
          ),
          value = 1000000,
          min = 0
        ),
        shinyBS::bsTooltip(id = "variance_offspring_phase2",
                           title = "This variable controls the variance of the negative binomial distribution that is used to determine the number of offspring that each mating pair produces. If the user requires that the Ne/Nc ratio to be equal to 1, a large enough value in the parameter variance_offspring (e.g., > 1,000) should be used. If the user requires that the Ne/Nc ratio to be different from 1, a calibration process can be performed for this end, see tutorial.")
      ),
      
      column(
        4,
        numericInput(
          "number_offspring_phase2",
          tags$div(tags$i(HTML(
            "number_offspring_phase2<br/>"
          )),
          "Mean number offspring per mating"),
          value = 10,
          min = 0
        ),
        shinyBS::bsTooltip(id = "number_offspring_phase2",
                           title = "This variable controls the mean of the negative binomial distribution. This variable allows to control the number of offspring per mating which is convenient when there is a need that each pair of parents produce enough offspring in each generation for the population not to become extinct. However, in the simulations the mean number of offspring per mating each generation is effectively equal to two for two reasons: a) the population size remains constant from generation to generation, this means that on average, across generations and replicates, two offspring per pair of parents are selected to become the parents of the next generation; and b) there is no variance in reproductive success (whether or not an individual gets to reproduce at all) because all individuals reproduce once")
      )

    ),
    
    hr(),
    
    fluidRow(
      
      column(
        4,
        shinyWidgets::awesomeCheckbox(inputId = "dispersal_phase2",
                        label = h4("Dispersal"),
                        value = FALSE,
                        status = "success"),
        shinyBS::bsTooltip(id = "dispersal_phase2",
                           title = "Dispersal between populations is symmetric and constant across generations. Dispersal between populations can be further paramaterised using the function gl.sim.create_dispersal")
      )
      
    ),
    
    fluidRow(
      
      column(
        4,
        radioButtons(
          "dispersal_type_phase2",
          tags$div(tags$i(HTML(
            "dispersal_type_phase2<br/>"
          )),
          "Type of dispersal"),
          choices = list(
            "All connected" = "all_connected",
            "Circle" = "circle",
            "Line" = "line"
          ),
          selected = "all_connected"
        ),
        shinyBS::bsTooltip(id = "dispersal_type_phase2",
                  title = "")
      ),
      
      column(
        4,
        numericInput(
          "number_transfers_phase2",
          tags$div(
            tags$i(HTML("number_transfers_phase2<br/>")),
            "Number of dispersing individuals in each dispersal event"),
          value = 1,
          min = 0
        ),
        shinyBS::bsTooltip(id = "number_transfers_phase2",
                  title = "")
      ),
    
      column(
        4,
        numericInput(
          "transfer_each_gen_phase2",
          tags$div(
            tags$i(HTML("transfer_each_gen_phase2<br/>")),
            "Interval of number of generations in which a dispersal event occurs"),
          value = 1,
          min = 0
        ),
        shinyBS::bsTooltip(id = "transfer_each_gen_phase2",
                  title = "")
      )
      
      ),
    
    hr(),
    
    fluidRow(
      
      column(
        4,
        shinyWidgets::awesomeCheckbox(inputId = "selection_phase2",
                        label = h4("Selection"),
                        value = FALSE,
                        status = "success"),
        shinyBS::bsTooltip(id = "selection_phase2",
                           title = "Selection is directional (replaces one allele by another) and multiplicative")
      )
      
    ),
    
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
        shinyBS::bsTooltip(id = "natural_selection_model",
                           title = "In the relative model, also called soft selection or density dependent selection, the fitness of each individual is dependent on the fitness of other individuals in the population. In the absolute model, also called hard selection or density independent selection is based on genetic load which measures the fraction of the population that fails to survive or reproduce. See tutorial.")
      ),
      
      column(
        4,
        numericInput(
          "genetic_load",
          tags$div(tags$i(HTML("genetic_load<br/>")),
                   "Approximation of the genetic load"),
          value = 0.8,
          min = 0
        ),    
        shinyBS::bsTooltip(id = "genetic_load",
                           title = "This variable is used in the absolute fitness model. See tutorial.")
      )
      
    ),
    
    fluidRow(
      
      column(
        4,
        textInput(
          "local_adap",
          tags$div(
            tags$i(HTML("local_adap<br/>")),
            "Populations in which local adaptation occurs (must be space delimited)"),
          value = NULL
        ),
        shinyBS::bsTooltip(id = "local_adap",
                           title = "Selection coefficients are set to 0 in loci and mutations that are advantageous in populations in which local adpatation does not occur")
      ),
      
      column(
        4,
        textInput(
          "clinal_adap",
          tags$div(
            tags$i(HTML("clinal_adap<br/>")),
            "Starting and ending populations in which clinal adaptation occurs (two values must be provided and be space delimited)"),
          value = NULL
        ),
        shinyBS::bsTooltip(id = "clinal_adap",
                           title = "Selection coefficient is decreased in loci and mutations that are advantageous by population along the cline")
      ),
      
      column(
        4,
        numericInput(
          "clinal_strength",
          tags$div(tags$i(HTML("clinal_strength<br/>")),
                   "Percentage of the decrease of the selection coefficient along the cline"),
          value = 10,
          min = 0
        ),    
        shinyBS::bsTooltip(id = "clinal_strength",
                           title = "Selection coefficient is decreased in loci and mutations that are advantageous by population along the cline")
      )
      
    ),
    
    hr(),
    
    fluidRow(
      
      column(
        4,
        shinyWidgets::awesomeCheckbox(inputId = "mutation",
                        label = h4("Mutation"),
                        value = FALSE,
                        status = "success"),
        shinyBS::bsTooltip(id = "mutation",
                           title = "Pending")
      ),
      
      column(
        4,
        numericInput(
          "mut_rate",
          tags$div(tags$i(HTML("mut_rate<br/>")),
                   "Mutation rate per genome per generation"),
          value = 0.01,
          min = 0
        ),    
        shinyBS::bsTooltip(id = "mut_rate",
                           title = "Real values range from 1-2 per genome per generation; see Keightley, 2012")
      )
      
    ),

    hr(),
    
    fluidRow(
      
      column(
        4,
        shinyWidgets::awesomeCheckbox(inputId = "recombination",
                        label = h4("Recombination"),
                        value = TRUE,
                        status = "success"),
        shinyBS::bsTooltip(id = "recombination",
                           title = "")
      ),
      
      column(
        4,
        shinyWidgets::awesomeCheckbox(inputId = "recombination_males",
                        label = tags$div(tags$i(HTML(
                          "recombination_males<br/>"
                        )),
                        "Recombination occurs in males"),
                        value = TRUE,
                        status = "success"),
        shinyBS::bsTooltip(id = "recombination_males",
                           title = "")
      )
    
    ),
    
    hr(),
    
    h5(
      em(
        "Simulations can have 2 phases (phase 1 and phase 2). Variable values are constant across generations in each phase but variable values can be different in each phase. The default is to run just phase 2, select phase1 to run it"
      )
    ),
    
    fluidRow(
      
      column(
        4,
        shinyWidgets::awesomeCheckbox(inputId = "phase1",
                        label = h4("Phase 1"),
                        value = FALSE,
                        status = "success"),
        shinyBS::bsTooltip(id = "phase1",
                           title = "")
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
          "Number of populations"),
          value = 2,
          min = 0
        ),
        shinyBS::bsTooltip(id = "number_pops_phase1",
                           title = "Must be the same as number of populations in phase 2")
      ),
      
      column(
        4,
        numericInput(
          "gen_number_phase1",
          tags$div(tags$i(HTML(
            "gen_number_phase1<br/>"
          )),
          "Number of generations"),
          value = 10,
          min = 0
        ),
        shinyBS::bsTooltip(id = "gen_number_phase1",
                           title = "Generations are not overlapping (i.e., parents and offspring do not coexist)")
      ),
      
      column(
        4,
        radioButtons(
          "same_line",
          tags$div(tags$i(HTML("same_line<br/>")),
                   "Sample phase 2 populations from the same phase 1 population or from different phase 1 populations"),
          choices = list(
            "Same population" = TRUE,
            "Different populations" = FALSE
          ),
          selected = FALSE
        ),    
        shinyBS::bsTooltip(id = "same_line",
                           title = "")
      )
      
      ),
    
    fluidRow(
      
      column(
        4,
        textInput(
          "population_size_phase1",
          tags$div(
            tags$i(HTML("population_size_phase1<br/>")),
            "Census population size of each population (must be even and space delimited)"
          ),
          value = c("20", "20")
        ),
        textOutput("equal_phase1"),
        tags$head(tags$style("#equal_phase1{color: red}")),
        shinyBS::bsTooltip(id = "population_size_phase1",
                           title = "")
      ),
      
      column(
        4,
        numericInput(
          "variance_offspring_phase1",
          tags$div(
            tags$i(HTML("variance_offspring_phase1<br/>")),
            "Coefficient that determines the variance in the number of offspring per mating"
          ),
          value = 1000000,
          min = 0
        ),
        shinyBS::bsTooltip(id = "variance_offspring_phase1",
                           title = "This variable controls the variance of the negative binomial distribution that is used to determine the number of offspring that each mating pair produces")
      ),
      
      column(
        4,
        numericInput(
          "number_offspring_phase1",
          tags$div(tags$i(HTML(
            "number_offspring_phase1<br/>"
          )),
          "Mean number offspring per mating"),
          value = 10,
          min = 0
        ),
        shinyBS::bsTooltip(id = "number_offspring_phase1",
                           title = "This variable controls the mean of the negative binomial distribution. This variable allows to control the number of offspring per mating which is convenient when there is a need that each pair of parents produce enough offspring in each generation for the population not to become extinct. However, in the simulations the mean number of offspring per mating each generation is effectively equal to two for two reasons: a) the population size remains constant from generation to generation, this means that  on average, across generations and replicates, two offspring per pair of parents are selected to become the parents of the next generation; and b) there is no variance in reproductive success (whether or not an individual gets to reproduce at all) because all individuals reproduce once")
      )
      
    ),
    
    fluidRow(
      
      column(
        4,
        shinyWidgets::awesomeCheckbox(inputId = "dispersal_phase1",
                        label = h5("Dispersal"),
                        value = FALSE,
                        status = "success"),
        shinyBS::bsTooltip(id = "dispersal_phase1",
                           title = "Pending")
      )
      
    ),
    
    fluidRow(
      
      column(
        4,
        radioButtons(
          "dispersal_type_phase1",
          tags$div(tags$i(HTML(
            "dispersal_type_phase1<br/>"
          )),
          "Type of dispersal"),
          choices = list(
            "All connected" = "all_connected",
            "Circle" = "circle",
            "Line" = "line"
          ),
          selected = "all_connected"
        ),
        shinyBS::bsTooltip(id = "dispersal_type_phase1",
                           title = "Information pending")
      ),
      
      column(
        4,
        numericInput(
          "number_transfers_phase1",
          tags$div(
            tags$i(HTML("number_transfers_phase1<br/>")),
            "Number of dispersing individuals in each dispersal event"
          ),
          value = 1,
          min = 0
        ),
        shinyBS::bsTooltip(id = "number_transfers_phase1",
                           title = "Information pending")
      ),
      
      column(
        4,
        numericInput(
          "transfer_each_gen_phase1",
          tags$div(
            tags$i(HTML("transfer_each_gen_phase1<br/>")),
            "Interval of number of generations in which a dispersal event occurs"
          ),
          value = 1,
          min = 0
        ),
        shinyBS::bsTooltip(id = "transfer_each_gen_phase1",
                           title = "")
      )
      
    ),
    
    fluidRow(
      
      column(
        4,
        shinyWidgets::awesomeCheckbox(inputId = "selection_phase1",
                        label = h5("Selection"),
                        value = FALSE,
                        status = "success"),
        shinyBS::bsTooltip(id = "selection_phase1",
                           title = "Selection is directional (replaces one allele by another) and multiplicative")
      )
      
    ),
    
    hr(),
    
    fluidRow(
      
      column(
        4,
        shinyWidgets::awesomeCheckbox(inputId = "real_dataset",
                        label = h4("Real information"),
                        value = FALSE,
                        status = "success")
      )
      
    ),
    
    fluidRow(
      
      column(
        4,
        radioButtons(
          "real_pops",
          tags$div(tags$i(HTML("real_pops<br/>")),
                   "Extract number of populations from genlight object"),
          choices = list("TRUE" = TRUE,
                         "FALSE" = FALSE),
          selected = FALSE
        ),    
        shinyBS::bsTooltip(id = "real_pops",
                           title = "")
      ),
      
      column(
        4,
        radioButtons(
          "real_pop_size",
          tags$div(tags$i(HTML("real_pop_size<br/>")),
                   "Extract census population sizes from genlight object"),
          choices = list("TRUE" = TRUE,
                         "FALSE" = FALSE),
          selected = FALSE
        ),    
        shinyBS::bsTooltip(id = "real_pop_size",
                           title = "Odd population sizes are converted to even")
      ),
      
      column(
        4,
        radioButtons(
          "real_loc",
          tags$div(tags$i(HTML("real_loc<br/>")),
                   "Extract location for neutral loci from genlight object"),
          choices = list("TRUE" = TRUE,
                         "FALSE" = FALSE),
          selected = FALSE
        ),    
        shinyBS::bsTooltip(id = "real_loc",
                           title = "")
      )
      
    ),
      
      fluidRow(
      
      column(
        4,
        radioButtons(
          "real_freq",
          tags$div(tags$i(HTML("real_freq<br/>")),
                   "Extract allele frequencies for neutral loci from genlight object"),
          choices = list("TRUE" = TRUE,
                         "FALSE" = FALSE),
          selected = FALSE
        ),    
        shinyBS::bsTooltip(id = "real_freq",
                           title = "")
      ),
      
      column(
        4,
        textInput(
          "chromosome_name",
          tags$div(tags$i(HTML("chromosome_name<br/>")),
                   "Chromosome name from where to extract allele frequencies"),
          value = "1"
        ),    
        shinyBS::bsTooltip(id = "chromosome_name",
                           title = "")
      )
      
    ),
    
    hr(),
    
    fluidRow(column(
      12,
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
    
    output$equal <- renderText({
      req(input$population_size_phase2)
      if (length(unlist(strsplit(input$population_size_phase2, " "))) != input$number_pops_phase2) {
        validate("Number of population sizes is not equal to number of populations")
      }
    })
    
    output$equal_phase1 <- renderText({
      req(input$population_size_phase1)
      if (length(unlist(strsplit(input$population_size_phase1, " "))) != input$number_pops_phase1) {
        validate("Number of population sizes is not equal to number of populations")
      }
    })
    
    observeEvent(input$phase1, {

      shinyjs::toggleElement(
        id = "same_line",
        condition = input$phase1 == TRUE
        )
      
      shinyjs::toggleElement(
        id = "number_pops_phase1",
        condition = input$phase1 == TRUE
      )
      
      shinyjs::toggleElement(
        id = "population_size_phase1",
        condition = input$phase1 == TRUE
      )

      shinyjs::toggleElement(
        id = "gen_number_phase1",
        condition = input$phase1 == TRUE
      )

      shinyjs::toggleElement(
        id = "dispersal_phase1",
        condition = input$phase1 == TRUE
      )

      shinyjs::toggleElement(
        id = "dispersal_type_phase1",
        condition = input$phase1 == TRUE
      )

      shinyjs::toggleElement(
        id = "number_transfers_phase1",
        condition = input$phase1 == TRUE
      )

      shinyjs::toggleElement(
        id = "transfer_each_gen_phase1",
        condition = input$phase1 == TRUE
      )

      shinyjs::toggleElement(
        id = "variance_offspring_phase1",
        condition = input$phase1 == TRUE
      )

      shinyjs::toggleElement(
        id = "number_offspring_phase1",
        condition = input$phase1 == TRUE
      )

      shinyjs::toggleElement(
        id = "selection_phase1",
        condition = input$phase1 == TRUE
      )
      
    })
    
    observeEvent(input$real_dataset, {
      
      shinyjs::toggleElement(
        id = "real_pops",
        condition = input$real_dataset == TRUE
      )
      
      shinyjs::toggleElement(
        id = "real_pop_size",
        condition = input$real_dataset == TRUE
      )
      
      shinyjs::toggleElement(
        id = "real_freq",
        condition = input$real_dataset == TRUE
      )
      
      shinyjs::toggleElement(
        id = "chromosome_name",
        condition = input$real_dataset == TRUE
      )
      
      shinyjs::toggleElement(
        id = "real_loc",
        condition = input$real_dataset == TRUE
      )
      
    })
    
    observeEvent(input$dispersal_phase2, {
      
      shinyjs::toggleElement(
        id = "dispersal_type_phase2",
        condition = input$dispersal_phase2 == TRUE
      )
      
      shinyjs::toggleElement(
        id = "number_transfers_phase2",
        condition = input$dispersal_phase2 == TRUE
      )
      
      shinyjs::toggleElement(
        id = "transfer_each_gen_phase2",
        condition = input$dispersal_phase2 == TRUE
      )
      
    })
    
    observeEvent(input$selection_phase2, {
      
      shinyjs::toggleElement(
        id = "natural_selection_model",
        condition = input$selection_phase2 == TRUE
      )
      
      shinyjs::toggleElement(
        id = "genetic_load",
        condition = input$selection_phase2 == TRUE
      )
      
      shinyjs::toggleElement(
        id = "local_adap",
        condition = input$selection_phase2 == TRUE
      )
      
      shinyjs::toggleElement(
        id = "clinal_adap",
        condition = input$selection_phase2 == TRUE
      )
      
      shinyjs::toggleElement(
        id = "clinal_strength",
        condition = input$selection_phase2 == TRUE
      )
      
    })
    
    observeEvent(input$mutation, {
      
      shinyjs::toggleElement(
        id = "mut_rate",
        condition = input$mutation == TRUE
      )
      
    })
    
    observeEvent(input$close, {
      
      sim_vars_temp <- as.data.frame(cbind(
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
          "phase1",
          "same_line",
          "real_dataset",
          "real_freq",
          "real_pops",
          "real_pop_size",
          "recombination",
          "recombination_males",
          "genetic_load",
          "natural_selection_model",
          "chromosome_name",
          "mutation",
          "mut_rate",
          "real_loc",
          "clinal_strength",
          "clinal_adap",
          "local_adap"
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
          input$phase1,
          input$same_line,
          input$real_dataset,
          input$real_freq,
          input$real_pops,
          input$real_pop_size,
          input$recombination,
          input$recombination_males,
          input$genetic_load,
          input$natural_selection_model,
          input$chromosome_name,
          input$mutation,
          input$mut_rate,
          input$real_loc,
          input$clinal_strength,
          input$clinal_adap,
          input$local_adap
        )))
      
      colnames(sim_vars_temp) <- c("variable","value")
      
      stopApp(sim_vars_temp)
      
    }
  )
    
  }

  runApp(shinyApp(ui, server))
  
}
