#' @name gl.diagnostics.sim
#' @title Comparing simulations against theoretical expectations
#' @param x [default ]
#' @param iteration [default ]
#' @param Ne [required]
#' @param pop_he [default ]
#' @param pops_fst [default ]
#' @param plot_colors Vector with two color names for the significant and
#' not-significant loci [default two_colors_contrast].
#' @param plot_theme User specified theme [default theme_dartR()].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session
#' temporary directory (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log ; 3, progress and results summary; 5, full report
#' [default NULL, unless specified using gl.set.verbosity].
#' @export

gl.diagnostics.sim <- function(x,
                               Ne,
                               iteration=1,
                               pop_he = 1 ,
                               pops_fst= c(1,2),
                               plot_colors = two_colors_contrast,
                               plot_theme = theme_dartR(),
                               save2tmp = FALSE,
                               verbose = NULL){
  
  Ne_hold <- Ne
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jody",
                   verbosity = verbose)
  
  # CHECK DATATYPE
  # datatype <- utils.check.datatype(x, verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
  # Check that call rate is up to date and recalculate if necessary
  
  # if (!x@other$loc.metrics.flags$CallRate) {
  # x <- utils.recalc.callrate(x, verbose = 0)
  # }
  
  # DO THE JOB
  
  lab<-gen<-He<-value<-variable<-fst_obs<- NULL
  
  x <- x[[iteration]]
  sep_pops <- lapply(x,seppop)

  ####################### He #######################
  
  he_pop_temp <- lapply(sep_pops,"[",pop_he)
  
  he_pop <- lapply(he_pop_temp,function(y){
    temp <- lapply(y,gl.He)
    temp <- unname(unlist(lapply(temp,mean)))
    temp <- as.data.frame(t(temp))
    return(temp)
  })
  
  he_pop <- rbindlist(he_pop)

  Ne <- gsub('\"', "", Ne , fixed = TRUE)
  Ne <- stringr::str_split(Ne,pattern=" ")[[1]]
  Ne <- as.numeric(Ne)[pop_he]

  Ne <- seq(Ne[1],Ne[1]*2, Ne[1]/4)
  
  rate_loss <-  1 - (1 / (2 * Ne))
  generations_sim <- unlist(lapply(x,function(y){y$other$sim.vars$generation}))
  
  expected_het <- lapply(1:length(Ne),function(x){sapply(mean(unlist(he_pop[1,])), "*", (rate_loss[x] ^ generations_sim))})
  
  expected_het_2 <- lapply(1:length(Ne),function(x){
    expected_het_temp <- as.data.frame(expected_het[[x]])
    expected_het_temp$Ne <- as.character(Ne[x])
    expected_het_temp$gen <- generations_sim
    colnames(expected_het_temp) <- c("He","Ne","gen")
    return(expected_het_temp)
    
  })
  
  colors_plot_ne <- discrete_palette(length(expected_het_2))
  
  for(i in 1:length(colors_plot_ne)){
    expected_het_2[[i]]$col <- colors_plot_ne[i]
  }
  
  expected_het_3 <- Reduce(rbind,expected_het_2)
  expected_het_3$lab <- paste("Ne ",expected_het_3$Ne)

  he_pop$gen <- generations_sim
  colnames(he_pop) <- c(paste0("pop",pop_he),"gen")
  he_pop <- reshape2::melt(he_pop,id="gen")
  
  expected_het_3$lab <- factor(expected_het_3$lab, levels = unique(expected_het_3$lab))
  
  p1 <- ggplot(data=expected_het_3, aes(x=gen,y=He,color=lab)) +
    geom_line(size=0.75,linetype = "dashed") +
    geom_line(data=he_pop,aes(x=gen,y=value,color=variable),size=1.5) +
    labs(x="Generations",
         y="He", 
         title=paste("Rate of loss of heterozygosity\nacross generations population",
                     paste(pop_he,collapse = " ") ))+ 
    plot_theme +
    theme(legend.title=element_blank())
  
  ####################### FST #######################
  
  fst_pop_temp <- lapply(sep_pops,"[",pops_fst)
  
  fst_pop <- lapply(fst_pop_temp,function(y){
    merge_pop <- Reduce(rbind,y)
    return(merge_pop)
  })
  
  fst_pop_hier <- lapply(fst_pop,function(y){
    temp <- hierfstat::genind2hierfstat(gl2gi(y, verbose = 0))
    return(temp)
  })
  
  fst_gen <- lapply(fst_pop_hier,hierfstat::pairwise.neifst)
  
  fst_gen <- unlist(unname(lapply(fst_gen,function(y){
    y[lower.tri(y)]
  })))
  
  Ne_fst <- gsub('\"', "", Ne_hold, fixed = TRUE)
  Ne_fst <- stringr::str_split(Ne_fst,pattern=" ")[[1]]
  Ne_fst <- as.numeric(Ne_fst)[pops_fst]
  
  number_transfers <- x[[1]]$other$sim.vars$number_transfers_phase2
  number_transfers <- gsub('\"', "", number_transfers , fixed = TRUE)
  number_transfers <- stringr::str_split(number_transfers,pattern=" ")[[1]]
  number_transfers <- as.numeric(number_transfers)[pops_fst]
  
  transfer_each_gen <- x[[1]]$other$sim.vars$transfer_each_gen_phase2
  transfer_each_gen <- gsub('\"', "", transfer_each_gen , fixed = TRUE)
  transfer_each_gen <- stringr::str_split(transfer_each_gen,pattern=" ")[[1]]
  transfer_each_gen <- as.numeric(transfer_each_gen)[pops_fst]
 
  population_size <- x[[1]]$other$sim.vars$population_size_phase2
  population_size <- gsub("'", "", population_size , fixed = TRUE)
  population_size <- as.numeric(unlist(strsplit(population_size, " ")))[pops_fst]
  
  dispersal_rate <- (number_transfers / transfer_each_gen) / (population_size)
  
  Fst_expected <- 1 / ((4 * Ne_fst * dispersal_rate) * ((2 / (2 - 1)) ^ 2) + 1)
  
  fst_equilibrium <- (log(1/2) / log( (1- dispersal_rate)^2 * (1-(1/(2*Ne[1]))) )) * 2
  
  generations_fst <- data.frame("gen"=generations_sim,"fst_obs"=fst_gen)

  p2 <- ggplot(generations_fst) +
    geom_line(aes(x=gen,y=fst_obs,color="brown"),size=1) +
    geom_hline(aes(yintercept = Fst_expected[1],color="Fst expected"),size=1)+
    geom_vline(aes(xintercept = fst_equilibrium[1],color="Fst equilibrium"),size=1)+
    labs(x="Generations", y="Fst", title=paste("Fst between populations:",paste(pops_fst,collapse = " ")))+ 
    scale_color_manual(values = c(plot_colors[1],"blue","chartreuse4"),labels=c("Fst observed", "Fst equilibrium" ,"Fst expected")) +
    plot_theme +
    theme(legend.title=element_blank())
  
  # PRINTING OUTPUTS
  # using package patchwork
  p3 <- (p1 / p2)
  print(p3)
  
  # print(hwe_summary, row.names = FALSE)
  
  # SAVE INTERMEDIATES TO TEMPDIR
  
  # # creating temp file names
  # if (save2tmp) {
  #   temp_plot <- tempfile(pattern = "Plot_")
  #   match_call <-
  #     paste0(names(match.call()),
  #            "_",
  #            as.character(match.call()),
  #            collapse = "_")
  #   # saving to tempdir
  #   saveRDS(list(match_call, p3), file = temp_plot)
  #   
  #   if (verbose >= 2) {
  #     cat(report("  Saving the ggplot to session tempfile\n"))
  #   }
  #   
  #   temp_table <- tempfile(pattern = "Table_")
  #   saveRDS(list(match_call, hwe_summary), file = temp_table)
  #   
  #   if (verbose >= 2) {
  #     cat(report("  Saving tabulation to session tempfile\n"))
  #     cat(
  #       report(
  #         "  NOTE: Retrieve output files from tempdir using gl.list.reports() and gl.print.reports()\n"
  #       )
  #     )
  #   }
  # }
  # 
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # RETURN
  
  # return(invisible(hwe_summary))
  

  
}
