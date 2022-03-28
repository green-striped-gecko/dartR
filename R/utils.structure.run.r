#' Utility function to run strcture
#' 
#' These functions were copied from package strataG, which is no longer on CRAN (maintained by Eric Archer)
#' @export
#' @author Bernd Gruber (bugs? Post to
#'  \url{https://groups.google.com/d/forum/dartr}); original implementation of
#'  Eric Archer \url{https://github.com/EricArcher/strataG}
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param k.range vector of values to for \code{maxpop} in multiple runs. If set
#'   to \code{NULL}, a single STRUCTURE run is conducted with \code{maxpops}
#'   groups. If specified, do not also specify \code{maxpops}.
#' @param num.k.rep number of replicates for each value in \code{k.range}.
#' @param label label to use for input and output files
#' @param delete.files logical. Delete all files when STRUCTURE is finished?
#' @param exec name of executable for STRUCTURE. Defaults to "structure".
#' @param ... arguments to be passed to \code{structureWrite}.
#'
#' @return \describe{ \item{\code{structureRun}}{a list where each element is a
#' list with results from \code{structureRead} and a vector of the filenames
#' used} \item{\code{structureWrite}}{a vector of the filenames used by
#' STRUCTURE} \item{\code{structureRead}}{a list containing: \describe{
#' \item{\code{summary}}{new locus name, which is a combination of loci in
#' group} \item{\code{q.mat}}{data.frame of assignment probabilities for each
#' id} \item{\code{prior.anc}}{list of prior ancestry estimates for each
#' individual where population priors were used} \item{\code{files}}{vector of
#' input and output files used by STRUCTURE} \item{\code{label}}{label for the
#' run} } } }

utils.structure.run <- 
  function (g, k.range = NULL, num.k.rep = 1, label = NULL, delete.files = TRUE, exec = "structure", ...) 
  {
################################################################
    .structureParseQmat <-    function (q.mat.txt, pops) 
    {
      q.mat.txt <- sub("[*]+", "", q.mat.txt)
      q.mat.txt <- sub("[(]", "", q.mat.txt)
      q.mat.txt <- sub("[)]", "", q.mat.txt)
      q.mat.txt <- sub("[|][ ]+$", "", q.mat.txt)
      cols1to4 <- c("row", "id", "pct.miss", 
                    "orig.pop")
      strsplit(q.mat.txt, " ") %>% purrr::map(function(q) {
        q <- q[!q %in% c("", " ", ":")] %>% 
          as.character() %>% rbind() %>% as.data.frame(stringsAsFactors = FALSE)
        stats::setNames(q, c(cols1to4, paste("Group", 1:(ncol(q) - 
                                                           4), sep = ".")))
      }) %>% dplyr::bind_rows() %>% dplyr::mutate_at(dplyr::vars("row", 
                                                                 "pct.miss", "orig.pop", dplyr::starts_with("Group.")), 
                                                     as.numeric) %>% dplyr::mutate(orig.pop = if (!is.null(pops)) {
                                                       pops[.data$orig.pop]
                                                     }
                                                     else {
                                                       .data$orig.pop
                                                     })
    }
    
    
################################################################
    .alleles2integer <- function (g, min.val = 0) 
    {
      g@data %>% dplyr::group_by(.data$locus) %>% dplyr::mutate(allele = min.val - 
                                                                  1 + as.integer(factor(.data$allele))) %>% dplyr::ungroup()
    }    
################################################################

################################################################
    .stackedAlleles <-    function (g, alleles2integer = FALSE, na.val = NULL, ...) 
    {
      x <- if (alleles2integer) 
        .alleles2integer(g, ...)
      else g@data
      if (!is.null(na.val)) 
        x$allele[is.na(x$allele)] <- na.val
      x %>% dplyr::arrange(.data$id, .data$locus) %>% dplyr::mutate(a = rep(1:g@ploidy, 
                                                                            dplyr::n()/g@ploidy)) %>% tidyr::spread(.data$locus, 
                                                                                                                        .data$allele) %>% dplyr::rename(allele = "a") %>% 
        dplyr::select(.data$id, .data$stratum, .data$allele, 
                      dplyr::everything())
    }    
        
####################################################    
    .getFileLabel <- function (g, label = NULL) 
    {
      desc <- g@description
      label <- if (!is.null(label)) {
        label
      }
      else if (!is.null(desc)) {
        desc
      }
      else "strataG.gtypes"
      gsub("[[:punct:]]", ".", label)
    }
####################################################
    structureWrite <-   function (g, label = NULL, maxpops = 1:(dplyr::n_distinct(g@data$stratum)), burnin = 1000, 
              numreps = 1000, noadmix = TRUE, freqscorr = FALSE, randomize = TRUE, 
              seed = 0, pop.prior = NULL, locpriorinit = 1, maxlocprior = 20, 
              gensback = 2, migrprior = 0.05, pfrompopflagonly = TRUE, 
              popflag = NULL, inferalpha = FALSE, alpha = 1, unifprioralpha = TRUE, 
              alphamax = 20, alphapriora = 0.05, alphapriorb = 0.001) 
    {
     
      if (!is.null(pop.prior)) {
        if (!pop.prior %in% c("locprior", "usepopinfo")) {
          stop("'pop.prior' must be 'locprior' or 'usepopinfo'.")
        }
      }
      if (is.null(popflag)) 
        popflag <- rep(1, length(unique(g@data$id)))
      if (length(popflag) != length(unique(g@data$id))) {
        stop("'popflag' should be the same length as the number of individuals in 'g'.")
      }
      if (!all(popflag %in% c(0, 1))) {
        stop("all values in 'popflag' must be 0 or 1.")
      }
      if (is.null(names(popflag))) 
        names(popflag) <- unique(g@data$id)
      in.file <- ifelse(is.null(label), "data", paste(label, 
                                                      "data", sep = "_"))
      out.file <- ifelse(is.null(label), "out", paste(label, 
                                                      "out", sep = "_"))
      main.file <- ifelse(is.null(label), "mainparams", paste(label, 
                                                              "mainparams", sep = "_"))
      extra.file <- ifelse(is.null(label), "extraparams", 
                           paste(label, "extraparams", sep = "_"))
      mat <- .stackedAlleles(g, alleles2integer = T, na.val = -9) %>% 
        dplyr::select(-.data$allele) %>% dplyr::mutate(id = gsub(" ", 
                                                                 "_", .data$id), stratum = as.numeric(factor(.data$stratum)), 
                                                       popflag = popflag[.data$id]) %>% dplyr::select(.data$id, 
                                                                                                      .data$stratum, .data$popflag, dplyr::everything()) %>% 
        as.matrix()
      write(paste(  sort(unique(g@data[["locus"]])), collapse = " "), file = in.file)
      for (i in 1:nrow(mat)) {
        write(paste(mat[i, ], collapse = " "), file = in.file, 
              append = TRUE)
      }
      main.params <- c(paste("MAXPOPS", as.integer(maxpops)), 
                       paste("BURNIN", as.integer(burnin)), paste("NUMREPS", 
                                                                  as.integer(numreps)), paste("INFILE", in.file), 
                       paste("OUTFILE", out.file), paste("NUMINDS", 
                                                         length(unique(g@data$id))), paste("NUMLOCI", length(unique(g@data$locus))), 
                       "MISSING -9", "LABEL 1", "POPDATA 1", 
                       "POPFLAG 1", "LOCDATA 0", "PHENOTYPE 0", 
                       "EXTRACOLS 0", "MARKERNAMES 1")
      main.params <- paste("#define", main.params)
      write(main.params, file = main.file)
      extra.params <- c(paste("NOADMIX", as.integer(noadmix)), 
                        paste("FREQSCORR", as.integer(freqscorr)), paste("INFERALPHA", 
                                                                         as.integer(inferalpha)), paste("ALPHA", as.numeric(alpha)), 
                        "FPRIORMEAN 0.01", "FPRIORSD 0.05", "LAMBDA 1.0", 
                        paste("UNIFPRIORALPHA", as.integer(unifprioralpha)), 
                        paste("ALPHAMAX", as.numeric(alphamax)), paste("ALPHAPRIORA", 
                                                                       as.numeric(alphapriora)), paste("ALPHAPRIORB", 
                                                                                                       as.numeric(alphapriorb)), "COMPUTEPROB 1", 
                        paste("ADMBURNIN", max(0, as.integer(burnin/2))), 
                        "ALPHAPROPSD 0.025", "STARTATPOPINFO 0", 
                        paste("RANDOMIZE", as.integer(randomize)), paste("SEED", 
                                                                         as.integer(seed)), "METROFREQ 10", "REPORTHITRATE 0")
      if (!is.null(pop.prior)) {
        pop.prior <- tolower(pop.prior)
        prior.params <- if (pop.prior == "locprior") {
          c("LOCPRIOR 1", "LOCISPOP 1", paste("LOCPRIORINIT", 
                                              locpriorinit), paste("MAXLOCPRIOR", maxlocprior))
        }
        else if (pop.prior == "usepopinfo") {
          c("USEPOPINFO 1", paste("GENSBACK", trunc(gensback)), 
            paste("MIGRPRIOR", migrprior), paste("PFROMPOPFLAGONLY", 
                                                 as.integer(pfrompopflagonly)))
        }
        extra.params <- c(extra.params, prior.params)
      }
      extra.params <- extra.params[!is.na(extra.params)]
      extra.params <- paste("#define", extra.params)
      write(extra.params, file = extra.file)
      invisible(list(files = c(data = in.file, mainparams = main.file, 
                               extraparams = extra.file, out = out.file), pops = sort(unique(g@data$stratum))))
    }
###########################################################
    structureRead <-    function (file, pops = NULL) 
    {
      if (!file.exists(file)) {
        stop(paste("the file '", file, "' can't be found.", 
                   sep = ""))
      }
      result <- scan(file, "character", quiet = TRUE)
      loc <- grep("Estimated", result, ignore.case = FALSE, 
                  value = FALSE)
      est.ln.prob <- as.numeric(result[loc[1] + 6])
      loc <- grep("likelihood", result, ignore.case = FALSE, 
                  value = FALSE)
      mean.lnL <- as.numeric(result[loc[1] + 2])
      var.lnL <- as.numeric(result[loc[2] + 2])
      loc <- grep("MAXPOPS", result, value = F)
      maxpops <- result[loc]
      maxpops <- sub("MAXPOPS=", "", maxpops)
      maxpops <- as.integer(sub(",", "", maxpops))
      loc <- grep("GENSBACK", result, value = F)
      gensback <- result[loc]
      gensback <- sub("GENSBACK=", "", gensback)
      gensback <- as.integer(sub(",", "", gensback))
      smry <- c(k = maxpops, est.ln.prob = est.ln.prob, mean.lnL = mean.lnL, 
                var.lnL = var.lnL)
      result <- scan(file, "character", sep = "\n", 
                     quiet = TRUE)
      first <- grep("(%Miss)", result, value = FALSE) + 1
      last <- grep("Estimated Allele", result, value = FALSE) - 
        1
      tbl.txt <- result[first:last]
      tbl.txt <- sub("[*]+", "", tbl.txt)
      tbl.txt <- sub("[(]", "", tbl.txt)
      tbl.txt <- sub("[)]", "", tbl.txt)
      tbl.txt <- sub("[|][ ]+$", "", tbl.txt)
      prior.lines <- grep("[|]", tbl.txt)
      no.prior <- if (length(prior.lines) < length(tbl.txt)) {
        no.prior.q.txt <- if (length(prior.lines) == 0) 
          tbl.txt
        else tbl.txt[-prior.lines]
        .structureParseQmat(no.prior.q.txt, pops)
      }
      else NULL
      if (maxpops == 1) {
        no.prior$row <- NULL
        return(list(summary = smry, q.mat = no.prior, prior.anc = NULL))
      }
      has.prior <- if (length(prior.lines) > 0) {
        prior.txt <- strsplit(tbl.txt[prior.lines], "[|]")
        prior.q.txt <- unlist(lapply(prior.txt, function(x) x[1]))
        df <- .structureParseQmat(prior.q.txt, pops)
        prior.anc <- purrr::map(prior.txt, function(x) {
          anc.mat <- matrix(NA, nrow = maxpops, ncol = gensback + 
                              1)
          rownames(anc.mat) <- paste("Pop", 1:nrow(anc.mat), 
                                     sep = ".")
          colnames(anc.mat) <- paste("Gen", 0:gensback, 
                                     sep = ".")
          x <- sapply(strsplit(x[-1], "\\s|[:]"), function(y) {
            y <- y[y != ""]
            y[-1]
          })
          for (i in 1:ncol(x)) {
            pop <- as.numeric(x[1, i])
            anc.mat[pop, ] <- as.numeric(x[-1, i])
          }
          anc.mat
        }) %>% stats::setNames(df$id)
        prob.mat <- t(sapply(1:nrow(df), function(i) {
          pop.probs <- rowSums(prior.anc[[i]])
          pop.probs[is.na(pop.probs)] <- df$Group.1[i]
          pop.probs
        }))
        colnames(prob.mat) <- paste("Group", 1:ncol(prob.mat), 
                                    sep = ".")
        df$Group.1 <- NULL
        df <- cbind(df, prob.mat)
        list(df = df, prior.anc = prior.anc)
      }
      else NULL
      has.prior.df <- if (is.null(has.prior)) 
        NULL
      else has.prior$df
      q.mat <- rbind(no.prior, has.prior.df)
      q.mat <- q.mat[order(q.mat$row), ]
      q.mat$row <- NULL
      rownames(q.mat) <- NULL
      q.mat[, -(1:3)] <- t(apply(q.mat[, -(1:3)], 1, function(i) i/sum(i)))
      prior.anc <- if (is.null(has.prior)) 
        NULL
      else has.prior$prior.anc
      list(summary = smry, q.mat = q.mat, prior.anc = prior.anc)
    }

###########################################################
    
    label <- .getFileLabel(g, label)
    label <- paste(label, "structureRun", sep = ".")
    label <- gsub("[[:space:]]", ".", label)
    unlink(label, recursive = TRUE, force = TRUE)
    dir.create(label)
    if (!utils::file_test("-d", label)) {
      stop(paste("'", label, "' is not a valid folder.", 
                 sep = ""))
    }
    label <- file.path(label, label)
    if (is.null(k.range)) 
      k.range <- 1:(dplyr::n_distinct(g@data$stratum))
    rep.df <- expand.grid(rep = 1:num.k.rep, k = k.range)
    rownames(rep.df) <- paste(label, ".k", rep.df$k, ".r", 
                              rep.df$rep, sep = "")
    out.files <- lapply(rownames(rep.df), function(x) {
      sw.out <- structureWrite(g, label = x, maxpops = rep.df[x,"k"])
      files <- sw.out$files
      cmd <- paste0(exec, " -m ", files["mainparams"], 
                    " -e ", files["extraparams"], " -i ", 
                    files["data"], " -o ", files["out"])
      err.code <- system(cmd)
      if (err.code == 127) {
        stop("You do not have STRUCTURE installed.")
      }
      else if (!err.code == 0) {
        stop(paste("Error running STRUCTURE. Error code", 
                   err.code, "returned."))
      }
      files["out"] <- paste(files["out"], "_f", 
                            sep = "")
      result <- structureRead(files["out"], sw.out$pops)
      if (file.exists("seed.txt")) 
        file.remove("seed.txt")
      files <- if (delete.files) 
        NULL
      else files
      result <- c(result, list(files = files, label = basename(x)))
      fname <- paste(x, ".ws.rdata", sep = "")
      save(result, file = fname)
      fname
    })
    run.result <- lapply(out.files, function(f) {
      result <- NULL
      load(f)
      result
    })
    names(run.result) <- sapply(run.result, function(x) x$label)
    class(run.result) <- c("structure.result", class(run.result))
    if (delete.files) 
      unlink(dirname(label), recursive = TRUE, force = TRUE)
    run.result
  }


