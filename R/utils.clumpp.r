# Functions from package starmie for merging Q matrices from Structure runs
# using the CLUMPP algorithms.
# @param Q_list A list of of Q matrices.
# @param method The algorithm to use to infer the correct permutations. One of
#  'greedy' or 'greedyLargeK' or 'stephens'
# @param iter The number of iterations to use if running either 'greedy' or
#  'greedyLargeK'
# @importFrom purrr map_dbl
#' @export

utils.clumpp <- function(Q_list,
                   method,
                   iter) {
  # i/o checks
  if (!(method %in% c("greedy", "greedyLargeK", "stephens"))) {
    stop("Not a valid CLUMPP method, please use on of: 'greedy', 'greedyLargeK' or 'stephens'")
  }
  
  if (!all.equal(iter, as.integer(iter)) || iter < 0) {
    stop("number of iterations must be a positive integer")
  }
  
  # if length of Q_list is 1, clummping is not necessary
  if (length(Q_list) == 1) {
    return(Q_list)
  }
  
  # check dims of input Q_list are all equal
  dim_Q <- dim(Q_list[[1]])
  if (!all(unlist(lapply(Q_list[-1], function(x)
    all(dim(
      x
    ) == dim_Q))) == TRUE)) {
    stop(error("  Size of all matrices in Q_list must be equal"))
  }
  
  K <- dim_Q[2]
  
  # gracefully return Q_list if K=1
  if (K == 1) {
    message(warn(
      "  Number of assumed populations = 1 for all Q matrices, clummping is unecessary, returning Q_list"
    ))
    return(Q_list)
  }
  
  if (method == "greedy") {
    #Greedy clumpp algorithm
    perms <-
      replicate(iter, sample(
        1:length(Q_list),
        size = length(Q_list),
        replace = FALSE
      ))
    if (K > 8) {
      permQs <- apply(perms, 2, function(p)
        iterativeGreedy(Q_list[p]))
      Hs <-
        lapply(permQs, function(x)
          averagePairWiseSimilarityH(x$Q_list))
      Q_list <- permQs[[which.max(Hs)]]
    } else{
      permQs <- apply(perms, 2, function(p)
        memoryGreedy(Q_list[p]))
      Hs <-
        lapply(permQs, function(x)
          averagePairWiseSimilarityH(x$Q_list))
      Q_list <- permQs[[which.max(Hs)]]
    }
    
  } else if (method == "greedyLargeK") {
    #Use LargeKGreedy algorithm
    perms <-
      replicate(iter, sample(
        1:length(Q_list),
        size = length(Q_list),
        replace = FALSE
      ))
    permQs <- apply(perms, 2, function(p){
      largeKGreedy(Q_list[p])
    })
    Hs <-
      lapply(permQs, function(x)
        averagePairWiseSimilarityH(x$Q_list))
    Q_list <- permQs[[which.max(Hs)]]
    
  } else if (method == "stephens") {
    Q_list <- getStephens(Q_list)
    
  }
  return(Q_list)
}

memoryGreedy <- function(Q_list) {
  #Create matrix to store permutations
  K <- ncol(Q_list[[1]])
  permutations <- matrix(0, nrow = length(Q_list), ncol = K)
  permutations[1, ] <- seq(1, K)
  
  #Faster but with a high memory footprint for large K
  for (i in 1:(length(Q_list) - 1)) {
    permuations <- permn(1:ncol(Q_list[[i + 1]]))
    perm_scores <- purrr::map_dbl(permuations, J_perm, Q_list[[i + 1]], Q_list[1:i])
    perm <- permuations[[which.max(perm_scores)]]
    Q_list[[i + 1]] <- Q_list[[i + 1]][, perm]
    permutations[i + 1, ] <- perm
  }
  
  #Check for bug in using the same column twice
  if (!all(unlist(lapply(Q_list, function(x) {
    length(unique(colnames(x))) == ncol(x)
  }))))
    stop("Duplicated column names in output Q matrices")
  
  #Rename columns
  column_names <- paste("Cluster ", seq(1, ncol(Q_list[[1]])))
  Q_list <- lapply(Q_list, function(x) {
    colnames(x) <- column_names
    return(x)
  })
  
  return(list(Q_list = Q_list, permutations = permutations))
}

J_perm <- function(perm, Q_x, Q_sub_list) {
  J(Q_x[, perm], Q_sub_list)
}

J <- function(Q_x, Q_sub_list) {
  sum(unlist(purrr::map_dbl(Q_sub_list, G, Q_x))) / length(Q_sub_list)
}

G <- function(Q_1, Q_2) {
  W <- matrix(1, nrow(Q_1), ncol(Q_1)) / ncol(Q_1)
  1 - norm(Q_1 - Q_2, type = "F") / sqrt(norm(Q_1 - W, type = "F") * norm(Q_2 -
                                                                            W, type = "F"))
}

# function from package combinat
# Generates all permutations of the elements of x
permn <- function (x, fun = NULL, ...) {
  if (is.numeric(x) && length(x) == 1 && x > 0 && trunc(x) == x) {
    x <- seq(x)
  }
  n <- length(x)
  nofun <- is.null(fun)
  out <- vector("list", gamma(n + 1))
  p <- ip <- seqn <- 1:n
  d <- rep(-1, n)
  d[1] <- 0
  m <- n + 1
  p <- c(m, p, m)
  i <- 1
  use <- -c(1, n + 2)
  while (m != 1) {
    out[[i]] <- if (nofun)
      x[p[use]]
    else
      fun(x[p[use]], ...)
    i <- i + 1
    m <- n
    chk <- (p[ip + d + 1] > seqn)
    m <- max(seqn[!chk])
    if (m < n)
      d[(m + 1):n] <- -d[(m + 1):n]
    index1 <- ip[m] + 1
    index2 <- p[index1] <- p[index1 + d[m]]
    p[index1 + d[m]] <- m
    tmp <- ip[index2]
    ip[index2] <- ip[m]
    ip[m] <- tmp
  }
  out
}

# Compute average pairwise similarity between Q matrices.
# @param Q_list A list of of Q matrices.
# @details Implementation of the pairwise similarity metric as
# defined in Jakobsson, M. and Rosenberg, N. A., 2007.
# @examples
# # Read in Structure files
# multiple_runs_k10 <- exampleStructure("mcmc_diagnostics")
# Q_list <- lapply(multiple_runs_k10, getQ)
# avgQ <- averagePairWiseSimilarityH(Q_list)
averagePairWiseSimilarityH <- function(Q_list) {
  #i/o checks
  if (!all(unlist(lapply(Q_list, inherits, "matrix"))))
    stop(error("cluster runs must be a list of Q matrices"))
  
  R <- length(Q_list)
  H <- 0
  for (i in 1:(R - 1)) {
    for (j in (i + 1):R) {
      H <- H + G(Q_list[[i]], Q_list[[j]])
    }
  }
  H <- 2 / (R * (R - 1)) * H
  return(H)
}

# Use the Stephen's method to permute sample labels
# @param Q_list A list of of Q matrices.
getStephens <- function(Q_list) {
  # Create 3-dimensional array for input into stephens function
  # dimensions are equal to R by n by K
  # R := number of runs
  # n := number of rows in Q matrix
  # K := number of columns in Q matrix
  # convert list to array then transpose columms
  p <- aperm(simplify2array(Q_list), c(3, 1, 2))
  
  pkg <- "label.switching"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }
  
  perm <- label.switching::stephens(p)

  # reorder columns in according to new permuations
  # Rename columns
  column_names <- paste("Cluster ", seq_len(dim(p)[3]))
  Q_update <- lapply(seq_len(dim(p)[1]),
                     function(i) {
                       q_perm <- Q_list[[i]][, perm$permutations[i,]]
                       colnames(q_perm) <- column_names
                       q_perm
                     })
  
  
  return(list(Q_list = Q_update, permutations = perm$permutations))
}

iterativeGreedy <- function(Q_list) {
  #Create matrix to store permutations
  K <- ncol(Q_list[[1]])
  permutations <- matrix(0, nrow = length(Q_list), ncol = K)
  permutations[1, ] <- seq(1, K)
  
  #slower but memory efficient
  for (i in 1:(length(Q_list) - 1)) {
    permuations <- iterpc::iterpc(ncol(Q_list[[i + 1]]), ordered = TRUE)
    
    max_perm <- 1:ncol(Q_list[[i + 1]])
    max <- -Inf
    j <- 0
    while (j < iterpc::getlength(permuations)) {
      value <- J_perm(iterpc::getnext(permuations), Q_list[[i + 1]], Q_list[1:i])
      if (value > max) {
        max <- value
        max_perm <- iterpc::getcurrent(permuations)
      }
      j <- j + 1
    }
    Q_list[[i + 1]] <- Q_list[[i + 1]][, max_perm]
    permutations[i + 1, ] <- max_perm
  }
  
  #Check for bug in using the same column twice
  if (!all(unlist(lapply(Q_list, function(x) {
    length(unique(colnames(x))) == ncol(x)
  }))))
    stop("Duplicated column names in output Q matrices")
  
  #Rename columns
  column_names <- paste("Cluster ", seq(1, ncol(Q_list[[1]])))
  Q_list <- lapply(Q_list, function(x) {
    colnames(x) <- column_names
    return(x)
  })
  
  return(list(Q_list = Q_list, permutations = permutations))
}

largeKGreedy <- function(Q_list) {
  #Initial iteration
  #calculate pairwise column comparisons
  column_pairs <-
    iterpc::getall(iterpc::iterpc(ncol(Q_list[[1]]), 2, ordered = TRUE, replace = TRUE))
  pair_comparisons <- apply(column_pairs, 1, function(x) {
    list(G = G(Q_list[[1]][, x[1], drop = FALSE], Q_list[[2]][, x[2], drop =
                                                                FALSE]),
         Qy = x[1],
         Qz = x[2])
  })
  pair_comparisons_df <-
    data.frame(matrix(
      unlist(pair_comparisons),
      ncol = 3,
      byrow = TRUE
    ))
  
  #Create matrix to store permutations
  K <- ncol(Q_list[[1]])
  permutations <- matrix(0, nrow = length(Q_list), ncol = K)
  permutations[1, ] <- seq(1, K)
  
  #permute the second Q matrix
  perm <- get_best_permutation(pair_comparisons_df)
  Q_list[[2]] <- Q_list[[2]][, perm]
  permutations[2, ] <- perm
  
  if (length(Q_list) > 2) {
    for (i in 3:length(Q_list)) {
      #remaining iterations
      pair_comparisons <- apply(column_pairs, 1, function(x) {
        list(
          J = J_largeK(Q_list[1:(i - 1)],
                       Q_list[[i]], x[1],
                       x[2]),
          Qy = x[1],
          Qz = x[2]
        )
      })
      pair_comparisons_df <-
        data.frame(matrix(
          unlist(pair_comparisons),
          ncol = 3,
          byrow = TRUE
        ))
      perm <- get_best_permutation(pair_comparisons_df)
      Q_list[[i]] <- Q_list[[i]][, perm]
      permutations[i, ] <- perm
    }
  }
  
  #Check for bug in using the same column twice
  if (!all(unlist(lapply(Q_list, function(x) {
    length(unique(colnames(x))) == ncol(x)
  })))) {
    stop(error("Duplicated column names in output Q matrices"))
  }
  
  #Rename columns
  column_names <- paste("Cluster ", seq(1, ncol(Q_list[[1]])))
  Q_list <- lapply(Q_list, function(x) {
    colnames(x) <- column_names
    return(x)
  })
  
  return(list(Q_list = Q_list, permutations = permutations))
}

get_best_permutation <- function(pair_df) {
  #returns the best permutation of columns based on a dataframe of pairwise comparisons
  K <- max(pair_df[, 2])
  best_pairs_df <- data.frame(x = rep(0.0, K),
                              y = rep(0, K),
                              z = rep(0, K))
  for (i in 1:K) {
    best_pairs_df[i, ] <- pair_df[which.max(pair_df[, 1]), ]
    dont_keep <- pair_df[, 2] != pair_df[which.max(pair_df[, 1]), 2]
    dont_keep <-
      dont_keep & (pair_df[, 3] != pair_df[which.max(pair_df[, 1]), 3])
    pair_df <- pair_df[dont_keep,]
  }
  best_pairs_df <- best_pairs_df[order(best_pairs_df[, 2]), ]
  return(best_pairs_df[, 3])
}

J_largeK <- function(Q_sub_list, Q_x, y, z) {
  1 / length(Q_sub_list) * sum(unlist(lapply(Q_sub_list, function(Q) {
    G(Q[, y, drop = FALSE], Q_x[, z, drop = FALSE])
  })))
}

calcThreshold <- function(simMatrix) {
  #We want to choose a threshold t such that the number of singletons
  # is less than 10% of nodes and the mean node degree is at least 50%
  # of the total number of nodes.
  thresholds <- as.vector(simMatrix)
  thresholds <- thresholds[order(thresholds)]
  
  is_valid <- unlist(lapply(thresholds, function(t) {
    degrees <- apply(simMatrix, 1, function(r)
      sum(r[r > t]))
    ((sum(degrees == 1) / length(degrees)) < 0.1) &
      (mean(degrees) >= 0.5 * nrow(simMatrix))
  }))
  threshold <- thresholds[is_valid][sum(is_valid)]
  return(threshold)
}

# function  from MCL package
mcl <- function(x,
                addLoops = NULL,
                expansion = 2,
                inflation = 2,
                allow1 = FALSE,
                max.iter = 100,
                ESM = FALSE) {
  if (is.null(addLoops)) {
    stop("addLoops has to be TRUE or FALSE")
  }
  
  if (addLoops)
    diag(x) <- 1
  
  adj.norm <- apply(
    x[, ],
    MARGIN = 2,
    FUN = function(Var) {
      Var / sum(Var)
    }
  )
  
  a <- 1
  
  repeat {
    expans <- expm::`%^%`(adj.norm,expansion) 
    infl <- expans ^ inflation
    infl.norm <- apply(
      infl[, ],
      MARGIN = 2,
      FUN = function(Var) {
        Var / sum(Var)
        
      }
    )
    
    if (identical(infl.norm, adj.norm)) {
      ident <- TRUE
      break
    }
    
    if (a == max.iter) {
      ident <- FALSE
      a <- a + 1
      break
    }
    
    adj.norm <- infl.norm
    a <- a + 1
  }
  
  if (!is.na(infl.norm[1, 1]) & ident) {
    count <- 0
    for (i in 1:ncol(infl.norm)) {
      if (sum(abs(infl.norm[i, ])) != 0) {
        count <- count + 1
      }
    }
    
    neu <- matrix(nrow = count, ncol = ncol(infl.norm))
    
    zeile <- 1
    for (i in 1:nrow(infl.norm)) {
      if (sum(infl.norm[i, ]) != 0) {
        for (j in 1:ncol(infl.norm)) {
          neu[zeile, j] <- infl.norm[i, j]
        }
        zeile <- zeile + 1
      }
    }
    
    for (i in 1:nrow(neu)) {
      for (j in 1:ncol(neu)) {
        if ((neu[i, j] < 1) & (neu[i, j] > 0)) {
          neu[, j] <- 0
          neu[i, j] <- 1
        }
      }
    }
    
    for (i in 1:nrow(neu)) {
      for (j in 1:ncol(neu)) {
        if (neu[i, j] != 0) {
          neu[i, j] <- i
        }
      }
    }
    
    ClusterNummern <- sum(neu[, 1])
    for (j in 2:ncol(neu)) {
      ClusterNummern <- c(ClusterNummern, sum(neu[, j]))
    }
    
  }
  
  ifelse(!(!is.na(infl.norm[1, 1]) &
             ident),
         output <- paste("An Error occurred at iteration", a - 1),
         {
           if (!allow1) {
             dub <-
               duplicated(ClusterNummern) + duplicated(ClusterNummern, fromLast = T)
             for (i in 1:length(dub)) {
               if (dub[[i]] == 0)
                 ClusterNummern[[i]] <- 0
             }
           }
           
           #### dimnames for infl.norm
           dimnames(infl.norm) <-
             list(1:nrow(infl.norm), 1:ncol(infl.norm))
           
           output <- list()
           output[[1]] <- length(table(ClusterNummern))
           output[[2]] <- a - 1
           output[[3]] <- ClusterNummern
           output[[4]] <- infl.norm
           
           names(output) <- c("K",
                              "n.iterations",
                              "Cluster",
                              "Equilibrium.state.matrix")
         })
  ifelse(ESM == TRUE, return(output), return(output[-4]))
}
