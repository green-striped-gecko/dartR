#' Calculates the genomic relatedness matrix
#' 
#' The G matrix is calculated by centering the allele frequency matrix  of the second allele by substracting 2 times the allefrequency
#'@param gl -- a genlight object 
#'@param plotheatmap -- a switch if a heatmap should be shown [Default:TRUE] 
#'@param return.imputed switch if loci with imputed data should be returned (see ?A.mat in package rrBLUP)
#'@param ... parameters passed to function A.mat from package rrBLUP
#'@return a genomic relatedness matrix 
#'@importFrom stats heatmap cov var
#'@export
#'  
#'@examples
#'gl.grm(foxes.gl[1:5,1:10])  


gl.grm <- function(gl, plotheatmap=TRUE, return.imputed=FALSE, ...)
{
 
  Amat <-  function (X, min.MAF = NULL, max.missing = NULL, impute.method = "mean",   tol = 0.02, n.core = 1, shrink = FALSE, return.imputed = FALSE) 
  {
    if (mode(shrink) == "list") {
      shrink.method <- shrink$method
      if (!is.element(shrink.method, c("EJ", "REG"))) {
        stop("Invalid shrinkage method.")
      }
      shrink.iter <- shrink$n.iter
      n.qtl <- shrink$n.qtl
      shrink <- TRUE
    }
    else {
      if (shrink) {
        shrink.method <- "EJ"
      }
    }
    shrink.coeff <- function(i, W, n.qtl, p) {
      m <- ncol(W)
      n <- nrow(W)
      qtl <- sample(1:m, n.qtl)
      A.mark <- tcrossprod(W[, -qtl])/sum(2 * p[-qtl] * (1 - 
                                                           p[-qtl]))
      A.qtl <- tcrossprod(W[, qtl])/sum(2 * p[qtl] * (1 - p[qtl]))
      x <- as.vector(A.mark - mean(diag(A.mark)) * diag(n))
      y <- as.vector(A.qtl - mean(diag(A.qtl)) * diag(n))
      return(1 - cov(y, x)/var(x))
    }
    impute.EM <- function(W, cov.mat, mean.vec) {
      n <- nrow(W)
      m <- ncol(W)
      S <- matrix(0, n, n)
      for (i in 1:m) {
        Wi <- matrix(W[, i], n, 1)
        missing <- which(is.na(Wi))
        if (length(missing) > 0) {
          not.NA <- setdiff(1:n, missing)
          Bt <- solve(cov.mat[not.NA, not.NA], cov.mat[not.NA, 
                                                       missing])
          Wi[missing] <- mean.vec[missing] + crossprod(Bt, 
                                                       Wi[not.NA] - mean.vec[not.NA])
          C <- cov.mat[missing, missing] - crossprod(cov.mat[not.NA, missing], Bt)
          D <- tcrossprod(Wi)
          D[missing, missing] <- D[missing, missing] + 
            C
          W[, i] <- Wi
        }
        else {
          D <- tcrossprod(Wi)
        }
        S <- S + D
      }
      return(list(S = S, W.imp = W))
    }
    cov.W.shrink <- function(W) {
      m <- ncol(W)
      n <- nrow(W)
      Z <- t(scale(t(W), scale = FALSE))
      Z2 <- Z^2
      S <- tcrossprod(Z)/m
      target <- mean(diag(S)) * diag(n)
      var.S <- tcrossprod(Z2)/m^2 - S^2/m
      b2 <- sum(var.S)
      d2 <- sum((S - target)^2)
      delta <- max(0, min(1, b2/d2))
      print(paste("Shrinkage intensity:", round(delta, 2)))
      return(target * delta + (1 - delta) * S)
    }
    X <- as.matrix(X)
    n <- nrow(X)
    frac.missing <- apply(X, 2, function(x) {
      length(which(is.na(x)))/n
    })
    missing <- max(frac.missing) > 0
    freq <- apply(X + 1, 2, function(x) {
      mean(x, na.rm = missing)
    })/2
    MAF <- apply(rbind(freq, 1 - freq), 2, min)
    if (is.null(min.MAF)) {
      min.MAF <- 1/(2 * n)
    }
    if (is.null(max.missing)) {
      max.missing <- 1 - 1/(2 * n)
    }
    markers <- which((MAF >= min.MAF) & (frac.missing <= max.missing))
    m <- length(markers)
    var.A <- 2 * mean(freq[markers] * (1 - freq[markers]))
    one <- matrix(1, n, 1)
    mono <- which(freq * (1 - freq) == 0)
    X[, mono] <- 2 * tcrossprod(one, matrix(freq[mono], length(mono), 
                                            1)) - 1
    freq.mat <- tcrossprod(one, matrix(freq[markers], m, 1))
    W <- X[, markers] + 1 - 2 * freq.mat
    if (!missing) {
      if (shrink) {
        if (shrink.method == "EJ") {
          W.mean <- rowMeans(W)
          cov.W <- cov.W.shrink(W)
          A <- (cov.W + tcrossprod(W.mean))/var.A
        }
        else {
          if ((n.core > 1) & requireNamespace("parallel", 
                                              quietly = TRUE)) {
            it <- split(1:shrink.iter, factor(cut(1:shrink.iter, 
                                                  n.core, labels = FALSE)))
            delta <- unlist(parallel::mclapply(it, function(ix, 
                                                            W, n.qtl, p) {
              apply(array(ix), 1, shrink.coeff, W = W, 
                    n.qtl = n.qtl, p = p)
            }, W = W, n.qtl = n.qtl, p = freq.mat[1, ], 
            mc.cores = n.core))
          }
          else {
            delta <- apply(array(1:shrink.iter), 1, shrink.coeff, 
                           W = W, n.qtl = n.qtl, p = freq.mat[1, ])
          }
          delta <- mean(delta, na.rm = T)
          print(paste("Shrinkage intensity:", round(delta, 
                                                    2)))
          A <- tcrossprod(W)/var.A/m
          A <- (1 - delta) * A + delta * mean(diag(A)) * 
            diag(n)
        }
      }
      else {
        A <- tcrossprod(W)/var.A/m
      }
      rownames(A) <- rownames(X)
      colnames(A) <- rownames(A)
      if (return.imputed) {
        return(list(A = A, imputed = X))
      }
      else {
        return(A)
      }
    }
    else {
      isna <- which(is.na(W))
      W[isna] <- 0
      if (toupper(impute.method) == "EM") {
        if (m < n) {
          print("Linear dependency among the lines: imputing with mean instead of EM algorithm.")
        }
        else {
          mean.vec.new <- matrix(rowMeans(W), n, 1)
          cov.mat.new <- cov(t(W))
          if (qr(cov.mat.new)$rank < nrow(cov.mat.new) - 
              1) {
            print("Linear dependency among the lines: imputing with mean instead of EM algorithm.")
          }
          else {
            W[isna] <- NA
            A.new <- (cov.mat.new + tcrossprod(mean.vec.new))/var.A
            err <- tol + 1
            print("A.mat converging:")
            while (err >= tol) {
              A.old <- A.new
              cov.mat.old <- cov.mat.new
              mean.vec.old <- mean.vec.new
              if ((n.core > 1) & requireNamespace("parallel", 
                                                  quietly = TRUE)) {
                it <- split(1:m, factor(cut(1:m, n.core, 
                                            labels = FALSE)))
                pieces <- parallel::mclapply(it, function(mark2) {
                  impute.EM(W[, mark2], cov.mat.old, mean.vec.old)
                }, mc.cores = n.core)
              }
              else {
                pieces <- list()
                pieces[[1]] <- impute.EM(W, cov.mat.old, 
                                         mean.vec.old)
              }
              n.pieces <- length(pieces)
              S <- matrix(0, n, n)
              W.imp <- numeric(0)
              for (i in 1:n.pieces) {
                S <- S + pieces[[i]]$S
                W.imp <- cbind(W.imp, pieces[[i]]$W.imp)
              }
              mean.vec.new <- matrix(rowMeans(W.imp), n, 
                                     1)
              cov.mat.new <- (S - tcrossprod(mean.vec.new) * 
                                m)/(m - 1)
              A.new <- (cov.mat.new + tcrossprod(mean.vec.new))/var.A
              err <- norm(A.old - A.new, type = "F")/n
              print(err, digits = 3)
            }
            rownames(A.new) <- rownames(X)
            colnames(A.new) <- rownames(A.new)
            if (return.imputed) {
              Ximp <- W.imp - 1 + 2 * freq.mat
              colnames(Ximp) <- colnames(X)[markers]
              rownames(Ximp) <- rownames(X)
              return(list(A = A.new, imputed = Ximp))
            }
            else {
              return(A.new)
            }
          }
        }
      }
      if (shrink) {
        if (shrink.method == "EJ") {
          W.mean <- rowMeans(W)
          cov.W <- cov.W.shrink(W)
          A <- (cov.W + tcrossprod(W.mean))/var.A
        }
        else {
          if ((n.core > 1) & requireNamespace("parallel", 
                                              quietly = TRUE)) {
            it <- split(1:shrink.iter, factor(cut(1:shrink.iter, 
                                                  n.core, labels = FALSE)))
            delta <- unlist(parallel::mclapply(it, function(ix, 
                                                            W, n.qtl) {
              apply(array(ix), 1, shrink.coeff, W = W, 
                    n.qtl = n.qtl)
            }, W = W, n.qtl = n.qtl, mc.cores = n.core))
          }
          else {
            delta <- apply(array(1:shrink.iter), 1, shrink.coeff, 
                           W = W, n.qtl = n.qtl)
          }
          delta <- mean(delta, na.rm = T)
          print(paste("Shrinkage intensity:", round(delta, 
                                                    2)))
          A <- tcrossprod(W)/var.A/m
          A <- (1 - delta) * A + delta * mean(diag(A)) * 
            diag(n)
        }
      }
      else {
        A <- tcrossprod(W)/var.A/m
      }
      rownames(A) <- rownames(X)
      colnames(A) <- rownames(A)
      if (return.imputed) {
        Ximp <- W - 1 + 2 * freq.mat
        colnames(Ximp) <- colnames(X)[markers]
        rownames(Ximp) <- rownames(X)
        return(list(A = A, imputed = Ximp))
      }
      else {
        return(A)
      }
    }
  }  
  
   
G <- Amat(as.matrix(gl)-1,return.imputed = return.imputed, ...)
if (plotheatmap & return.imputed==FALSE) heatmap(G) else heatmap(G$A)

# ff <- as.matrix(gl)
# alf <- colMeans(ff, na.rm = T)/2
# pjm <- matrix(rep(alf,nInd(gl)), ncol=nLoc(gl), nrow=nInd(gl), byrow=T)
# W <- ff  - (2*pjm)
# het <- 2*sum(alf *(1-alf))
#             
# G <- (W %*% t(W) )/het
# GG <-G
# ii<- !(colMeans(is.na(G))==1)
# GG <- GG[, ii ]
# ii<- !(rowMeans(is.na(G))==1)
# GG <- GG[ii,  ]
# 
# if (plotheatmap & nrow(GG)>0) heatmap(GG)
return (G)
}