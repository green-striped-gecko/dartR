#' Util function for evanno plots
#' 
#' These functions were copied from package strataG, which is no longer on CRAN (maintained by Eric Archer)
#' @export
#' @author Bernd Gruber (bugs? Post to
#'  \url{https://groups.google.com/d/forum/dartr}); original implementation of
#'   Eric Archer \url{https://github.com/EricArcher/strataG}
#' @param sr structure run object
#' @param plot should the plots be returned


utils.structure.evanno <- function (sr, plot = TRUE) 
{
  if (!"structure.result" %in% class(sr)) {
    stop("'sr' is not a result from 'structure.run'.")
  }
  k.tbl <- table(sapply(sr, function(x) x$summary["k"]))
  if (length(k.tbl) < 3) 
    stop("must have at least two values of k.")
  sr.smry <- t(sapply(sr, function(x) x$summary))
  ln.k <- tapply(sr.smry[, "est.ln.prob"], sr.smry[, 
                                                   "k"], mean)
  sd.ln.k <- tapply(sr.smry[, "est.ln.prob"], sr.smry[, 
                                                      "k"], stats::sd)
  ln.pk <- diff(ln.k)
  ln.ppk <- abs(diff(ln.pk))
  delta.k <- sapply(2:(length(ln.k) - 1), function(i) {
    abs(ln.k[i + 1] - (2 * ln.k[i]) + ln.k[i - 1])/sd.ln.k[i]
  })
  df <- data.frame(k = as.numeric(names(ln.k)), reps = as.numeric(table(sr.smry[, 
                                                                                "k"])), mean.ln.k = as.numeric(ln.k), sd.ln.k = as.numeric(sd.ln.k), 
                   ln.pk = c(NA, ln.pk), ln.ppk = c(NA, ln.ppk, NA), delta.k = c(NA, 
                                                                                 delta.k, NA))
  rownames(df) <- NULL
  df$sd.min <- df$mean.ln.k - df$sd.ln.k
  df$sd.max <- df$mean.ln.k + df$sd.ln.k
  plot.list <- list(mean.ln.k = ggplot2::ggplot(df, ggplot2::aes_string(x = "k", 
                                                                        y = "mean.ln.k")) + ggplot2::ylab("mean LnP(K)") + 
                      ggplot2::geom_segment(ggplot2::aes_string(x = "k", 
                                                                xend = "k", y = "sd.min", yend = "sd.max")), 
                    ln.pk = ggplot2::ggplot(df[!is.na(df$ln.pk), ], ggplot2::aes_string(x = "k", 
                                                                                        y = "ln.pk")) + ggplot2::ylab("LnP'(K)"), 
                    ln.ppk = ggplot2::ggplot(df[!is.na(df$ln.ppk), ], ggplot2::aes_string(x = "k", 
                                                                                          y = "ln.ppk")) + ggplot2::ylab("LnP''(K)"))
  if (!all(is.na(df$delta.k))) {
    plot.list$delta.k <- ggplot2::ggplot(df[!is.na(df$delta.k), 
    ], ggplot2::aes_string(x = "k", y = "delta.k")) + 
      ggplot2::ylab(expression(Delta(K)))
  }
  for (i in 1:length(plot.list)) {
    plot.list[[i]] <- plot.list[[i]] + ggplot2::geom_line() + 
      ggplot2::geom_point(fill = "white", shape = 21, 
                          size = 3) + ggplot2::xlim(c(1, max(df$k))) + 
      ggplot2::theme(axis.title.x = ggplot2::element_blank())
  }
  if (plot) {
    p <- plot.list %>% purrr::map(function(x) {
      ggplot2::ggplot_gtable(ggplot2::ggplot_build(x))
    })
    maxWidth <- do.call(grid::unit.pmax, purrr::map(p, function(x) x$widths[2:3]))
    for (i in 1:length(p)) p[[i]]$widths[2:3] <- maxWidth
    p$bottom <- "K"
    p$ncol <- 2
    do.call(gridExtra::grid.arrange, p)
  }
  df$sd.min <- df$sd.max <- NULL
  invisible(list(df = df, plots = plot.list))
}
