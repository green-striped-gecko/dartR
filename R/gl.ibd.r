#' Performs isolation by distance analysis
#'
#' This function performs an isolation by distance analysis based on a Mantel
#' test and also produces an isolation by distance plot. If a genlight object
#' with coordinates is provided, then an Euclidean and genetic distance matrices
#' are calculated.'
#' @importFrom vegan mantel
#' @importFrom MASS kde2d
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics abline title points
#' @importFrom stats as.dist lm
#' @importFrom StAMPP stamppFst stamppNeisD
#' @importFrom stats coef
#' @param x Genlight object. If provided a standard analysis on Fst/1-Fst and
#' log(distance) is performed [required].
#' @param distance Type of distance that is calculated and used for the
#' analysis. Can be either population based 'Fst' [\link[StAMPP]{stamppFst}],
#' 'D' [\link[StAMPP]{stamppNeisD}] or individual based 'propShared',
#'  [gl.propShared], 'euclidean' [gl.dist.ind, method='Euclidean']
#'  [default "Fst"].
#' @param coordinates Can be either 'latlon', 'xy' or a two column data.frame
#' with column names 'lat','lon', 'x', 'y'). Coordinates are provided via
#' \code{gl@other$latlon} ['latlon'] or via \code{gl@other$xy} ['xy']. If latlon
#' data will be projected to meters using Mercator system [google maps] or if
#' xy then distance is directly calculated on the coordinates.
#' @param Dgen Genetic distance matrix if no genlight object is provided
#' [default NULL].
#' @param Dgeo Euclidean distance matrix if no genlight object is provided
#' [default NULL].
#' @param Dgeo_trans Transformation to be used on the Euclidean distances. See
#' Dgen_trans [default "Dgeo"].
#' @param Dgen_trans You can provide a formula to transform the genetic
#' distance. The transformation can be applied as a formula using Dgen as the
#'  variable to be transformed. For example: \code{Dgen_trans = 'Dgen/(1-Dgen)'.
#'   Any valid R expression can be used here
#'    [default 'Dgen', which is the identity function.]}
#' @param permutations Number of permutations in the Mantel test [default 999].
#' @param plot.out Should an isolation by distance plot be returned
#' [default TRUE].
#' @param paircols Should pairwise dots colored by 'pop'ulation/'ind'ividual
#' pairs [default 'pop']. You can color pairwise individuals by pairwise
#'  population colors.
#' @param plot_theme Theme for the plot. See details for options
#' [default theme_dartR()].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session
#' temporary directory (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log ; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @details
#' Currently pairwise Fst and D between populations and
#' 1-propShared and Euclidean distance between individuals are
#' implemented. Coordinates are expected as lat long and converted to Google
#' Earth Mercator projection. If coordinates are already projected, provide them
#' at the x@other$xy slot.
#'
#' You can provide also your own genetic and Euclidean distance matrices. The
#' function is based on the code provided by the adegenet tutorial
#' (\url{http://adegenet.r-forge.r-project.org/files/tutorial-basics.pdf}),
#' using the functions  \link[vegan]{mantel} (package vegan),
#' \link[StAMPP]{stamppFst}, \link[StAMPP]{stamppNeisD} (package StAMPP) and
#' gl.propShared or gl.dist.ind. For transformation you need to have the dismo
#' package installed. As a new feature you can plot pairwise relationship using
#' double colored points (paircols=TRUE). Pairwise relationship can be
#' visualised via populations or individuals, depending which distance is
#' calculated. Please note: Often a problem arises, if an individual based 
#' distance is calculated (e.g. propShared) and some individuals have identical
#'  coordinates as this results in distances of zero between those pairs of 
#'  individuals.
#'  
#' If the standard transformation [log(Dgeo)] is used, this results in an 
#' infinite value, because of trying to calculate'log(0)'. To avoid this, the 
#' easiest fix is to change the transformation from log(Dgeo) to log(Dgeo+1) or 
#' you could add some "noise" to the coordinates of the individuals (e.g. +- 1m,
#'  but be aware if you use lat lon then you rather want to add +0.00001 degrees
#'   or so).
#' @return Returns a list of the following components: Dgen (the genetic
#' distance matrix), Dgeo (the Euclidean distance matrix), Mantel (the
#' statistics of the Mantel test).
#' @export
#' @author Bernd Gruber (bugs? Post to 
#' \url{https://groups.google.com/d/forum/dartr})
#' @seealso \link[vegan]{mantel}, \link[StAMPP]{stamppFst}
#' @references
#' Rousset, F. (1997). Genetic differentiation and estimation of gene flow from
#' F-statistics under isolation by distance. Genetics, 145(4), 1219-1228.
#' @examples
#'  \donttest{
#' #because of speed only the first 100 loci
#' ibd <- gl.ibd(bandicoot.gl[,1:100], Dgeo_trans='log(Dgeo)' ,
#' Dgen_trans='Dgen/(1-Dgen)')
#' #because of speed only the first 10 individuals)
#' ibd <- gl.ibd(bandicoot.gl[1:10,], distance='euclidean', paircols='pop', 
#' Dgeo_trans='Dgeo')
#' }
#' #only first 100 loci
#' ibd <- gl.ibd(bandicoot.gl[,1:100])

gl.ibd <- function(x = NULL,
                   distance = "Fst",
                   coordinates = "latlon",
                   Dgen = NULL,
                   Dgeo = NULL,
                   Dgeo_trans = "Dgeo",
                   Dgen_trans = "Dgen",
                   permutations = 999,
                   plot.out = TRUE,
                   paircols = NULL,
                   plot_theme = theme_dartR(),
                   save2tmp = FALSE,
                   verbose = NULL) {
  
    # CHECK IF PACKAGES ARE INSTALLED
  pkg <- "dismo"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  } else {
        
        # TRAP COMMAND
        funname <- match.call()[[1]]
        
        # GENERAL ERROR CHECKING
        
        verbose <- gl.check.verbosity(verbose)
        
        if (!is.null(x)){
            dt <- utils.check.datatype(x, verbose = 0)
        }
        
        # specific error checks
        
        if (!is.null(Dgen) & !is.null(Dgeo)) {
            if (verbose > 0)
                cat(
                    report(
                        "Analysis performed using provided genetic and Euclidean distance matrices. If a genlight object is provided, it is ignored.\n"
                    )
                )
            ta <-"dgendgeo"
            # make sure both matrices are distance objects if provided via Dgen 
            # and Dgeo directly
            Dgen <- as.dist(Dgen)
            Dgeo <- as.dist(Dgeo)
        }
        
        if (is(x, "genlight")) {
            if (verbose > 0)
                cat(report("Analysis performed on the genlight object.\n"))
            ta <-"genlight"
        }
        
        # check coordinates (if no Dgen and Dgeo is provided)
        if (ta == "genlight") {
            coords <- NULL
            if (is(coordinates, "character")) {
                if (coordinates == "latlon") {
                    if (is.null(x@other$latlon))
                        stop(error(
                            "Cannot find coordinates in x@other$latlon.\n"
                        ))
                    coords <-
                        dismo::Mercator(x@other$latlon[, c("lon", "lat")])
                    if (verbose > 0) {
                        cat(
                            report(
                                "Coordinates transformed to Mercator (google) projection to calculate distances in meters.\n"
                            )
                        )
                    }
                    coordstring <-"x@other$latlon (Mercator transformed)"
                }
                
                if (coordinates == "xy") {
                    if (is.null(x@other$xy))
                        stop(error("Cannot find coordinates in x@other$xy.\n"))
                    coords <- x@other$xy
                    coordstring <-"x@other$xy"
                }
                
            }
            
            if (is(coordinates, "data.frame")) {
                if (length(setdiff(colnames(coordinates), c("lat", "lon"))) == 0) {
                    coords <- dismo::Mercator(coordinates[, c("lon", "lat")])
                    coordstring <- "data.frame lat/lon (Mercator transformed)"
                }
                
                if (length(setdiff(colnames(coordinates), c("x", "y"))) == 0) {
                    coords <- coordinates[, c("x", "y")]
                    coordstring <-"data.frame x/y"
                }
                
                if (is.null(coords)){
                    stop(
                        error(
                            "No valid coordinates provided. check the provided data.frame and its format.\n"
                        )
                    )
            }
            }
            
            if (is.null(coords)){
                stop(error("No valid coordinates provided!\n"))
            }
            
            # make sure coordinates have the correct length
            if (nrow(coords) != nInd(x) & ta == "genlight"){
                stop(error(
                    "Cannot find coordinates for each individual in slot @other$latlon.\n"
                ))
            }
            
            typedis <-NULL
            if (distance == "Fst" | distance == "D") {
                typedis <-"pop"
            }
            
            if (distance == "propShared" |
                distance == "euclidean") {
                typedis <-"ind"
            }
            
            if (typedis == "pop" & nPop(x) < 2){
                stop(
                    error(
                        "You specified a population based distance, but there is either no population or only a single population specified within your genlight object. Check via table(pop(genlight)).\n"
                    )
                )
            }
            
            if (is.null(Dgeo) & typedis == "pop") {
                if (nPop(x) > 1) {
                    pop.xy <-
                        apply(coords, 2, function(a)
                            tapply(a, pop(x), mean, na.rm = TRUE))
                    Dgeo <- dist(pop.xy)
                } else {
                    stop(
                        error(
                            "Less than 2 populations provided, therefore no pairwise distances can be calculated.\n"
                        )
                    )
                }
            }
            
            if (is.null(Dgeo) & typedis == "ind") {
                if (nInd(x) > 1) {
                    Dgeo <- dist(coords)
                } else {
                    stop(
                        error(
                            "Less than 2 individuals provided, therefore no pairwise distances can be calculated.\n"
                        )
                    )
                }
            }
            
            # apply logarithm to distance
            
            if (is.null(Dgen) & distance == "Fst") {
                class(x)<- "genlight" #stampp issue
                Dgen <-
                    as.dist(StAMPP::stamppFst(x, nboots = 1))
            }
            
            if (is.null(Dgen) & distance == "D") {
              class(x)<- "genlight" #stampp issue
                Dgen <-
                    as.dist(StAMPP::stamppNeisD(x, pop = TRUE))
            }
            
            if (is.null(Dgen) & distance == "propShared") {
                Dgen <- as.dist(1 - gl.propShared(x))
            }
            
            if (is.null(Dgen) & distance == "euclidean") {
                Dgen <- as.dist(dist(as.matrix(x)))
            }
            
            ### order both matrices to be alphabetically as levels in genlight (ind or pop)
            if (is(x, "genlight")) {
                if (typedis == "pop") {
                    oo <- order(colnames(as.matrix(Dgen)))
                    Dgen <- as.dist(as.matrix(Dgen)[oo, oo])
                    oo <- order(colnames(as.matrix(Dgeo)))
                    Dgeo <- as.dist(as.matrix(Dgeo)[oo, oo])
                }
            }
        } else {
            # end of ta=='genlight' ta='dgendgeo
            coordstring <- "Dgeo provided."
            distance <- "Dgen provided"
            typedis <- "ind"
        }
        
        # use tranformations
        Dgen <- eval(parse(text = Dgen_trans))
        Dgeo <- eval(parse(text = Dgeo_trans))
        
        if (sum(is.infinite(Dgeo)) > 0) {
            stop(
                error(
                    "Most likely some pairwise individual distances were zero and the transformation created missing values [e.g. log(Dgeo)]. This affects the Mantel test and points are omitted from the plot. Consider adding a suitable tranformation e.g. an offset to your Dgeo transformation if using a log transformation [e.g. Dgeo_trans='log(Dgeo+1)'] or adding some 'noise' to the coordinates.\n"
                )
            )
            
        }
        
        if (is.null(Dgeo))
            stop(error(
                "Cannot calculate distance matrix or no distance matrix provided\n!"
            ))
        if (is.null(Dgen))
            stop(
                error(
                    "Cannot calculate genetic distance matrix or no genetic distance matrix provided!\n"
                )
            )
        
        manteltest <-
            vegan::mantel(Dgen, Dgeo, na.rm = TRUE, permutations = permutations)
        
        lm_eqn <-
            function(df,
                     r = manteltest$statistic,
                     pp = manteltest$signif) {
                m <- lm(Dgen ~ Dgeo, df)
                eq <-
                    substitute(
                        italic(y) == a + b %.% italic(x) * "," ~ ~ italic(R) ^ 2 ~ "=" ~ r2 * "," ~ ~
                            italic(p) ~ "=" ~ pp,
                        list(
                            a = format(unname(coef(m)[1]),
                                       digits = 2),
                            b = format(unname(coef(m)[2]), digits = 2),
                            r2 = format(summary(m)$r.squared, digits = 3),
                            pp = format(pp, digits = 3)
                        )
                    )
                as.character(as.expression(eq))
            }
        
        ####### Printing outputs, using package patchwork
        
        res <-
            data.frame(Dgen = as.numeric(Dgen), Dgeo = as.numeric(Dgeo))
        if (is.null(paircols)) {
          
            p3 <-
                ggplot(res, aes(x = Dgeo, y = Dgen)) + 
              geom_point() + 
              geom_smooth(method = "lm", se = TRUE) + 
              ylab(Dgen_trans) + 
              xlab(Dgeo_trans) +
                annotate(
                    "text",
                    label = lm_eqn(res),
                    x = Inf,
                    y = -Inf,
                    parse = TRUE,
                    hjust = 1.05,
                    vjust = 0) + 
              plot_theme
            
        } else {
            Legend <- col2 <- NA  #ggplot bug
            cols <- which(lower.tri(as.matrix(Dgen)), arr.ind = T)
            c1 <- cols[, 2]
            c2 <- cols[, 1]
            cn <- colnames(as.matrix(Dgen))
            # if someone wants to color pairwise individuals by pairwise colors
            if (typedis == "ind" & paircols == "pop") {
                if (is(x, "genlight"))
                    cn <- pop(x)
                else
                    cn <-rownames(as.matrix(Dgen))
            }
            res <-
                data.frame(
                    Dgen = as.numeric(Dgen),
                    Dgeo = as.numeric(Dgeo),
                    Legend = cn[c1],
                    col2 = cn[c2]
                )
            p3 <-
                ggplot(res) + 
              geom_point(aes(Dgeo, Dgen, col = Legend), size = 5) + 
              geom_point(aes(Dgeo, Dgen, col = col2), size = 2) + 
              geom_point(aes(Dgeo, Dgen),size = 2,shape = 1) + 
              guides(size = "none",color = guide_legend(title = "Populations")) + 
              geom_smooth(aes(x = Dgeo, y = Dgen),method = "lm", se = TRUE) + 
              ylab(Dgen_trans) + 
              annotate("text",label = lm_eqn(res), 
                       x = Inf,
                       y = -Inf, 
                       parse = TRUE,
                       hjust = 1.05,
                       vjust = 0) +
              xlab(Dgeo_trans) + plot_theme
            
        }
        
        if (plot.out) {
            suppressMessages(print(p3))
        }
        
        if (verbose > 0) {
            cat(report("  Coordinates used from:", coordstring, "\n"))
            cat(report("  Transformation of Dgeo:", Dgeo_trans, "\n"))
            cat(report("  Genetic distance:", distance, "\n"))
            cat(report("  Tranformation of Dgen: ", Dgen_trans, "\n"))
            print(manteltest)
        }
        
        # SAVE INTERMEDIATES TO TEMPDIR
        if (save2tmp & plot.out) {
            # creating temp file names
            match_call <-
                paste0(names(match.call()),
                       "_",
                       as.character(match.call()),
                       collapse = "_")
            temp_plot <- tempfile(pattern = "Plot_")
            
            # saving to tempdir
            saveRDS(list(match_call, p3), file = temp_plot)
            if (verbose >= 2) {
                cat(report("  Saving the ggplot to session tempfile\n"))
            }
        }
        
        if (save2tmp) {
            temp_table <- tempfile(pattern = "Table_")
            saveRDS(list(match_call, manteltest), file = temp_table)
            if (verbose >= 2) {
                cat(report("  Saving the report to the session tempfile\n"))
            }
        }
        
        # FLAG SCRIPT END
        
        if (verbose >= 1) {
            cat(report("\nCompleted:", funname, "\n\n"))
        }
        
        out <- list(Dgen = Dgen,
                    Dgeo = Dgeo,
                    mantel = manteltest)
        return(out)
    }
}
