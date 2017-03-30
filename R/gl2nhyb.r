#' Create an input file for the program NewHybrids
#'
#' This function compares two sets of parental populations to identify loci
#' that exhibit a fixed difference, returns an genlight object with the reduced data,
#' and creates an input file for the program NewHybrids using the top 200 loci. In
#' the absence of two identified parental populations, the script will
#' select a random set 200 loci only (method=random) or the first 200 loci ranked
#' on information content (AvgPIC).
#'
#' A fixed difference occurs when a SNP allele is present in all individuals
#' of one population and absent in the other. There is provision for setting
#' a level of tollerance, e.g. t = 0.05 which considers alleles present
#' at greater than 95% in one population and less than 5% in the other to be
#' a fixed difference. Only the 200 loci are retained, because of limitations
#' of NewHybids.
#'
#' It is important to stringently filter the
#' data on RepAvg and CallRate if using the random option. One might elect to repeat
#' the analysis using the random option and combine the resultant posterior probabilites
#' should 200 loci be considered insufficient.
#'
#' @param gl -- name of the genlight object containing the SNP data [required]
#' @param outfile -- name of the file to become the input file for NewHybrids [default nhyb.txt]
#' @param p0 -- list of populations to be regarded as parental population 0 [default NULL]
#' @param p1 -- list of populations to be regarded as parental population 1 [default NULL]
#' @param t -- sets the level at which a gene frequency difference is considered to be fixed [default 0]
#' @param m -- specifies the method to use to select 200 loci for NewHybrids [default random]
#' @return The reduced genlight object
#' @export
#' @importFrom MASS write.matrix
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' \dontrun{
#' m <- gl2nhyb(gl, c("Pop1", "Pop4"), c("Pop7", "Pop9"), t=0, m="random")
#' }

gl2nhyb <- function(gl, outfile="nhyb.txt", p0=NULL, p1=NULL, t=0, m="random") {

  # Author: Arthur Georges
  #
  # Debug
  #  outfile <- "nhyb.txt"
  #  t <- 0
  #  p0 <- NULL
  #  p1 <- NULL
  #

  if(class(gl)!="genlight") {
    cat("Fatal Error: genlight object required for gl2nhyb.r!\n"); stop()
  }

  # EXTRACT THE SNP DATA
    cat("Extracting the SNP data\n")
    gl.tmp <- gl
    thold<-t

  # PROCESS AS FOLLOWS IF BOTH PARENTAL POPULATIONS ARE SPECIFIED
    if (!is.null(p0) & !is.null(p1)) {
      cat("  Both parental populations have been specified \n")

      # Replace those names with Z0 for Parental Population 0, and z1 for Parental Population 1
      indNames(gl.tmp) <- replace(indNames(gl.tmp), is.element(pop(gl.tmp), p0), "z0")
      indNames(gl.tmp) <- replace(indNames(gl.tmp), is.element(pop(gl.tmp), p1), "z1")

      # Error checks
      if (length(indNames(gl.tmp)[indNames(gl.tmp)=="z0"]) < 1  | length(indNames(gl.tmp)[indNames(gl.tmp)=="z1"]) < 1) {
        cat("Fatal Error [gl2nhyb]: One of the two parental populations contains no individuals\n"); stop()
      }

      # Create a vector containing the flags for the parental population
      par.names <- indNames(gl.tmp)
      par.names <- replace(par.names, (par.names != "z0" & par.names != "z1"), " ")

      # Discard non-parental populations
      gl.tmp <- gl.tmp[(indNames(gl.tmp)=="z0" | indNames(gl.tmp)=="z1"),]
      pop(gl.tmp) <- indNames(gl.tmp)

      # Reformat the data broken down by population and locus, and calculate allele frequencies
      gl2 <- gl.percent.freq(gl.tmp)

      # IDENTIFY LOCI WITH FIXED DIFFERENCES BETWEEN P0 AND P1
      cat("  Identifying loci with fixed difference between parental stocks\n")

      # Establish a vector to hold the loci
      npops<-2
      nloci<-nlevels(gl2$locus)
      fixed.loci <- NA
      length(fixed.loci) <- nloci*2

      # Cycle through the data to identify the fixed loci
      for (i in seq(1, nloci*2, 2)) {
        if (as.character(gl2$locus[i]) != as.character(gl2$locus[i+1])) {
          cat("Warning: Loci do not agree for the is.fixed comparison\n")
        }
        if (!is.na(is.fixed(gl2$frequency[i],gl2$frequency[i+1], t=thold))) {
          if (is.fixed(gl2$frequency[i],gl2$frequency[i+1], t=thold)) {
            fixed.loci[i] <- as.character(gl2$locus[i])
          }
        }
      }

      # Remove the NAs
      fixed.loci <- fixed.loci[!is.na(fixed.loci)]

      # If no loci remain, set the data matrix to be the original matrix
      if (length(fixed.loci) == 0) {
        cat("  No fixed differences between parental populations \n")
        if(m=="random") {cat("    Selecting ca 200 random loci\n")}
        if (m=="avgpic") {cat("    Selecting 200 loci with most information content (AvgPIC)\n")}
        tmp <- gl.subsample.loci(gl, 200, method=m)
        gl2 <- as.matrix(tmp)
        flag <- "bothparnonefixed"
      } else {
      # Set the data matrix to contain only the fixed differences
        gl.reduced <- gl[, (locNames(gl) %in% fixed.loci)]
        gl.reduced@other$loc.metrics <- gl@other$loc.metrics[(locNames(gl) %in% fixed.loci),]
        k <- min(nLoc(gl.reduced),200)
        gl2 <- as.matrix(gl.reduced[,1:k])
        cat("  Only",k,"loci with fixed differences between parental populations have been retained\n")
        flag <- "bothpar"
      }

    }

    # PROCESS AS FOLLOWS IF ONLY ONE PARENTAL POPULATION IS SPECIFIED
    if ((!is.null(p0) & is.null(p1)) || (is.null(p0) & !is.null(p1))) {
      cat("  Only one parental population specified \n")
      if(m=="random") {cat("    Selecting ca 200 random loci\n")}
      if (m=="avgpic") {cat("    Selecting 200 loci with most information content (AvgPIC)\n")}
      tmp <- gl.subsample.loci(gl, 200, method=m)
      gl2 <- as.matrix(tmp)
      gl.reduced <- NULL
      flag <- "onepar"
    }

    # PROCESS AS FOLLOWS IF NO PARENTAL POPULATION IS SPECIFIED
    if (is.null(p0) & is.null(p1)) {
      cat("  No parental population specified \n")
      if(m=="random") {cat("    Selecting ca 200 random loci\n")}
      if (m=="avgpic") {cat("    Selecting 200 loci with most information content (AvgPIC)\n")}
      tmp <- gl.subsample.loci(gl, 200, method=m)
      gl2 <- as.matrix(tmp)
      gl.reduced <- NULL
      flag <- "nopar"
    }
    # Repecify the parental names (only done above for case where both pops were specified)
    # Replace those names with Z0 for Parental Population 0, and z1 for Parental Population 1
    gl.tmp <- gl
    indNames(gl.tmp) <- replace(indNames(gl.tmp), is.element(pop(gl.tmp), p0), "z0")
    indNames(gl.tmp) <- replace(indNames(gl.tmp), is.element(pop(gl.tmp), p1), "z1")
    # Create a vector containing the flags for the parental population
    par.names <- indNames(gl.tmp)
    par.names <- replace(par.names, (par.names != "z0" & par.names != "z1"), " ")

    # CREATE THE NEWHYBRIDS INPUT FILE
    # Convert to NewHrbrids lumped format
    # Recode values
    #  0 to 11
    #  1 to 12
    #  2 to 22
    #  NA to 0
    cat("Converting data to NewHybrids format\n")
    gl2[gl2 == 2] <- 22
    gl2[gl2 == 1] <- 12
    gl2[gl2 == 0] <- 11
    gl2[is.na(gl2)] <- 0
    n.loci <- ncol(gl2)

    # Create sequential row number
    rownum <- seq(1:nrow(gl2))

    # Bind to the matrix
    gl2 <- data.frame(gl2)

    #NOTE: NewHybrids does not seem to work with the addition of specimen names.
    # ind.names <- indNames(gl)
    #if (flag=="bothpar" || flag == "onepar") {
    #  gl2 <- cbind(rownum, " n ", ind.names, par.names, gl2)
    #  metarows <- 4
    #  cat("  Adding sequential number, individual names, and flagging parental stock\n")
    #}
    #if (flag=="nopar") {
    #  gl2 <- cbind(rownum, " n ", ind.names, gl2)
    #  metarows <- 3
    #  cat("  Adding sequential number and individual names, parental stock not identified\n")
    #}

    if (flag=="bothpar" || flag == "onepar" || flag == "bothparnonefixed") {
      gl2 <- cbind(rownum, par.names, gl2)
      metarows <- 2
      cat("  Adding sequential number and flagging parental stock\n")
    }
    if (flag=="nopar") {
      gl2 <- cbind(rownum, gl2)
      metarows <- 1
      cat("  Adding sequential number\n")
    }

    # Clear row and column names
    rownames(gl2) <- NULL
    colnames(gl2) <- NULL

    # Output data file
    cat("Writing the NewHybrids input file", outfile, "\n")
    cat(c("  NumIndivs ", nrow(gl2), "\n"))
    cat(c("  NumLoci ", n.loci, " \n"))
    cat(c("  Digits 1\n  Format Lumped \n"))
    sink(outfile)
    cat(c("NumIndivs ", nrow(gl2), "\n"))
    cat(c("NumLoci ", n.loci, " \n"))
    cat(c("Digits 1\nFormat Lumped \n"))
    write.matrix(gl2[,1:(ncol(gl2))], sep=" ")
    sink()

    return(gl.reduced)

}
