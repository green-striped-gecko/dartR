#' Create an input file for the program NewHybrids and run it if NewHybrids is installed
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
#' If you specify a directory for the NewHybrids executable file, then the script
#' will create the input file from the snp data then run NewHybrids. If the directory
#' is set to NULL, the exectution will stop once the input file (nhyb.txt) has been 
#' written to disk.
#' 
#' Refer to the New Hybrids manual for further information on the parameters to set
#' -- http://ib.berkeley.edu/labs/slatkin/eriq/software/new_hybs_doc1_1Beta3.pdf
#'
#' It is important to stringently filter the
#' data on RepAvg and CallRate if using the random option. One might elect to repeat
#' the analysis (method=random) and combine the resultant posterior probabilites
#' should 200 loci be considered insufficient.
#'
#' @param gl -- name of the genlight object containing the SNP data [required]
#' @param outfile -- name of the file that will be the input file for NewHybrids [default nhyb.txt]
#' @param outpath -- path where to save the output file (set to tempdir by default)
#' @param p0 -- list of populations to be regarded as parental population 0 [default NULL]
#' @param p1 -- list of populations to be regarded as parental population 1 [default NULL]
#' @param t -- sets the level at which a gene frequency difference is considered to be fixed [default 0]
#' @param method -- specifies the method (random or AvgPIC) to select 200 loci for NewHybrids [default random]
#' @param nhyb.directory -- directory that holds the NewHybrids executable file e.g. C:/NewHybsPC [default NULL]
#' @param BurnIn -- number of sweeps to use in the burn in [default 10000]
#' @param sweeps -- number  of  sweeps  to  use  in  computing  the  actual Monte Carlo averages [default 10000]
#' @param GtypFile -- name of a file containing the genotype frequency classes [default TwoGensGtypFreq.txt]
#' @param AFPriorFile -- name of the file containing prior allele frequency information [default NULL]
#' @param PiPrior -- Jeffreys-like priors or Uniform priors for the parameter pi [default Jeffreys]
#' @param ThetaPrior -- Jeffreys-like priors or Uniform priors for the parameter theta [default Jeffreys]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return The reduced genlight object, if parentals are provided; output of NewHybrids to disk
#' @export
#' @importFrom MASS write.matrix
#' @author Arthur Georges (bugs? Post to https://groups.google.com/d/forum/dartr)
#' @examples
#' \donttest{
#' m <- gl.nhybrids(testset.gl, c("Pop1", "Pop4"), c("Pop7", "Pop9"), t=0, method="random")
#' 
#' m <- gl.nhybrids(testset.gl, outfile="nhyb.txt", 
#' p0=NULL, p1=NULL, 
#' nhyb.directory="C:/workspace/R_analysis/NewHybsPC",
#' BurnIn=100,
#' sweeps=1000,
#' v=3)
#' }

gl.nhybrids <- function(gl, outfile="nhyb.txt", outpath=tempdir(),
                    p0=NULL, p1=NULL, 
                    t=0, 
                    method="random",
                    nhyb.directory=NULL,
                    BurnIn=10000,
                    sweeps=10000,
                    GtypFile = "TwoGensGtypFreq.txt",
                    AFPriorFile = NULL,
                    PiPrior = "Jeffreys",
                    ThetaPrior = "Jeffreys",
                    v=2) {

  outfile <- file.path(outpath, outfile)
  if(class(gl)!="genlight") {
    cat("Fatal Error: genlight object required for gl2nhyb.r!\n"); stop()
  }
  if (v > 0) {cat("Starting gl.nhybrids: Assigning individual to hybrid categories\n")}

    gl.tmp <- gl
    thold<-t
    
  # Housekeeping on the paths
    if (!is.null(nhyb.directory)) {
       nhyb.directory.win <- gsub("/","\\\\",nhyb.directory)
       wd.hold <- getwd()
       wd.hold.win <- gsub("/","\\\\",wd.hold)
    }

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
        if (!is.na(is.fixed(gl2$frequency[i],gl2$frequency[i+1], tloc=thold))) {
          if (is.fixed(gl2$frequency[i],gl2$frequency[i+1], tloc=thold)) {
            fixed.loci[i] <- as.character(gl2$locus[i])
          }
        }
      }

      # Remove the NAs
      fixed.loci <- fixed.loci[!is.na(fixed.loci)]

      # If no loci remain, set the data matrix to be the original matrix
      if (length(fixed.loci) == 0) {
        cat("  No fixed differences between parental populations \n")
        if(method=="random") {cat("    Selecting ca 200 random loci\n")}
        if (method=="avgpic") {cat("    Selecting 200 loci with most information content (AvgPIC)\n")}
        tmp <- gl.subsample.loci(gl, 200, method=method)
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
      if(method=="random") {cat("    Selecting ca 200 random loci\n")}
      if (method=="avgpic") {cat("    Selecting 200 loci with most information content (AvgPIC)\n")}
      tmp <- gl.subsample.loci(gl, 200, method=method)
      gl2 <- as.matrix(tmp)
      gl.reduced <- NULL
      flag <- "onepar"
    }

    # PROCESS AS FOLLOWS IF NO PARENTAL POPULATION IS SPECIFIED
    if (is.null(p0) & is.null(p1)) {
      cat("  No parental population specified \n")
      if(method=="random") {cat("    Selecting ca 200 random loci\n")}
      if (method=="avgpic") {cat("    Selecting 200 loci with most information content (AvgPIC)\n")}
      tmp <- gl.subsample.loci(gl, 200, method=method)
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


# Run New Hybrids   
    
    if (!is.null(nhyb.directory)) {
      
      #nhyb.directory.win<- gsub("/","\\\\",nhyb.directory)
      cp <- paste("copy nhyb.txt", nhyb.directory.win)
      shell(cp)
      
      # Shift default directory
      #wd.hold <- getwd()
      #wd.hold.win <- gsub("/","\\\\",wd.hold)
      
      setwd(nhyb.directory)
      aaa <- getwd()
      cat(aaa)
      
      # Set the parameters to conform with NewHybrids input
      if (GtypFile == "TwoGensGtypFreq.txt") {GtypFile <- "0"}
      if (is.null(AFPriorFile)) {AFPriorFile <- "0"}
      
      if (PiPrior == "Jeffreys" || PiPrior == "jeffreys") {
        PiPrior <- "0"
      } else if (PiPrior == "Uniform" || PiPrior == "uniform") {
        PiPrior <- "1"
      }  else {
        cat("Fatal Error: PiPrior parameter must be Jeffreys or Uniform\n"); stop()
      }
      if (ThetaPrior == "Jeffreys" || ThetaPrior == "jeffreys" ) {
        ThetaPrior <- "0"
      } else if (ThetaPrior == "Uniform" || ThetaPrior == "uniform") {
        ThetaPrior <- "1"
      }  else {
        cat("Fatal Error: ThetaPrior parameter must be Jeffreys or Uniform\n"); stop()
      }
      rand1 <- sample(1:10,1)
      rand2 <- sample(11:20,1)
      
      # Windows batch file to run new hybrids
      #  (
      #    echo nhyb.txt # Name of the data file
      #    echo 0        # Name of the file containing the genotype frequency classes
      #    echo 0        # Name of the file containing the prior allele frequency information
      #    echo 3 8      # two small numbers as seeds
      #    echo 0        # prior type for pi
      #    echo 0        # prior type for theta
      #    echo 10000    # burnin sweeps
      #    echo 100000   # sweeps
      #  ) | NewHybrids_PC_1_1_WOG.exe
      
      # Create the batch file
      sink("nhyb.cmd")
      cat("(\n")
      cat("echo nhyb.txt\n")
      cat("echo",GtypFile,"\n")
      cat("echo",AFPriorFile,"\n")
      cat("echo",rand1, rand2,"\n")
      cat("echo",PiPrior,"\n")
      cat("echo",ThetaPrior,"\n")
      cat("echo",BurnIn,"\n")
      cat("echo",sweeps,"\n")
      cat(") | NewHybrids_PC_1_1_WOG.exe")
      sink()

  # Run New Hybrids
      
      shell("nhyb.cmd")
      
  # Add in individual labels
      tbl <- read.table("aa-PofZ.txt", stringsAsFactors = FALSE)
      names(tbl) <- tbl[1,]
      tbl <- tbl[-1,-1]
      tbl <- cbind(indNames(gl),pop(gl),tbl)
      names(tbl)[1] <- "id"
      names(tbl)[2] <- "pop"
      write.csv(tbl,file="aa-pofZ.csv",row.names=FALSE)

  # Transfer files to default directory and housekeeping
      cp <- paste("copy aa-LociAndAlleles.txt", wd.hold.win); shell(cp)
      cp <- paste("del aa-LociAndAlleles.txt"); shell(cp) 
      cp <- paste("copy aa-ProcessedPriors.txt", wd.hold.win); shell(cp)
      cp <- paste("del aa-ProcessedPriors.txt"); shell(cp) 
      cp <- paste("copy aa-Pi.hist", wd.hold.win); shell(cp)
      cp <- paste("del aa-Pi.hist"); shell(cp)
      cp <- paste("copy aa-PofZ.csv", wd.hold.win); shell(cp)
      cp <- paste("del aa-PofZ.csv"); shell(cp)
      cp <- paste("del aa-PofZ.txt"); shell(cp)
      cp <- paste("copy aa-Theta.hist", wd.hold.win); shell(cp)
      cp <- paste("del aa-Theta.hist"); shell(cp)
      
      cp <- paste("del nhyb.txt"); shell(cp)
      cp <- paste("del nhyb.cmd"); shell(cp)

      setwd(wd.hold)
      
    }
    
    if (v > 0) {cat("gl.nhybrids Completed\n")}

    return(gl.reduced)

}

