#' @name gl.nhybrids
#' @title Creates an input file for the program NewHybrids and runs it if
#'  NewHybrids is installed
#' @description
#' This function compares two sets of parental populations to identify loci that
#' exhibit a fixed difference, returns an genlight object with the reduced
#' data, and creates an input file for the program NewHybrids using the top 200
#' (or hard specified loc.limit) loci. In the absence of two identified
#' parental populations, the script will select a random set 200 loci only
#' (method='random') or the first 200 loci ranked on information content
#' (method='AvgPIC').
#'
#' A fixed difference occurs when a SNP allele is present in all individuals
#' of one population and absent in the other. There is provision for setting
#' a level of tolerance, e.g. threshold = 0.05 which considers alleles present
#' at greater than 95% in one population and less than 5% in the other to be
#' a fixed difference. Only the 200 loci are retained, because of limitations
#' of NewHybids.
#'
#' If you specify a directory for the NewHybrids executable file, then the
#' script will create the input file from the SNP data then run NewHybrids. If
#' the directory is set to NULL, the execution will stop once the input file
#' (default='nhyb.txt') has been written to disk. Note: the executable option
#' will not work on a Mac; Mac users should generate the NewHybrids input file
#' and run this on their local installation of NewHybrids.
#'
#' Refer to the New Hybrids manual for further information on the parameters to
#' set
#' -- http://ib.berkeley.edu/labs/slatkin/eriq/software/new_hybs_doc1_1Beta3.pdf
#'
#' It is important to stringently filter the data on RepAvg and CallRate if
#' using the random option. One might elect to repeat the analysis
#' (method='random') and combine the resultant posterior probabilities should
#' 200 loci be considered insufficient.
#'
#' The F1 individuals should be homozygous at all loci for which the parental
#' populations are fixed and different, assuming parental populations have been
#' specified. Sampling errors can result in this not being the case, especially
#' where the sample sizes for the parental populations are small. Alternatively,
#' the threshold for posterior probabilities used to determine assignment
#' (pprob) or the definition of a fixed difference (threshold) may be too lax.
#' To assess the error rate in the determination of assignment of F1
#' individuals, a plot of the frequency of homozygous reference, heterozygotes
#' and homozygous alternate (SNP) can be produced by setting plot=TRUE (the
#' default).
#'
#' @param gl Name of the genlight object containing the SNP data [required].
# @param outfile Name of the file that will be the input file for NewHybrids
# [default nhyb.txt].
#' @param outpath Path where to save the output file [default tempdir()].
#' @param p0 List of populations to be regarded as parental population 0
#' [default NULL].
#' @param p1 List of populations to be regarded as parental population 1
#' [default NULL].
#' @param threshold Sets the level at which a gene frequency difference is
#' considered to be fixed [default 0].
#' @param plot If TRUE, a plot of the frequency of homozygous reference,
#' heterozygotes and homozygous alternate (SNP) is produced for the F1
#' individuals
#' [default TRUE, applies only if both parental populations are specified].
#' @param plot_theme User specified theme [default theme_dartR()].
#' @param plot_colors Vector with two color names for the borders and fill
#' [default two_colors].
#' @param pprob Threshold level for assignment to likelihood bins
#' [default 0.95, used only if plot=TRUE].
#' @param method Specifies the method (random or AvgPIC) to select 200 loci for
#'  NewHybrids [default random].
#' @param nhyb.directory Directory that holds the NewHybrids executable file
#' e.g. C:/NewHybsPC [default NULL].
#' @param BurnIn Number of sweeps to use in the burn in [default 10000].
#' @param sweeps Number  of  sweeps  to  use  in  computing  the  actual Monte
#' Carlo averages [default 10000].
#' @param GtypFile Name of a file containing the genotype frequency classes
#' [default TwoGensGtypFreq.txt].
#' @param AFPriorFile Name of the file containing prior allele frequency
#' information [default NULL].
#' @param PiPrior Jeffreys-like priors or Uniform priors for the parameter pi
#' [default Jeffreys].
#' @param ThetaPrior Jeffreys-like priors or Uniform priors for the parameter
#' theta [default Jeffreys].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @return The reduced genlight object, if parentals are provided; output of
#'  NewHybrids is saved to the working directory.
#' @export
#' @importFrom MASS write.matrix
#' @references Anderson, E.C. and Thompson, E.A.(2002). A model-based method for identifying 
#' species hybrids using multilocus genetic data. Genetics. 160:1217-1229. 
#' @author Custodian: Arthur Georges -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \dontrun{
#' m <- gl.nhybrids(testset.gl, 
#' p0=NULL, p1=NULL,
#' nhyb.directory='D:/workspace/R/NewHybsPC', # Specify as necessary
#' outpath="D:/workspace",  # Specify as necessary, usually getwd() [= workspace]
#' BurnIn=100,
#' sweeps=100,
#' verbose=3)
#' }

gl.nhybrids <- function(gl,
#                       outfile = "nhyb.txt",
                        outpath = tempdir(),
                        p0 = NULL,
                        p1 = NULL,
                        threshold = 0,
                        method = "random",
                        plot = TRUE,
                        plot_theme = theme_dartR(),
                        plot_colors = two_colors,
                        pprob = 0.95,
                        nhyb.directory = NULL,
                        BurnIn = 10000,
                        sweeps = 10000,
                        GtypFile = "TwoGensGtypFreq.txt",
                        AFPriorFile = NULL,
                        PiPrior = "Jeffreys",
                        ThetaPrior = "Jeffreys",
                        verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbosity = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(gl, verbose = verbose)
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    if (!(method == "random" |
          method == "AvgPic" |
          method == "avgPic" | method == "AvgPIC")) {
        cat(warn(
            "  Warning: method must be either 'random' or AvgPic, set to 'random'\n"
        ))
        method <- "random"
    }
    
    if (pprob < 0 | pprob > 1) {
        cat(
            warn(
                "  Warning: threshold posterior probability for assignment, pprob, must be between 0 and 1, typically close to 1, set to 0.99\n"
            )
        )
        pprob <- 0.99
    }
    
    # Housekeeping on the outfile specifications
    outfile <- "nhyb.txt"
    outfile <- file.path(outpath, outfile)
    outfile.win <- gsub("/", "\\\\", outfile)
    if (!is.null(nhyb.directory)) {
      nhyb.directory.win <- gsub("/", "\\\\", nhyb.directory)
      wd.hold <- getwd()
      wd.hold.win <- gsub("/", "\\\\", wd.hold)
    }
    
    # DO THE JOB
    
    gl.tmp <- gl
    thold <- threshold
    loc.limit <- 200
    
    # PROCESS AS FOLLOWS IF BOTH PARENTAL POPULATIONS ARE SPECIFIED
    if (!is.null(p0) & !is.null(p1)){
        if (verbose >= 3) {
            cat(report("  Both parental populations have been specified \n"))
        }
        
        # Replace those names with Z0 for Parental Population 0, and z1 for Parental Population 1
        indNames(gl.tmp) <-
            replace(indNames(gl.tmp), is.element(pop(gl.tmp), p0), "z0")
        indNames(gl.tmp) <-
            replace(indNames(gl.tmp), is.element(pop(gl.tmp), p1), "z1")
        
        # Error checks
        if (length(indNames(gl.tmp)[indNames(gl.tmp) == "z0"]) < 1 |
            length(indNames(gl.tmp)[indNames(gl.tmp) == "z1"]) < 1) {
            stop(
                error(
                    "Fatal Error [gl2nhyb]: One or both of two specified parental populations contains no individuals\n"
                )
            )
        }
        
        # Create a vector containing the flags for the parental population
        par.names <- indNames(gl.tmp)
        par.names <-
            replace(par.names, (par.names != "z0" &
                                    par.names != "z1"), " ")
        
        # Discard non-parental populations
        gl.tmp <-
            gl.tmp[(indNames(gl.tmp) == "z0" |
                        indNames(gl.tmp) == "z1"), ]
        pop(gl.tmp) <- indNames(gl.tmp)
        
        # Reformat the data broken down by population and locus, and calculate allele frequencies
        gl2 <- gl.percent.freq(gl.tmp, verbose = verbose)
        
        # IDENTIFY LOCI WITH FIXED DIFFERENCES BETWEEN P0 AND P1
        if (verbose >= 3) {
            cat(
                report(
                    "  Identifying loci with fixed difference between parental stocks\n"
                )
            )
        }
        
        # Establish a vector to hold the loci
        npops <- 2
        nloci <- nlevels(gl2$locus)
        fixed.loci <- NA
        length(fixed.loci) <- nloci * 2
        
        # Cycle through the data to identify the fixed loci
        for (i in seq(1, nloci * 2, 2)) {
            if (as.character(gl2$locus[i]) != as.character(gl2$locus[i + 1])) {
                cat(warn(
                    "  Warning: Loci do not agree for the is.fixed comparison\n"
                ))
            }
            if (!is.na(is.fixed(gl2$frequency[i], gl2$frequency[i + 1], tloc = thold))) {
                if (is.fixed(gl2$frequency[i], gl2$frequency[i + 1], tloc = thold)) {
                    fixed.loci[i] <- as.character(gl2$locus[i])
                }
            }
        }
        
        # Remove the NAs
        fixed.loci <- fixed.loci[!is.na(fixed.loci)]
        nloci <- length(fixed.loci)
        if (nloci == 0) {
            flag <- "bothparnonefixed"
        } else {
            flag <- "bothpar"
        }
        
        # Report on the analysis
        
        if (verbose >= 3) {
            cat(report("  No. of fixed loci identified:", nloci, "\n"))
        }
        if (nloci > loc.limit) {
            if (verbose >= 3) {
                cat(
                    report(
                        "  Selecting",
                        loc.limit,
                        "loci showing fixed differences between parentals at random\n"
                    )
                )
            }
            gl.fixed.used <-
                gl.subsample.loci(gl,
                                  loc.limit,
                                  method = "random",
                                  verbose = 0)
            gl.fixed.all <- gl[, (locNames(gl) %in% fixed.loci)]
            gl.fixed.all@other$loc.metrics <-
                gl@other$loc.metrics[(locNames(gl) %in% fixed.loci), ]
            gl2nhyb <- gl.fixed.used
        } else {
            if (method == "random") {
                if (verbose >= 3) {
                    cat(
                        report(
                            "  Selecting",
                            nloci,
                            "loci showing fixed differences between parentals, supplementing with",
                            loc.limit - nloci,
                            "other loci selected at random\n"
                        )
                    )
                }
                gl.fixed.all <-
                    gl[, (locNames(gl) %in% fixed.loci)]
                gl.fixed.used <- gl.fixed.all
                tmp <-
                    gl.subsample.loci(gl,
                                      200 - nloci,
                                      method = "random",
                                      verbose = 0)
                gl2nhyb <- cbind(gl.fixed.used, tmp)
                gl2nhyb@other$loc.metrics <-
                    gl@other$loc.metrics[locNames(gl) %in% locNames(gl2nhyb), ]
            } else {
                if (verbose >= 3) {
                    cat(
                        report(
                            "  Selecting",
                            nloci,
                            "loci showing fixed differences between parentals, supplementing with",
                            loc.limit - nloci,
                            "other loci selected based on AvgPic\n"
                        )
                    )
                }
                gl.fixed.all <-
                    gl[, (locNames(gl) %in% fixed.loci)]
                gl.fixed.used <- gl.fixed.all
                tmp <-
                    gl.subsample.loci(gl,
                                      loc.limit - nloci,
                                      method = "AvgPic",
                                      verbose = 0)
                gl2nhyb <- cbind(gl.fixed.used, tmp)
                gl2nhyb@other$loc.metrics <-
                    gl@other$loc.metrics[locNames(gl) %in% locNames(gl2nhyb), ]
            }
        }
    }  # Finished -- loc.limit selected loci in gl2nhyb
    
    # PROCESS AS FOLLOWS IF ONLY ONE PARENTAL POPULATION IS SPECIFIED
    if ((!is.null(p0) &
         is.null(p1)) || (is.null(p0) & !is.null(p1))) {
        if (verbose >= 3) {
            cat(report("  Only one parental population specified \n"))
        }
        
        # Replace those names with Z0 for Parental Population 0, and z1 for Parental Population 1
        if (!is.null(p0)) {
            indNames(gl.tmp) <-
                replace(indNames(gl.tmp), is.element(pop(gl.tmp), p0), "z0")
        }
        if (!is.null(p1)) {
            indNames(gl.tmp) <-
                replace(indNames(gl.tmp), is.element(pop(gl.tmp), p1), "z1")
        }
        
        # Error checks
        if (length(indNames(gl.tmp)[indNames(gl.tmp) == "z0"]) < 1 &
            length(indNames(gl.tmp)[indNames(gl.tmp) == "z1"]) < 1) {
            stop(
                error(
                    "Fatal Error [gl2nhyb]: Specified parental population contains no individuals\n"
                )
            )
        }
        
        # Create a vector containing the flags for the parental population
        par.names <- indNames(gl.tmp)
        par.names <-
            replace(par.names, (par.names != "z0" &
                                    par.names != "z1"), " ")
        
        if (method == "random" & verbose >= 3) {
            cat("    Selecting", loc.limit, "random loci\n")
        }
        if (method == "avgpic" & verbose >= 3) {
            cat("    Selecting",
                loc.limit,
                "loci with most information content (AvgPIC)\n")
        }
        gl2nhyb <-
            gl.subsample.loci(gl, loc.limit, method = method, verbose = 0)
        gl2nhyb@other$loc.metrics <-
            gl@other$loc.metrics[locNames(gl) %in% locNames(gl2nhyb), ]
        flag <- "onepar"
    }
    
    # PROCESS AS FOLLOWS IF NO PARENTAL POPULATION IS SPECIFIED
    if (is.null(p0) & is.null(p1)) {
        if (verbose >= 3) {
            cat("  No parental population specified \n")
        }
        
        if (method == "random" & verbose >= 3) {
            cat(report("  Selecting", loc.limit, "random loci\n"))
        }
        if (method == "avgpic" & verbose >= 3) {
            cat(
                report(
                    "  Selecting",
                    loc.limit,
                    "loci with most information content (AvgPIC)\n"
                )
            )
        }
        gl2nhyb <- gl.subsample.loci(gl, loc.limit, method = method, verbose = 0)
        gl2nhyb@other$loc.metrics <-
            gl@other$loc.metrics[locNames(gl) %in% locNames(gl2nhyb), ]
        flag <- "nopar"
    }
    
    # CREATE THE NEWHYBRIDS INPUT FILE
    if (verbose >= 3) {
        cat(report("\n  Converting data to NewHybrids format\n"))
    }
    gl2 <- as.matrix(gl2nhyb)
    gl2[gl2 == 2] <- 22
    gl2[gl2 == 1] <- 12
    gl2[gl2 == 0] <- 11
    gl2[is.na(gl2)] <- 0
    n.loci <- ncol(gl2)
    
    # Create sequential row number
    rownum <- seq(1:nrow(gl2))
    
    # Bind to the matrix
    gl2 <- data.frame(gl2)
    
    if (flag == "bothpar" ||
        flag == "onepar" || flag == "bothparnonefixed") {
        gl2 <- cbind(rownum, par.names, gl2)
        metarows <- 2
        if (verbose >= 3) {
            cat(report(
                "  Adding sequential number and flagging parental stock\n"
            ))
        }
    }
    if (flag == "nopar") {
        gl2 <- cbind(rownum, gl2)
        metarows <- 1
        if (verbose >= 3) {
            cat(report("  Adding sequential number\n"))
        }
    }
    
    # Clear row and column names
    rownames(gl2) <- NULL
    colnames(gl2) <- NULL
    
    # Output data file
    if (verbose >= 3) {
        cat("  Writing the NewHybrids input file", outfile, "\n")
        cat(c("    NumIndivs ", nrow(gl2), "\n"))
        cat(c("    NumLoci ", n.loci, " \n"))
        cat(c("    Digits 1\n  Format Lumped \n"))
    }
    sink(outfile)
    cat(c("NumIndivs ", nrow(gl2), "\n"))
    cat(c("NumLoci ", n.loci, " \n"))
    cat(c("Digits 1\nFormat Lumped \n"))
    MASS::write.matrix(gl2[, 1:(ncol(gl2))], sep = " ")
    sink()
    
##### IF AN EXECUTABLE DIRECTORY IS SPECIFIED
    OS <- Sys.info()[1]
    
    ##### IF WINDOWS ON PC
    
    if (OS == "Windows"){
      if (verbose >= 2){
        cat(report("  Windows operating system\n"))
      }  
    
    # Checking for directories and files
    if(nchar(outfile.win) > nchar(gsub(" ","",outfile.win))){
      stop(error("Fatal Error: NewHybrids will not accept file or directory names that contain spaces\n"))
    }
    if(!is.null(nhyb.directory)){
      if(nchar(nhyb.directory) > nchar(gsub(" ","",nhyb.directory))){
        stop(error("Fatal Error: NewHybrids will not accept file or directory names that contain spaces\n"))
      }
      if(!dir.exists(nhyb.directory)){
        stop(error("Fatal Error: Directory for the NewHybrids executable does not exist\n"))
      }
      if (outpath == nhyb.directory){
          stop(error("Fatal Error: Directory for the NewHybrids executable cannot be the same as",
            "the directory to receive the output\n"))
      }
    }
    
    # Check the installation of New Hybrids
    tmp1 <- file.exists(paste0(nhyb.directory.win,"/NewHybrids_PC_1_1_WOG.exe"))
    if(!tmp1){
        stop(error("Fatal Error: New Hybrids executable not found in",nhyb.directory.win,"; required\n"))
    }
    tmp2 <- file.exists(paste0(nhyb.directory.win,"/TwoGensGtypFreq.txt"))
    if(!tmp2){
        stop(error("Fatal Error: New Hybrids Genotype Frequency file not found in",nhyb.directory.win,"; required\n"))
    }
    if (verbose >= 2){
        if(tmp1 & tmp2){
            cat(report("  New Hybrids executable files found\n"))
        }
    }
    
    # Run New Hybrids
    
    if (!is.null(nhyb.directory)) {
        if (verbose >= 2) {
            cat(
                report(
                    "  Copying New Hybrids input file",
                    outfile,
                    "to",
                    nhyb.directory.win,
                    "\n"
                )
            )
            cat(report("  Passing control to New Hybrids executable\n"))
        }
        tmp <- file.copy(from=outfile.win, to=nhyb.directory.win, overwrite = TRUE)
        if (verbose >= 2) {
            if(tmp){
              cat(report("  .... success\n"))
            } else {
                cat(stop("  .... failed to copy nhyb.txt to",nhyb.directory.win,"-- check permissions\n"))
            }
        }
        
        setwd(nhyb.directory)
        
        # Set the parameters to conform with NewHybrids input
        if (GtypFile == "TwoGensGtypFreq.txt") {
            GtypFile <- "0"
        }
        if (is.null(AFPriorFile)) {
            AFPriorFile <- "0"
        }
        
        if (PiPrior == "Jeffreys" || PiPrior == "jeffreys") {
            PiPrior <- "0"
        } else if (PiPrior == "Uniform" ||
                   PiPrior == "uniform") {
            PiPrior <- "1"
        } else {
            stop(error(
                "Fatal Error: PiPrior parameter must be Jeffreys or Uniform\n"
            ))
        }
        if (ThetaPrior == "Jeffreys" ||
            ThetaPrior == "jeffreys") {
            ThetaPrior <- "0"
        } else if (ThetaPrior == "Uniform" ||
                   ThetaPrior == "uniform") {
            ThetaPrior <- "1"
        } else {
            stop(error(
                "Fatal Error: ThetaPrior parameter must be Jeffreys or Uniform\n"
            ))
        }
        rand1 <- sample(1:10, 1)
        rand2 <- sample(11:20, 1)
        
        # Create the batch file
        sink("nhyb.cmd")
        cat("(\n")
        cat("echo", outfile, "\n")
        cat("echo", GtypFile, "\n")
        cat("echo", AFPriorFile, "\n")
        cat("echo", rand1, rand2, "\n")
        cat("echo", PiPrior, "\n")
        cat("echo", ThetaPrior, "\n")
        cat("echo", BurnIn, "\n")
        cat("echo", sweeps, "\n")
        cat(") | NewHybrids_PC_1_1_WOG.exe")
        sink()
        
        # Run New Hybrids
        
        system("nhyb.cmd")
        
        # Add in individual labels
        tbl <-
            read.table("aa-PofZ.txt", stringsAsFactors = FALSE)
        names(tbl) <- tbl[1, ]
        tbl <- tbl[-1, -1]
        tbl <- cbind(indNames(gl), pop(gl), tbl)
        names(tbl) <-
            c("id", "pop", "P0", "P1", "F1", "F2", "F1xP0", "F1xP1")

        write.csv(tbl, file = "aa-PofZ.csv", row.names = FALSE)
        
        # Transfer files to default directory and housekeeping
        if (verbose == 2){
          cat(report("  Transferring output files to output directory",outpath,"\n"))
        }
        tmp <- file.copy(from="aa-LociAndAlleles.txt", to=outpath, overwrite = TRUE)
        tmp <- file.remove("aa-LociAndAlleles.txt")

        tmp <- file.copy(from="aa-ProcessedPriors.txt", to=outpath, overwrite = TRUE)
        tmp <- file.remove("aa-ProcessedPriors.txt")

        tmp <- file.copy(from="aa-Pi.hist", to=nhyb.directory.win, overwrite = TRUE)
        tmp <- file.remove("aa-Pi.hist")

        tmp <- file.copy(from="aa-PofZ.csv", to=outpath, overwrite = TRUE)
        tmp <- file.remove("aa-PofZ.csv")
        tmp <- file.remove("aa-PofZ.txt")

        tmp <- file.copy(from="aa-Theta.hist", to=outpath, overwrite = TRUE)
        tmp <- file.remove("aa-Theta.hist")
        
        tmp <- file.remove(basename(outfile))
        tmp <- file.remove("nhyb.cmd")

        setwd(wd.hold)
        
    }
    } ##### END WINDOWS BLOCK
    
    ##### FOR THE MAC -- Luis to add code
    
    # Use a sys.info to determine operating system?
    
      # 1. check that the necessary files exist in the executable directory
      # 2. Make sure there are no spaces in the user specified file or directory names -- stop, error
      # 3. Trap the user specifying the same directory for the executables as for the output
      # 2. transfer the nhyb.txt file to the executable directory
      # 3. Run the executable (no gui)
      # 4. Transfer the relevant output files to the user specified directory
      # 5. Delete those files from the executable directory, including any cmd file
    
    #### MAC CODE #####
    
    ##### IF MAC
    
    if (OS != "Windows"){
      if (verbose >= 2){
        cat(report("  Unix operating system\n"))
      }  
      
      outfile.mac  <- outfile
      nhyb.directory.mac <- nhyb.directory
      
      # Checking for directories and files
      if(!is.null(nhyb.directory)){
        if(!dir.exists(nhyb.directory)){
          stop(error("Fatal Error: Directory for the NewHybrids executable does not exist\n"))
        }
        if (outpath == nhyb.directory){
          stop(error("Fatal Error: Directory for the NewHybrids executable cannot be the same as",
                     "the directory to receive the output\n"))
        }
      }
      
      if (grepl("\\s", outfile.mac) | grepl("\\s", nhyb.directory)) {
        stop(
          error(
            "Fatal Error: The path to the executable for NewHybrids or the outfile name has spaces. Please move it to a path without spaces or choose a file name without spaces.\n"
          )
        )
      }
      
      
      if(nchar(outfile.mac) > nchar(gsub(" ","",outfile.mac))){
        stop(error("Fatal Error: NewHybrids will not accept filenames that contain spaces\n"))
      }
      
      # Check the installation of New Hybrids
      tmp1 <- file.exists(paste0(nhyb.directory.mac,"/newhybs"))
      if(!tmp1){
        stop(error("Fatal Error: New Hybrids executable not found in",nhyb.directory.mac,"; required\n"))
      }
      # tmp2 <- file.exists(paste0(nhyb.directory.mac,"/TwoGensGtypFreq.txt"))
      # if(!tmp2){
      #   stop(error("Fatal Error: New Hybrids Genotype Frequency file not found in",nhyb.directory.win,"; required\n"))
      # }
      if (verbose >= 2){
        # if(tmp1 & tmp2){
        if(tmp1){
          
          cat(report("  New Hybrids executable files found\n"))
        }
      }
      
      # Run New Hybrids
      
      if (!is.null(nhyb.directory)) {
        if (verbose >= 2) {
          cat(
            report(
              "  Copying New Hybrids input file",
              outfile,
              "to",
              nhyb.directory.mac,
              "\n"
            )
          )
          cat(report("  Passing control to New Hybrids executable\n"))
        }
        
        tmp <- file.copy(from=outfile.mac, to=nhyb.directory.mac, overwrite = TRUE)
        if (verbose >= 2) {
          if(tmp){
            cat(report("  .... success\n"))
          } else {
            cat(stop("  .... failed to copy nhyb.txt to",nhyb.directory.mac,"-- check permissions\n"))
          }
        }
        
        setwd(nhyb.directory)
        
        # Set the parameters to conform with NewHybrids input
        if (GtypFile == "TwoGensGtypFreq.txt") {
          GtypFile <- "0"
        }
        if (is.null(AFPriorFile)) {
          AFPriorFile <- "0"
        }
        
        rand1 <- sample(1:10, 1)
        rand2 <- sample(11:20, 1)
        
        system(paste(paste0(nhyb.directory.mac,"/newhybs"),
                     "--no-gui",
                     "--data-file",paste0(nhyb.directory.mac,"/",basename(outfile)),
                     "--seeds",rand1, rand2,
                     "--pi-prior", PiPrior,
                     "--theta-prior", ThetaPrior,
                     "--burn-in", BurnIn,
                     "--num-sweeps", sweeps
                     # "--gtyp-cat-file",GtypFile,
        ))
        
        # Add in individual labels
        tbl <-
          read.table("aa-PofZ.txt", stringsAsFactors = FALSE)
        names(tbl) <- tbl[1, ]
        tbl <- tbl[-1, -1]
        tbl <- cbind(indNames(gl), pop(gl), tbl)
        names(tbl) <-
          c("id", "pop", "P0", "P1", "F1", "F2", "F1xP0", "F1xP1")
        
        write.csv(tbl, file = "aa-PofZ.csv", row.names = FALSE)
        
        # Transfer files to default directory and housekeeping
        if (verbose == 2){
          cat(report("  Transferring output files to output directory",outpath,"\n"))
        }
        tmp <- file.copy(from="aa-LociAndAlleles.txt", to=outpath, overwrite = TRUE)
        tmp <- file.remove("aa-LociAndAlleles.txt")
        
        # tmp <- file.copy(from="aa-ProcessedPriors.txt", to=outpath, overwrite = TRUE)
        # tmp <- file.remove("aa-ProcessedPriors.txt")
        
        tmp <- file.copy(from="aa-Pi.hist", to=nhyb.directory.mac, overwrite = TRUE)
        tmp <- file.remove("aa-Pi.hist")
        
        tmp <- file.copy(from="aa-PofZ.csv", to=outpath, overwrite = TRUE)
        tmp <- file.remove("aa-PofZ.csv")
        tmp <- file.remove("aa-PofZ.txt")
        
        tmp <- file.copy(from="aa-Theta.hist", to=outpath, overwrite = TRUE)
        tmp <- file.remove("aa-Theta.hist")
        
        tmp <- file.remove(basename(outfile))
        
        tmp <- file.remove("aa-ThetaAverages.txt")
        tmp <- file.remove("aa-EchoedGtypFreqCats.txt")
        tmp <- file.remove("EchoedGtypData.txt")
        
        # tmp <- file.remove("nhyb.cmd")
        
        setwd(wd.hold)
        
      }
    } ##### END MAC BLOCK
    
    
    ##### Analyse the F1 genotypes
    
    if (flag == "bothpar" & plot == TRUE) {
        # Read in the results of the New Hybrids analysis
        F1.test <- read.csv(file = "aa-Pofz.csv")
        # Pull out results for F1 hybrids only, defined by posterior probability >= pprob
        F1.test <- F1.test[(F1.test$F1 >= pprob), ]
        # Use the id of the F1 hybrids to subset the genlight object containing the loci with fixed differences used in the analysis
        F1.only <-
            gl.keep.ind(
                gl.fixed.used,
                as.character(F1.test$id),
                mono.rm = FALSE,
                verbose = 0
            )
        # Invert the matrix of data for F1 individuals only
        mat <- t(as.matrix(F1.only))
        # Tally and plot the counts of homozygous reference, heterozygotes, and homozygous alternate (SNP)
        if (verbose >= 3) {
            cat("  Plotting genotypes of",
                length(as.character(F1.test$id)),
                "F1 individuals\n")
        }
        mat_plot <- as.data.frame(table(mat))
        colnames(mat_plot) <- c("Genotype","Count")
        
        print(
          ggplot(mat_plot, aes(x= Genotype,y = Count)) +
            geom_col(color = plot_colors[1], fill = plot_colors[2]) +
            plot_theme +
            theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()) +
            ggtitle("F1 Genotypes")
        )
        # Report the results
        if (verbose >= 3) {
            cat("  No. of F1 individuals:",
                nInd(F1.only),
                "[",
                indNames(F1.only),
                "]\n")
            cat(
                "  No. of loci with fixed differences used in the analysis:",
                nLoc(F1.only),
                "\n"
            )
            cat(
                "  No. of heterozygous loci for the F1s:",
                table(mat)["1"],
                "(",
                round(table(mat)["1"] * 100 / sum(table(mat)), 2),
                "%)\n"
            )
            cat(
                "  No. of homozygous reference loci (errors) for the F1s:",
                table(mat)["0"],
                "(",
                round(table(mat)["0"] * 100 / sum(table(mat)),
                      2),
                "%)\n"
            )
            cat(
                "  No. of homozygous alternate loci (errors) for the F1s:",
                table(mat)["2"],
                "(",
                round(table(mat)["2"] * 100 / sum(table(mat)),
                      2),
                "%)\n"
            )
        }
    }
    
    # if (verbose < 2) {sink()}
    if (verbose >= 3) {
        cat(report("  Results are stored in file aa-PofZ.csv\n"))
        cat(
            report(
                "  Returning data used by New Hybrids in generating the results, as a genlight object\n"
            )
        )
    }
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(gl2nhyb)
    
}
