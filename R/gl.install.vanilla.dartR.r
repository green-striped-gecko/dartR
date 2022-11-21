#' Installs all required packages for using all functions
#' available in dartR
#'
#' The function compares the installed packages with the the currently available
#' ones on CRAN. Be aware this function only works if a version of dartR is
#' already installed on your system. You can choose if you also want to have a
#' specific version of dartR installed ('CRAN', 'master', 'beta' or 'dev' ). 'master'
#' , 'beta' and 'dev' are installed from Github. Be aware that the dev version from github is
#'  not fully tested and most certainly will contain untested functions.
#' @param flavour The version of R you want to install. If NULL
#' then only packages needed for the current version will be installed. If
#' 'CRAN' current CRAN version will be installed. 'master' installs the GitHub
#' master branch, 'beta' installs the latest stable version, and 'dev' installs 
#' the experimental development branch from
#' GitHub [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log ; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @return Returns a message if the installation was successful/required.
#' @export
#' @importFrom stringr str_trim
#' @rawNamespace import(crayon, except = c(chr, '%+%'))
#' @importFrom utils installed.packages install.packages available.packages
#' @author Custodian: Bernd Gruber (Post to
#' \url{https://groups.google.com/d/forum/dartr})

gl.install.vanilla.dartR <- function(flavour = NULL,
                                     verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Josh",
                     verbosity = verbose)
    
    # ERROR CHECKING
    
    if(is.null(flavour)){
        flavour <- "CRAN"
    }
    
    pkg <- "devtools"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      cat(error(
        "Package",
        pkg,
        " needed for this function to work. Please install it.\n"
      ))
      return(-1)
    }
    
    pkg <- "stringr"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      cat(error(
        "Package",
        pkg,
        " needed for this function to work. Please install it.\n"
      ))
      return(-1)
    }
    
    # DO THE JOB
    
    err <- NULL
    if (!is.null(flavour)) {
        if (flavour == "CRAN") {
            detach("package:dartR", unload = TRUE)
            if (verbose >= 2) {
                cat(report("  Installing dartR from CRAN (latest version)\n"))
            }
            install.packages("dartR")
        }
        if (flavour == "master") {
            if (verbose >= 2) {
                cat(report("  Installing dartR from github (master)\n"))
            }
            detach("package:dartR", unload = TRUE)
            devtools::install_github("green-striped-gecko/dartR",
                                     ref = "master",
                                     dependencies = TRUE)
        }
        if (substr(flavour, 1, 3) == "dev") {
            detach("package:dartR", unload = TRUE)
            if (verbose >= 2) {
                cat(report("  Installing dartR from github (dev)\n"))
            }
            devtools::install_github("green-striped-gecko/dartR",
                                     ref = flavour,
                                     dependencies = TRUE)
        }
        if (substr(flavour, 1, 4) == "beta") {
            detach("package:dartR", unload = TRUE)
            if (verbose >= 2) {
                cat(report("  Installing dartR from github (beta)\n"))
            }
            devtools::install_github("green-striped-gecko/dartR",
                                     ref = flavour,
                                     dependencies = TRUE)
        }
    }
    
    # check all depends #check all imports check all suggests
     ip <- installed.packages()[, "Package"]
     suggests <-
        stringr::str_trim(unlist(strsplit(
            packageDescription("dartR", fields = "Suggests"), ","
        )))
     toinstall <- suggests[!suggests %in% ip]
     available <- available.packages()[, "Package"]
     toinsav <- toinstall[toinstall %in% available]
     
     if(verbose>=2){
       cat(report("  Installing package ggtern from provisional GitHub Repository\n"))
     }
       
     #Installing package ggtern from personal github repository because it seems
     #that this package is not maintained anymore
     devtools::install_github("mijangos81/ggtern")
     
    
    if (length(toinsav) > 0) {
        if (verbose >= 2) {
            cat(report(
                    "  The following packages will be installed (also available on CRAN):\n"
                )
            )
        }
        
        if (verbose >= 2) {
            cat("  ",paste0("  ",toinsav, "\n"))
        }
        for (ii in 1:length(toinsav)) {
            install.packages(toinsav[ii])
            
            if (verbose >= 2) {
                cat(report(paste(
                    "  Package:", toinsav[ii], "installed.\n"
                )))
            }
        }
        
        if (verbose >= 2) {
            cat(report(
                    "  All required packages are now installed. If still errors, update them using "
                )
            )
            cat(code("update.packages()"))
        }
        
        if (verbose >= 2) {
            cat(report("  ",paste(
                    "\n  You have installed dartR",
                    packageVersion("dartR")
                )
            ))
        }
        
        if (verbose >= 2) {
            cat(report("\nHave fun using Vanilla dartR!\n"))
        }
        # package installed but not available
    } else if (length(toinstall) > 0) {
        if (verbose >= 2) {
            cat(warn(
                    "  Warning: The following packages need to be installed, but are not available via CRAN.\n  Look to github or other repositories such as bioconductor.\n"
                )
            )            
            cat("  ",paste0(toinstall,"\n"))
        }
        
        if (verbose >= 2) {
            cat(warn(
                    "  Warning: If the packages qvalue and/or SNPrelate are missing run the code and then gl.install.vanilla.dartR() again\n"
                )
            )
        }
        
        if (verbose >= 2) {
            cat(warn(
                paste(
                    "    install.packages('devtools') \n    library(devtools)\n    install.packages('BiocManager')\n    BiocManager::install(c('SNPRelate', 'qvalue'))\n"
                )
            ))
        }
        
        err <- TRUE
    } else {
        if (verbose >= 2) {
            cat(
                report(
                    "  All required packages are already installed.\n  If there are still errors you might need to update them using "
                )
            )
        }
        
        if (verbose >= 2) {
            cat(code("update.packages()"))
            cat("\n")
        }
        
        # if (verbose) {
        #     cat(report(
        #         paste(
        #             "  You have installed dartR:",
        #             packageVersion("dartR"),"\n"
        #         )
        #     ))
        # }
    }
    
    if (verbose>= 2 & is.null(err) & !is.null(flavour)) {
        if (flavour != "CRAN") {
            fl <-paste0("Github [", flavour, "]")
        } else {
            fl <-"CRAN"
            cat(report(
                paste(
                    "  You have installed dartR",
                    packageVersion("dartR"),
                    "from",
                    fl,
                    "\n"
                )
            ))
        }
        
        
        if (verbose>= 2) {
            message(report("\n  Have fun using Vanilla dartR!\n"))
        }
        
        # FLAG SCRIPT END
        
        if (verbose > 0) {
            cat(report("Completed:", funname, "\n"))
        }
        
    }
    library(dartR) 
    
}
