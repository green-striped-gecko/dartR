#' This functions installs all required packages for using all functions available in dartR
#'  
#' The function compares the installed packages with the the currently available ones on cran. Be aware this function only works if a version of dartR is already installed on your system. You can choose if you also want to have a specific version of dartR installed ("CRAN", "master" or "dev" ). "master" and "dev" are installed from Github. Be aware the dev version from github is not fully tested and most certainly will contain untested functions.
#' @param flavour -- If and which version of R you want to install. If NULL then only needed packages for the current version will be installed. If "CRAN" current CRAN version will be installed. "master" installs the GitHub master branch and "dev" installs the experimental development branch from GitHub.
#' @param verbose returns information on packages and dartR versions
#'  @return returns a message if the installation was successful/required
#' @export
#' @importFrom stringr str_trim
#' @importFrom devtools install_github
#' @rawNamespace import (crayon, except = c(chr, '%+%') )
#' @importFrom utils installed.packages install.packages available.packages
#' @author Bernd Gruber (bugs? Post to \url{https://groups.google.com/d/forum/dartr})

gl.install.vanilla.dartR <- function(flavour=NULL, verbose=TRUE)
{
  pkg <- "devtools"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop("Package ",pkg," needed for this function to work. Please install it using install.packages('devtools').") }   
  
  pkg <- "stringr"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop("Package ",pkg," needed for this function to work. Please install it.") }   
  
  pkg <- "crayon"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop("Package ",pkg," needed for this function to work. Please install it.") }   
  err <- NULL
  report <- crayon::green
  warn <- crayon::red
  note <- crayon::cyan
  code <- crayon::blue
  important <- crayon::cyan$bold
  
  if (!is.null(flavour))
  {
    if (flavour=="CRAN"){
      detach("package:dartR", unload=TRUE)
      if (verbose) message(report("Installing dartR from CRAN (latest version)"))
      install.packages("dartR")
    }
    
    if (flavour=="master"){
      if (verbose) message(report("Installing dartR from github (master)"))
      detach("package:dartR", unload=TRUE)
      devtools::install_github("green-striped-gecko/dartR", ref="master", dependencies = TRUE)

    }
    if (substr(flavour,1,3)=="dev"){
      detach("package:dartR", unload=TRUE)
      if (verbose) message(report("Installing dartR from github (dev)"))
      devtools::install_github("green-striped-gecko/dartR", ref=flavour, dependencies = TRUE)
      
    }
  }
  
  
  
  
  ip <- installed.packages()[,"Package"]
  #check all depends #check all imports
  #check all suggests
  suggests <- stringr::str_trim(unlist(strsplit(packageDescription("dartR", fields = "Suggests"),",") ))
  toinstall <- suggests[!suggests %in% ip]
  available <- available.packages()[,"Package"]
  toinsav <- toinstall[toinstall %in% available]
  
  
  
  if (length(toinsav)>0)  {
    if (verbose) message(report("The following packages will be installed (and are available on Cran):\n"))
    if (verbose) message(report(paste0(toinsav,"\n")))
    for (ii in 1:length(toinsav))
    {
      
      install.packages(toinsav[ii])
      if (verbose) message(report(paste("Package:",toinsav[ii],"installed.\n")))
    }
    if (verbose) message(report("All required packages are now installed. If there are still errors you might need to update them using"))
    if (verbose) message(code("update.packages()"))  
    if (verbose) message(report(paste("\nYou have installed dartR",packageVersion("dartR"))))
    if (verbose) message(note("\nHave fun using Vanilla dartR!\n"))    
    
    
  } else  if(length(toinstall)>0)  #package installed but not available
    { if (verbose) message(warn("The following packages need to be installed, but are not available via Cran. You need to try and find them on github or other repositories such as bioconductor."))
    if (verbose) message(warn(paste0(toinstall,"\n")))
    if (verbose) message(important("If the packages qvalue and/or SNPrelate are missing run the following lines of code and then gl.install.vanilla.dartR() agian\n"))
    if (verbose) message(important(paste("install.packages('devtools') \nlibrary(devtools)\ninstall.packages('BiocManager')\nBiocManager::install(c('SNPRelate', 'qvalue'))")))
  err <- TRUE
  } else {
    if (verbose) message(report("All required packages are already installed.\n If there are still errors you might need to update them using"))
    if (verbose) message(code("update.packages()"))  
    
    if (verbose) message(report(paste("You have installed dartR:",packageVersion("dartR"))))
    if (verbose) message(note("\nHave fun using Vanilla dartR!\n"))    
  }
  
  if (verbose & is.null(err) & !is.null(flavour)) {
    if (flavour!="CRAN") fl=paste0("Github [",flavour,"]") else fl="CRAN"
      message(report(paste("You have installed dartR",packageVersion("dartR"),"from", fl,".\n")))
      
}

}