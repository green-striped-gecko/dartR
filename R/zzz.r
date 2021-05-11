.onAttach <- function(...) {
  error <- crayon::red 
  warn <- crayon::yellow
  report <- crayon::green
  important <- crayon::cyan 
  wes_palette_Zissou1 <- c("#3B9AB2" ,"#78B7C5" ,"#EBCC2A" ,"#E1AF00" ,"#F21A00")
  packageStartupMessage(important("**** Welcome to dartR ****\n"))
  packageStartupMessage(report("Be aware that owing to CRAN requirements and compatibility reasons not all functions of the packages may run yet, as some dependencies could be missing. Hence for a most enjoyable experience we recommend to run the function "))
  packageStartupMessage(report("gl.install.vanilla.dartR()"))                 
  
  packageStartupMessage(report("This installs all missing and required packages for your version of dartR. \nFor citation information please use:"))
  packageStartupMessage(report("citation('dartR')"))
  packageStartupMessage(important("\n**** Have fun using dartR! ****"))
}

