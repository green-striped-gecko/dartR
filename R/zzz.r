.onAttach <- function(...) {

  report <- crayon::green
  note <- crayon::cyan
  code <- crayon::blue
  
  packageStartupMessage(note("**** Welcome to dartR ****\n"))
  packageStartupMessage(report("Be aware that owing to CRAN requirements and compatibility reasons not all functions of the packages may run yet, as some dependencies could be missing. Hence for a most enjoyable experience we recommend to run the function "))
  packageStartupMessage(code("gl.install.vanilla.dartR()"))                 
                 
  packageStartupMessage(report("This installs all missing and required packages for your version of dartR. \nFor citation information please use:"))
  packageStartupMessage(code("citation('dartR')"))
  packageStartupMessage(note("\n**** Have fun using dartR! ****"))
}
