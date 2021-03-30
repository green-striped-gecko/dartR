.onLoad <- function(libname, pkgname) {
  
 
  report <- crayon::green
  note <- crayon::cyan
  code <- crayon::blue
  
  message(note("**** Welcome to dartR ****\n"))
  message(report("Be aware that owing to CRAN requirements and compatibility reasons not all functions of the packages may run yet, as some dependencies could be missing. Hence for a most enjoyable experience we recommend to run the function "))
  message(code("gl.install.vanilla.dartR()"))                 
                 
  message(report("This installs all missing and required packages for your version of dartR. \nFor citation information please use:"))
  message(code("citation('dartR')"))
   message(note("\n**** Have fun using dartR! ****"))
}
