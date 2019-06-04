
gl.play.history <- function(x, history=NULL)
{
  
  if (is.null(history)) hist2 <- x@other$history
  if (is.numeric(history)) hist2 <- x@other$history[history]
  if (is.list(history))  hist2 <- history
  for (i in 1:length(hist2))
  {
    glhist <- hist2[[i]]
    narg = length(glhist)
    ll <- list()
    ll[1:(narg - 1)] <- glhist[2:narg]
    names(ll) <- names(glhist[2:narg])
    #check if gl.read.dart or history on different gl
    if (as.character(glhist[[1]])!="gl.read.dart") {
      if (i==1) ll[[1]] <- x else ll[[1]]<- gout
    }
    #run history one by one
    gout <- do.call(as.character(glhist)[1], ll)
  }
  
  if (as.character(hist2[[1]][[1]]) == "gl.read.dart")
    gout@other$history <- hist2
  else
    gout@other$history <- c(x@other$history, hist2)
  gout
}
