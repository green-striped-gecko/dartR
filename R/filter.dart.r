#' Filter function to facilitate analysing of dart data [deprecated use explicit filter functions gl.filter...]
#' 
#' @param ... vectors of true false created via comparisons between metrics (see examples)
#' @return returns the combined index of filters. this function is unusual as it has no predefined parameters. The idea is here that you can provide any type of filter (and index of true false that is used to filter the loci of your genlight object afterwards. This functions combines your filter and provides a plot on their single and combined effect. see example
#' @export
#' @examples
#' \dontrun{
#' #gl is a genlight object created with the read.dart and dart2genlight functions 
#' index.repro <- gl@other$loc.metrics[,"RepAvg"] > 0.98
#' index.callrate <-  gl@other$loc.metrics[,"CallRate"] > 0.90
#' index.highhet <- fox.gl.keep@other$loc.metrics[,"FreqHets"] <0.75
#' index.comb <- filter.dart(index.repro, index.callrate, index.coverage)
#' }


filter.dart<- function(...)
{
  #library(knitr, quietly = T)
  filnames <-   deparse(substitute(list(...)))
  
  filnames <-  gsub("[()]","",filnames)
  filnames <-  gsub(" ","",filnames)
  filnames <- gsub("list","", filnames)
  filters <- eval(eval(quote(substitute(list(...)))))
  nfilters <- length(filters)
  if (nfilters<2) stop("Only one filter provided, so combination does not make sense!")
  #filters <- gsub("list( | ")
  names(filters) <- unlist(strsplit( filnames,","))
  lenfil <- unlist(lapply(filters, length))
  tabfilters <- sapply(filters, table)
  tabfilters <- do.call(rbind, list(tabfilters))
  tabfilters <- cbind(tabfilters, sum=lenfil, freq=round(tabfilters[,2]/lenfil,3))
  #print(kable(tabfilters))
  nsnp <- min(lenfil)
  dffil <- do.call(rbind, filters)
  index.all <- colSums(dffil)==nfilters
  index.comb <- rbind(index.all,dffil)
  index.combout <- index.comb
  index.combout[1,] <- index.combout[1,]*2
  rs <-rowSums(index.comb)
  image(1:ncol(index.comb), 1:nrow(index.comb),t(index.combout), col=c("white", "red", "orange"), axes=F, xlab="loci", ylab="",main="Filters for each loci")
  axis(2, at= 1:nrow(index.comb),labels= gsub("index.","",row.names(index.comb)), las=2)
  text(nsnp/2,1:nrow(index.comb), rs)
  axis(1)#, at=c(1,round(nsnp/2), nsnp), label=c(1,round(nsnp/2), nsnp))
  lines(c(1,nsnp),c(rep(nrow(index.comb),2))-0.5)
  box()
  #filters
  data.frame(t(index.comb))
}
