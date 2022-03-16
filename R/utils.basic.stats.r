#' Calculates mean observed heterozygosity, mean expected heterozygosity and Fis
#'  per population and Fst across all populations 
#' @param x A genlight object containing the SNP genotypes [required].
#' @return A list with with the statistics for each population
#' @export
#' @author Luis Mijangos (bugs? Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' out <- utils.basic.stats(platypus.gl)
#' @export

utils.basic.stats <- function(x) {
  
  pop.names <- popNames(x)
  sgl_mat <- lapply(seppop(x),as.matrix)
  pop.vec <- 1:nPop(x)
  n.ind <- table(pop(x))
  if(any(n.ind <= 1)){
    cat(error(" There are populations with one individual. Please remove populations with one individual or merged them with other populations fot his function to work"))
  }
  h_mean <- 1/mean(1/n.ind)
  n.pop <- length(pop.names)
  
  Ho <- lapply(pop.vec,function(y){
    colMeans(sgl_mat[[y]] == 1, na.rm = T)
  })
  Ho <- Reduce(cbind,Ho)
  colnames(Ho) <- pop.names
  mHo <- rowMeans(Ho,na.rm=TRUE) 
  
  Hs <- lapply(pop.vec,function(y){
    q <- colMeans(sgl_mat[[y]], na.rm = T)/2
    He <- 2 * (1-q) * q 
    n <- nrow(sgl_mat[[y]]) - colSums(apply(sgl_mat[[y]],2,is.na))
    n_factor <- (2 * n) / ((2 * n ) - 1 ) 
    return(n_factor * He)
  })
  Hs <- Reduce(cbind,Hs)
  colnames(Hs) <- pop.names
  mHs <- rowMeans(Hs,na.rm=TRUE)

  q_tmp <- colMeans(as.matrix(x), na.rm = T)/2
  Ht <- 2 * (1-q_tmp) * q_tmp
  
  Ht <- Ht + (mHs/(h_mean*n.pop)) - (mHo/(2*h_mean*n.pop))
  
  Dst <- Ht - mHs
  
  Fst <- Dst / Ht
  
  Fis <- 1 - (Ho/Hs)

  mFis <- 1 - (mHo/mHs)
    
  res <- cbind(mHo,mHs,Ht,Dst,Fst,mFis)
  
  colnames(res) <- c("Ho","Hs","Ht","Dst","Fst","Fis")
  
  overall <- colMeans(res, na.rm=TRUE)
  overall["Fst"] <- overall["Dst"] / overall["Ht"]
  overall["Fis"] <- 1 - (overall["Ho"] / overall["Hs"])
  
  all.res <- list("Ho" = as.data.frame(round(Ho, 4)), 
                  "Hs" = as.data.frame(round(Hs, 4)), 
                  "Fis" = as.data.frame(round(Fis, 4)), 
                  perloc = as.data.frame(round(res, 4)), 
                  overall = round(overall, 4))
  
  return(all.res)
  
}
