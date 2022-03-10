#' Calculates mean observed heterozygosity, mean expected heterozygosity and Fis
#'  per population and Fst across all populations 
#' @param x A genlight object containing the SNP genotypes [required].
#' @return A list with with the statistics for each population
#' @export
#' @author Luis Mijangos (bugs? Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' out <- utils.basic.stats(testset.gl)
#' @export

utils.basic.stats <- function(x) {
  
  pop.names <- popNames(x)
  sgl_mat <- lapply(seppop(x),as.matrix)
  pop.vec <- 1:nPop(x)
  
  Ho <- lapply(pop.vec,function(y){
    colMeans(sgl_mat[[y]] == 1, na.rm = T)
  })
  Ho <- Reduce(cbind,Ho)
  colnames(Ho) <- pop.names
  mHo <- rowMeans(Ho,na.rm=TRUE) 
  
  Hs <- lapply(pop.vec,function(y){
    q <- colMeans(sgl_mat[[y]], na.rm = T)/2
    return( 2 * (1-q) * q )
  })
  Hs <- Reduce(cbind,Hs)
  colnames(Hs) <- pop.names
  mHs <- rowMeans(Hs,na.rm=TRUE)

  q_tmp <- colMeans(as.matrix(x), na.rm = T)/2
  Ht <- 2 * (1-q_tmp) * q_tmp

  Fst <- (Ht-mHs)/Ht
  
  Fis <- 1 - (Ho/Hs)

  mFis <- 1 - (mHo/mHs)
    
  res <- cbind(mHo,mHs,Ht,Fst,mFis)
  
  colnames(res) <- c("Ho","Hs","Ht","Fst","Fis")
  
  overall <- colMeans(res, na.rm=TRUE)
  overall["Fst"] <- (overall["Ht"] - overall["Hs"]) / overall["Ht"]
  overall["Fis"] <- 1 - (overall["Ho"]/overall["Hs"])
  
  all.res <- list("Ho" = as.data.frame(round(Ho, 4)), 
                  "Hs" = as.data.frame(round(Hs, 4)), 
                  "Fis" = as.data.frame(round(Fis, 4)), 
                  perloc = as.data.frame(round(res, 4)), 
                  overall = round(overall, 4))
  
  return(all.res)
  
}
