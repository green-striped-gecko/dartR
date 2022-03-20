#' Calculates mean observed heterozygosity, mean expected heterozygosity and Fis
#'  per locus, per population and Fst across all populations 
#'  
#'  This is a re-implementation of \code{hierfstat::basics.stats} specifically 
#'  for genlight objects. Formula (and hence results) match exactly the original 
#' version of \code{hierfstat::basics.stats} (although Dstp and Fstp are not 
#' currently computed) but it is much faster.
#' 
#' @param x A genlight object containing the SNP genotypes [required].
#' @return A list with with the statistics for each population
#' @export
#' @author Luis Mijangos and Carlo Pacioni (bugs? Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' out <- utils.basic.stats(platypus.gl)
#' @export

utils.basic.stats <- function(x) {
  
  n.ind <- table(pop(x))
  
    if(any(n.ind <= 1)){
    cat(error(" There are populations with one individual. Please remove populations with one individual or merged them with other populations for his function to work\n"))
    }
  
  pop.names <- popNames(x)
  sgl_mat <- lapply(seppop(x),as.matrix)
  pop.vec <- 1:nPop(x)
  
  n.pop <- lapply(pop.vec,function(y){
    apply(sgl_mat[[y]], 2, function(x){
      all(is.na(x))
    })
  })    
  
  n.pop <- Reduce("+",n.pop)
  n.pop <- nPop(x) - n.pop
  
  np <- lapply(pop.vec,function(y){
    colSums(!is.na(sgl_mat[[y]])) 
  })
  
  np <- Reduce(cbind,np)
  
  mn <- apply(np,1,function(y){
    1/mean(1/y)
  })

  Ho <- lapply(pop.vec,function(y){
    colMeans(sgl_mat[[y]] == 1, na.rm = TRUE)
  })
  Ho <- Reduce(cbind,Ho)
  colnames(Ho) <- pop.names
  mHo <- rowMeans(Ho, na.rm=TRUE) 
  
  q <- lapply(pop.vec,function(y){
    colMeans(sgl_mat[[y]], na.rm = TRUE)/2
  })

  Hs <- lapply(pop.vec,function(y){
    n <- nrow(sgl_mat[[y]]) - colSums(apply(sgl_mat[[y]],2,is.na))
    Hs <- 2 * (1-q[[y]]) * q[[y]] - Ho[,y]/2/n
    #n_factor <- (2 * n) / ((2 * n ) - 1 ) 
    return(n/(n - 1) * Hs)
  })
  Hs <- Reduce(cbind,Hs)
  colnames(Hs) <- pop.names
  #mHs <- rowMeans(Hs,na.rm=TRUE)
  q_m <- Reduce(cbind, q)
  sd2 <- q_m^2 + (1 - q_m)^2
  msp2 <- rowMeans(sd2, na.rm = TRUE)
  mHs <- mn/(mn - 1) * (1 - msp2 - mHo/2/mn)
  
  q_mean <- rowMeans(Reduce(cbind,q),na.rm = T)
  Ht <- 2 * (1-q_mean) * q_mean
  Ht <- Ht + mHs/mn/n.pop - mHo/2/mn/n.pop
  
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
