############################################ 
### pairwise LD function across all loci
############################################
#' Calculates pairwise population based Linkage Disequilibirum across all loci using the specifyied number of cores
#' @description this function is implemented in a parallel fashion to speed up the process. There is also the ability to restart the function if crashed by specifying the chunkfile names or restarting the function exactly in the same way as in the first run. This is implemented as sometimes due to connectivity loss between cores the function my crash half way. 
#' 
#' @param gi a genind object created via \code{\link{gl2gi}}
#' @param name character string for rdata file. If not given genind object name is used
#' @param save switch if results are saved in a file
#' @param nchunks how many subchunks will be used (the less the faster, but if the routine crashes more bits are lost
#' @param ncores how many cores should be used
#' @param chunkname the name of the chunks for saving, default is NULL 
#' @return returns calculation of pairwise LD across all loci between subpopulation. This functions uses if specified many cores on your computer to speed up. And if save is used can restart (if save=TRUE is used) with the same command starting where it crashed.
#' @export
#' @import data.table 
#' @import parallel 
#' @import foreach
#' @import data.table
#' @importFrom doParallel registerDoParallel
#' @author Bernd Gruber (glbugs@@aerg.canberra.edu.au)


gi.report.ld <- function(gi, name=NULL, save=TRUE,  nchunks=2, ncores=1, chunkname=NULL)
{
  #library(doParallel)
  #library(adegenet)
  #library(data.table)
  cat(paste("Start to calculate LD for all pairs of loci...\n"))
  cat(paste("Using", ncores,"cores in", nchunks," chunks.\n"))
  cat("Depending on the number of loci this may take a while...\n")
  cat("nchunks specifies the number of steps in the progress bar and the number of intermediate saves, but slows the computation a bit. nchunks = 1 is fastest.\n")
  cat(paste("Seperate all",length(indNames(gi)),"loci...\n"))
  flush.console()
  
  #convert into list of 
  slg <- seploc(gi)
    for (i in 1:length(slg)) slg[[i]] <- slg[[i]]@tab
  
  
  cat(paste("Generate all possible pairs:", length(slg)*(length(slg)-1)/2,"...\n"))
  flush.console()
  allp <- combn(length(slg),2)
  resnames <- c("loc1" ,"loc2","D", "Dprime", "r", "R2", "n", "X2", "p")
  lddone <- NULL
  chunknr <- 0  #to make sure old chunks are not overridden
  if (!is.null(chunkname)) 
  {
    cat("You specified results from a previous runs ...\n")
    cat(paste("Loooking for LD_chunks_", chunkname, "files.\n"))
    chunkfiles <- list.files(pattern=paste0("LD_chunks_",chunkname))
    if (length(chunkfiles>0))
      {
      cat(paste("Found", length(chunkfiles), "file(s).\n"))
      for (i in 1:length(chunkfiles))
        {
        load(paste0("LD_chunks_", chunkname,"_",i,".rdata"))
        lddone[[i]] <- ldc
        }
      chunknr<-length(chunkfiles)
      lddone <- rbindlist(lddone)
      
      cat(paste("Found", nrow(lddone),"pairs...\n"))
      setnames(lddone,resnames)
      done <- nrow(lddone)
      if (done==ncol(allp)) {cat("Already everyting is calculated. If you want to recalculate please delete al LD_chunk files or specify a different chunkname.\n Aborting function...\n");return(lddone)}
      allp <- allp[,-c(1:done)]  
      cat(paste("...only", ncol(allp), "pairs left to be done...\n"))
      } else cat(paste("No chunkfiles with LD_chunks_",chunkname,"_x.rdata found. \nTherefore I restart to calculate all pairs.\n"))
  }
  cat(paste("Calculate LD for all pairs...\n"))
  flush.console()
  n<- ncol(allp)
  ptm <- proc.time()[3]
  #ld functions here
  
  LD.fast <- function(g1,g2,...)
  {
    prop.A <- colMeans(g1, na.rm=T)/2
    if(length(prop.A)!=2) return(list(NA,NA,NA,NA,NA,NA,NA,NA))
    names(prop.A) <- c("A","B")
    prop.B <- colMeans(g2, na.rm=T)/2
    if(length(prop.B)!=2) return(list(NA,NA,NA,NA,NA,NA,NA,NA))
    names(prop.B) <- c("A","B")
    
    major.A <- names(prop.A)[which.max(prop.A)]
    major.B <- names(prop.B)[which.max(prop.B)]
    pA <- max(prop.A, na.rm=TRUE)
    pB <- max(prop.B, na.rm=TRUE)
    pa <- 1-pA
    pb <- 1-pB
    
    Dmin <- max(-pA*pB, -pa*pb)
    pmin <- pA*pB + Dmin;
    
    Dmax <- min(pA*pb, pB*pa);
    pmax <- pA*pB + Dmax;
    
    mja <- which.max(colSums(g1, na.rm=T))
    mjb <- which.max(colSums(g2, na.rm=T))
    counts <- table(g1[,mja], g2[,mjb], useNA="no")
    
    #counts <- counts[c("1","2"),c("1","2")]
    
    #     counts <- table(
    #       allele.count(g1, major.A),
    #       allele.count(g2, major.B) )
    #rownames(counts)<-1:2
    #colnames(counts)<-1:2
    
    n3x3 <- matrix(0, nrow=3, ncol=3)
    colnames(n3x3) <- rownames(n3x3) <- 0:2
    
    # ensure the matrix is 3x3, with highest frequency values in upper left
    for(i in rownames(counts))
      for(j in colnames(counts))
        n3x3[3-as.numeric(i),3-as.numeric(j)] <- counts[i,j]
    
    
    loglik <- function(pAB,...)
    {
      (2*n3x3[1,1]+n3x3[1,2]+n3x3[2,1])*log(pAB) +
        (2*n3x3[1,3]+n3x3[1,2]+n3x3[2,3])*log(pA-pAB) +
        (2*n3x3[3,1]+n3x3[2,1]+n3x3[3,2])*log(pB-pAB) +
        (2*n3x3[3,3]+n3x3[3,2]+n3x3[2,3])*log(1-pA-pB+pAB) + 
        n3x3[2,2]*log(pAB*(1-pA-pB+pAB) + (pA-pAB)*(pB-pAB))
    }
    
    # SAS code uses:
    #
    #s <- seq(pmin+0.0001,pmax-0.0001,by=0.0001)
    #lldmx <- loglik(s)
    #maxi <- which.max(lldmx)
    #pAB <- s[maxi]
    
    # but this should be faster:
    solution <- optimize(
      loglik,
      lower=pmin+.Machine$double.eps,
      upper=pmax-.Machine$double.eps,
      maximum=TRUE
    )
    pAB <- solution$maximum
    
    estD <- pAB - pA*pB
    if (estD>0)  estDp <- estD / Dmax else    estDp <- estD / Dmin
    
    n <-  sum(n3x3)
    
    corr <- estD / sqrt( pA * pB * pa * pb )
    
    dchi <- (2*n*estD^2)/(pA * pa * pB* pb)
    dpval <- 1 - pchisq(dchi,1)
    
    retval <- list(
      call=match.call(),
      "D"=estD,
      "D'"=estDp,
      "r" = corr,
      "R^2" = corr^2,
      "n"=n,
      "X^2"=dchi,
      "P-value"=dpval
    )
    class(retval) <- "LD"
    retval
  }
  
  #split into nchunks steps for the progress bar...
  runs <- 1:n
  if (nchunks>length(runs)) nchunks <- length(runs)
  chunks <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
  if (nchunks<2 ) splitruns <- list( runs) else splitruns <- chunks(runs, nchunks)
  ldchunks <- list()
  cl<-makeCluster(ncores) #adjust the number of cores of your computer!!!!
  registerDoParallel(cl)
  pb <- txtProgressBar(min=0, max=nchunks, style=3, initial=NA)
  for (i in 1:nchunks)
  {
    iter <- splitruns[[i]]
    if (ncores <2) it2 <- list(iter) else it2 <- chunks(iter, ncores)
    ip <- NA
    ll <- foreach (ip=1:length(it2), .combine=rbind) %dopar% {
        res <- matrix(NA, ncol=9, nrow=length(it2[[ip]]))
        for (ii in 1:length(it2[[ip]]))
          {
          l1 <- allp[1,it2[[ip]][ii]]
          l2 <- allp[2,it2[[ip]][ii]]
          s1 <- slg[[l1]]
          s2 <- slg[[l2]]
          if (  (length(colSums(s1, na.rm=T)) +  length(colSums(s2, na.rm=T)) ) ==4)
          {
          r <- LD.fast(s1,s2)
          res[ii,] <-c(l1,l2, do.call(cbind, r[2:8]))
          } else res[ii,] <- c(l1,l2, rep(NA,7))
          }
        res
        }
  ldchunks[[i]] <-as.data.frame(ll)
  setTxtProgressBar(pb, i)
  ldc <- ldchunks[[i]]
  save(ldc, file=paste0("LD_chunks_",chunkname,"_",i+chunknr,".rdata"))
  }
stopCluster(cl)
LDres2 <- rbindlist(ldchunks)
#LDres2 <- t(do.call(cbind, ldchunks))
setnames(LDres2,resnames)
if (!is.null(lddone)) LDres2 <- rbindlist(list(lddone, LDres2))
#colnames(LDres2)<- resnames

cat(paste("\n# Simulations:", n,". Took", round(proc.time()[3]-ptm),"seconds.\n"))
if (save) 
{
  if (!is.null(name)) 
  { 
    nobj <- name
    filename <- paste0("LD_",name,".rdata")
  } else
  {
    fx <- deparse(substitute(gi))
    #cat(paste('name:',fx,"\n"))
    fx <- gsub("\\[|\\]|[.]|[:]|[,]|[&]|[%]","_",fx) 
    fx <- gsub("__","_", fx)
    fx <- gsub(" ","", fx)
    fx <- gsub("_$","",fx)
    if (length(fx)>0) 
    {
      nobj <- paste0("LD_",fx)
      filename=paste0(fx,".rdata") 
    } else
    {
      filename="LDallp.rdata"
      nobj <- "LDallp"
    }
  }
  cat(paste0("\n Results are saved as object ", nobj," under ", filename,".\n"))
  (cat(paste("Once you have checked you can delete your LD_chunks_",chunkname,"files.\n")))
  assign(nobj, LDres2)
  save(list=nobj, file=filename)
}
LDres2
}

