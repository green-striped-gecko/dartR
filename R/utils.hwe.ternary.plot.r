###########################################################
######################### GenerateSamples #######################
###########################################################
`GenerateSamples` <- function(n=5) {
  # generates all possible samples of size n.
  Res <- NULL
  for (i in 0:n) {
    AA <- i
    for (j in 0:(n-i)) {
      AB <- j
      BB <- (n-(AA+AB))
      sam <- c(AA,AB,BB)
      Res <- rbind(Res,sam)
    }
  }
  rownames(Res) <- 1:nrow(Res)
  colnames(Res) <- c("AA","AB","BB")
  return(Res)
}

###########################################################
######################### CritSam #######################
###########################################################
`CritSam` <- function(n=5,Dpos=TRUE,alphalimit=0.05,pvaluetype="dost") {
  X <- GenerateSamples(n)
  ncomp <- nrow(X)
  Res <- NULL
  Ds <- NULL
  pval <- NULL
  fA <- NULL
  for (i in 1:nrow(X)) {
    fA <- c(fA,(2*X[i,1]+X[i,2])/(2*n))
    Ds <- c(Ds,HardyWeinberg::HWChisq(X[i,],verbose =FALSE)$D)
    pval <- c(pval,HardyWeinberg::HWExact(X[i,],alternative="two.sided",pvaluetype=pvaluetype,verbose =FALSE)$pval)
  }
  
  Y <- data.frame(X[,1],X[,2],X[,3],fA,Ds,pval)
  colnames(Y) <- c("AA","AB","BB","fA","Ds","pval")
  if(Dpos) Y <- Y[Y$Ds>0,] else Y <- Y[Y$Ds<0,]
  Y <- Y[Y$pval<alphalimit,]
  fre <- unique(fA)
  for (i in 1:length(fre)) {
    Ys <- Y[Y$fA == fre[i],]
    if (nrow(Ys) > 0) {
      indi <- which.max(Ys$pval)
      Ys <- Ys[indi,]
      Res <- rbind(Res,c(Ys$AA,Ys$AB,Ys$BB))
    }
  }
  Xn <- Res/n
  return(list(Xn=Xn,Ds=Ds,fA=fA))
}
###########################################################
######################### alowcurve #######################
###########################################################
`alowcurve` <- function(p,q,chiquant,n,cc=0.5) {
  y <- (p*q-cc*(1-p*q)/n)^2 - (cc^2*p*q*(p*q-0.5)/n^2 + p^2*q^2*chiquant/n)
  return(y)
}


###########################################################
######################### HWChisqccCurve ##################
###########################################################

`HWChisqccCurve` <- function(p,q,chiquant,n,cc=0.5,curvetype="DposUp") {
  switch(curvetype,
         DposUp  = 2*p*q+2*cc*(1-p*q)/n+2*sqrt(cc^2*p*q*(p*q-0.5)/n^2 + p^2*q^2*chiquant/n),
         DposLow = 2*p*q+2*cc*(1-p*q)/n-2*sqrt(cc^2*p*q*(p*q-0.5)/n^2 + p^2*q^2*chiquant/n),
         DnegUp  = 2*p*q-2*cc*(1-p*q)/n+2*sqrt(cc^2*p*q*(p*q-0.5)/n^2 + p^2*q^2*chiquant/n),
         DnegLow = 2*p*q-2*cc*(1-p*q)/n-2*sqrt(cc^2*p*q*(p*q-0.5)/n^2 + p^2*q^2*chiquant/n))
}

###########################################################
######################### CheckRoots ######################
###########################################################

`CheckRoots` <- function(roots,limits=c(0,1),mini=TRUE,verbose=FALSE) {
  if(verbose) {  
    cat("checking roots\n")
    print(roots)
  }
  roots <- roots[round(roots$rAB,digits=5)==round(roots$check,digits=5),]
  roots <- roots[roots$"Im<1e-10"==TRUE,]
  roots <- roots[(roots$Re < limits[2] & roots$Re > limits[1]),]
  roots <- roots$Re
  if(length(roots)>1)
    if(verbose) cat("multiple roots\n")
  if(length(roots)!= 0) {
    if(mini)
      y <- min(roots)
    else y <- max(roots)
  }
  else y<-NULL
  if(verbose) print(y)
  return(y)
}

###########################################################
######################### DnegLow #########################
###########################################################

`DnegLow` <- function(r,k,n,cc,chiquant,verbose=FALSE,cex=1,curcol="black",curtyp="solid") {
  # draw the lower cl for D<0, chisquare with cc
  ll <- r/(1+r)
  ul <- 1/(1+r)
  
  kla <- k^2
  klb <- -2*k - 1.5*k^2
  klc <- 1 + 4*k + 1.5*k^2 - chiquant/n
  kld <- -2 - 4*k + 2*chiquant/n
  kle <- 1 + 2*k - chiquant/n
  
  lcoef <- c(kla,klb,klc,kld,kle)
  
  lout <- polyroot(lcoef)
  Imlout <- round(abs(Im(lout)),digits=6)
  Relout <- round(Re(lout),digits=6)
  nlroot <- length(lout)
  
  lroots <- data.frame(1:nlroot,lout,Imlout,Relout,Imlout<1e-10,
                       round(HWChisqccCurve(1-Relout,Relout,chiquant,n,curvetype="DnegLow"),digits=6),
                       round(alowcurve(1-Relout,Relout,chiquant,n),digits=6))
  
  colnames(lroots) <- c("Root","Value","Im","Re","Im<1e-10","rAB","check")
  
  if(verbose) cat("Roots lower curve\n")
  
  ql <- CheckRoots(lroots,verbose=verbose)
  
  DnegLL <- ql
  
  if (!is.null(ql)) {
    q <- seq(ql,1-ql,by=0.005)
    p <- 1-q
    pt <-  2*(p-0.5)/sqrt(3)
    llc <- 2*p*q-2*cc*(1-p*q)/n-2*sqrt(cc^2*p*q*(p*q-0.5)/n^2 + p^2*q^2*chiquant/n)
    points(pt,llc,type="l",col=curcol,cex=cex,lty=curtyp)
  }
  return(DnegLL=DnegLL)
}


###########################################################
######################### DnegUp ##########################
###########################################################

`DnegUp` <- function(r,k,n,cc,chiquant,verbose=FALSE,cex=1,curcol="black",curtyp="dotted")
{
  # draw the upper cl for D<0, chisquare with cc
  
  kra <- k^2
  krb <- -1.5*k^2
  krc <- 1.5*k^2 + 2*k - chiquant/n
  krd <- 2*chiquant/n - 2*k
  kre <- 1 + 2*k - chiquant/n
  
  coefr <- c(kra,krb,krc,krd,kre)
  
  out <- polyroot(coefr)
  Imout <- round(abs(Im(out)),digits=6)
  Reout <- round(Re(out),digits=8)
  nroot <- length(out)
  
  lmidroots <- data.frame(1:nroot,out,Imout,Reout,Imout<1e-10,HWChisqccCurve(1-Reout,Reout,chiquant,n,curvetype="DnegUp"),round(2*Reout,digits=8))
  
  colnames(lmidroots) <- c("Root","Value","Im","Re","Im<1e-10","rAB","check")
  
  if(verbose) cat("D<0; Upper Curve\n")
  
  # paint the upper part of the curve
  ql <- CheckRoots(lmidroots,verbose=verbose,mini=FALSE)
  
  DnegUL <- ql
  
  if (!is.null(ql)) {
    q <- seq(ql,1-ql,by=0.005)
    p <- 1-q
    qt <-  2*(q-0.5)/sqrt(3)
    ul <- 2*p*q-2*cc*(1-p*q)/n+2*sqrt(cc^2*p*q*(p*q-0.5)/n^2 + p^2*q^2*chiquant/n)
    points(qt,ul,col=curcol,type="l",cex=cex,lty=curtyp)
  }
  
  # paint a second lower part of the curve.
  
  d1 <- (-3*cc^2 - 4*n*cc)/(2*n^2)
  d2 <- (2*cc-chiquant)/n + 1
  
  kla <- (cc/n)^2
  klb <- d1
  klc <- d2-d1
  kld <- -2*d2
  kle <- d2
  
  lcoef <- c(kla,klb,klc,kld,kle)
  
  lout <- polyroot(lcoef)
  Imlout <- round(abs(Im(lout)),digits=6)
  Relout <- Re(lout)
  nlroot <- length(lout)
  
  lroots <- data.frame(1:nlroot,lout,Imlout,Relout,Imlout<1e-10,round(HWChisqccCurve(1-Relout,Relout,chiquant,n,curvetype="DnegUp"),digits=10),0)
  
  colnames(lroots) <- c("Root","Value","Im","Re","Im<1e-10","rAB","check")
  
  
  if(verbose) {
    cat("Zeros Upper curve::\n")
    print(lroots)
  }
  
  ql <- CheckRoots(lroots,verbose=verbose)
  
  qlu <- CheckRoots(lmidroots,verbose=verbose,mini=TRUE)
  
  ql <- min(ql,qlu)
  qlu <- max(ql,qlu)
  if(verbose) print(ql)
  
  DnegUL <- ql
  
  if (!is.null(ql)) {
    q <- seq(ql,qlu,by=0.005)
    p <- 1-q
    qt <-  2*(q-0.5)/sqrt(3)
    ul <- 2*p*q-2*cc*(1-p*q)/n+2*sqrt(cc^2*p*q*(p*q-0.5)/n^2 + p^2*q^2*chiquant/n)
    points(qt,ul,col=curcol,type="l",cex=cex,lty=curtyp)
    
    q <- seq(1-qlu,1-ql,by=0.005)
    p <- 1-q
    qt <-  2*(q-0.5)/sqrt(3)
    ul <- 2*p*q-2*cc*(1-p*q)/n+2*sqrt(cc^2*p*q*(p*q-0.5)/n^2 + p^2*q^2*chiquant/n)
    points(qt,ul,col=curcol,type="l",cex=cex,lty=curtyp)
    
  }
  return(DnegUL=DnegUL)     
}

###########################################################
######################### DposLow #########################
###########################################################

`DposLow` <-
  function(r,k,n,cc,chiquant,verbose=FALSE,cex=1,curcol="black",curtyp="dotted") {
    # draw the lower cl for D>0, chisquare with cc
    ka <- k^2
    kb <- -1.5*k^2
    kc <- 1.5*k^2 - 2*k - chiquant/n
    kd <- 2*chiquant/n + 2*k
    ke <- 1 - 2*k - chiquant/n
    
    coef <- c(ka,kb,kc,kd,ke)
    
    out <- polyroot(coef)
    Imout <- round(abs(Im(out)),digits=6)
    Reout <- round(Re(out),digits=8)
    nroot <- length(out)
    
    umidroots <- data.frame(1:nroot,out,Imout,Reout,Imout<1e-10,HWChisqccCurve(1-Reout,Reout,chiquant,n,curvetype="DposLow"),round(2*Reout,digits=8))
    
    
    
    colnames(umidroots) <- c("Root","Value","Im","Re","Im<1e-10","rAB","check")
    
    if(verbose) cat("upper (D>0) mid curve\n")
    
    ql <- CheckRoots(umidroots,verbose=verbose)
    
    DposLL <- ql
    
    if (!is.null(ql)) {
      q <- seq(ql,1-ql,by=0.005)
      p <- 1-q
      qt <-  2*(q-0.5)/sqrt(3)
      ul <- 2*p*q+2*cc*(1-p*q)/n-2*sqrt(cc^2*p*q*(p*q-0.5)/n^2 + p^2*q^2*chiquant/n)
      points(qt,ul,col=curcol,type="l",cex=cex,lty=curtyp)
    }
    return(DposLL=DposLL)
  }

###########################################################
######################### DposUp ##########################
###########################################################

`DposUp` <-
  function(r,k,n,cc,chiquant,verbose=FALSE,cex=1,curcol="black",curtyp="solid") {
    # draw the upper cl for D>0, chisquare with cc
    ka <- -k^2
    kb <- 1.5*k^2
    kc <- chiquant/n - 1.5*k^2 + 2*k
    kd <- - 2*chiquant/n - 2*k
    ke <- chiquant/n - 1 + 2*k
    
    coef <- c(ka,kb,kc,kd,ke)
    
    out <- polyroot(coef)
    Imout <- round(abs(Im(out)),digits=6)
    Reout <- round(Re(out),digits=8)
    nroot <- length(out)
    
    uroots <- data.frame(1:nroot,out,Imout,Reout,Imout<1e-10,HWChisqccCurve(1-Reout,Reout,chiquant,n,
                                                                            curvetype="DposUp"),round(2*Reout,digits=8))
    
    colnames(uroots) <- c("Root","Value","Im","Re","Im<1e-10","rAB","check")
    
    if(verbose) cat("Roots upper curve\n")
    
    ql <- CheckRoots(uroots,verbose=verbose)
    
    DposUL <- ql
    
    if (!is.null(ql)) {
      q <- seq(ql,1-ql,by=0.005)
      p <- 1-q
      qt <-  2*(q-0.5)/sqrt(3)
      ul <- 2*p*q+2*cc*(1-p*q)/n+2*sqrt(cc^2*p*q*(p*q-0.5)/n^2 + p^2*q^2*chiquant/n)
      points(qt,ul,type="l",cex=cex,col=curcol,lty=curtyp)
    }
    return(DposUL=DposUL)
  }

###########################################################
######################### HWTernaryPlot ####################
###########################################################
`HWTernaryPlot_correction` <- function(X, n=NA, addmarkers=TRUE, newframe=TRUE, hwcurve=TRUE, vbounds=TRUE, mafbounds=FALSE, mafvalue=0.05, axis=0, region=1, vertexlab=colnames(X), alpha = 0.05, vertex.cex = 1, pch = 19, cc = 0.5, markercol = "black", markerbgcol= "black", cex=0.75, axislab ="", verbose=FALSE, markerlab=NULL, markerpos=NULL, mcex=1, connect = FALSE, curvecols=rep("black",5) , signifcolour=TRUE, curtyp = "solid", ssf = "max", pvaluetype = "selome", multi_comp=multi_comp, multi_comp_method=multi_comp_method, ...)
{
  # plot a ternary diagram that represents all rows of X as points. Acceptance regions for various tests for HWE
  #    can be added.
  if(is.vector(X)) {
    if(length(X)!=3) {
      stop("X must have three elements")
    }
    else {
      X <- matrix(X,ncol=3,dimnames=list(c("1"),names(X)))
    }
  }
  X <- as.matrix(X)
  nr <- nrow(X)
  nc <- ncol(X)
  if(any(X<0)) stop("X must be non-negative")
  if(nc != 3) stop("X must have three columns")
  if(is.na(n)) {
    if((sum(apply(X,1,sum))==nr))  # data are compositions
      stop("argument n (the sample size) should be supplied")
    else { # raw counts
      ssf <- match.fun(ssf)
      vsums <- as.matrix(apply(X,1,sum),ncol=1)
      n <- apply(vsums,2,ssf)
      Xr <- X
      if (nrow(X) == 1) 
        Xcom <- X/sum(X)
      else {
        Xcom <- HardyWeinberg::HWClo(X)
      }
    }
  }
  else {
    if((sum(apply(X,1,sum))==nr)) {  # data are compositions
      Xr <- round(n*X)
      Xcom <- X
    }
    else { # raw counts
      Xr <- X
      if (nrow(X) == 1) 
        Xcom <- X/sum(X)
      else {
        Xcom <- HardyWeinberg::HWClo(X)
      }
    }
  }
  chiquant <- qchisq(1-alpha,1)
  r <- sqrt(chiquant/n)
  k <- cc/n
  M <- matrix(c(-1/sqrt(3),0,0,1,1/sqrt(3),0),ncol=2,byrow=T)
  nsignif <- NA
  
  markerq <- (Xcom[,2]+2*Xcom[,3])/2
  
  if(newframe) {
    opar <- par(pty="m",xpd=TRUE)
    on.exit(par(opar))
    
    plot(M[,1],M[,2], type="n", axes=FALSE, xlab="", ylab="", pch=19, asp=1, cex.main=2, ... )
    polygon(M)
    
    eps <- 0.04 * vertex.cex
    Mlab <- M + matrix(c(-eps,0,0,eps,eps,0),ncol=2,byrow=T)
    text(Mlab[,1],Mlab[,2], vertexlab, cex=vertex.cex)
    
    text(0,-0.1,axislab,cex=vertex.cex)
  }
  
  if (axis==1) {
    AXA <- rbind(c(0,0.5,0.5),c(1,0,0))
    AXA <- AXA%*%M
    lines(AXA[,1],AXA[,2],...)      
  }
  
  if (axis==2) {
    AXAB <- rbind(c(0.5,0,0.5),c(0,1,0))
    AXAB <- AXAB%*%M
    lines(AXAB[,1],AXAB[,2],...)
  }
  
  if (axis==3) {
    AXB <- rbind(c(0.5,0.5,0),c(0,0,1))
    AXB <- AXB%*%M
    lines(AXB[,1],AXB[,2],...)
  }
  
  
  if(hwcurve) {
    p <- seq(0,1,by=0.005)
    HW <- cbind(p^2,2*p*(1-p),(1-p)^2)
    HWc <- HW%*%M
    points(HWc[,1],HWc[,2],type="l",col=curvecols[1])
  }
  
  minp  <- sqrt(5/n)
  minpt <- 2*(minp-0.5)/sqrt(3)
  maxp  <- 1-sqrt(5/n)
  maxpt <- 2*(maxp-0.5)/sqrt(3)
  
  inrange <- sum((markerq >= minp & markerq <= maxp))
  percinrange <- round(100*inrange/nr,digits=2)
  
  ind1 <- markerq==1
  ind0 <- markerq==0
  nfixed <- sum(ind1) + sum(ind0)
  
  D <- 0.5*(Xcom[,2] - 2*(1-markerq)*markerq)
  Dpos <- sum(Xcom[,2] > 2*(1-markerq)*markerq)
  Dneg <- sum(Xcom[,2] < 2*(1-markerq)*markerq)
  Dzer <- sum(Xcom[,2] == 2*(1-markerq)*markerq)
  Dtot <- Dpos+Dneg
  
  
  #       cat("D>0:",Dpos,"D<0:",Dneg,"D=0:",Dzer,"nfix:",nfixed,"Tot:",Dtot,"mD:",mean(D),"medD:",median(D),"\n")
  #       cat("D>0:",round(100*Dpos/Dtot,digits=5),"D<0:",round(100*Dneg/Dtot,digits=5),"D=0:",
  #           round(100*Dzer/Dtot,digits=5),"\n")
  
  if(vbounds) {
    if (n >= 20) {
      lines(c(minpt,minpt),c(0,2*minp),lty="dashed")
      lines(c(maxpt,maxpt),c(0,2-2*maxp),lty="dashed")
    }
  }
  
  if(mafbounds) {
    minaf <- mafvalue
    minaft <- 2*(minaf-0.5)/sqrt(3)
    maxaf <- 0.95
    maxaft <- 2*(maxaf-0.5)/sqrt(3)
    lines(c(minaft,minaft),c(0,2*minaf),lty="dashed")
    lines(c(maxaft,maxaft),c(0,2-2*maxaf),lty="dashed")
  }    
  

  if(region==2) { # all curves for hw with cc
    
    # lower curve for D<0
    
    DnegLL <- DnegLow(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[4],curtyp=curtyp)
    
    # upper curve for D<0
    
    DnegUL <- DnegUp(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[4],curtyp=curtyp)
    
    # lower curve for D>0
    
    DposLL <- DposLow(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[3],curtyp=curtyp)
    
    # upper curve for D>0
    
    DposUL <- DposUp(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[3],curtyp=curtyp)
    
    if(verbose) {
      cat("D<0 LL",round(DnegLL,digits=2),"\n")
      cat("D<0 UL",round(DnegUL,digits=2),"\n")
      cat("D>0 LL",round(DposLL,digits=2),"\n")
      cat("D>0 UL",round(DposUL,digits=2),"\n")
    }
  }
  
  if(region==7) { # For Haldane's Exact test
    
    Crit <- CritSam(n,alphalimit=alpha,pvaluetype=pvaluetype)$Xn
    Critcar <- Crit%*%M        # cartesian coordinates
    points(Critcar[,1],Critcar[,2],pch=pch,col=curvecols[5],cex=cex,type="l",lty=curtyp)
    
    Crit <- CritSam(n,Dpos=FALSE,alphalimit=alpha,pvaluetype=pvaluetype)$Xn
    Critcar <- Crit%*%M        # cartesian coordinates
    points(Critcar[,1],Critcar[,2],pch=pch,col=curvecols[5],cex=cex,type="l",lty=curtyp)
    
  }
  
  if(addmarkers) {
    Xc <- Xcom%*%M # cartesian coordinates
    
    if (signifcolour==TRUE) { 
      
       if (region == 2) {
        pvals <- numeric(nr)
        
        pvals <- HardyWeinberg::HWChisqStats(Xr)
        
        if(multi_comp == TRUE){
          pvals <- stats::p.adjust(pvals, method = multi_comp_method)
        }
        
        markerbgcol <- rep("green",nr)             
        markerbgcol[pvals<alpha] <- "red"
        markercol <- rep("green",nr)             
        markercol[pvals<alpha] <- "red"
        nsignif <- sum(pvals<alpha)
       }
      
      if (region == 7) {
        pvals <- numeric(nr)
        
        pvals <- HardyWeinberg::HWExactStats(Xr)
        
        if(multi_comp == TRUE){
          pvals <- stats::p.adjust(pvals, method = multi_comp_method)
        }
        
        markerbgcol <- rep("green",nr)             
        markerbgcol[pvals<alpha] <- "red"
        markercol <- rep("green",nr)             
        markercol[pvals<alpha] <- "red"
        nsignif <- sum(pvals<alpha)
      }
    }
    
    if (connect) {
      points(Xc[,1],Xc[,2],pch=pch,col=curvecols[1],cex=cex,type="l",lty=curtyp)
      points(Xc[,1],Xc[,2],pch=pch,bg=markerbgcol,col=markercol,cex=cex)
    }
    else
      points(Xc[,1],Xc[,2],pch=pch,bg=markerbgcol,col=markercol,cex=cex)
    text(Xc[,1],Xc[,2],markerlab,cex=mcex,pos=markerpos)
  }
  results <- list(minp=minp,maxp=maxp,inrange=inrange,percinrange=percinrange,nsignif=nsignif)
}