gl.report.pa2 <- function(gl1, gl2=NULL, verbose=NULL){
	
	# TRAP COMMAND, SET VERSION
	
	funname <- match.call()[[1]]
	build <- "Jacob"
	
	# SET VERBOSITY
	
	if (is.null(verbose)){ 
		if(!is.null(gl1@other$verbose)){ 
			verbose <- gl1@other$verbose
		} else { 
			verbose <- 2
		}
	} 
	
	if (verbose < 0 | verbose > 5){
		cat(paste("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n"))
		verbose <- 2
	}
	
	# FLAG SCRIPT START
	
	if (verbose >= 1){
		if(verbose==5){
			cat("Starting",funname,"[ Build =",build,"]\n")
		} else {
			cat("Starting",funname,"\n")
		}
	}
	
	# STANDARD ERROR CHECKING
	
	if(class(gl1)!="genlight") {
		stop("Fatal Error: genlight object required!\n")
	}
	if (all(gl1@ploidy == 1)){
		stop("Cannot calculate minor allele frequences for Tag presence/absence data. Please provide a SNP dataset.\n")
	} else if (all(gl1@ploidy == 2)){
		if(verbose>=2){cat("  Processing a SNP dataset\n")}
	} else {
		stop("Fatal Error: Ploidy must be universally 1 (Tag P/A data) or 2 (SNP data)!\n")
	}
	
	if(!is.null(gl2)){
		if(class(gl2)!="genlight") {
			stop("Fatal Error: genlight object required for gl2!\n")
		}
		if (all(gl2@ploidy == 1)){
			stop("Fatal Error: Private alleles can only be calculated for SNP data. Please provide a SNP dataset for gl2\n")
		} else if (all(gl2@ploidy == 2)){
			if (verbose >= 2){cat("  Processing a SNP dataset",gl2,"\n")}
		} else {
			stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
		}
	}
	
	# FUNCTION SPECIFIC ERROR CHECKING
	
	if (!is.null(gl2)) pops <- list(pop1=gl1, pop2=gl2) else 
	{
		if (length(unique(pop(gl1)))>1) pops <- seppop(gl1) else stop("Only one population provided. Check the @pop slot in your genlight object.\n ")
	}
	
	# DO THE JOB
	
	pc <- t(combn(length(pops),2))
	palist <- list()
	pall <- data.frame(p1=pc[,1], p2=pc[,2], pop1=names(pops)[pc[,1]], pop2=names(pops)[pc[,2]], N1=NA, N2=NA,fixed=NA, priv1=NA, priv2=NA, totalpriv=NA, mdf=NA)
	
	for (i in 1:nrow(pc))
	{
		
		i1 =pall[i,1]
		i2 =pall[i,2]
		
		palist[[i]] <- list(fixed=NA,priv1=NA, priv2=NA)
		names(palist)[i]<- paste(names(pops)[i1], names(pops)[i2], sep="-"
		)		
		p1 <- as.matrix(pops[[i1]])
		p2 <- as.matrix(pops[[i2]])
		p1alf <- colMeans(p1, na.rm = T)/2
		p2alf <- colMeans(p2, na.rm = T)/2
		
		pall[i,5:6] <- c(nrow(p1), nrow(p2))
		pall[i,7] = sum(abs(p1alf-p2alf)==1, na.rm=T)
		dummy <- locNames(pops[[i1]])[abs(p1alf-p2alf)==1]
		palist[[i]]$fixed <- dummy[!is.na(dummy)]
		pall[i,8] =  sum(p2alf==0 & p1alf!=0, na.rm=T) + sum(p2alf==1 & p1alf!=1, na.rm = T) 
		dummy <-  locNames(pops[[i1]])[p2alf==0 & p1alf!=0 | p2alf==1 & p1alf!=1]
		palist[[i]]$priv1<- dummy[!is.na(dummy)]
	
		pall[i,9] =  sum(p1alf==0 & p2alf!=0, na.rm=T) + sum(p1alf==1 & p2alf!=1, na.rm = T)  
		dummy <- locNames(pops[[i1]])[p1alf==0 & p2alf!=0 | p1alf==1 & p2alf!=1]
		palist[[i]]$priv2<- dummy[!is.na(dummy)]
		pall[i,10] = pall[i,8]+pall[i,9]
		pall[i,11] = round(mean(abs(p1alf-p2alf), na.rm=T),3)
	}
	
	if(verbose >= 3){
		print(pall)
	}
	
	if(verbose >= 2){cat("  Table of private alleles and fixed differences returned\n")}
	
	# FLAG SCRIPT END
	
	if (verbose >= 1) {
		cat("Completed:",funname,"\n")
	}
	
	return(list(pa.table=pall, pa.list=palist))
	
}
