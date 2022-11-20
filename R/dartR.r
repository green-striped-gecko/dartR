
#seppop needs to be imported to work for dartR
#also internal functions for "[" methods
seppop <- getFromNamespace("seppop", "adegenet")
.seppop_internal <- getFromNamespace(".seppop_internal", "adegenet")
.get_pop_inds <- getFromNamespace(".get_pop_inds","adegenet")

setClass("dartR", contains="genlight",)


### adding new slots...
#setClass("dartR",slots=list(what="integer"), contains="genlight")
###

########################
## show dartR ##
########################
setMethod ("show", "dartR", function(object){
  ## HEADER
  cat(" ********************\n")
  cat(" *** DARTR OBJECT ***\n")
  cat(" ********************")
  marker <- "mixed markers"
  if (all(!is.na(ploidy(object)))) {
  if (all(ploidy(object)==2)) marker <-"SNPs" 
  if (all(ploidy(object)==1)) marker <- "silicoDarts (P/A) "
  } 
  
  cat("\n\n **", format(nInd(object), big.mark=","), "genotypes, ",
      format(nLoc(object), big.mark=","), marker,", size:", format(object.size(object), units="auto"))
  
  temp <- sapply(object@gen, function(e) length(e@NA.posi))
  if(length(temp>1)){
    cat("\n\n    missing data: ", sum(temp), " (=", round((sum(temp)/(nInd(object)*nLoc(object))) *100,2)," %) scored as NA", sep="")
  }
  
  ## BASIC CONTENT
  cat("\n\n ** Genetic data")
  cat("\n   @gen: list of", length(object@gen), "SNPbin")

  if(!is.null(object@ploidy)){
    ploidytxt <- paste("(range: ", paste(range(object@ploidy), collapse="-"), ")", sep="")
    cat("\n   @ploidy: ploidy of each individual ", ploidytxt)
  }
  
  ## Additional data
  cat("\n\n ** Additional data")
  optional <- FALSE
  
  if(!is.null(object@ind.names)){
    optional <- TRUE
    cat("\n   @ind.names: ", length(object@ind.names), "individual labels")
  } else cat("\n   @ind.names: ", "no individual labels")
  
  if(!is.null(object@loc.names)){
    optional <- TRUE
    cat("\n   @loc.names: ", length(object@loc.names), "locus labels")
  } else cat("\n   @loc.names: ", "no locus labels")
  
  if(!is.null(object@loc.all)){
    optional <- TRUE
    cat("\n   @loc.all: ", length(object@loc.all), "allele labels")
  } else cat("\n   @loc.all: "," no allele labels") 
  
  if(!is.null(object@chromosome)){
    optional <- TRUE
    cat("\n   @chromosome: factor storing chromosomes of the", marker)
  } 
  
  if(!is.null(object@position)){
    optional <- TRUE
    cat("\n   @position: integer storing positions of the",marker,"[within 69 base sequence]")
  } 
  if(!is.null(object@pop)){
    optional <- TRUE
    poptxt <- paste("(group size range: ", paste(range(table(object@pop)), collapse="-"), ")", sep="")
    cat("\n   @pop:", paste("population of each individual", poptxt))
  } else cat("\n   @pop:", "no population lables for individuals")
  
  if (!is.null(object@strata)){
    optional <- TRUE
    cat("\n   @strata: ")
    levs <- names(object@strata)
    if (length(levs) > 6){
      levs <- paste(paste(head(levs), collapse = ", "), "...", sep = ", ")
    } else {
      levs <- paste(levs, collapse = ", ")
    }
    cat("a data frame with", length(object@strata), "columns (", levs, ")")
  }
  
  if (!is.null(object@hierarchy)){
    optional <- TRUE
    cat("\n   @hierarchy:", paste(object@hierarchy, collapse = ""))
  }
  
  if(!is.null(object@other)){
    optional <- TRUE
    cat("\n   @other: ")
    cat("a list containing: ")
    cat(ifelse(is.null(names(object@other)), "elements without names", paste(names(object@other), collapse= ", ")), "\n")
  }
  
  if(!is.null(object@other$ind.metrics)){
    optional <- TRUE
    cat("    @other$ind.metrics: ")
    cat(ifelse(is.null(names(object@other$ind.metrics)), "elements without names", paste(names(object@other$ind.metrics), collapse= ", ")), "\n")
  }

  if(!is.null(object@other$ind.metrics)){
    optional <- TRUE
    cat("    @other$loc.metrics: ")
    cat(ifelse(is.null(names(object@other$loc.metrics)), "elements without names", paste(names(object@other$loc.metrics), collapse= ", ")), "\n")
  }
    
  if(!optional) cat("\n   - empty -")
  
  
    cat("   @other$latlon[g]:")
    if(!is.null(object@other$latlon)){
      if (nrow(object@other$latlon)==nInd(object))
    cat(" coordinates for all individuals are attached") else cat(" number of coordinates does not match number of individuals")
  } else cat(" no coordinates attached")
  cat("\n")
 
}) # end show method


#################
## subset dartR
#################

#' indexing dartR objects correctly...
#' 
#' @param x dartR object
#' @param i index for individuals
#' @param j index for loci
#' @param ... other parameters
#' @param pop list of populations to be kept
#' @param treatOther elements in other (and ind.metrics & loci.metrics) as indexed as well. default: TRUE
#' @param quiet warnings are suppressed. default: TRUE
#' @param drop reduced to a vector if a single individual/loci is selected. default: FALSE [should never set to TRUE]


## dartR
setMethod("[", signature(x = "dartR", i = "ANY", j = "ANY", drop = "ANY"),
          function(x, i, j, ..., pop=NULL, treatOther=TRUE, quiet=TRUE, drop=FALSE) {
            if (missing(i)) i <- TRUE
            if (missing(j)) j <- TRUE
            
            ori.n <- nInd(x)
            ori.p <- nLoc(x)
            
            ## recycle logicals if needed
            if(!is.null(i) && is.logical(i)) i <- rep(i, length=ori.n)
            if(!is.null(j) && is.logical(j)) j <- rep(j, length=ori.p)
            
            if (!is.null(pop) && !is.null(pop(x))){
              i <- .get_pop_inds(x, pop)
            }
            
            ## SUBSET INDIVIDUALS ##
            ## genotypes
            x@gen <- x@gen[i]
            
            ## ind names
            x@ind.names <- x@ind.names[i]
            
            ## ploidy
            if(!is.null(x@ploidy)) {
              ori.ploidy <- ploidy(x) <- ploidy(x)[i]
            } else {
              ori.ploidy <- NULL
            }
            
            ## pop
            if(!is.null(pop(x))) {
              ori.pop <- pop(x) <- factor(pop(x)[i])
            } else {
              ori.pop <- NULL
            }
            
            ## strata
            if(!is.null(x@strata)) {
              ori.strata <- x@strata <- x@strata[i, , drop = FALSE]
            } else {
              ori.strata <- NULL
            }
            
            ## HANDLE 'OTHER' SLOT ##
            nOther <- length(other(x))
            namesOther <- names(other(x))
            counter <- 0
            if(treatOther & !(is.logical(i) && all(i))){
              f1 <- function(obj,n=ori.n){
                counter <<- counter+1
                if(!is.null(dim(obj)) && nrow(obj)==ori.n) { # if the element is a matrix-like obj
                  obj <- obj[i,,drop=FALSE]
                } else if(length(obj) == ori.n) { # if the element is not a matrix but has a length == n
                  obj <- obj[i]
                  if(is.factor(obj)) {obj <- factor(obj)}
                } else {if(!quiet) warning(paste("cannot treat the object",namesOther[counter]))}
                
                return(obj)
              } # end f1
              
              other(x) <- lapply(x@other, f1) # treat all elements
              
            } # end treatOther
            
            ## SUBSET LOCI ##
            
            ## handle ind.names, loc.names, chromosome, position, and alleles
            if (is.character(j)){
              j <- match(j, x@loc.names, nomatch = 0)
            }
            x@loc.names   <- x@loc.names[j]
            x@chromosome  <- chr(x)[j]
            x@position    <- position(x)[j]
            x@loc.all     <- alleles(x)[j]
            x@gen         <- lapply(x@gen, function(e) e[j])
            x@n.loc       <- x@gen[[1]]@n.loc
            #subset also loc.metrics (if this data.frame exists)
            if (!is.null(x@other$loc.metrics)) x@other$loc.metrics <- x@other$loc.metrics[j,]
            
            return(x)
          }) # end [] for genlight

###############################################################
#' adjust cbind for dartR
#' 
#' cbind is a bit lazy and does not take care for the metadata (so data in the 
#' other slot is lost). You can get most of the loci metadata back using 
#' gl.compliance.check.
#' @param ... list of dartR objects
#' @export 
cbind.dartR <- function(...){
  ## store arguments
  dots <- list(...)
  
  ## extract arguments which are genlight objects
  myList <- dots[sapply(dots, inherits, "genlight")]
  
  ## keep the rest in 'dots'
  dots <- dots[!sapply(dots, inherits, "genlight")]
  if(length(myList)==1 && is.list(myList[[1]])) myList <- myList[[1]]
  if(!all(sapply(myList, function(x) inherits(x,"genlight")))) stop("some objects are not genlight objects")
  ## remove empty objects
  myList <- myList[sapply(myList,nLoc)>0 & sapply(myList,nInd)>0]
  if(length(myList)==0) {
    warning("All objects are empty")
    return(NULL)
  }
  
  ## different checks
  if(length(unique(sapply(myList, nInd))) > 1 ) stop("objects have different numbers of individuals")
  n.obj <- length(myList)
  n.ind <- nInd(myList[[1]])
  if(n.ind==0){
    warning("All objects are empty")
    return(NULL)
  }
  temp <- as.matrix(as.data.frame(lapply(myList, ploidy)))
  if(any(apply(temp,1,function(r) length(unique(r)))>1)) stop("non-consistent ploidy across datasets")
  ori.ploidy <- ploidy(myList[[1]])
  
  ## merge one individual at a time ##
  res <- list()
  for(i in 1:n.ind){
    res[[i]] <- Reduce(function(a,b) {cbind(a,b,checkPloidy=FALSE)}, lapply(myList, function(e) e@gen[[i]]) )
  }
  
  dots$gen <- res
  dots$Class <- "dartR"
  res <- do.call(new, dots)
  
  ## handle loc.names, alleles, etc. ##
  indNames(res) <- indNames(myList[[1]])
  locNames(res) <- unlist(lapply(myList, locNames))
  alleles(res) <- unlist(lapply(myList, alleles))
  pop(res) <- pop(myList[[1]])
  res@strata <- myList[[1]]@strata
  ploidy(res) <- ori.ploidy
  
  ## return object ##
  return(res)
} # end cbind.dartR

#' adjust rbind for dartR
#' 
#' rbind is a bit lazy and does not take care for the metadata (so data in the 
#' other slot is lost). You can get most of the loci metadata back using
#'  gl.compliance.check.
#' @param ... list of dartR objects
#' @export 
rbind.dartR <- function(...){
  ## store arguments
  dots <- list(...)
  
  ## extract arguments which are genlight objects
  myList <- dots[sapply(dots, inherits, "genlight")]
  
  ## keep the rest in 'dots'
  dots <- dots[!sapply(dots, inherits, "genlight")]
  
  if(!all(sapply(myList, function(x) inherits(x,"genlight")))) stop("some objects are not genlight objects")
  
  ## remove empty objects
  myList <- myList[sapply(myList,nLoc)>0 & sapply(myList,nInd)>0]
  if(length(myList)==0) {
    warning("All objects are empty")
    return(NULL)
  }
  
  if(length(unique(sapply(myList, nLoc))) !=1 ) stop("objects have different numbers of SNPs")
  
  ## build output
  dots$Class <- "dartR"
  dots$gen <- Reduce(c, lapply(myList, function(e) e@gen))
  res <- do.call(new, dots)
  locNames(res) <- locNames(myList[[1]])
  alleles(res)  <- alleles(myList[[1]])
  indNames(res) <- unlist(lapply(myList, indNames))
  pop(res)      <- factor(unlist(lapply(myList, pop)))
  
  
  #hierachies are ignored in dart objects here
  # Hierarchies are tricky. Using dplyr's bind_rows.
  #res <- .rbind_strata(myList, res)
  
  ## return object ##
  return(res)
  
} # end rbind.genlight



