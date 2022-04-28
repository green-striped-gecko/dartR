#' structure util functions
#' 
#' These functions were copied from package strataG, which is no longer on CRAN (maintained by Eric Archer)
#' @export
#' @author Bernd Gruber (bugs? Post to
#'  \url{https://groups.google.com/d/forum/dartr}); original implementation of
#'   Eric Archer \url{https://github.com/EricArcher/strataG}
#' @param x a genind object
#' @return a gtypes object 


utils.structure.genind2gtypes <- 
function (x) 
{
  .parseLocusNames <- function(locus.names, ploidy) {
    if(ploidy == 1) return(locus.names)
    loc.i <- matrix(1:length(locus.names), nrow = ploidy)
    apply(loc.i, 2, function(i) {
      this.loc <- locus.names[i]
      max.length <- max(sapply(this.loc, nchar))
      # find location of first difference
      ptr <- 1
      all.same <- length(unique(substr(this.loc, 1, ptr))) == 1
      while(all.same & ptr < max.length) {
        ptr <- ptr + 1
        all.same <- length(unique(substr(this.loc, 1, ptr))) == 1
      }
      ptr <- ptr - 1
      if(ptr == 0) {
        this.loc[1] # return first name if all names are different
      } else {
        # remove last same character if it is not alphanumeric
        if(!substr(this.loc[1], ptr, ptr) %in% c(LETTERS, letters, 0:9)) {
          ptr <- ptr - 1
        }
        substr(this.loc[1], 1, ptr)
      }
    })
  }
  
  
############################################################  
  togtypes <- function( gen.data, ploidy, ind.names = NULL,
                        sequences = NULL, strata = NULL, schemes = NULL,
                        description = NULL, other = NULL,
                        remove.sequences = FALSE) {
    
    if(is.null(gen.data) | is.null(ploidy)) return(NULL)
    
    # check gen.data
    if(is.vector(gen.data)) {
      gen.data <- cbind(gen.data)
      colnames(gen.data) <- "Haplotype"
    }
    if(!(is.matrix(gen.data) | is.data.frame(gen.data))) {
      stop("'gen.data' is not a vector, matrix, or data.frame", call. = FALSE)
    }
    gen.data <- as.matrix(gen.data)
    
    # check ploidy
    ploidy <- as.integer(ploidy)
    if(ncol(gen.data) %% ploidy != 0) {
      stop(
        "the number of columns in 'gen.data' is not a multiple of 'ploidy'",
        call. = FALSE
      )
    }
    if(ploidy > 1 & !is.null(sequences)) {
      stop(
        "'sequences' can't be present if 'ploidy' is greater than 1", 
        call. = FALSE
      )
    }
    
    # check ind.names
    ind.names <- if(!is.null(ind.names)) {
      if(!is.vector(ind.names)) stop("'ind.names' must be a vector")
      if(length(ind.names) != nrow(gen.data)) {
        stop(
          "the length of 'ind.names' must equal the number of rows in 'gen.data'",
          call. = FALSE
        )
      }
      as.character(ind.names)
    } else {
      if(!is.null(rownames(gen.data))) rownames(gen.data) else 1:nrow(gen.data)
    }
    if(any(duplicated(ind.names))) {
      dup.names <- unique(ind.names[duplicated(ind.names)])
      dup.names <- paste(dup.names, collapse = ", ")
      stop("there are duplicated individual names: ", dup.names, call. = FALSE)
    }
    rownames(gen.data) <- ind.names
    
    # check strata
    if(!is.null(strata)) {
      if(is.null(names(strata))) {
        if(length(strata) == length(ind.names)) names(strata) <- ind.names
      } 
    } else strata <- rep("Default", nrow(gen.data))
    if(length(strata) != nrow(gen.data)) {
      warning("the length of 'strata' is not the same as the number of individuals. strata will be recycled.")
    }
    
    # check schemes
    if(!is.null(schemes)) {
      # check that schemes is a data.frame
      if(!(is.matrix(schemes) | is.data.frame(schemes))) {
        stop("'schemes' is not a matrix or data.frame", call. = FALSE)
      } 
      schemes <- as.data.frame(schemes)
      
      # check that 'id' column in schemes exists
      if(!"id" %in% colnames(schemes)) {
        if(is.null(rownames(schemes))) {
          if(nrow(schemes) != nrow(gen.data)) {
            stop(
              "'schemes' doesn't have an 'id' column or rownames and is not as long as 'gen.data'",
              call. = FALSE
            )
          } else {
            schemes$id <- 1:nrow(schemes)
          }
        } else {
          schemes$id <- rownames(schemes)
        }
      }
      rownames(schemes) <- NULL
      schemes <- dplyr::select(schemes, .data$id, dplyr::everything())
      
      # check that ids in schemes can be found
      if(length(intersect(schemes$id, rownames(gen.data))) == 0) {
        stop(
          "no ids in 'schemes' are in 'gen.data' or 'ind.names'", 
          call. = FALSE
        )
      }
      
      schemes$id <- as.character(schemes$id)
    }
    
    # check description
    if(is.null(description)) {
      description <- paste(
        "gtypes created on", format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      )
    }
    
    # create locus names
    nloc <- ncol(gen.data) / ploidy
    if(is.null(colnames(gen.data))) {
      # return generic names if no colnames assigned
      nums <- formatC(1:nloc, digits = floor(log10(nloc)), flag = "0")
      generic.locus.names <- paste0("Locus", "_", nums)
      colnames(gen.data) <- .expandLocusNames(generic.locus.names, ploidy)
    } 
    locus.names.lookup <- stats::setNames(
      rep(.parseLocusNames(colnames(gen.data), ploidy), each = ploidy),
      colnames(gen.data)
    )
    
    # check sequences
    if(!is.null(sequences)) {
      sequences <- as.multidna(sequences)
      if(getNumLoci(sequences) != ncol(gen.data)) {
        stop(
          "the number of genes in 'sequences' is not equal to the number of loci",
          call. = FALSE
        )
      }
      setLocusNames(sequences) <- colnames(gen.data)
      for(loc in colnames(gen.data)) {
        haps <- stats::na.omit(unique(as.character(gen.data[, loc])))
        seq.names <- apex::getSequenceNames(sequences)[[loc]]
        missing <- setdiff(haps, seq.names)
        if(length(missing) > 0) {
          stop(
            "the following haplotypes can't be found in sequences for locus '", 
            loc, "': ", paste(missing, collapse = ", "), 
            call. = FALSE
          )
        }
      }
    }
    
    gen.data <- cbind(
      id = rownames(gen.data), 
      stratum = as.character(strata), 
      gen.data
    ) %>% 
      as.data.frame(stringsAsFactors = FALSE) %>% 
      tidyr::gather("locus", "allele", -.data$id, -.data$stratum) %>% 
      dplyr::mutate(locus = locus.names.lookup[as.character(.data$locus)])
    
    data.table::setDT(gen.data, key = c("id", "stratum", "locus"))
    
    # create and return gtypes object
    #g <- .Object
    g<-list()
    g$data <- gen.data
    g$ploidy <- ploidy
    g$sequences <- sequences
    g$schemes <- schemes
    g$description <- description
    g$other <- if(is.null(other)) list() else other
    
    # Check for samples missing data for all loci
    #g <- .removeIdsMissingAllLoci(g)
    
    # Remove unreferenced sequences
    #if(remove.sequences) g <- removeUnusedSequences(g)
    
    g
  }    
  
  
  
##################################################
  df2gtypes <- function (x, ploidy, id.col = 1, strata.col = 2, loc.col = 3, 
            sequences = NULL, schemes = NULL, description = NULL, other = NULL) 
  {
    if (!(is.matrix(x) | is.data.frame(x))) {
      stop("'x' must be a matrix or data.frame")
    }
    x <- as.data.frame(x)
    if (!is.null(id.col)) {
      if (is.character(id.col)) 
        id.col <- match(id.col, colnames(x))
      id.col <- as.numeric(id.col)
    }
    if (!is.null(strata.col)) {
      if (is.character(strata.col)) 
        strata.col <- match(strata.col, colnames(x))
      strata.col <- as.numeric(strata.col)
    }
    loc.col <- as.numeric(loc.col)
    if (loc.col < max(id.col, strata.col, loc.col)) {
      stop("'loc.col' must be greater than 'id.col' and 'strata.col'")
    }
    ind.names <- if (is.null(id.col)) {
      if (is.null(rownames(x))) 
        1:nrow(x)
      else rownames(x)
    }    else {
      x[, id.col]
    }
    ind.names <- as.character(ind.names)
    strata <- if (is.null(strata.col))  NULL     else x[, strata.col]
    gen.data <- x[, loc.col:ncol(x), drop = FALSE]
    togtypes( gen.data = gen.data, ploidy = ploidy, 
                 ind.names = ind.names, strata = strata, schemes = schemes, 
                 sequences = sequences, description = description, other = if (is.null(other)) 
                   list()
                 else other)
  }
##################################################  
  
  
  
  gen.mat <- adegenet::genind2df(x, usepop = TRUE, oneColPerAll = TRUE)
  gen.mat[gen.mat == "NA"] <- NA
  has.pop <- !is.null(x@pop)
  df2gtypes(x = gen.mat, ploidy = x@ploidy[1], id.col = NULL, 
            strata.col = if (has.pop) 1 else NULL, loc.col = if (has.pop) 2  else 1, schemes = x@strata, other = list(genind = adegenet::other(x)))
}
