gl2genepop <- function (xx, output = "data.frame") 
{
  if (!inherits(xx, "genlight")) {
    stop("Input 'x' must be an object of class 'genind'.")
  }
  x <- gl2gi(xx, verbose = 0,probar = FALSE)
  data <- as.matrix(x[order(pop(x)),])
  
  
  pop_names <- x@pop
  
  
  
  if (all(unlist(unique(x@all.names)) %in% c("A", "T", 
                                             "C", "G"))) {
    m_type <- "snp"
    #message("Your dataset is treated as a SNP dataset.\n            Alleles initially coded A, T, C, G were respectively coded\n            01, 02, 03 and 04")
    colnames(data) <- gsub(colnames(data), pattern = "\\.A", 
                           replacement = ".01")
    colnames(data) <- gsub(colnames(data), pattern = "\\.T", 
                           replacement = ".02")
    colnames(data) <- gsub(colnames(data), pattern = "\\.C", 
                           replacement = ".03")
    colnames(data) <- gsub(colnames(data), pattern = "\\.G", 
                           replacement = ".04")
  }
  
  loci_names_l <- x@loc.fac
  loc_all <- data.frame(col = colnames(data))
  if (all(stringr::str_count(colnames(data), "\\.") == 
          1) != TRUE) {
    stop("The columns' names of x@tab must be of form 'locus.allele' with only 1\n         '.' between locus and allele")
  }
  loc_all <- tidyr::separate(loc_all, col = 1, sep = "\\.", 
                             into = c("locus", "allele"))
  loci_names <- as.character(loci_names_l[-which(duplicated(loci_names_l))])
  n.loci <- length(loci_names_l[-which(duplicated(loci_names_l))])
  data_gpop <- data.frame(id = paste(pop_names, "_", 
                                     row.names(data), ",", sep = ""))
  for (i in 1:n.loci) {
    loc <- loci_names[i]
    a <- c()
    for (j in 1:nrow(data)) {
      col_loc <- which(loc_all[, "locus"] == loc)
      hom <- which(data[j, col_loc] == 2)
      het <- which(data[j, col_loc] == 1)
      if (length(hom) != 0) {
        a[j] <- paste(loc_all[col_loc[hom], "allele"], 
                      loc_all[col_loc[hom], "allele"], sep = "")
      }
      else if (length(het) != 0) {
        if (as.character(loc_all[col_loc[het[1]], "allele"]) < 
            as.character(loc_all[col_loc[het[2]], "allele"])) {
          a[j] <- paste(loc_all[col_loc[het[1]], "allele"], 
                        loc_all[col_loc[het[2]], "allele"], 
                        sep = "")
        }
        else {
          a[j] <- paste(loc_all[col_loc[het[2]], "allele"], 
                        loc_all[col_loc[het[1]], "allele"], 
                        sep = "")
        }
      }
      else {
        if (nchar(loc_all[1, "allele"] == 6)) {
          a[j] <- "000000"
        }
        else {
          a[j] <- "0000"
        }
      }
    }
    data_gpop <- cbind(data_gpop, a)
  }
  colnames(data_gpop) <- c("ID", as.character(loci_names))
  data_gpop[, ] <- apply(data_gpop,c(1,2), as.character)
  
  
  dummy<- paste("Genepop output. Loci:", nLoc(x), "Populations:",nPop(x))
  dummy[2] <- paste(locNames(x),collapse=",")
  cs <- c(cumsum(table(pop(x))))
  from=c(1,(cs[-length(cs)]+1))
  to= cs
  for (i in 1:nPop(x))
  {
    da <- apply(data_gpop[from[i]:to[i],], 1, function(y) paste(y, collapse = " "))  
    dummy <- c(dummy,"Pop",da)
  }
  
  data_gpop2<- dummy
  
  if (output == "data.frame") {
    return(data.frame(lines=data_gpop2))
  }
  else if (stringr::str_sub(output, -4, -1) == ".gen") {
    utils::write.table(data_gpop2, file = output, quote = FALSE, 
                       row.names = FALSE, col.names = FALSE)
  }
  else {x
    stop("'output' parameter must be 'data.frame' (then the function returns an\nobject of class data.frame) or a path and file name ending in '.txt'. In the\n          latter case, the functions creates a text file in the directory\n         specified by the path or in the current working directory when a\n         path is not specified.")
  }
}

