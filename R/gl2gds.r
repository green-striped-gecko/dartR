#' @name gl2gds
#' @title Convert a genlight object to gds format
#' @description 
#' Package SNPRelate relies on a bit-level representation of a SNP dataset that 
#' competes with \{adegenet\} genlight objects and associated files. This 
#' function saves a genlight object to a gds format file.
#' @param x Name of the genlight object containing the SNP data [required].
#' @param outfile File name of the output file (including extension) 
#' [default gl2gds.gds].
#' @param outpath Path where to save the output file 
#' [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath="."
#'  when calling this function to direct output files to your working directory.
#' @param snp_pos Field name from the slot loc.metrics where the SNP position is 
#' stored [default NULL].
#' @param snp_chr Field name from the slot loc.metrics where the chromosome of 
#' each is is stored [default NULL].
#' @param chr_format Whether chromosome information is stored as "numeric" or as 
#' "character", see details [default "character"].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity].
#' @details 
#' This function orders the SNPS by chromosome and by position before converting 
#' to SNPRelate format, as required by this package. 
#' 
#' The chromosome of each SNP can be a character or numeric, as described in the 
#' vignette of SNPRelate:
#' "snp.chromosome, an integer or character mapping for each chromosome. 
#' Integer: numeric values 1-26, mapped in order from 1-22, 23=X, 24=XY 
#' (the pseudoautosomal region), 25=Y, 26=M (the mitochondrial probes), and 0 
#' for probes with unknown positions; it does not allow NA. Character: “X”,
#'  “XY”, “Y” and “M” can be used here, and a blank string indicating unknown 
#'  position."
#'  
#' When using some functions from package SNPRelate with datasets other than 
#' humans it might be necessary to use the option autosome.only=TRUE to avoid 
#' detecting chromosome coding, so we recommended to read the documentation of 
#' each SNPRelate function. 
#' 
#' The chromosome information for unmapped SNPS is coded as 0, as required by 
#' SNPRelate. 
#' 
#' Remember to close the GDS file before working in a different GDS object with 
#' the function \link[SNPRelate]{snpgdsClose} (package SNPRelate).
#' @return NULL
#' @export
#' @importFrom SNPRelate snpgdsCreateGeno snpgdsOpen snpgdsSummary snpgdsClose
#' @author Custodian: Luis Mijangos (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' \donttest{
#' gl2gds(platypus.gl,snp_pos="ChromPos_Platypus_Chrom_NCBIv1", snp_chr = "Chrom_Platypus_Chrom_NCBIv1")
#' }

gl2gds <- function(x, 
                   outfile = "gl2gds.gds",
                   outpath = tempdir(), 
                   snp_pos = NULL,
                   snp_chr = NULL,
                   chr_format = "character",
                   verbose = NULL) {
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func=funname,build="Jody", verbosity =verbose)
  
  # CHECK DATATYPE 
  datatype <- utils.check.datatype(x, verbose=verbose)
  
# DO THE JOB
  
  #ordering loc.metrics by chromosome and snp position 
  snp_order_temp <- x$other$loc.metrics
  # using CloneID as SNP ID
  snp_order_temp$snp_id <- as.numeric(snp_order_temp$CloneID)
  #adding snps order to be used to order snp matrix 
  snp_order_temp$order <- 1:nLoc(x)
  

  if(is.null(snp_chr)){
    snp_order_temp$chrom <- 0
  }else{
    if(chr_format=="numeric"){
      snp_order_temp$chrom <- as.numeric(unname(unlist(snp_order_temp[snp_chr])))
    }
    if(chr_format=="character"){
      snp_order_temp$chrom <- as.character(unname(unlist(snp_order_temp[snp_chr])))
    }
    
  }
  
  if(is.null(snp_pos)){
    snp_order_temp$snp_pos <- 0
  }else{
    snp_order_temp$snp_pos <- as.numeric(unname(unlist(snp_order_temp[snp_pos])))
  }

  # Convert any NA values to 0 (genlight objects have NA for missing; SNPRelate has 0 in this instance)  
  snp_order_temp[is.na(snp_order_temp$snp_pos),"snp_pos"] <- 0
  # Convert any NA values to 0 (genlight objects have NA for missing; SNPRelate has 0 in this instance)  
  snp_order_temp[snp_order_temp$snp_pos == 0,"chrom"] <- 0
  snp_order_temp <- snp_order_temp[with(snp_order_temp, order(chrom,snp_pos)),]

  #ordering snp matrix based on 
  genmat_temp <- t(as.matrix(x))
  genmat_temp <- genmat_temp[order(snp_order_temp$order),]
  genmat_temp[is.na(genmat_temp)] <- 3

  snp.id_temp <- locNames(x)
  snp.id_temp <- snp.id_temp[order(snp_order_temp$order)]
  
  snp.allele_temp <- x@loc.all
  snp.allele_temp <-   snp.allele_temp[order(snp_order_temp$order)]
  
  sample.id_temp <- indNames(x)
  sample.id_temp <- gsub(" ",replacement = "_",sample.id_temp)
  
  geno_list <- list(sample.id =sample.id_temp ,snp.id=snp.id_temp, snp.position=snp_order_temp$snp_pos, snp.chromosome=snp_order_temp$chrom ,snp.allele = snp.allele_temp ,genotype=genmat_temp)
  
# Create the gds file
  if (verbose >= 2){
    cat(report("  Converting SNP data to gds format\n"))
    }
  
  # create a gds file
  with(geno_list, SNPRelate::snpgdsCreateGeno(outfile, genmat=genotype,
                                     sample.id=sample.id, snp.id=snp.id, snp.chromosome=snp.chromosome,
                                     snp.position=snp.position, snp.allele=snp.allele, snpfirstdim=TRUE))
  
# Open the GDS file, which will print out a summary of contents
  if (verbose >= 2 ) {
    cat(report("  Writing data to file",outfilespec,"\n"))
    }
  genofile <- SNPRelate::snpgdsOpen(outfilespec)
  cat(important("Structure of gds file\n\n"))
  SNPRelate::snpgdsSummary(genofile)
  print(genofile)

# Close the GDS file
  if (verbose >= 2 ) {
    cat(report("  Closing file",outfilespec,"\n"))
    }
  SNPRelate::snpgdsClose(genofile)

# FLAG SCRIPT END

  if (verbose > 0) {
    cat(report("Completed:",funname,"\n"))
  }

  invisible(NULL)

}

