#' @name gl2vcf
#' @title Converts a genlight object into vcf format
#' @description
#' This function exports a genlight object into VCF format and save it into a
#' file.
#' @param x Name of the genlight object containing the SNP data [required].
#' @param plink_path Path of PLINK binary file [default getwd())].
#' @param outfile File name of the output file [default 'gl_vcf'].
#' @param outpath Path where to save the output file
#' [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath='.'
#'  when calling this function to direct output files to your working directory.
#' @param snp_pos Field name from the slot loc.metrics where the SNP position is
#' stored [default '0'].
#' @param snp_chr Field name from the slot loc.metrics where the chromosome of
#' each is stored [default '0'].
#' @param chr_format Whether chromosome information is stored as 'numeric' or as
#' 'character', see details [default 'character'].
#' @param pos_cM A vector, with as many elements as there are loci, containing
#' the SNP position in morgans or centimorgans [default '0'].
#' @param ID_dad A vector, with as many elements as there are individuals,
#' containing the ID of the father, '0' if father isn't in dataset [default '0'].
#' @param ID_mom A vector, with as many elements as there are individuals,
#' containing the ID of the mother, '0' if mother isn't in dataset [default '0'].
#' @param sex_code A vector, with as many elements as there are individuals,
#' containing the sex code ('male', 'female', 'unknown') [default  'unknown'].
#' @param phen_value A vector, with as many elements as there are individuals,
#' containing the phenotype value. '1' = control, '2' = case, '0' = unknown
#' [default '0'].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @details
#' This function requires to download the binary file of PLINK 1.9 and provide
#' its path (plink_path).
#' The binary file can be downloaded from:
#' \url{https://www.cog-genomics.org/plink/}
#'
#' The chromosome information for unmapped SNPS is coded as 0.
#' Family ID is taken from  x$pop
#' Within-family ID (cannot be '0') is taken from indNames(x)
#' Variant identifier is taken from locNames(x)
#' @return  returns no value (i.e. NULL)
#' @references
#' Danecek, P., Auton, A., Abecasis, G., Albers, C. A., Banks, E., DePristo, M.
#'  A., ... & 1000 Genomes Project Analysis Group. (2011). The variant call
#'  format and VCFtools. Bioinformatics, 27(15), 2156-2158.
#' @export
#' @author Custodian: Luis Mijangos (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' \dontrun{
#' require("dartR.data")
#' gl2vcf(platypus.gl,snp_pos='ChromPos_Platypus_Chrom_NCBIv1',
#'  snp_chr = 'Chrom_Platypus_Chrom_NCBIv1')
#' }

gl2vcf <- function(x,
                   plink_path = getwd(),
                   outfile = "gl_vcf",
                   outpath = tempdir(),
                   snp_pos = "0",
                   snp_chr = "0",
                   chr_format = "character",
                   pos_cM = "0",
                   ID_dad = "0",
                   ID_mom = "0",
                   sex_code = "unknown",
                   phen_value = "0",
                   verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbosity = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # DO THE JOB
    
    # assigning SNP position information 
    if(snp_pos == "0"){
      x$position <- rep(as.integer(0),nLoc(x))
      
    }else{
      
      if(snp_pos %in% names(x$other$loc.metrics)){
        
        if(verbose>=2){
          cat(report("  Using the SNP position information in the field",snp_chr, "from loc.metrics.\n"))
        }
        
        x$position <- unname(unlist(x$other$loc.metrics[snp_pos]))
        
      }else{
        
        stop(error("  The field",snp_pos, "with the SNP position information is not present in loc.metrics.\n"))
        
      }
    }
    
    # assigning chromosome information 
    if(is.null(x$chromosome)){
   
      if(snp_chr == "0" ){
        x$chromosome <- rep(as.factor("0"),nLoc(x))
        
        if(verbose>=2){
          
          cat(report("  Chromosome information is not present in the slot 'chromosome'. Setting '0' as the name chromosome for all the SNPs.\n"))
          
        }
        
      }else{
        
        if(snp_chr %in% names(x$other$loc.metrics)){
          
          if(verbose>=2){
            cat(report("  Using the chromosome information in the field",snp_chr, "from loc.metrics.\n"))
          }
          
          x$chromosome <- as.factor(unname(unlist(x$other$loc.metrics[snp_chr])))
          
        }else{
          
          stop(error("  The field",snp_chr, "with the chromosome information is not present in loc.metrics.\n"))
          
        }
      }
    }
   
    
    gl2plink(
        x = x,
        outfile = "gl_plink_temp",
        outpath = outpath,
        chr_format = chr_format,
        pos_cM = pos_cM,
        ID_dad = ID_dad,
        ID_mom = ID_mom,
        sex_code = sex_code,
        phen_value = phen_value,
        verbose = NULL
    )
    
    prefix.in_temp <- paste0(outpath, "/gl_plink_temp")
    prefix.out_temp <- file.path(outpath, outfile)
    
    allele_tmp <- gsub("/"," ", x$loc.all)
    allele_tmp <- strsplit(allele_tmp,split = " ")
    allele_tmp <- Reduce(rbind,allele_tmp)[,2]
    allele_tmp <- cbind(locNames(x), allele_tmp)
    write.table(allele_tmp,
                file = file.path(tempdir(),"mylist.txt"),
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE
    )
    
    make_plink <-
        function(plink.path,
                 prefix.in = prefix.in_temp,
                 prefix.out = prefix.out_temp,
                 autosome.only = FALSE,
                 extra.options = "") {
            bedfile.out <- paste0(prefix.out, ".bed")
            system_verbose(
                paste(
                    plink.path,
                    "--file",
                    prefix.in,
                    "--recode",
                    "vcf",
                    if (autosome.only)
                        "--autosome"
                    else
                        "",
                    "--allow-no-sex",
                    paste("--reference-allele",file.path(tempdir(),'mylist.txt')),
                    # "--keep-allele-order",
                    # "--real-ref-alleles",
                    # paste("--a1-allele", file.path(outpath,'alleles.csv'),"1"),
                    # paste("--a2-allele", file.path(outpath,'alleles.csv'),"2"),
                    "--out",
                    prefix.out,
                    extra.options
                )
            )
        }
    
    system_verbose <-function(...) {
        report <-system(..., intern = T)
        message(
            paste0(
                "\n\n----------Output of function start:\n\n",
                paste(report, collapse = "\n"),
                "\n\n----------Output of function finished...\n\n"
            )
        )
    }
    
    make_plink(plink.path = paste0(plink_path, "/plink"),
               extra.options = "--aec")
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    invisible(NULL)
}
