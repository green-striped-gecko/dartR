#' @name gl2plink
#' @title Convert a genlight object into PLINK format
#' @description
#' This function exports a genlight object into PLINK format and save it into a
#' file.
#' This function produces the following PLINK files: bed, bim, fam, ped and map.
#' @param x Name of the genlight object containing the SNP data [required].
#' @param plink_path Path of PLINK binary file [default getwd())].
#' @param bed_file Whether create PLINK files .bed, .bim and .fam
#' [default FALSE].
#' @param outfile File name of the output file [default 'gl_plink'].
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
#' To create PLINK files .bed, .bim and .fam (bed_file = TRUE), it is necessary
#' to download the binary file of PLINK 1.9 and provide its path (plink_path).
#' The binary file can be downloaded from:
#' \url{https://www.cog-genomics.org/plink/}
#'
#' The chromosome of each SNP can be a character or numeric. The chromosome
#' information for unmapped SNPS is coded as 0.
#' Family ID is taken from  x$pop.
#' Within-family ID (cannot be '0') is taken from indNames(x).
#' Variant identifier is taken from locNames(x).
#' @return NULL
#' @references
#' Purcell, Shaun, et al. 'PLINK: a tool set for whole-genome association and
#' population-based linkage analyses.' The American journal of human genetics
#' 81.3 (2007): 559-575.
#' @export
#' @author Custodian: Luis Mijangos (Post to
#'  \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' \donttest{
#' gl2plink(platypus.gl,snp_pos='ChromPos_Platypus_Chrom_NCBIv1',
#' snp_chr = 'Chrom_Platypus_Chrom_NCBIv1')
#' }

gl2plink <- function(x,
                     plink_path = getwd(),
                     bed_file = FALSE,
                     outfile = "gl_plink",
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
    
    outfilespec <- file.path(outpath, outfile)
    
    snp_temp <- x$other$loc.metrics
    
    if (snp_chr == "0") {
        snp_temp$chrom <- 0
    } else {
        if (chr_format == "numeric") {
            snp_temp$chrom <- as.numeric(unname(unlist(snp_temp[snp_chr])))
        }
        if (chr_format == "character") {
            snp_temp$chrom <- as.character(unname(unlist(snp_temp[snp_chr])))
        }
    }
    
    if (snp_pos == "0") {
        snp_temp$snp_pos <- 0
    } else {
        snp_temp$snp_pos <- as.numeric(unname(unlist(snp_temp[snp_pos])))
    }
    
    # Convert any NA values to 0
    snp_temp[is.na(snp_temp$snp_pos), "snp_pos"] <- 0
    # Convert any NA values to 0
    snp_temp[snp_temp$snp_pos == 0, "chrom"] <- 0
    
    # Chromosome code
    snp_chr <- snp_temp$chrom
    # Variant identifier
    var_id <- locNames(x)
    
    # Base-pair coordinate
    pos_bp <- snp_temp$snp_pos
    
    gl_map <- cbind(snp_chr, var_id, pos_cM, pos_bp)
    
    write.table(
        gl_map,
        file = paste0(outfilespec, ".map"),
        quote = F,
        row.names = F,
        col.names = F
    )
    
    ########## .fam (PLINK sample information file)
    
    sample.id_temp <- indNames(x)
    sample.id_temp <-
        gsub(" ", replacement = "_", sample.id_temp)
    
    # Family ID ('FID')
    FID <- as.character(x$pop)
    FID <- gsub(" ","_",FID)
    # Within-family ID ('IID'; cannot be '0')
    IID <- sample.id_temp
    IID <- gsub(" ","_",IID)

    # Sex code ('1' = male, '2' = female, '0' = unknown)
    if (length(sex_code) > 1) {
        sex_code <- as.character(sex_code)
        sex_code[startsWith(sex_code, "f") |
                     startsWith(sex_code, "F")] <- "2"
        sex_code[startsWith(sex_code, "m") |
                     startsWith(sex_code, "M")] <- "1"
        sex_code[startsWith(sex_code, "u") |
                     startsWith(sex_code, "U")] <- "0"
        sex_code[nchar(sex_code) == 0] <- "0"
    }
    
    gl_fam <-
        cbind(FID, IID, ID_dad, ID_mom, sex_code, phen_value)
    
    x_mat <- as.matrix(x[, ])
    homs1 <-
        paste(substr(x@loc.all, 1, 1), "/", substr(x@loc.all, 1, 1), sep = "")
    hets <- x@loc.all
    homs2 <-
        paste(substr(x@loc.all, 3, 3), "/", substr(x@loc.all, 3, 3), sep = "")
    xx <- matrix(NA, ncol = ncol(x_mat), nrow = nrow(x_mat))
    for (i in 1:nrow(x_mat)) {
        for (ii in 1:ncol(x_mat)) {
            inp <- x_mat[i, ii]
            if (!is.na(inp)) {
                if (inp == 0)
                    xx[i, ii] <- homs1[ii]
                else if (inp == 1)
                    xx[i, ii] <- hets[ii]
                else if (inp == 2)
                    xx[i, ii] <- homs2[ii]
            } else
                xx[i, ii] = "0/0"
        }
    }
    xx <- gsub("/", " ", xx)
    xx <- apply(xx, 1, paste, collapse = " ")
    xx <- cbind(gl_fam, xx)
    
    write.table(
        xx,
        file = paste0(outfilespec, ".ped"),
        quote = F,
        row.names = F,
        col.names = F
    )
    
    if (bed_file) {
        prefix.in_temp <- outfilespec
        prefix.out_temp <- outfilespec
        
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
                        if (autosome.only)
                            "--autosome"
                        else
                            "",
                        "--allow-no-sex",
                        "--keep-allele-order",
                        "--out",
                        prefix.out,
                        extra.options
                    )
                )
                bedfile.out
            }
        
        system_verbose = function(...) {
            report = system(..., intern = T)
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
    }
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    invisible(NULL)
}
