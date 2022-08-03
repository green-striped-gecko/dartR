WC_FST_Diploids_2Alleles <- function(Sample_Mat) {
    ## Calculate both Fst and Fst NoCorr at the same time, from WC84 Sample Mat
    # has three columns (homo_p,m heterozygotes, and homo_q) and
    ## a row for each population
    
    sample_sizes <-rowSums(Sample_Mat)
    n_ave <-mean(sample_sizes)
    n_pops <-nrow(Sample_Mat)  #r
    r <-n_pops
    n_c <-(n_pops * n_ave - sum(sample_sizes ^ 2) / (n_pops * n_ave)) /
        (n_pops - 1)
    p_freqs <-(Sample_Mat[, 1] + Sample_Mat[, 2] / 2) / sample_sizes
    p_ave <-sum(sample_sizes * p_freqs) / (n_ave * n_pops)
    
    s2 <-sum(sample_sizes * (p_freqs - p_ave) ^ 2) / ((n_pops - 1) * n_ave)
    if (s2 == 0) {
        return(0)
    }
    
    h_freqs <-Sample_Mat[, 2] / sample_sizes
    h_ave <-sum(sample_sizes * h_freqs) / (n_ave * n_pops)
    
    a <-
        n_ave / n_c * (s2 - 1 / (n_ave - 1) * (p_ave * (1 - p_ave) - ((r - 1) /
                                                                          r) * s2 - (1 / 4) * h_ave))
    
    b <-
        n_ave / (n_ave - 1) * (p_ave * (1 - p_ave) - (r - 1) / r * s2 - (2 * n_ave - 1) /
                                   (4 * n_ave) * h_ave)
    
    c <- 1 / 2 * h_ave
    
    aNoCorr <- n_ave / n_c * (s2)
    
    bNoCorr <-
        (p_ave * (1 - p_ave) - (r - 1) / r * s2 - (2 * n_ave) / (4 * n_ave) * h_ave)
    
    cNoCorr <- 1 / 2 * h_ave
    
    He <- 1 - sum(p_ave ^ 2, (1 - p_ave) ^ 2)
    
    FST <- a / (a + b + c)
    FSTNoCorr <-aNoCorr / (aNoCorr + bNoCorr + cNoCorr)
    
    return(
        list(
            He = He,
            FST = FST,
            T1 = a,
            T2 = (a + b + c),
            STNoCorr = FSTNoCorr,
            T1NoCorr = aNoCorr,
            T2NoCorr = (aNoCorr + bNoCorr + cNoCorr),
            meanAlleleFreq = p_ave
        )
    )
}

#' Creates OutFLANK input file from individual genotype info.
#' @param SNPmat This is an array of genotypes with a row for each individual.
#' There should be a column for each SNP, with the number of copies of the focal
#'  allele (0, 1, or 2) for that individual. If that individual is missing data
#'  for that SNP, there should be a 9, instead.
#' @param locusNames A list of names for each SNP locus. There should be the
#' same number of locus names as there are columns in SNPmat.
#' @param popNames A list of population names to give location for each
#' individual. Typically multiple individuals will have the same popName. The
#' list popNames should have the same length as the number of rows in SNPmat.
#' @return Returns a data frame in the form needed for the main OutFLANK
#' function.
#' @export

utils.outflank.MakeDiploidFSTMat <- function(SNPmat, 
                                             locusNames,
                                             popNames) {
        # SNPmat is a matrix with individuals in rows and snps in columns 0, 1,
    #or 2 represent the number of copies of the focal allele, and 9
        # is for missing data locusNames is a character vector of names of each 
    #SNP popNames is a character vector with the population
        # identifier for each individual
        
        locusname <- unlist(locusNames)
        popname <- unlist(popNames)
        
        ### Check that SNPmat has appropriate values (0, 1, 2, or 9, only)
        snplevs <- levels(as.factor(unlist(SNPmat)))
        ls <- paste(snplevs, collapse = "")
        if (ls != "012" & ls != "0129") {
            print(error("Error: Your snp matrix does not have 0,1, and 2"))
        }
        
        ### Checking that locusNames and popNames have the same lengths as the
        #columns and rows of SNPmat
        if (dim(SNPmat)[1] != length(popname)) {
            print(error("Error: your population names do not match your SNP matrix"))
        }
        
        if (dim(SNPmat)[2] != length(locusname)) {
            print(error("Error:  your locus names do not match your SNP matrix"))
        }
        
        writeLines(report("Calculating FSTs, may take a few minutes..."))
        
        nloci <- length(locusname)
        FSTmat <- matrix(NA, nrow = nloci, ncol = 8)
        for (i in 1:nloci) {
            FSTmat[i,] <-unlist(getFSTs_diploids(popname, SNPmat[, i]))
            if (i %% 10000 == 0) {
                print(paste(i, "done of", nloci))
            }
        }
        outTemp <-as.data.frame(FSTmat)
        outTemp <-cbind(locusname, outTemp)
        
        colnames(outTemp) <-c(
            "LocusName",
            "He",
            "FST",
            "T1",
            "T2",
            "FSTNoCorr",
            "T1NoCorr",
            "T2NoCorr",
            "meanAlleleFreq"
        )
        return(outTemp)
        
    }

### Calculates FST etc. from a single locus from a column of individual data
getFSTs_diploids <- function(popNameList,
                            SNPDataColumn) {
    # eliminating the missing data for this locus
    popnames <-unlist(as.character(popNameList))
    popNameTemp <-popnames[which(SNPDataColumn != 9)]
    snpDataTemp <-SNPDataColumn[SNPDataColumn != 9]
    
    HetCounts <-
        tapply(snpDataTemp, list(popNameTemp, snpDataTemp), length)
    HetCounts[is.na(HetCounts)] <-0
    
    # Case: all individuals are genetically identical at this locus
    if (dim(HetCounts)[2] == 1) {
        return(
            list(
                He = NA,
                FST = NA,
                T1 = NA,
                T2 = NA,
                FSTNoCorr = NA,
                T1NoCorr = NA,
                T2NoCorr = NA,
                meanAlleleFreq = NA
            )
        )
    }
    
    if (dim(HetCounts)[2] == 2) {
        if (paste(colnames(HetCounts), collapse = "") == "01") {
            HetCounts <-cbind(HetCounts, `2` = 0)
        }
        if (paste(colnames(HetCounts), collapse = "") == "12") {
            HetCounts <-cbind(`0` = 0, HetCounts)
        }
        if (paste(colnames(HetCounts), collapse = "") == "02") {
            HetCounts <-cbind(HetCounts[, 1], `1` = 0, HetCounts[, 2])
        }
    }
    
    out <-WC_FST_Diploids_2Alleles(HetCounts)
    return(out)
}
