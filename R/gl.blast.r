#'@name gl.blast
#'
#'@title Aligns nucleotides sequences against those present in a target database 
#'using blastn
#'
#'@description Basic Local Alignment Search Tool (BLAST; Altschul et al., 1990 &
#'  1997) is a sequence comparison algorithm optimized for speed used to search
#'  sequence databases for optimal local alignments to a query. This function
#'  creates fasta files, creates databases to run BLAST, runs blastn and filters
#'  these results to obtain the best hit per sequence.
#'
#'  This function can be used to run BLAST alignment of short-read (DArTseq
#'  data) and long-read sequences (Illumina, PacBio... etc). You can use
#'  reference genomes from NCBI, genomes from your private collection, contigs,
#'  scaffolds or any other genetic sequence that you would like to use as
#'  reference.
#'
#'@param x Either a genlight object containing a column named
#'  'TrimmedSequence' containing the sequence of the SNPs (the sequence tag)
#'  trimmed of adapters as provided by DArT; or a path to a fasta file with the
#'  query sequences [required].
#'@param ref_genome Path to a reference genome in fasta of fna format
#'  [required].
#'@param task Four different tasks are supported: 1) “megablast”, for very
#'  similar sequences (e.g, sequencing errors), 2) “dc-megablast”, typically
#'  used for inter-species comparisons, 3) “blastn”, the traditional program
#'  used for inter-species comparisons, 4) “blastn-short”, optimized for
#'  sequences less than 30 nucleotides [default 'megablast'].
#'@param Percentage_identity Not a very sensitive or reliable measure of
#'  sequence similarity, however it is a reasonable proxy for evolutionary
#'  distance. The evolutionary distance associated with a 10 percent change in
#'  Percentage_identity is much greater at longer distances. Thus, a change from
#'  80 – 70 percent identity might reflect divergence 200 million years earlier
#'  in time, but the change from 30 percent to 20 percent might correspond to a
#'  billion year divergence time change [default 70].
#'@param Percentage_overlap Calculated as alignment length divided by the
#'  query length or subject length (whichever is shortest of the two lengths,
#'  i.e.  length / min(qlen,slen) ) [default 0.8].
#'@param bitscore A rule-of-thumb for inferring homology, a bit score of 50
#'  is almost always significant [default 50].
#'@param number_of_threads Number of threads (CPUs) to use in blastn search
#'  [default 2].
#'@param verbose verbose= 0, silent or fatal errors; 1, begin and end; 2,
#'progress log ; 3, progress and results summary; 5, full report
#'[default 2 or as specified using gl.set.verbosity]
#'
#'@details \strong{Installing BLAST}
#'
#'  You can download the BLAST installs from:
#'  \url{https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/}
#'
#'  It is important to install BLAST in a path that does not contain spaces for
#'  this function to work.
#'
#'  \strong{Running BLAST}
#'
#'  Four different tasks are supported: 
#'  \itemize{ 
#'  \item “megablast”, for very
#'  similar sequences (e.g, sequencing errors) 
#'  \item “dc-megablast”, typically
#'  used for inter-species comparisons 
#'  \item “blastn”, the traditional program
#'  used for inter-species comparisons
#'  \item “blastn-short”, optimized for
#'  sequences less than 30 nucleotides }
#'
#'  If  you  are  running  a  BLAST alignment of  similar  sequences,  for
#'  example  Turtle  Genome  Vs Turtle Sequences, the recommended parameters
#'  are: task = “megablast”, Percentage_identity = 70, Percentage_overlap =  0.8
#'  and bitscore = 50.
#'
#'  If you are running a BLAST alignment of highly dissimilar sequences because
#'  you are probably looking for sex linked  hits in  a distantly  related
#'  species,  and  you  are aligning for example sequences of Chicken Genome Vs
#'  Bassiana, the recommended parameters are: task = “dc-megablast”,
#'  Percentage_identity = 50, Percentage_overlap =  0.01 and bitscore = 30.
#'
#'  Be aware that running BLAST might take a long time (i.e. days) depending of
#'  the size of your query, the size of your database and the number of threads
#'  selected for your computer.
#'
#'  \strong{BLAST output}
#'
#'  The BLAST output is formatted as a table using output format 6, with columns
#'  defined in the following order: \itemize{ \item qseqid - Query Seq-id \item
#'  sacc - Subject accession \item stitle - Subject Title \item qseq - Aligned
#'  part of query sequence \item sseq - Aligned part of subject sequence \item
#'  nident - Number of identical matches \item mismatch - Number of mismatches
#'  \item pident - Percentage of identical matches \item length - Alignment
#'  length \item evalue - Expect value \item bitscore - Bit score \item qstart -
#'  Start of alignment in query \item qend - End of alignment in query \item
#'  sstart - Start of alignment in subject \item send - End of alignment in
#'  subject \item gapopen - Number of gap openings \item gaps - Total number of
#'  gaps \item qlen - Query sequence length \item slen - Subject sequence length
#'  \item PercentageOverlap - length / min(qlen,slen) }
#'
#'  Databases containing unfiltered aligned sequences, filtered aligned
#'  sequences and one hit per sequence are saved to the temporal directory
#'  (tempdir) and can be accessed with the function
#'  \code{\link{gl.print.reports}} and listed with the function
#'  \code{\link{gl.list.reports}}. Note that they can be accessed only in the
#'  current R session because tempdir is cleared each time that the R session is
#'  closed.
#'
#'  \strong{BLAST filtering}
#'
#'  BLAST output is filtered by ordering the hits of each sequence first by the
#'  highest percentage identity, then the highest percentage overlap and then
#'  the highest bitscore. Only one hit per sequence is kept based on these
#'  selection criteria.
#'
#'@return If the input is a genlight object: returns a genlight object with one
#'  hit per sequence merged to the slot $other$loc.metrics. If the input is a
#'  fasta file: returns a dataframe with one hit per sequence.
#'
#'@author Berenice Talamantes Becerra & Luis Mijangos (Post to
#'  \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' \dontrun{
#' res <- gl.blast(x= testset.gl,ref_genome = 'sequence.fasta')
#' # display of reports saved in the temporal directory
#' gl.list.reports()
#' # open the reports saved in the temporal directory
#' blast_databases <- gl.print.reports(1)
#' }
#'
#'@references
#'\itemize{
#'\item Altschul, S. F., Gish, W., Miller, W., Myers, E. W., & Lipman, D.
#'  J. (1990). Basic local alignment search tool. Journal of molecular biology,
#'  215(3), 403-410.
#'\item Altschul, S. F., Madden, T. L., Schäffer, A. A., Zhang, J., Zhang,
#'  Z., Miller, W., & Lipman, D. J. (1997). Gapped BLAST and PSI-BLAST: a new
#'  generation of protein database search programs. Nucleic acids research,
#'  25(17), 3389-3402.
#'\item Pearson, W. R. (2013). An introduction to sequence similarity
#'  (“homology”) searching. Current protocols in bioinformatics, 42(1), 3-1.
#'  }
#'
#'@seealso \code{\link{gl.print.history}}
#'
#'@family reference genomes
#'
#'@export
#'

gl.blast <- function(x,
                     ref_genome,
                     task = "megablast",
                     Percentage_identity = 70,
                     Percentage_overlap = 0.8,
                     bitscore = 50,
                     number_of_threads = 2,
                     verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbose = verbose)
    
    # Check if the x@other$loc.metrics$TrimmedSequence slot exists
    if (class(x)[1] == "genlight") {
        if (is.null(x$other$loc.metrics$TrimmedSequence)) {
            stop(error(
                "\n\nFatal Error: TrimmedSequence column is required!.\n\n"
            ))
        }
    }
    
    # DO THE JOB
    
    # getting the query fasta files
    if (class(x)[1] == "genlight" | class(x)[1] == "dartR") {
        fasta.input <-
            c(rbind(
                paste("> ", 1:nLoc(x)),
                as.character(x$other$loc.metrics$TrimmedSequence)
            ))
        write.table(
            fasta.input,
            file = paste0(tempdir(), "/fasta.input"),
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE
        )
    } else {
        file.copy(
            from = x,
            to = paste0(tempdir(), "/fasta.input"),
            overwrite = TRUE
        )
    }
    
    # Find executable makeblastdb if unix
    if (grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
        path_makeblastdb <- Sys.which("makeblastdb")
    }
    ## if windows
    if (!grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
        path_makeblastdb <-
            tryCatch(
                system(sprintf("where %s", "makeblastdb"), intern = TRUE)[1],
                warning = function(w)
                    "",
                error = function(e)
                    ""
            )
        if (grepl("\\s", path_makeblastdb)) {
        stop(
            error(
              "  The path to the executable for makeblastdb has spaces. 
              Please move it\n to a path without spaces so BLAST can work.\n\n"
                )
            )
        }
    }
    
    ## if not found
    if (all(path_makeblastdb == "")) {
        stop(
            error(
                "  Executable for makeblastdb not found! Please make sure that 
                the software\n is correctly installed.\n\n"
            )
        )
    }
    
    # creating BLAST databases
    system(paste(
        path_makeblastdb,
        "-in",
        ref_genome,
        "-dbtype",
        "nucl",
        "-out",
        paste0(tempdir(), "/db_blast")
    ))
    
    # Find executable blastn if unix
    if (grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
        path_blastn <- Sys.which("blastn")
    }
    ## if windows
    if (!grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
        path_blastn <-
            tryCatch(
                system(sprintf("where %s", "blastn"), intern = TRUE)[1],
                warning = function(w)
                    "",
                error = function(e)
                    ""
            )
        if (grepl("\\s", path_blastn)) {
            stop(
                error(
                    "The path to the executable for blastn has spaces. Please 
                    move it to a\n path without spaces so BLAST can work.\n\n"
                )
            )
        }
    }
    ## if not found
    if (all(path_blastn == "")) {
        stop(
            error(
                "Executable for blastn not found! Please make sure that the 
                software is\n correctly installed.\n\n"
            )
        )
    }
    
    if (verbose >= 1) {
        cat(report("Starting BLASTing\n\n"))
    }
    
    # BLASTing
    system(
        paste(
            path_blastn,
            " -task ",
            task,
            " -db ",
            paste0(tempdir(), "/db_blast"),
            " -query ",
            paste0(tempdir(), "/fasta.input"),
            " -out ",
            paste0(tempdir(), "/output_blast.txt"),
            " -perc_identity ",
            Percentage_identity,
            " -num_threads ",
            number_of_threads,
            " -outfmt ",
            "\"",
            6,
            " qseqid sacc stitle qseq sseq nident mismatch pident length evalue 
            bitscore qstart qend sstart send gapopen gaps qlen slen\""
        )
    )
    
    file_size <- file.info(paste0(tempdir(), "/output_blast.txt"))$size
    
    # reading file for filtering
    if (file_size > 0) {
        blast_res_unfiltered <-
            read.table(
                file = paste0(tempdir(), "/output_blast.txt"),
                header = FALSE,
                sep = "\t",
                quote = "\"",
                dec = ".",
                fill = TRUE,
                comment.char = "",
                stringsAsFactors = FALSE
            )
    } else {
        cat(report("No sequences were aligned\n\n"))
        return(x)
    }
    
    if (verbose >= 1) {
        cat(report("Starting filtering\n\n"))
    }
    
    colnames(blast_res_unfiltered) <-
        c(
            "qseqid",
            "sacc",
            "stitle",
            "qseq",
            "sseq",
            "nident",
            "mismatch",
            "pident",
            "length",
            "evalue",
            "bitscore",
            "qstart",
            "qend",
            "sstart",
            "send",
            "gapopen",
            "gaps",
            "qlen",
            "slen"
        )
    
    # calculate percentage overlap ratio of the alignment length divided by the
    # query length or subject length (whichever is shortest of
    # the two lengths)
    min_length_BF <-
        apply(blast_res_unfiltered[, c("qlen", "slen")], 1, min)
    blast_res_unfiltered$PercentageOverlap <-
        blast_res_unfiltered$length / min_length_BF
    # filtering first by percentage overlap and bitscore
    blast_res_filtered <-
        blast_res_unfiltered[which(
            blast_res_unfiltered$PercentageOverlap >= Percentage_overlap &
                blast_res_unfiltered$bitscore >=
                bitscore
        ),]
    # splitting hits by sequence
    all_hits <-
        split(x = blast_res_filtered, f = blast_res_filtered$qseqid)
    # ordering by first considering the highest percentage identity, then the 
    # highest percentage overlap, then the highest bitscore. Only
    # one hit per sequence is kept based on these selection criteria.
    one_hit_temp <- lapply(all_hits, function(x) {
        x[order(x$pident,
                x$PercentageOverlap,
                x$bitscore,
                decreasing = T),][1,]
    })
    
    one_hit <- plyr::rbind.fill(one_hit_temp)
    # merging one hit per sequence with genlight object
    if (class(x)[1] == "genlight") {
        one_hit_temp <- x$other$loc.metrics
        one_hit_temp$qseqid <- 1:nLoc(x)
        if (!is.null(one_hit)) {
            x$other$loc.metrics <-
                merge(one_hit_temp,
                      one_hit,
                      by = "qseqid",
                      all = T)
        }
    }
    
    if (verbose >= 1) {
        cat(report(paste(
            nrow(one_hit), " sequences were aligned after filtering"
        )))
    }
    
    match_call <-
        paste0(names(match.call()),
               "_",
               as.character(match.call()),
               collapse = "_")
    
    # creating file names
    temp_blast_unfiltered <- tempfile(pattern = "Blast_unfiltered_")
    temp_blast_filtered <- tempfile(pattern = "Blast_filtered_")
    temp_one_hit <- tempfile(pattern = "Blast_one_hit_")
    
    # saving to tempdir
    saveRDS(list(match_call, blast_res_unfiltered), 
            file = temp_blast_unfiltered)
    saveRDS(list(match_call, blast_res_filtered), file = temp_blast_filtered)
    saveRDS(list(match_call, one_hit), file = temp_one_hit)
    
    if (verbose >= 2) {
        cat(
            report(
                "  NOTE: Retrieve output files from tempdir using 
                gl.list.reports() and gl.print.reports()\n"
            )
        )
    }
    
    # ADD TO HISTORY
    if (class(x)[1] == "genlight") {
        nh <- length(x@other$history)
        x@other$history[[nh + 1]] <- match.call()
    }
    
    # FLAG SCRIPT END
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n\n"))
    }
    
    if (class(x)[1] == "genlight") {
        return(x)
    } else {
        return(one_hit)
    }
    
}
