library(dartR)
library(data.table)
library(ggplot2)



gl.hwe.diagnostics <- function(x, alpha_val=0.05, brk = seq(0, 1, 1 / 20)) {
  
  #### Function helpers ####
  hwe.dist.plot <- function(hweout, brk) {
    
    ggplot(hweout, aes(Prob)) + geom_histogram(breaks=brk) +
      geom_hline(yintercept = nrow(hweout)/(length(brk)-1), col="red")
  }
  
  # abandoned this for now
  # FstvsFis.plot <- function(Fstats) {
  #   ggplot(Fstats, aes(Fst, Fis)) + geom_point() + geom_smooth(method = "lm") #+ facet_grid(~.id)
  # }
  # 
  nLoBynPop_plot <- function(nTimes) {
    ggplot(nTimes, aes(nPop, Freq, col=Data, fill=Data)) + 
      geom_histogram(stat="identity", position = "dodge") + 
      scale_y_log10()
  }
  
  
  # Distribution of pvlues by qual bins
  suppressWarnings(
  hweout <- gl.report.hwe(x, sig_only = F, verbose = 0)
  )
  print(hwe.dist.plot(hweout, brk))
  
  # Fst vs Fis scatter plot with linear regression
  # lpops <- seppop(x)
  # lFstats <- lapply(lpops, gl.basic.stats, verbose = 0)
  # lFstats <- lapply(lFstats, "[[", "perloc")
  # Fstats <- rbindlist(l = lFstats, use.names = TRUE, idcol = TRUE)
  Fstats <- gl.basic.stats(x, verbose = 0)
  # FstatsLoc <- Fstats$perloc
  # print(FstvsFis.plot(FstatsLoc))
  
  # Number of loci out of HWE as a function of a population
  hweout.dt <- data.table(hweout)
  nTimesBypop <- hweout.dt[, .N, by=c("Locus", "Sig")]
  setkey(nTimesBypop, Sig)

  nTimesBypop.df <- as.data.frame(table(nTimesBypop["sig", N]))
  
  # Include the non-sig tests
  nTimesBypop.df <- rbind(data.frame(Var1=0, Freq=nTimesBypop["no_sig", sum(N)]), nTimesBypop.df)
  nTimesBypop.df$Data <- "Observed"

  # Generate the null distribution
  nullDist <- as.data.frame(table(rbinom(length(hweout.dt[, unique(Locus)]), 
                                         size = length(hweout.dt[, unique(Population)]), 
                                         prob = alpha_val)))
  nullDist$Data <- "Null expectation"
  
  # Compile the data for the plot
  nTimesBypop.fin <- rbind(nTimesBypop.df, nullDist)
  names(nTimesBypop.fin)[1] <- "nPop"
  print(nLoBynPop_plot(nTimesBypop.fin))
  
  # Collate HWE tests and Fis per locus and pop
  FisPops <- data.table(Fstats$Fis, keep.rownames = TRUE)
  
  # fix the ehadings when there is only one pop
  if(length(levels(pop(x))) == 1) {
    FisPops[, dumpop := NULL]
    setnames(FisPops, "1", levels(pop(x)))
  }
  setnames(FisPops, "rn", "Locus")
  FisPopsLong <- melt(FisPops, id.vars = "Locus", variable.name = "Population", value.name = "Fis")
  FisPopsLong[, Locus := sub("^X", "", Locus)]
  hweout.dt[, Locus := gsub("-|/", replacement = ".", x = Locus)]
  hwe_Fis <- merge(hweout.dt, FisPopsLong, by = c("Locus", "Population"))
  hwe_Fis[, Deficiency := Fis > 0]
  hwe_Fis[,  Excess:= Fis < 0]
  setkey(hwe_Fis, Sig)

  hwe_summary <- hwe_Fis["sig", .(nSig=.N, nExpected = alpha_val * nLoc(x),
                                  Deficiency=sum(Deficiency, na.rm = TRUE), 
                                  Excess=sum(Excess, na.rm = TRUE),
                                  PropDeficiency=sum(Deficiency, na.rm = TRUE)/.N),
                         by=Population]
  
  chsq <- hwe_Fis[, .(ChiSquare=-2*(sum(log(Prob)))), by= Population]
  chsq[, pvalue:=pchisq(ChiSquare, 2*nLoc(x), lower.tail = FALSE)]
  hwe_summary <- merge(hwe_summary, chsq, by = "Population")
  return(hwe_summary)
}

