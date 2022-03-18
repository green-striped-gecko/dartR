### Calculate Allele Frequencies from Genotype Counts
allele.freq <- function(genotypeCounts) {
  n = sum(genotypeCounts) - genotypeCounts["NN"]
  p = ((2*genotypeCounts["AA"]) + genotypeCounts["Aa"])/(2*n)
  q = 1-p
  freqs = c(p,q)
  names(freqs) = c("p", "q")
  return(freqs)
}
### Calculate the LD parameter R-squared from 2 rows of a VCF file
calc_r2 <- function(row1, row2) {
  g1 = get.field(row1[10:length(row1)], row1[9], "GT")
  g2 = get.field(row2[10:length(row2)], row2[9], "GT")
  pA = unname(allele.freq(count.genotypes(g1))["p"])
  pB = unname(allele.freq(count.genotypes(g2))["p"])
  h = get.haplotypes(g1, g2)
  pAB = (length(h[h=="00"]))/(length(h))
  D = pAB - (pA * pB)
  rsq = (D**2)/(pA*(1-pA)*pB*(1-pB))
  return(rsq)
}
### Count the AA, Aa, and aa genotypes in a sample
count.genotypes <- function(genotypes) {
  genotypes = gsub("(\\||/)", "", genotypes) 
  gen.patterns = c("00", "01", "10", "11", "..") 
  my.counts=table(factor(genotypes, levels=gen.patterns)) 
  final.counts = c(my.counts[1], (my.counts[2]+my.counts[3]), my.counts[4:5]) 
  names(final.counts) = c("AA", "Aa", "aa", "NN") 
  return(final.counts)
}
### Count the number of derived alleles for a VCF row
derivedCount <- function(row) {
  row=as.vector(row, mode="character")
  x=count.genotypes(get.field(row[10:length(row)], row[9], "GT"))
  dc=(2*x["aa"])+x["Aa"]
  return(unname(dc))
}
### Calculate Nucleotide Divergence (Dxy)
dxy <- function(vcf1, vcf2, perBP=TRUE) {
  g1=t(apply(vcf1, 1, function(x) get.field(x[10:length(x)], x[9], "GT")))
  g2=t(apply(vcf2, 1, function(x) get.field(x[10:length(x)], x[9], "GT")))
  af1=apply(g1, 1, function(x) allele.freq(count.genotypes(x))["p"])
  af2=apply(g2, 1, function(x) allele.freq(count.genotypes(x))["p"])
  ## Let x be the allele frequency (p) in pop1 * (1-p) in Pop2
  x = af1 * (1-af2)
  ## Let y be the allele frequency (p) in pop2 * (1-p) in Pop1
  y = af2 * (1-af1)
  dxy=sum((x+y))
  if (perBP) {
    c = unique(vcf1$CHROM)
    s = sapply(c, function(x,y) min(y[which(y$CHROM==x),2]), y=vcf1)
    e = sapply(c, function(x,y) max(y[which(y$CHROM==x),2]), y=vcf1)
    bp=sum(e-s)
    return(dxy/bp)
  } else { 
    return(dxy) 
  }
}
### Calculate Expected Het. and return, p, n, and H
expected.het <- function(genotypes) {
  obs.counts = count.genotypes(genotypes)
  n = sum(obs.counts) - obs.counts["NN"]
  freqs = allele.freq(obs.counts)
  Hexp = 2 * freqs[1] * freqs[2]
  res = c(freqs["p"], n, Hexp)
  res=as.numeric(unname(res))
  return(res)
}
### Determine which SNPs are Polymorphic vs Fixed in 2 species
fixed.poly <- function(vcf1, vcf2) {
  g1=t(apply(vcf1, 1, function(x) get.field(x[10:length(x)], x[9], "GT")))
  g2=t(apply(vcf2, 1, function(x) get.field(x[10:length(x)], x[9], "GT")))
  af1=apply(g1, 1, function(x) allele.freq(count.genotypes(x))["p"])
  af2=apply(g2, 1, function(x) allele.freq(count.genotypes(x))["p"])
  res = rep("Polymorphic", nrow(g1))
  res[which(abs(af1-af2)==1)] = "Fixed"
  return(res)
}
### Get a Specified Field From a VCF Sample/Genotype String
get.field <- function(samples, format, fieldName) {
  x=strsplit(samples, split=":")
  fields=unlist(strsplit(format, split=":")) 
  i=which(fields==fieldName)
  if (!(fieldName %in% fields)) stop('fieldName not found in format fields') 
  return(sapply(x, `[[`, i)) 
}
### Get the HAPLOTYPES for a pair of genotype strings
get.haplotypes <- function(genotypes1, genotypes2) {
  a1 = gsub("\\|", "", genotypes1) 
  a2 = gsub("\\|", "", genotypes2)
  a1=unlist(strsplit(paste0(a1, collapse=""), split="")) 
  a2=unlist(strsplit(paste0(a2, collapse=""), split=""))
  haps = paste0(a1,a2)
  return(haps)
}
### Calculate minor allele frequency
maf <- function(vcf.row) {
  temp=as.vector(vcf.row, mode="character")
  af=allele.freq(count.genotypes(get.field(temp[10:length(temp)], temp[9], "GT")))
  maf=min(unname(af))
  return(maf)
}
### Calculate Nucleotide Diversity (Pi)
pi.diversity <- function(vcf, perBP=TRUE) {
  J=apply(vcf, 1, derivedCount)
  N=apply(vcf, 1, function(x) sum(count.genotypes(get.field(x[10:length(x)], x[9], "GT"))))
  C=2*N
  pi = sum((2*J*(C-J))/(C*(C-1)))
  if (perBP) {
    c = unique(vcf$CHROM)
    s = sapply(c, function(x,y) min(y[which(y$CHROM==x),2]), y=vcf)
    e = sapply(c, function(x,y) max(y[which(y$CHROM==x),2]), y=vcf)
    bp=sum(e-s)
    return(pi/bp)
  } else { return(pi) }
}	      
### Read in a VCF file as a table
read.vcf <- function(file, special.char="##", ...) {
  my.search.term=paste0(special.char, ".*")
  all.lines=readLines(file)
  clean.lines=gsub(my.search.term, "",  all.lines)
  clean.lines=gsub("#CHROM", "CHROM", clean.lines)
  read.table(..., text=paste(clean.lines, collapse="\n"))
}
### Calculate the variance for Tajima's D
variance.d <- function(n,S) {
  a1=sum(1/(seq(from=1, to=(n-1), by=1)))
  a2=sum(1/((seq(from=1, to=(n-1), by=1))**2))
  b1=(n+1)/(3*(n-1))
  b2=(2*((n**2)+n+3))/((9*n)*(n-1))
  c1=b1 - (1/a1)
  c2=b2-((n+2)/(a1*n)) + (a2/(a1**2))
  e1=c1/a1
  e2=c2/((a1**2)+a2)
  var=(e1*S) + (e2*S*(S-1))
  return(var)
}
### Calculate Waterson's theta
waterson.theta <- function(data, perBP=TRUE) {
  num.bp=nrow(data)
  check=apply(data, 1, FUN=maf)
  filter=data[which(check>0),]
  Sn=nrow(filter)
  n.ind=ncol(filter)-9
  i=seq(1, ((2*n.ind)-1))
  theta=Sn/sum(1/i)
  if (perBP) { return(theta/num.bp) }
  else { return(theta) }
}
