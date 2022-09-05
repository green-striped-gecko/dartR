
utils.pa.Chao <- function(x,pop1_m,pop2_m){

p1_m <- as.matrix(pop1_m)
p2_m <- as.matrix(pop2_m)
p1alf_m <- colMeans(p1_m, na.rm = T) / 2
p2alf_m <- colMeans(p2_m, na.rm = T) / 2

#identifying private alleles
pop1_pop2_m <- unique(unname(unlist(c(which(p2alf_m == 0 & p1alf_m != 0),
                                      which(p2alf_m == 1 & p1alf_m != 1)))))

pop2_pop1_m <- unique(unname(unlist(c(which(p1alf_m == 0 & p2alf_m != 0),
                                      which(p1alf_m == 1 & p2alf_m != 1)))))

# keeping only the loci with private alleles
pa_pop1_pop2_m <- gl.keep.loc(x,loc.list = locNames(x)[pop1_pop2_m],verbose =0)
pa_pop2_pop1_m <- gl.keep.loc(x,loc.list = locNames(x)[pop2_pop1_m],verbose =0)
# recalculating allele frequencies 
pa_pop1_pop2_m <- gl.recalc.metrics(pa_pop1_pop2_m,verbose =0)
pa_pop2_pop1_m <- gl.recalc.metrics(pa_pop2_pop1_m,verbose =0)

# calculating the number of singletons 
table_pop1_pop2_m <- as.data.frame(table(pa_pop1_pop2_m$other$loc.metrics$maf*(nInd(pa_pop1_pop2_m)*2)))
colnames(table_pop1_pop2_m) <- c("number_alleles","count")
table_pop1_pop2_m$number_alleles <- as.numeric(table_pop1_pop2_m$number_alleles)
table_pop1_pop2_m$number_alleles <- ceiling(table_pop1_pop2_m$number_alleles)

table_pop2_pop1_m <- as.data.frame(table(pa_pop2_pop1_m$other$loc.metrics$maf*(nInd(pa_pop2_pop1_m)*2)))
colnames(table_pop2_pop1_m) <- c("number_alleles","count")
table_pop2_pop1_m$number_alleles <- as.numeric(table_pop2_pop1_m$number_alleles )
table_pop2_pop1_m$number_alleles <- ceiling(table_pop2_pop1_m$number_alleles)

# estimating the number of un detected private alleles 
# equation 2c 
# Deciphering the enigma of undetected species, phylogenetic, and functional 
# diversity based on Good-Turing theory

f1_pop1_pop2 <- table_pop1_pop2_m[1,"count"]
f2_pop1_pop2 <- table_pop1_pop2_m[2,"count"]
n_pop1_pop2 <- sum(table_pop1_pop2_m$count)
pa_Chao_pop1_pop2 <- ((n_pop1_pop2-1)/n_pop1_pop2) * ((f1_pop1_pop2^2) / (2*f2_pop1_pop2))

f1_pop2_pop1 <- table_pop2_pop1_m[1,"count"]
f2_pop2_pop1 <- table_pop2_pop1_m[2,"count"]
n_pop2_pop1 <- sum(table_pop2_pop1_m$count)
pa_Chao_pop2_pop1 <- ((n_pop2_pop1-1)/n_pop2_pop1) * ((f1_pop2_pop1^2) / (2*f2_pop2_pop1))

return(list(pa_Chao_pop1_pop2,pa_Chao_pop2_pop1))

}