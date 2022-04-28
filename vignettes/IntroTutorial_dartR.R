## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("dartR")

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("devtools")
#  library(devtools)
#  install.packages("BiocManager")
#  BiocManager::install(c("SNPRelate", "qvalue"))
#  install_github("green-striped-gecko/dartR")
#  library(dartR)
#  

## -----------------------------------------------------------------------------
 library(dartR)

## ---- eval=FALSE--------------------------------------------------------------
#  BiocManager::install(c("SNPRelate", "qvalue"))

## ---- eval=FALSE--------------------------------------------------------------
#  # Set the default working directory (change this to suit)
#  setwd("c:/your.working.directory/")

## -----------------------------------------------------------------------------
# Rename the test genlight object to gl, something simple
gl <- testset.gl

## -----------------------------------------------------------------------------
m <- as.matrix(gl)

## -----------------------------------------------------------------------------
as.matrix(gl)[1:5,1:3]

## ---- eval=FALSE--------------------------------------------------------------
#  gl <- gl.read.dart(filename = "testset.csv", ind.metafile = " ind_metrics.csv")

## -----------------------------------------------------------------------------
dartfile <- system.file("extdata","testset_SNPs_2Row.csv", package="dartR")

metafile <- system.file("extdata","testset_metadata.csv", package="dartR")
gl <- gl.read.dart(filename=dartfile, ind.metafile = metafile, probar=FALSE)

## -----------------------------------------------------------------------------
gl

## -----------------------------------------------------------------------------
as.matrix(gl)[1:3,1:3]

## -----------------------------------------------------------------------------
nLoc(gl)
nInd(gl)
nPop(gl)

## -----------------------------------------------------------------------------
levels(pop(gl))[1:5]

## ---- eval=FALSE--------------------------------------------------------------
#  
#  gl <- gl.read.dart(filename="mydata.csv", ind.metafile = "my.metadata.csv")

## -----------------------------------------------------------------------------
#Only the entries for the first ten individuals are shown
gl@other$loc.metrics$RepAvg[1:10]

## -----------------------------------------------------------------------------
names(gl@other$loc.metrics)

## -----------------------------------------------------------------------------
#only first 10 entries showns
gl@other$ind.metrics$sex[1:10]

## -----------------------------------------------------------------------------
read.csv( paste(.libPaths()[1],"/dartR/extdata/platy.csv",sep="" ))

## -----------------------------------------------------------------------------
#you might need to install PopGenReport via
#install.packages("PopGenReport")
library(PopGenReport) 
platy <- read.genetable( paste(.libPaths()[1],"/dartR/extdata/platy.csv",
sep="" ), ind=1, pop=2, lat=3, long=4, other.min=5, other.max=6, oneColPerAll=FALSE,sep="/")

platy


## ---- eval=FALSE--------------------------------------------------------------
#  platy.gl <- (gi2gl(platy))

## ---- echo=FALSE--------------------------------------------------------------
platy.gl <- (gi2gl(platy, parallel = FALSE))

## -----------------------------------------------------------------------------
platy.gl@other$ind.metrics <- platy.gl@other$data

## -----------------------------------------------------------------------------
ts <- sapply(1:nLoc(platy.gl), function(x) paste(sample(c("A","T","G","C"), 50, replace = T),
                                                 collapse = ""))
df.loc <- data.frame(RepAvg = runif(nLoc(platy.gl)),  TrimmedSequence=ts)

platy.gl@other$loc.metrics <- df.loc



## -----------------------------------------------------------------------------
gl.report.callrate(testset.gl)
gl2 <- gl.filter.reproducibility(platy.gl, t=0.5)

## -----------------------------------------------------------------------------
platy.gl@position <- as.integer(runif(nLoc(platy.gl),2,49))
platy.gl@loc.all <- testset.gl@loc.all[1:6]

## ---- eval=FALSE--------------------------------------------------------------
#  gl2fasta(platy.gl)

## -----------------------------------------------------------------------------
gl2 <- gl.filter.callrate(gl, method = "loc", threshold = 0.95)

## -----------------------------------------------------------------------------
gl2 <- gl.filter.callrate(gl, method="ind", threshold = 0.90)


## -----------------------------------------------------------------------------
gl2 <- gl.filter.reproducibility(gl, t=1)

## ---- eval=F, echo=F----------------------------------------------------------
#  gl2 <- gl.filter.secondary(gl)
#  

## -----------------------------------------------------------------------------
gl2 <- gl.filter.monomorphs(gl, v=0)

## -----------------------------------------------------------------------------
gl2 <- gl.filter.hamming(testset.gl, threshold = 0.25, pb = F)

## ---- eval=F------------------------------------------------------------------
#  gl2 <- gl.filter.callrate(gl, method = "loc", threshold = 0.95)
#  gl3 <-  gl.filter.callrate(gl2, method="ind", threshold = 0.90)
#  gl4 <- gl.filter.repavg(gl3, t=1)

## -----------------------------------------------------------------------------
#population names (#30 populations)
levels(pop(gl))
#table on individuals per population
table(pop(gl))


## ---- fig.height=5------------------------------------------------------------
barplot(table(pop(gl)), las=2)

## ---- eval=FALSE--------------------------------------------------------------
#  glnew3 <- gl.keep.pop(gl, pop.list=c("EmsubRopeMata","EmvicVictJasp"))

## ---- eval=FALSE--------------------------------------------------------------
#  glnew3 <- gl.drop.pop(gl, pop.list=c("EmsubRopeMata","EmvicVictJasp"))

## ---- eval=FALSE--------------------------------------------------------------
#  glnew3 <- gl.merge.pop(gl, old=c("EmsubRopeMata","EmvicVictJasp"), new="outgroup")

## ---- eval=FALSE--------------------------------------------------------------
#  glnew3 <- gl.merge.pop(gl, old="EmsubRopeMata", new="Emydura_victoriae")

## -----------------------------------------------------------------------------
#individual names
indNames(gl)


## ---- eval=FALSE--------------------------------------------------------------
#  glnew3 <- gl.keep.ind(gl, ind.list=c("AA019073","AA004859"))

## ---- eval=FALSE--------------------------------------------------------------
#  glnew3 <- gl.drop.pop(gl, ind.list=c("AA019073","AA004859"))

## ---- eval=FALSE--------------------------------------------------------------
#  gl.make.recode.pop(gl, outfile = "new_pop_assignments.csv")
#  

## ---- eval=FALSE--------------------------------------------------------------
#  
#  glnew <- gl.recode.pop(gl, pop.recode="new_pop_assignments.csv")
#  

## -----------------------------------------------------------------------------
levels(pop(gl))

## ---- eval=FALSE--------------------------------------------------------------
#  
#  glnew2 <- gl.edit.recode.pop(gl)
#  

## -----------------------------------------------------------------------------
#only first 10 entries are shown
indNames(gl)[1:10]


## ---- eval=FALSE--------------------------------------------------------------
#  gl.make.recode.ind(gl, outfile="new_ind_assignments.csv")

## ---- eval=FALSE--------------------------------------------------------------
#  glnew3 <- gl.recode.ind(gl, ind.recode="new_ind_assignments.csv")
#  

## ---- eval=F------------------------------------------------------------------
#  gl <- gl.edit.recode.ind(gl, ind.recode="new_ind_assignments.csv")

## -----------------------------------------------------------------------------
gl_new <- gl[gl$pop!="EmmacBrisWive", ]

## -----------------------------------------------------------------------------
glsub <- gl[1:7, 1:3]
glsub

## -----------------------------------------------------------------------------
dim(glsub@other$ind.metrics)

dim(glsub@other$loc.metrics)

## -----------------------------------------------------------------------------

index.ind <- pop(gl)=="EmmacRussEube" | pop(gl)=="EmvicVictJasp"
#check if the index worked
table( pop(gl), index.ind)


index.loc <- sample(nLoc(gl), 30, replace = F)
index.loc


## -----------------------------------------------------------------------------
glsub2 <- gl[index.ind, index.loc]
glsub2@other$ind.metrics <- gl@other$ind.metrics[index.ind,] #not necessary
glsub2@other$loc.metrics <- gl@other$loc.metrics[index.loc,] #necessary

## -----------------------------------------------------------------------------
glsub2
dim(glsub2@other$ind.metrics)
dim(glsub2@other$loc.metrics)

## -----------------------------------------------------------------------------
gg <- gl[1:7,1:100]
#need to subset loc.metrics as well!!!
gg@other$loc.metrics <- gg@other$loc.metrics[1:100,]

gl.dist.pop(gg, method="euclidean")
gl.dist.pop(gg, method="reynolds")

## -----------------------------------------------------------------------------
glind7 <- gl[1:7,]  #copy and store the original dataset in glind
pop(glind7) <- indNames(glind7)# redefine the population information

gl.dist.pop(glind7[1:7,], method="euclidean")
data.frame(ind=1:7, indNames=indNames(gl)[1:7], pop=pop(gl)[1:7])

## ---- eval=FALSE--------------------------------------------------------------
#  library(StAMPP) #you may need to install the package
#  pwfst <-stamppFst(gl[1:20,], nboots=1, percent=95, nclusters=1)
#  round(pwfst,3)

## ---- eval=FALSE--------------------------------------------------------------
#  pwGst <-stamppNeisD(gl[1:20,]) #no parallel version :-()
#  round(pwGst,3)

## ---- warning=FALSE-----------------------------------------------------------
library(mmod) #you may need to install the package first
#for performance reason use only a subset (and recode the populations)
recpops<- factor(rep(LETTERS[1:5],50))
glsub <- testset.gl
pop(glsub)<-recpops

gi <- gl2gi(glsub, v=0) #v=0 suppresses output
round(pairwise_D(gi),4)
round(pairwise_Gst_Hedrick(gi),4)
round(pairwise_Gst_Nei(gi),4)

## -----------------------------------------------------------------------------
pc <- gl.pcoa(gl, nfactors=5)

## -----------------------------------------------------------------------------
names(pc)

## -----------------------------------------------------------------------------
barplot(pc$eig/sum(pc$eig)*100, )

## ---- fig.height=5------------------------------------------------------------
gl.pcoa.plot(pc, gl, pop.labels="pop", xaxis=1, yaxis=2)


## ----eval=FALSE---------------------------------------------------------------
#  glnew <- gl.edit.recode.pop(gl)

## -----------------------------------------------------------------------------
glnew <- testset.gl
levels(pop(glnew)) <- c(rep("Coast",5),rep("Cooper",3),rep("Coast",5),
rep("MDB",8),rep("Coast",7),"Em.subglobosa","Em.victoriae")
gl.pcoa.plot(pc, glnew, pop.labels="pop", xaxis=1, yaxis=2)

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("devtools")
#  library(devtools)
#  install_github("hadley/ggplot2")
#  library(ggplot2)

## ---- eval=FALSE--------------------------------------------------------------
#  gl.pcoa.plot(pc, glnew, interactive = TRUE, xaxis=1, yaxis=2)

## ---- eval=FALSE--------------------------------------------------------------
#  gl.pcoa.plot.3d(pc, glnew)
#  

## -----------------------------------------------------------------------------
gl.tree.nj(glnew, type="fan")


## ---- eval=FALSE--------------------------------------------------------------
#  fd <- gl.collapse.recursive(gl, t=0)

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  gl <- testset.gl
#  fd <- gl.collapse.recursive(gl, t=0)
#  

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  gl <- testset.gl
#  fd <- gl.collapse.recursive(gl, test=TRUE, delta=0.02, reps=1000, t=0, v=3)
#  fd.sig <- gl.collapse.pval(fd, prefix="fd_sig", delta=0.02, alpha=0.05, v=3)
#  

## ---- eval=FALSE--------------------------------------------------------------
#  gl <- fd.sig$gl
#  phy <- gl2phylip(gl, outfile="turtle.phy", bstrap=1000)

## ---- fig.height=4------------------------------------------------------------
gl <- gl.ibd(x=testset.gl[1:180,])

## ---- eval=FALSE--------------------------------------------------------------
#  gl <- testset.gl
#  phy <- gl2phylip(gl, outfile="turtle.phy", bstrap=1000)

## ---- eval=FALSE--------------------------------------------------------------
#  gl2fasta(gl, method=2, outfile="nohets.fasta")

## -----------------------------------------------------------------------------
gl.report.bases(testset.gl)

## ---- eval=FALSE--------------------------------------------------------------
#  gl2fasta(gl, method=4, outfile="nohets.fasta")
#  

## ---- eval=FALSE--------------------------------------------------------------
#  gl2fasta(gl, method=1, outfile="ambcodes.fasta")
#  

## ---- eval=FALSE--------------------------------------------------------------
#  gl2fasta(gl, method=3, outfile="ambcodes.fasta")
#  

## -----------------------------------------------------------------------------
gl.report.bases(testset.gl)

## -----------------------------------------------------------------------------
x <- gl.assign.pa(testset.gl, unknown = "UC_00146", nmin=10, t=0)

## ---- fig.height=4------------------------------------------------------------
#x <- gl.assign.pca(x, unknown ="UC_00146", plevel=0.95)

## -----------------------------------------------------------------------------
gl <- testset.gl
gi <- gl2gi(gl, v=0)


## ---- eval=FALSE--------------------------------------------------------------
#  gl2 <- gi2gl(gi)

## ---- eval=F------------------------------------------------------------------
#  glnew <- gl.nhybrids(testset.gl)

## ---- eval=FALSE--------------------------------------------------------------
#  glnew <- gl2phylip(outfile = "phyinput.txt")
#  

## ---- eval=FALSE--------------------------------------------------------------
#  gl.new <- gl2phylip(outfile = "phyinput.txt", bstrap = 1000)
#  

## ---- eval=FALSE--------------------------------------------------------------
#  gl2gds(gl, outfile="test.gds")

## ---- eval=FALSE--------------------------------------------------------------
#  gds <- snpgdsOpen("gl2gds.gds")

## ---- eval=T------------------------------------------------------------------
gl2faststructure(gl, outfile="myfile.fs", probar = FALSE)

