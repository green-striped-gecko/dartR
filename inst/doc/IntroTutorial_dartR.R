## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("dartR")

## ---- eval=FALSE---------------------------------------------------------
#  # Install and attach library dartR
#    install.packages("devtools")
#    library(devtools)
#    source("http://bioconductor.org/biocLite.R")
#    biocLite("qvalue", suppressUpdates=T)
#    biocLite("SNPRelate", suppressUpdates=T)
#    install_github("green-striped-gecko/dartR")

## ------------------------------------------------------------------------
 library(dartR)

## ---- eval=FALSE---------------------------------------------------------
#  # Set the default working directory (change this to suit)
#  setwd("c:/your.working.directory/")

## ------------------------------------------------------------------------
# Rename the test genlight object to gl, something simple
gl <- testset.gl

## ------------------------------------------------------------------------
m <- as.matrix(gl)

## ------------------------------------------------------------------------
as.matrix(gl)[1:5,1:3]

## ---- eval=FALSE---------------------------------------------------------
#  gl <- gl.read.dart(filename = "testset.csv", covfilename = " ind_metrics.csv")
#  

## ------------------------------------------------------------------------
dartfile <- system.file("extdata","testset_SNPs_2Row.csv", package="dartR")
covfilename <- system.file("extdata","testset_metadata.csv", package="dartR")
gl <- gl.read.dart(filename=dartfile, covfilename = covfilename, probar=FALSE)

## ------------------------------------------------------------------------
gl

## ------------------------------------------------------------------------
as.matrix(gl)[1:3,1:3]

## ------------------------------------------------------------------------
nLoc(gl)
nInd(gl)
nPop(gl)

## ------------------------------------------------------------------------
levels(pop(gl))[1:5]

## ---- eval=FALSE---------------------------------------------------------
#  gl <- gl.read.dart(filename="mydata.csv", covfilename = "my.metadata.csv")

## ------------------------------------------------------------------------
#Only the entries for the first ten individuals are shown
gl@other$loc.metrics$RepAvg[1:10]

## ------------------------------------------------------------------------
names(gl@other$loc.metrics)

## ------------------------------------------------------------------------
#only first 10 entries showns
gl@other$ind.metrics$sex[1:10]

## ------------------------------------------------------------------------
read.csv( paste(.libPaths()[1],"/dartR/extdata/platy.csv",sep="" ))

## ------------------------------------------------------------------------
#you might need to install PopGenReport via
#install.packages("PopGenReport")
library(PopGenReport) 
platy <- read.genetable( paste(.libPaths()[1],"/dartR/extdata/platy.csv",
sep="" ), ind=1, pop=2, lat=3, long=4, other.min=5, other.max=6, oneColPerAll=FALSE,sep="/")

platy


## ---- eval=FALSE---------------------------------------------------------
#  platy.gl <- (gi2gl(platy)

## ---- echo=FALSE---------------------------------------------------------
platy.gl <- (gi2gl(platy, parallel = FALSE))

## ------------------------------------------------------------------------
platy.gl@other$ind.metrics <- platy.gl@other$data

## ------------------------------------------------------------------------
ts <- sapply(1:nLoc(platy.gl), function(x) paste(sample(c("A","T","G","C"), 50, replace = T),
                                                 collapse = ""))
df.loc <- data.frame(RepAvg = runif(nLoc(platy.gl)),  TrimmedSequence=ts)

platy.gl@other$loc.metrics <- df.loc


## ------------------------------------------------------------------------
gl.report.callrate(platy.gl)
gl2 <- gl.filter.repavg(platy.gl, t=0.5)

## ------------------------------------------------------------------------
platy.gl@position <- as.integer(runif(nLoc(platy.gl),2,49))
platy.gl@loc.all <- testset.gl@loc.all[1:6]

## ---- eval=FALSE---------------------------------------------------------
#  gl2fasta(platy.gl)

## ------------------------------------------------------------------------
gl2 <- gl.filter.callrate(gl, method = "loc", threshold = 0.95)

## ------------------------------------------------------------------------
gl2 <- gl.filter.callrate(gl, method="ind", threshold = 0.90)


## ------------------------------------------------------------------------
gl2 <- gl.filter.repavg(gl, t=1)

## ---- eval=F, echo=F-----------------------------------------------------
#  gl2 <- gl.filter.secondary(gl)
#  

## ------------------------------------------------------------------------
gl2 <- gl.filter.monomorphs(gl, v=0)

## ------------------------------------------------------------------------
gl2 <- gl.filter.hamming(testset.gl, t=0.25, probar = F)

## ---- eval=F-------------------------------------------------------------
#  gl2 <- gl.filter.callrate(gl, method = "loc", threshold = 0.95)
#  gl3 <-  gl.filter.callrate(gl2, method="ind", threshold = 0.90)
#  gl4 <- gl.filter.repavg(gl3, t=1)

## ------------------------------------------------------------------------
#population names (#30 populations)
levels(pop(gl))
#table on individuals per population
table(pop(gl))


## ---- fig.height=5-------------------------------------------------------
barplot(table(pop(gl)), las=2)

## ---- eval=T-------------------------------------------------------------
gl.make.recode.pop(gl, outfile = file.path(tempdir(),"new_pop_assignments.csv"))

## ------------------------------------------------------------------------
glnew <- gl.recode.pop(gl, pop.recode=file.path(tempdir(),"new_pop_assignments.csv"))

## ------------------------------------------------------------------------
levels(pop(gl))

## ---- eval=FALSE---------------------------------------------------------
#  glnew2 <- gl.edit.recode.pop(gl, pop.recode = file.path(tempdir(),"new_pop_assingments.csv"))

## ------------------------------------------------------------------------
#only first 10 entries are shown
indNames(gl)[1:10]



## ------------------------------------------------------------------------
gl.make.recode.ind(gl, outfile=file.path(tempdir(),"new_ind_assignments.csv"))

## ---- eval=FALSE---------------------------------------------------------
#  glnew3 <- gl.recode.ind(gl, ind.recode=file.path(tempdir(),"new_ind_assignments.csv"))
#  

## ---- eval=F-------------------------------------------------------------
#  gl <- gl.edit.recode.ind(gl, ind.recode=file.path(tempdir(),"new_ind_assignments.csv"))

## ------------------------------------------------------------------------
gl_new <- gl[gl$pop!="EmmacBrisWive", ]

## ------------------------------------------------------------------------
glsub <- gl[1:7, 1:3]
glsub

## ------------------------------------------------------------------------
dim(glsub@other$ind.metrics)

dim(glsub@other$loc.metrics)

## ------------------------------------------------------------------------

index.ind <- pop(gl)=="EmmacRussEube" | pop(gl)=="EmvicVictJasp"
#check if the index worked
table( pop(gl), index.ind)


index.loc <- sample(nLoc(gl), 30, replace = F)
index.loc


## ------------------------------------------------------------------------
glsub2 <- gl[index.ind, index.loc]
glsub2@other$ind.metrics <- gl@other$ind.metrics[index.ind,] #not necessary
glsub2@other$loc.metrics <- gl@other$loc.metrics[index.loc,] #necessary

## ------------------------------------------------------------------------
glsub2
dim(glsub2@other$ind.metrics)
dim(glsub2@other$loc.metrics)

## ------------------------------------------------------------------------
d <- gl.dist(gl[1:5,1:10])
d

## ---- eval=FALSE---------------------------------------------------------
#  library(StAMPP) #you may need to install the package
#  pwfst <-stamppFst(gl[1:20,], nboots=1, percent=95, nclusters=1)
#  round(pwfst,3)

## ---- eval=FALSE---------------------------------------------------------
#  pwGst <-stamppNeisD(gl[1:20,]) #no parallel version :-()
#  round(pwGst,3)

## ---- warning=FALSE------------------------------------------------------
library(mmod) #you may need to install the package first
#for performance reason use only a subset (and recode the populations)
recpops<- factor(rep(LETTERS[1:5],50))
glsub <- gl
pop(glsub)<-recpops

gi <- gl2gi(glsub, probar = FALSE)
round(pairwise_D(gi),4)
round(pairwise_Gst_Hedrick(gi),4)
round(pairwise_Gst_Nei(gi),4)

## ------------------------------------------------------------------------
pc <- gl.pcoa(gl, nfactors=5)

## ------------------------------------------------------------------------
names(pc)

## ------------------------------------------------------------------------
barplot(pc$eig/sum(pc$eig)*100, )

## ---- fig.height=5-------------------------------------------------------
gl.pcoa.plot(pc, gl, labels="pop", xaxis=1, yaxis=2)


## ----eval=FALSE----------------------------------------------------------
#  glnew <- gl.edit.recode.pop(gl)

## ------------------------------------------------------------------------
glnew <- gl
levels(pop(glnew)) <- c(rep("Cooper",13), rep("MDB", 8 ), rep("Emmac_Coast",7),"EmsubRopeMata" ,  "EmvicVictJasp")
gl.pcoa.plot(pc, glnew, labels="pop", xaxis=1, yaxis=2)

## ---- eval=FALSE---------------------------------------------------------
#  gl.pcoa.plot(pc, glnew, labels="interactive", xaxis=1, yaxis=2)
#  ggplotly()
#  

## ---- fig.height=4-------------------------------------------------------
gl.pcoa.scree(pc)

## ---- eval=FALSE---------------------------------------------------------
#  gl.pcoa.plot.3d(pc, glnew)

## ------------------------------------------------------------------------
gl.tree.nj(glnew, type="fan")


## ---- eval=FALSE---------------------------------------------------------
#  gl.collapse.recursive(gl, t=0)

## ---- echo=FALSE, eval=FALSE---------------------------------------------
#  gl <- testset.gl
#  gl.collapse.recursive(gl, t=0)
#  

## ---- fig.height=4-------------------------------------------------------
gl <- gl.ibd(gl=testset.gl[1:180,])

## ---- eval=FALSE---------------------------------------------------------
#  gl <- testset.gl
#  phy <- gl2phylip(gl, outfile="turtle.phy", bstrap=1000)

## ---- eval=FALSE---------------------------------------------------------
#  gl2fasta(gl, method=2, outfile="nohets.fasta")

## ------------------------------------------------------------------------
gl.report.bases(testset.gl)

## ---- eval=FALSE---------------------------------------------------------
#  gl2fasta(gl, method=4, outfile="nohets.fasta")
#  
#  

## ---- eval=FALSE---------------------------------------------------------
#  gl2fasta(gl, method=1, outfile="ambcodes.fasta")
#  

## ---- eval=FALSE---------------------------------------------------------
#  gl2fasta(gl, method=3, outfile="ambcodes.fasta")
#  

## ------------------------------------------------------------------------
gl.report.bases(testset.gl)

## ------------------------------------------------------------------------
x <- gl.report.pa(testset.gl, id="UC_00146", nmin=10, t=0)

## ---- fig.height=4-------------------------------------------------------
x <- gl.assign(testset.gl, id="UC_00146", nmin=10, alpha=0.95, t=1)

## ------------------------------------------------------------------------
gl <- testset.gl
gi <- gl2gi(gl, probar=FALSE)


## ---- eval=FALSE---------------------------------------------------------
#  gl2 <- gi2gl(gi)

## ---- eval=T-------------------------------------------------------------
glnew <- gl2nhyb(gl, outfile = file.path(tempdir(),"nhyb.txt"))

## ---- eval=FALSE---------------------------------------------------------
#  gl.new <- gl2nhyb(gl, outfile = "nhyb.txt", p0 = NULL,p1 = NULL, t = 0,   m = "random")
#  

## ---- eval=FALSE---------------------------------------------------------
#  glnew <- gl2phylip(outfile = "phyinput.txt")
#  

## ---- eval=FALSE---------------------------------------------------------
#  gl.new <- gl2phylip(outfile = "phyinput.txt", bstrap = 1000)

## ---- eval=FALSE---------------------------------------------------------
#  gl2gds(gl, outfile="test.gds")

## ---- eval=FALSE---------------------------------------------------------
#  gds <- snpgdsOpen("gl2gds.gds")

## ---- eval=T-------------------------------------------------------------
gl2faststructure(gl, outfile=file.path(tempdir(),"myfile.fs"), probar = FALSE)

