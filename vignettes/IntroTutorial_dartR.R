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
#  platy.gl <- gi2gl(platy)

## ---- echo=FALSE--------------------------------------------------------------
platy.gl <- gi2gl(platy, parallel = FALSE)

## -----------------------------------------------------------------------------
platy.gl@other$ind.metrics <- platy.gl@other$data

## -----------------------------------------------------------------------------
ts <- sapply(1:nLoc(platy.gl), function(x) paste(sample(c("A","T","G","C"), 50, replace = T),
                                                 collapse = ""))
df.loc <- data.frame(RepAvg = runif(nLoc(platy.gl)),  TrimmedSequence=ts)

platy.gl@other$loc.metrics <- df.loc



