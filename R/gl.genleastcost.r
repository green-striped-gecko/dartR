#' Least-cost path analysis based on a friction matrix
#' 
#' This function calculates the pairwise distances (Euclidean, cost path distances and genetic distances) of populations using a friction matrix and a spatial genind object. The genind object needs to have coordinates in the same projected coordinate system as the friction matrix. The friction matrixcan be either a single raster of a stack of several layers. If a stack is provided the specified cost distance is calculated for each layer in the stack. The output of this function can be used with the functions \code{\link[PopGenReport]{wassermann}} or \code{\link[PopGenReport]{lgrMMRR}} to test for the significance of a layer on the genetic structure.
#' @param x a spatial genind object. see ?popgenreport how to provide coordinates in genind objects
#' @param fric.raster a friction matrix
#' @param gen.distance specification which genetic distance method should be used to calculate pairwise genetic distances between populations ( "D", "Gst.Nei", "Gst.Hedrick") or individuals ("Smouse", "Kosman", "propShared")
#' @param NN Number of neighbours used when calculating the cost distance (possible values 4,8 or 16). As the default is NULL a value has to be provided if pathtype='leastcost'. NN=8 is most commonly used. Be aware that linear structures may cause artefacts in the least-cost paths, therefore inspect the actual least-cost paths in the provided output.
#' @param pathtype Type of cost distance to be calculated (based on function in the \code{\link{gdistance}} package. Available distances are 'leastcost', 'commute' or 'rSPDistance'. See functions in the gdistance package for futher explanations. If the path type is set to 'leastcost' then paths and also pathlength are returned.
#' @param plotpath switch if least cost paths should be plotted (works only if pathtype='leastcost'. Be aware this slows down the computation, but it is recommended to do this to check least cost paths visually.
#' @param theta value needed for rSPDistance function. see \code{\link{rSPDistance}} in package \code{gdistance}.
#' @importFrom mmod pairwise_D pairwise_Gst_Hedrick pairwise_Gst_Nei
#' @importFrom PopGenReport lgrMMRR wassermann gd.smouse gd.kosman
#' @importFrom gdistance transition costDistance shortestPath geoCorrection
#' @importFrom stats step
#' @importFrom sp Line Lines SpatialLines SpatialLinesLengths
#' @return returns a list that consists of four pairwise distance matrixes (Euclidean, Cost, length of path and genetic) and the actual paths as spatial line objects.
#' @author Bernd Gruber (bugs? Post to \url{https://groups.google.com/d/forum/dartr})
#' @seealso \code{\link{landgenreport}}, \code{\link{popgenreport}}, \code{\link{wassermann}}, \code{\link{lgrMMRR}}
#' @references Cushman, S., Wasserman, T., Landguth, E. and Shirk, A. (2013). Re-Evaluating Causal Modeling with Mantel Tests in Landscape Genetics. Diversity, 5(1), 51-72.
#' Landguth, E. L., Cushman, S. A., Schwartz, M. K., McKelvey, K. S., Murphy, M. and Luikart, G. (2010). Quantifying the lag time to detect barriers in landscape genetics. Molecular ecology, 4179-4191.
#' Wasserman, T. N., Cushman, S. A., Schwartz, M. K. and Wallin, D. O. (2010). Spatial scaling and multi-model inference in landscape genetics: Martes americana in northern Idaho. Landscape Ecology, 25(10), 1601-1612.
#' @examples
#' \donttest{
#' data(possums.gl)
#' library(raster)  #needed for that example
#' landscape.sim <- readRDS(system.file("extdata","landscape.sim.rdata", package="dartR"))
#' glc <- gl.genleastcost(x=possums.gl,fric.raster=landscape.sim , 
#' gen.distance = "D", NN=8, pathtype = "leastcost",plotpath = TRUE)
#' library(PopGenReport)
#' wassermann(eucl.mat = glc$eucl.mat, cost.mat = glc$cost.mats,  gen.mat = glc$gen.mat)
#' lgrMMRR(gen.mat = glc$gen.mat, cost.mats = glc$cost.mats,  eucl.mat = glc$eucl.mat)
#' }
#'
#' @export
gl.genleastcost <- function(x, fric.raster, gen.distance, NN=NULL, pathtype="leastcost", plotpath=TRUE, theta=1)
{
# 
# pairwise_D2 <- function (x, linearized = FALSE) 
# {
#       pops <- seppop(x)
#       n.pops <- length(pops)
#       allP <- utils::combn(1:n.pops, 2)
#       pair <- function(index.a, index.b) {
#         a <- pops[[index.a]]
#         b <- pops[[index.b]]
#         temp <- repool(a, b)
#         return(D_Jost(temp)$global.het)
#       }
#       res <- sapply(1:dim(allP)[2], function(i) pair(allP[, i][1], 
#                                                      allP[, i][2]))
#       attributes(res) <- list(class = "dist", Diag = FALSE, Upper = FALSE, 
#                               Labels = popNames(x), Size = n.pops)
# 
#       if (linearized) res <- res/(1 - res)
#       return(res)
# }
#   
#   pairwise_Gst_Hedrick2 <-
#     function (x, linearized = FALSE) 
#     {
#       pops <- seppop(x)
#       n.pops <- length(pops)
#       allP <- utils::combn(1:n.pops, 2)
#       pair <- function(index.a, index.b) {
#         a <- pops[[index.a]]
#         b <- pops[[index.b]]
#         temp <- repool(a, b)
#         return(Gst_Hedrick(temp)$global)
#       }
#       res <- sapply(1:dim(allP)[2], function(i) pair(allP[, i][1], 
#                                                      allP[, i][2]))
#       attributes(res) <- list(class = "dist", Diag = FALSE, Upper = FALSE, 
#                               Labels = popNames(x), Size = n.pops)
#        if (linearized) res <- res/(1 - res)
#       return(res)
#      }
#   
#   pairwise_Gst_Nei2 <-
#     function (x, linearized = FALSE) 
#     {
#       pops <- seppop(x)
#       n.pops <- length(pops)
#       allP <- utils::combn(1:n.pops, 2)
#       pair <- function(index.a, index.b) {
#         a <- pops[[index.a]]
#         b <- pops[[index.b]]
#         temp <- repool(a, b)
#         return(Gst_Nei(temp)$global)
#       }
#       res <- sapply(1:dim(allP)[2], function(i) pair(allP[, i][1], 
#                                                      allP[, i][2]))
#       attributes(res) <- list(class = "dist", Diag = FALSE, Upper = FALSE, 
#                               Labels = popNames(x), Size = n.pops)
#       if (linearized) res <- res/(1 - res)
#       return(res)
#     }
#   
# 
#   
  
  
  if (is.null(NN) & pathtype=="leastcost") 
  {
    stop("NN is not specified!\nPlease specify the number of nearest neighbour to use for the least-cost path calculations (NN=4 or NN=8). If linear features are tested you may want to consider NN=4 otherwise NN=8 is the most commonly used and prefered option. In any case check the actual least-cost paths for artefacts by inspecting the plot on least-cost paths.\n")

  }
  
dist.type<-NA
if (gen.distance=="D" || gen.distance=="Gst.Hedrick" || gen.distance=="Gst.Nei") dist.type<- "pop" 

if (gen.distance=="Kosman" || gen.distance=="Smouse" || gen.distance=="propShared") dist.type<- "ind" 

if (is.na(dist.type)) 
  {stop("No valid genetic distance type was provided. Please check ?landgenreport for valid options\n")

}

if (is.null(x@other$xy)) 
  
{
  cat("No projected coordinates in @other$xy found. Hence will use latlongs (if provided), which are not projected, hence there might be distortions if the area covered is large or close to the poles. Be aware your resistance layer and coordinates in the genlight object need to have the same coordinate system.\n")
  x@other$xy <- x@other$latlong[,c("lon","lat")]
  if (is.null(x@other$xy)) step("No coordinates found in the genlight object!!\n")
}


if (dist.type=="pop")
{
#calculate the centers if population meassurment is wanted
c.x <- tapply(x@other$xy[,1],x@pop, mean)
c.y <- tapply(x@other$xy[,2],x@pop, mean)
cp<-cbind(c.x, c.y)
eucl.mat <- as.matrix(dist(cp))
dimnames(eucl.mat) <- list(popNames(x), popNames(x)) 
npop <- length(levels(x@pop))

} else 
{
cp <- cbind(x@other$xy[,1], x@other$xy[,2])
eucl.mat <- as.matrix(dist(cp))
dimnames(eucl.mat) <- list(indNames(x), indNames(x)) 
npop <- length(indNames(x))
}





#check if fric.raster is a stack or not...
mats <- list()
mats.names<- NA
mats.pathlength<- list()
mats.paths<- list()

pathlength.mat <- NULL
paths <- NULL



n.mats <- dim(fric.raster)[3] #number of rasters in the stack

for (ci in 1:n.mats)
{

plot(fric.raster[[ci]], main=paste(names(fric.raster)[ci],":",pathtype,", NN=",NN,sep=""))
 #image(fric.raster, col=fric.raster@legend@colortable, asp=1)

points(x@other$xy,cex=1, pch=16, col=rainbow(nPop(x))[as.numeric(pop(x))] )
if (dist.type=="pop")  points(cp,cex=1.5 , pch= 15, col="black")


#create friction matrix
fric.mat <- transition(fric.raster[[ci]],function(x) 1/x[2],NN)

#set distances to meters  if not projected already
fric.mat@crs@projargs<- "+proj=merc +units=m"
fric.mat.cor <- geoCorrection(fric.mat)

if (pathtype=="leastcost") cd.mat <-costDistance(fric.mat.cor, cp, cp)

if (pathtype=="rSPDistance") cd.mat <- rSPDistance(fric.mat.cor, cp, cp, theta=1)

if (pathtype=="commute") cd.mat <-as.matrix(commuteDistance(fric.mat.cor, cp))
dimnames(cd.mat) <- dimnames(eucl.mat) 




if (pathtype=="leastcost" & plotpath==TRUE)   #only show paths if leastcost otherwise not possible
{
comb <- t(combn(1:npop,2))


#pathlength matrix
pathlength.mat <- cd.mat
pathlength.mat[,] <- 0
paths<- list()

cols <- rainbow(dim(comb)[1], alpha=0.5)
for (i in 1:dim(comb)[1])
{

if (dist(rbind(cp[comb[i,1],], cp[comb[i,2],]))==0)
{
 ll <- Line(rbind( cp[comb[i,1],], cp[comb[i,2],]))
 S1 <- Lines(list(ll),ID="Null")
 sPath <- SpatialLines(list(S1))
} else 
{
sPath <- shortestPath(fric.mat.cor, cp[comb[i,1],], cp[comb[i,2],], output="SpatialLines")
}

lines(sPath, lwd=1.5, col=cols[i])
paths[[i]] <- sPath
ll <-  round(SpatialLinesLengths(sPath),3)
pathlength.mat[comb[i,1],comb[i,2]] <- ll
pathlength.mat[comb[i,2],comb[i,1]] <- ll
}

}

mats[[ci]] <- cd.mat
mats.names[[ci]] <- names(fric.raster)[ci]
mats.pathlength[[ci]] <- pathlength.mat
mats.paths[[ci]] <- paths



} #end of ci loop

#mats[[n.mats+1]] <- eucl.mat
#mats.names[n.mats+1]<- "Euclidean"
#
names(mats)  <- names(fric.raster)
#put other calculations here....
# Calculate genetic distances across subpopulations

xx <- gl2gi(x, v=0)

if (gen.distance=="Gst.Nei")
{
gendist.mat<-as.matrix(pairwise_Gst_Nei(xx))
}
if (gen.distance=="Gst.Hedrick")
{
gendist.mat<-as.matrix(pairwise_Gst_Hedrick(xx))
}
if (gen.distance=="D")
{
gendist.mat<-as.matrix(pairwise_D(xx))
}

if (gen.distance=="Smouse")
{
gendist.mat <- as.matrix(gd.smouse(xx,verbose=FALSE))
}
if (gen.distance=="Kosman")
{
gendist.mat <-as.matrix(as.dist(gd.kosman(xx)$geneticdist))
}
if (gen.distance=="propShared")
{
  gendist.mat <-as.matrix(as.dist(propShared(xx)))
}

  dimnames(gendist.mat)<-dimnames(eucl.mat)


return(list( gen.mat=gendist.mat, eucl.mat=eucl.mat,cost.matnames=mats.names, cost.mats=mats,pathlength.mats= mats.pathlength,  paths=mats.paths))
}


