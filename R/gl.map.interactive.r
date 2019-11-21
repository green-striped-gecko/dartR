#' Creates an interactive map (based on latlong) from a genlight object
#' 
#' @param x -- a genlight object [including coordinates within the latlong slot]
#' @param matrix -- a distance matrix between populations or individuals. The matrix is visualised as lines between individuals/populations. If matrix is asymmetric two lines with arrows are plotted.
#' @param standard -- if a matrix is provided line width will be standardised to be between 1 to 10, if set to true, otherwise taken as given.
#' @param symmetric -- if a symmetric matrix is provided only one line is drawn based on the lower triangle of the matrix. If set to false arrows indicating the direction are used instead.
#' @param ind.circles -- should individuals plotted as circles, default is TRUE
#' @param pop.labels -- population labels at the center of the individuals of populations, default is TRUE
#' @param pop.labels.cex -- size of population labels, default is 20.
#' @param provider -- passed to leaflet
#' @return plots a map
#' @importFrom leaflet leaflet addTiles addProviderTiles addCircleMarkers labelOptions addLabelOnlyMarkers addPolylines
#' @importFrom leaflet.minicharts addFlows
#' @export
#' @details A wrapper around the \pkg{leaflet} package. For possible background maps check as specified via the provider: \url{http://leaflet-extras.github.io/leaflet-providers/preview/index.html}
#' @author Bernd Gruber (glbugs@@aerg.canberra.edu.au)
#' @examples 
#' gl.map.interactive(bandicoot.gl)


gl.map.interactive <- function(x, matrix=NULL, standard=TRUE,symmetric=TRUE, ind.circles=TRUE, pop.labels=TRUE, pop.labels.cex=12,  provider="Esri.NatGeoWorldMap")
{

cols <- rainbow(nPop(x))
cols <- substr(cols, 1,7)
df <- x@other$latlong
centers <- apply(df, 2, function(xx) tapply(xx, pop(x), mean, na.rm=TRUE))
# Add default OpenStreetMap map tiles
m <- leaflet() %>%
  addTiles()
if (ind.circles) m <- m %>%  addCircleMarkers(lng=df$lon, lat=df$lat, popup=indNames(x),  color = cols[as.numeric(pop(x))], opacity = 0.8)

if (pop.labels)   m <- m %>% addLabelOnlyMarkers(lng=centers[,"lon"], lat=centers[,"lat"], label = popNames(x),labelOptions = labelOptions(noHide = T, direction = 'top', textOnly = T, textsize = paste0(pop.labels.cex,"px")))


if (!is.null(matrix)){
  
  #standarddise
  if (standard) matrix[,] <- ((matrix[,]-min(matrix, na.rm = T))/(max(matrix, na.rm = T)-min(matrix, na.rm = T)))*9+1
  
  
if (nrow(matrix)==nPop(x)) xys <- centers else xys <- df
if (symmetric) {
for (ii in 1:nrow(matrix)){
  for (i in ii:nrow(matrix)){
    if(!is.null(matrix[i,ii])) m <- m %>% addPolylines(lng=c(xys[i,"lon"],xys[ii,"lon"]), lat=c(xys[i,"lat"], xys[ii,"lat"]), weight = matrix[i,ii], color = "#0000FF" , opacity = 1)
  }
}
}
if (!symmetric) {
  for (i in 1:nrow(matrix))
  {
    for (ii in 1:nrow(matrix))
    {
      if (abs((i-ii))!=0) {
        from <- xys[i,]
        to <- xys[ii,]
        if (!is.null(matrix[i,ii]) & !is.null(matrix[ii,i])) {
        if (matrix[i,ii]>matrix[ii,i])  lcols="#FFAA00" else lcols="#00AAFF"
        if (matrix[i,ii]==matrix[ii,i]) lcols="#00AA00"
        } else lcols ="#333333"
        m <- m %>% addFlows(lng0 = as.numeric(from["lon"]), lng1 =  as.numeric(to["lon"]), lat0=as.numeric( from["lat"]),lat1 =  as.numeric(to["lat"]),flow = matrix[i,ii], color = lcols, maxThickness = 10, minThickness=1, maxFlow = max(matrix, na.rm=T), opacity = 0.8)
      }
    }
  }
  }
}

m %>% addProviderTiles(provider)
}
