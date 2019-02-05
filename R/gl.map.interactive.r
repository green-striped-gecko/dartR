#' Creates an interactive map (based on latlong) from a genlight object
#' 
#' @param x -- a genlight object [including coordinates within the latlong slot]
#' @param provider -- passed to leaflet
#' @return plots a map
#' @importFrom leaflet leaflet addTiles addProviderTiles addCircleMarkers
#' @export
#' @details A wrapper around the \pkg{leaflet} package. For possible background maps check as specified via the provider: \url{http://leaflet-extras.github.io/leaflet-providers/preview/index.html}
#' @author Bernd Gruber (glbugs@@aerg.canberra.edu.au)
#' @examples 
#' #gl.map(bandicoot.gl)


gl.map.interactive <- function(x, provider="Esri.NatGeoWorldMap")
{

cols <- rainbow(nPop(x))
cols <- substr(cols, 1,7)
df <- x@other$latlong
m <- leaflet() %>%
  addTiles() %>%  # Add default OpenStreetMap map tiles
  addCircleMarkers(lng=df$lon, lat=df$lat, popup=indNames(x),  color = cols[as.numeric(pop(x))], opacity = 0.8)
m %>% addProviderTiles(provider)
}
