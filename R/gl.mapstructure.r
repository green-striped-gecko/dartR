# cex =1.5
# 
# sx <- abs(diff(range(centers[,"lon"])))/(10*10)*cex
# sy <- sx/max(ff)* cex
# 
# ii=1
# m1 <- m
# for (ii in 1:5) {
#   for ( i in 1:5) {
#     oo <- (i-1)*sx
#     
#     m1 <- m1   %>% addRectangles(cx[ii]+oo, cy[ii], cx[ii]+oo+sx, cy[ii]+ff[ii,i]*sy, opacity = 0, color =  rainbow(5)[i], fillOpacity = 0.5)
#     
#   }
# }
# 
# m1 %>% leaflet::addProviderTiles(provider) %>% addLegend(labels=colnames(ff), colors=rainbow(5),position ="topright" )
