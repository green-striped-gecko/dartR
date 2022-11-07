# Ternary Coordinate System
# 
# \code{coord_tern} is a function which creates a transformation mechanism between the ternary system, and, the cartesian system.
# It inherits from the fixed coordinate system, employing fixed ratio between x and y axes once transformed.
# 
# @section Aesthetics (Required in Each Layer):
# \Sexpr[results=rd,stage=build]{ggtern:::rd_aesthetics("coord", "tern")}
# 
# Abovementioned limitations include the types of geometries which can be used (ie approved geometries), 
# or modifications to required aesthetic mappings. One such essential patch is, for approved geometries previously 
# requiring \code{x} and \code{y} coordinates, now require an additional \code{z} coordinate, and, 
# \code{\link{geom_segment}} goes one step further in that it requires both an additional 
# \code{z} and \code{zend} coordinate mappings. 
# 
# In essence, the required aesthetics are the product between what
# is required of each 'layer' and what is required of the 'coordinate system'.
# @param Tlim the range of T in the ternary space
# @param Llim the range of L in the ternary space
# @param Rlim the range of R in the ternary space
# @inheritParams ggplot2:::coord_cartesian
# @return \code{coord_tern} returns a CoordTern ggproto
# @rdname coord_tern
# @author Nicholas Hamilton
# @export
coord_tern <- function(Tlim = NULL, Llim = NULL, Rlim = NULL, expand = TRUE){
  rs      = CoordTern$required_scales
  ra      = CoordTern$required_aes
  mapping = sapply(rs,function(x) getOption(sprintf('tern.default.%s',x)) )
  if(!all(ra %in% as.character(mapping)))
    stop(sprintf("Options for %s are %s, must be assigned and NOT duplicated, alter the 'tern.default.X' (X=%s) global option",
                 joinCharacterSeries(rs,'and'),
                 joinCharacterSeries(ra,'and'),
                 joinCharacterSeries(rs,'or')),call.=FALSE)
  ggproto(NULL, CoordTern,
          mapping         = as.list(mapping),
          limits          = list(x = c(0,1), y = c(0,1)*.ratio(), T=Tlim, L=Llim, R=Rlim),
          ratio           = 1,
          expand          = expand,
          scales          = list(),
          labels_coord    = list(),
          theme           = theme_get()
  )
}

# @rdname coord_tern
# @name coord_tern
# @usage NULL
# @format NULL
# @export
CoordTern <- ggproto("CoordTern", CoordCartesian,
                     required_aes    = c("x","y","z"),
                     required_scales = c("T","L","R"),
                     aspect          = function(self, ranges) { 
                       (diff(ranges$y.range) / diff(ranges$x.range)) * self$ratio 
                     },
                     transform       = function(self, data, panel_params){
                       
                       #Variables and functions.
                       ix               = c('x','y')
                       angle            = .theme.get.rotation(self)
                       addOrigin        = function(d,ix,o){  d[,ix] = t(t(d[,ix,drop=FALSE]) + o); d }
                       data             = tlr2xy(data,self)
                       coordShift       = sapply(c('tern.axis.hshift','tern.axis.vshift'),function(x){ calc_element(x,self$theme) %||% 0})
                       
                       #For each coordinate group, conduct the re-centering, translation / rotation process
                       for(group in .get.sets(ix,names(data)) ){
                         ix.comb        = .combos(ix,group)
                         
                         #Run check to see if all the required columns are present within the group.
                         #A required column is the cartesian product of the group and the required index, so it could be (xend,yend), (xstart,ystart), etc...
                         #This shouldn't be a problem, since the same check is executed during the master transformation in the tlr2xy(...) operation above
                         if(!all(ix.comb %in% names(data))) stop(sprintf("Problem with proposed tranlation/rotation, require columns (%s) wich are not present.",
                                                                         joinCharacterSeries(setdiff(ix.comb,names(data)),lastWord='and') ))
                         
                         #Determine Origin
                         xtrm           = tlr2xy(.get.tern.extremes(self,panel_params,FALSE),self)
                         origin         = apply(xtrm[,ix],2,function(x){mean(x)}) ##ORIGIN TO BE CENTROID
                         
                         #Any multiple of 360 degrees can be disregarded
                         if(angle %% 360 != 0){
                           #Translation   
                           data           = addOrigin(data,ix.comb,-origin); 
                           xtrm           = addOrigin(xtrm,ix,     -origin)
                           #Rotation  
                           data[,ix.comb] = .rotation(data[,ix.comb],angle) 
                           xtrm[,ix]      = .rotation(xtrm[,ix],     angle)
                           #Translate Back
                           data           = addOrigin(data,ix.comb,+origin)
                           xtrm           = addOrigin(xtrm,ix,     +origin)
                         }
                         
                         #Re-Center
                         origin         = apply(xtrm[,ix],2,function(x){ mean(range(x))} ) ##ORIGIN TO BE MEDIAN
                         target         = c(mean(panel_params$x.range),mean(panel_params$y.range)) + as.numeric(coordShift)
                         data           = addOrigin(data,ix.comb,(target - origin))
                       }
                       
                       #self$super$super$transform(data,panel_params)
                       self$super()$super()$transform(data,panel_params)
                     },
                     render_axis_h   = function(panel_params, theme){
                       list(
                         top = zeroGrob(),
                         bottom = zeroGrob()
                       )
                     },
                     render_axis_v   = function(panel_params, theme){
                       list(
                         left = zeroGrob(),
                         right = zeroGrob()
                       )
                     },
                     
                     render_bg       = function(self,panel_params, theme){
                       items = list() 
                       extrm = .get.tern.extremes(self,panel_params)
                       items = .render.background(self,extrm,theme,items)
                       if(!.theme.get.gridsontop(theme)){
                         items = .render.grid(self,extrm,panel_params,theme,items)
                         if(!.theme.get.showmask(theme))
                           items = .render.fgset(self,extrm,panel_params,theme,items)
                       }
                       # grid::gTree(children = do.call("gList",items))
                       grid::gTree(children = do.call(grid::gList,items))
                       
                     },
                     
                     render_fg     = function(self,panel_params, theme){
                       items = list() 
                       extrm = .get.tern.extremes(self,panel_params)
                       if(.theme.get.gridsontop(theme)){
                         items = .render.grid( self,extrm,panel_params,theme,items)
                         items = .render.fgset(self,extrm,panel_params,theme,items)
                       }
                       items = .render.titles(self,extrm,panel_params,theme,items)
                       # grid::gTree(children = do.call("gList", items))
                       grid::gTree(children = do.call(grid::gList, items))
                       
                     },
                     
                     remove_labels = function(self,table){
                       #Determine the Layout
                       layout <- table$layout
                       #Remove Y-axis columns
                       ycols  <- layout[grepl("^ylab-l", layout$name), , drop = FALSE]
                       table  <- table[,-ycols$l]
                       #Remove X-axis Rows
                       xrows  <- layout[grepl("^xlab-b", layout$name), , drop = FALSE]
                       table  <- table[-xrows$t,]
                       table
                     },
                     
                     setup_panel_params = function(self, scale_x, scale_y, params = list() ) {
                       
                       train_cartesian <- function(panel_params,limits,name,continuousAmount) {
                         if (self$expand) {
                           # expand <- ggint$expand_default(panel_params,continuous=continuousAmount)
                           # expand = ggplot2:::expand_limits_continuous(limits, continuousAmount, panel_params)
                           # expand <- c(continuousAmount, continuousAmount)
                           expand <- c(0, 0)
                         } else {
                           expand <- c(0, 0)
                         }
                         print(limits)
                         if (is.null(limits)) {
                           range <- panel_params$dimension(expand)
                         } else {
                           range <- range(panel_params$transform(limits))
                           range <- expand_range(range, expand[1], expand[2])
                         }
                         
                         out <- panel_params$break_info(range)
                         names(out) <- paste(name, names(out), sep = ".")
                         out
                       }
                       
                       #Determine the epansiion amount
                       expand.amount = calc_element('tern.panel.expand',theme=self$theme)
                       
                       #Adjust for rotation
                       extremes        = .get.tern.extremes(self,list(x.range=self$limits$x, y.range=self$limits$y,x.rescale=FALSE,y.rescale=FALSE))
                       extremes        = extremes[,c('x','y')]
                       currentMidpoint = c(mean(self$limits$x),mean(self$limits$y))
                       ternaryMidpoint = apply(extremes,2,function(x){mean(range(x))})
                       shift           = 0.5*(currentMidpoint - ternaryMidpoint)
                       
                       #Execute the Training
                       c(
                         # train_cartesian(scale_x, self$limits$x - shift[1],"x",c(expand.amount,0) ),
                         # train_cartesian(scale_y, self$limits$y - shift[2],"y",c(expand.amount,0) )
                         view_scales_from_scale(scale_x, self$limits$x, TRUE, c(expand.amount,0)),
                         view_scales_from_scale(scale_y, self$limits$y, TRUE, c(expand.amount,0))
                       )
                     }
)

# Inspired by ggplot2
view_scales_from_scale <- function(scale, coord_limits = NULL, expand = TRUE, expansion) {
  # expansion <- default_expansion(scale, expand = expand)
  limits <- scale$get_limits()
  continuous_range <- ggint$expand_limits_scale(scale, expansion, limits, coord_limits = coord_limits)
  aesthetic <- scale$aesthetics[1]
  
  view_scales <- list(
    ggint$view_scale_primary(scale, limits, continuous_range),
    sec = ggint$view_scale_secondary(scale, limits, continuous_range),
    arrange = scale$axis_order(),
    range = continuous_range
  )
  names(view_scales) <- c(aesthetic, paste0(aesthetic, ".", names(view_scales)[-1]))
  
  view_scales
}

#----------------------------------------------------------------------------------
#Internals >>>> Ratio
#----------------------------------------------------------------------------------
.ratio = function(){ 0.5*tan(60*pi/180) }

#----------------------------------------------------------------------------------
#Internals >>>> Ternary Extremes
#----------------------------------------------------------------------------------
.get.tern.extremes <- function(self,panel_params,transform=TRUE){
  ex = .check.tern.extremes(self)
  if(transform) ex = self$transform(ex,panel_params)
  ex
}

.check.tern.extremes <- function(self,precision=7){
  ix           = self$required_scales
  scales       = self$scales[ix]
  limitsfunc   = Vectorize(function(r,c){ do.call(if(r == c)'max'else'min',list( scales[[ ix[c] ]]$limits %||% self$limits[[ ix[c] ]] %||% c(0,1) ))  })
  ex           = as.data.frame(outer(1:3,1:3,limitsfunc))
  colnames(ex) = as.character(self$mapping[ix])
  sums         = round(apply(ex,1,sum),precision)
  if(!all(sums == 1.0)){
    colnames(ex) = ix; ex$Sum = sums; print(ex)
    stop("Invalid Ternary Limits, Each Point Must Sum to Unity...",call.=FALSE) 
  }
  invisible(ex)
}

#----------------------------------------------------------------------------------
#Internals >>>> ANGLES
#Functions to determine the rotation angles for the various components
#----------------------------------------------------------------------------------
.get.angles             <- function(clockwise){ c(0,120,240) - as.logical(clockwise)*180 }
.get.angles.arrows      <- function(clockwise){ .get.angles(clockwise) + if(clockwise){-30}else{30} }
.get.angles.arrowmarker <- function(clockwise){ x = c(60,0,-60); if(clockwise){x}else{x[c(3,1:2)]}}
.get.angles.ticklabels  <- function(clockwise){ c(0,-60,60) }


#----------------------------------------------------------------------------------
#Internals >>>> Theme flags.
#----------------------------------------------------------------------------------
.theme.get.clockwise <- function(theme){
  clockwise = calc_element('tern.axis.clockwise',theme)
  ifthenelse(is.logical(clockwise),clockwise[1],getOption("tern.clockwise"))
}
.theme.get.gridsontop <- function(theme){
  ret = calc_element('tern.panel.grid.ontop',theme)
  ifthenelse(is.logical(ret),ret[1],FALSE)
}

.theme.get.bordersontop <- function(theme){
  ret = calc_element('tern.axis.line.ontop',theme)
  ifthenelse(is.logical(ret),ret[1],FALSE)
}

.theme.get.showtitles <- function(theme){
  showtitles = calc_element("tern.axis.title.show",theme=theme)
  ifthenelse(is.logical(showtitles),showtitles[1],getOption("tern.title.show"))
}
.theme.get.showlabels <- function(theme){
  showlabels = calc_element("tern.axis.text.show",theme=theme)
  ifthenelse(is.logical(showlabels),showlabels[1],getOption("tern.text.show"))
}
.theme.get.showgrid.major <- function(theme){
  showgrid   = calc_element("tern.panel.grid.major.show",theme=theme)
  ifthenelse(is.logical(showgrid),showgrid[1],getOption("tern.grid.major.show"))
}
.theme.get.showgrid.minor <- function(theme){
  showgrid   = calc_element("tern.panel.grid.minor.show",theme=theme)
  ifthenelse(is.logical(showgrid),showgrid[1],getOption("tern.grid.minor.show"))
}
.theme.get.outside       <- function(theme){
  outside     = calc_element("tern.axis.ticks.outside",theme=theme)
  ifthenelse(is.logical(outside),outside[1],getOption("tern.ticks.outside"))
}
.theme.get.showprimary   <- function(theme){
  showprimary = calc_element("tern.axis.ticks.primary.show",theme=theme)
  ifthenelse(is.logical(showprimary), showprimary[1],getOption("tern.ticks.primary.show"))
}
.theme.get.showsecondary <- function(theme){
  showsecondary = calc_element("tern.axis.ticks.secondary.show",theme=theme)
  ifthenelse(is.logical(showsecondary),showsecondary[1],getOption("tern.ticks.secondary.show"))
}
.theme.get.showarrows <- function(theme){
  showarrows   = calc_element('tern.axis.arrow.show',theme=theme)
  ifthenelse(is.logical(showarrows),showarrows[1],getOption("tern.arrow.show"))
}
.theme.get.showmask <- function(theme){
  showmask = calc_element('tern.panel.mask.show',theme=theme)
  ifthenelse(is.logical(showmask),showmask[1],getOption('tern.mask.show'))
} 

.theme.get.label <- function(self,n,d=n,suffix=''){
  x      = function(n,s) sprintf("%s%s",n,s)
  if(!is.character(n)) return('')
  labels = self$labels_coord
  ix     = unique(c(x(n,suffix),x(self$mapping[[n]],suffix),n,self$mapping[[n]])) 
  id     = which(ix %in% names(labels))
  if(length(id) == 0) d else labels[[ix[id[1]]]]
}
.theme.get.rotation <- function(self){
  tryCatch({
    angle = calc_element('tern.panel.rotate',self$theme)
    return(ifthenelse(is.finite(angle),angle,0)[1])
  },error=function(e){ 
    warning(e) 
  })
  return(0)
}

#----------------------------------------------------------------------------------
#Internals >>>> Render Components
# -Backgrounds
# -Borders
# -Grids
# -Precession Arrows & Markers
# _Apex Titles
#----------------------------------------------------------------------------------
.get.grid.data <- function(self,theme,data.extreme,X,major=TRUE,angle=0,angle.text=0){
  clockwise   = .theme.get.clockwise(theme)
  seq.tlr     = self$required_scales
  ix          = which(X == seq.tlr)
  existing    = data.frame()
  
  tryCatch({
    scale = self$scales[[X]]
    
    #Determine the limits
    limits = (if(inherits(scale$limits,'waiver')) NULL else scale$limits) %||% c(0,1)
    
    #Determine the Breaks
    breaks <- if(major) scale$breaks else scale$minor_breaks
    if(inherits(breaks,'waiver')){
      breaks = breaks_tern(limits = limits, isMajor = major)
    }
    
    #Bypass if necessary
    if(length(breaks) == 0) 
      return(existing)
    
    #Determine the labels
    labels <- if(major) (scale$labels %||% waiver()) else ""
    if(inherits(labels,'waiver'))
      labels = labels_tern(limits=limits,breaks=breaks)
    
    #major ticklength
    tl.major <- 0
    tryCatch({
      tl.major <- convertUnit(theme$tern.axis.ticks.length.major,"npc",valueOnly=T)
    },error=function(e){ warning(e) })
    
    #minor ticklength
    tl.minor <- 0
    tryCatch({
      tl.minor <- convertUnit(theme$tern.axis.ticks.length.minor,"npc",valueOnly=T)
    },error=function(e){  warning(e) })
    
    #Assign new id.
    id     <- (max(existing$ID,0) + 1)
    ix     <- min(ix,ifthenelse(major,length(tl.major),length(tl.minor)))
    majmin <- if(major) 'major' else 'minor' #Major or Minor Element Name part.
    
    #The new dataframe
    new            <- data.frame(ID = id,Scale = X, breaks, Labels = labels, Major = major)
    new            <- subset(new, breaks >= min(limits) & breaks <= max(limits))
    new$Prop       <- (new$breaks - min(limits)) / abs(diff(limits))
    new$TickLength <- ifthenelse(major,tl.major[ix],tl.minor[ix])
    new$NameText   <- sprintf("tern.axis.text.%s",X)
    new$NameTicks  <- sprintf("tern.axis.ticks.%s.%s",majmin,X)
    new$NameGrid   <- sprintf("tern.panel.grid.%s.%s",majmin,X)
    new$Major      <- major
    
    ##Start and finish positions of scale.
    out            <- c("x","y")
    
    #Start indexes.
    ix.s           <- which(seq.tlr == X);
    finish         <- as.numeric(data.extreme[ix.s,out])
    
    #For Ticks
    ix.f           <- ifthenelse(clockwise,if(ix.s == 3){1}else{ix.s+1},if(ix.s == 1){3}else{ix.s-1})
    start          <- as.numeric(data.extreme[ix.f,out])
    for(i in 1:length(out))
      new[,out[i]] <- new$Prop*(finish[i]-start[i]) + start[i]
    
    #FOR GRID
    ix.f           <- ifthenelse(clockwise,if(ix.s == 1){3}else{ix.s-1},if(ix.s == 3){1}else{ix.s+1})
    start          <- as.numeric(data.extreme[ix.f,out])
    for(i in 1:length(out))
      new[,paste0(out[i],"end.grid")] <- new$Prop*(finish[i]-start[i]) + start[i]
    
    #The tick angles.
    new$Angle      <- angle + .theme.get.rotation(self)
    new$Angle.Text <- .valid.angle(angle.text)
    
    #Determine the tick finish positions for segments.
    new$xend       <- cos(new$Angle*pi/180)*new$TickLength + new$x
    new$yend       <- sin(new$Angle*pi/180)*new$TickLength + new$y
    
    #Determine the secondary tick start and finish positions.
    new$x.sec      <- new$xend.grid
    new$y.sec      <- new$yend.grid
    new$xend.sec   <- cos((new$Angle+180)*pi/180)*new$TickLength + new$x.sec
    new$yend.sec   <- sin((new$Angle+180)*pi/180)*new$TickLength + new$y.sec
    
    return(new)
    
  },error=function(e){
    warning(e)
  })
  
  return(existing)
}

.render.fgset <- function(self,data.extreme,panel_params,theme,items){
  items = .render.ticks(      self,data.extreme,panel_params,theme,items)
  items = .render.border.main(self,data.extreme,              theme,items)
  items = .render.border.axis(self,data.extreme,              theme,items)
  items = .render.labels(     self,data.extreme,panel_params,theme,items)
  items = .render.arrows(     self,data.extreme,panel_params,theme,items)
  items
}

.render.background <- function(self,data.extreme,theme,items){
  tryCatch({
    e <- calc_element('tern.panel.background',theme=theme,verbose=F)
    if(inherits(e, "element_blank"))
      return(items)
    
    g <- polygonGrob( data.extreme$x, 
                      data.extreme$y, 
                      default.units = "npc",
                      id   = rep(1,nrow(data.extreme)),
                      gp   = gpar(  col  = e$colour %||% 'transparent',
                                    fill = alpha(e$fill %||% 'transparent',ifthenelse(!is.numeric(e$alpha),1,e$alpha)),
                                    lwd  = 0, #ifthenelse(!is.numeric(e$size),0,e$size)*find_global_tern(".pt"),
                                    lty  = e$linetype
                      )
    )
    items[[length(items) + 1]] <- g
  },error = function(e){ 
    warning(e) 
  })
  items
}

.render.border.main <- function(self,data.extreme,theme,items){
  tryCatch({
    #e      = calc_element('tern.panel.background',theme,verbose=F)
    e      = calc_element('panel.border',theme,verbose=F)
    if(inherits(e, "element_blank"))
      return(items)
    ex     = rbind(data.extreme,data.extreme[1,])
    grob   = ggint$element_grob.element_line(e,size=e$size,x=ex$x,y=ex$y)
    items[[length(items) + 1]] = grob
  },error=function(e){ 
    warning(e)
  })
  items
}

.render.border.axis <- function(self,data.extreme,theme,items,X=self$required_scales){
  
  #Only Unique Entries
  X = unique(X)
  
  #Checks
  seq.tlr = self$required_scales
  if(any(!{X %in% seq.tlr})) stop('Invalid X')
  
  grobs = function(name,s,f,items){
    tryCatch({
      grob = element_render(theme,name, x = data.extreme$x[c(s,f)], y = data.extreme$y[c(s,f)])
      items[[length(items) + 1]] <- grob
    },error=function(e){ 
      #warning(e) 
    })
    items
  }
  
  f = if(.theme.get.clockwise(theme)) c(2,3,1) else c(3,1,2)
  for(x in X){
    s      = which(seq.tlr == x)
    name   = sprintf("tern.axis.line.%s",seq.tlr[s])
    items  = grobs(name,s,f[s],items)
  }
  items
}

.render.ticks  <- function(self,data.extreme,panel_params,theme,items,X=self$required_scales){
  
  #Only Unique Entries
  X = unique(X)
  
  #Checks
  seq.tlr = self$required_scales
  if(any(!{X %in% seq.tlr})) stop('Invalid X')
  
  outside     = .theme.get.outside(theme)
  clockwise   = .theme.get.clockwise(theme)
  primary     = .theme.get.showprimary(theme)
  secondary   = .theme.get.showsecondary(theme)
  angle       = .get.angles(clockwise) + (!outside)*180
  angle.text  = .get.angles.ticklabels(clockwise) + .theme.get.rotation(self)
  
  #Function to render the grobs
  grobs  <- function(name,items,df,primary=TRUE){
    tryCatch({  
      append = function(x) paste0(x,if(primary) NULL else '.sec')
      ix.x   = append(c('x','xend'))
      ix.y   = append(c('y','yend'))
      
      e = calc_element(name,theme)
      if(inherits(e, "element_blank"))
        return(items)
      
      
      g = lapply(1:nrow(df),function(ix){ 
        segmentsGrob(x0 = df[ix,ix.x[1]], 
                     x1 = df[ix,ix.x[2]], 
                     y0 = df[ix,ix.y[1]], 
                     y1 = df[ix,ix.y[2]],
                     default.units ="npc",
                     arrow         = NULL,
                     gp            = gpar(col     = e$colour %||% 'transparent', 
                                          lty     = e$linetype,
                                          lineend = 'butt',
                                          lwd     = e$size %||% 0)
        )
      })
      
      items = c(items,g)
    },error = function(e){ 
      warning(e) 
    })
    items
  }
  
  #Iterate over the values of X
  for(x in X){
    
    #Determine the index
    ix = which(seq.tlr == x)
    
    #Generate the data
    df = ldply(c(T,F),function(major){
      .get.grid.data(self,theme,data.extreme,X = x, major = major, angle = angle[ix], angle.text = angle.text[ix])
    })
    
    #If Primary Ticks
    if(primary | secondary){
      for(name in unique(df$NameTicks)){
        df.sub = df[which(df$NameTicks == name),,drop = F]
        if(primary)   items = grobs(name = name, items = items, df = df.sub, primary = T)
        if(secondary) items = grobs(name = name, items = items, df = df.sub, primary = F)
      }
    }
  }
  
  #Done
  items
}

.render.labels <- function(self,data.extreme,panel_params,theme,items,X=self$required_scales){
  #Only Unique Entries
  X = unique(X)
  
  #Checks
  seq.tlr = self$required_scales
  if(any({!X %in% seq.tlr})) 
    stop('Invalid X')
  
  #Axis labels, ie labels next to ticks
  grobs <- function(name,items,df,outside,showprimary){ 
    tryCatch({
      df        = df[which(df$Labels != ''),]
      e         = calc_element(name,theme,verbose=FALSE)
      
      if(plyr::empty(df) || inherits(e, "element_blank"))
        return(items)
      
      latex     = calc_element('tern.plot.latex',theme)
      xts       = if(outside) df$x else df$xend
      xtf       = if(outside) df$xend else df$x
      yts       = if(outside) df$y else df$yend
      ytf       = if(outside) df$yend else df$y 
      angle     = is.numericor(e$angle,0) + is.numericor(unique(df$Angle.Text)[1],0)
      dA        = angle - atan2(ytf - yts, xtf - xts)*180/pi #DEGREES, Angle Difference between Ticks and Labels
      adj       = as.numeric(!outside)
      
      g = textGrob(
        label   = label_formatter(df$Labels,latex=latex),
        x       = ifthenelse(showprimary || !outside,xtf,xts) + convertX(cos(pi*(df$Angle/180 + adj))*unit(2,'pt'),'npc',valueOnly = T),
        y       = ifthenelse(showprimary || !outside,ytf,yts) + convertY(sin(pi*(df$Angle/180 + adj))*unit(2,'pt'),'npc',valueOnly = T),
        hjust   = +0.5*cos(pi*(dA/180 - 1)) + is.numericor(e$hjust,0), #BACK TO RADIANS
        vjust   = -0.5*sin(pi*(dA/180 - 1)) + is.numericor(e$vjust,0), #BACK TO RADIANS
        rot     = angle,
        default.units = 'npc',
        gp  = gpar(
          col        = e$colour %||% 'transparent',
          fontsize   = e$size %||% 0,
          fontfamily = e$family,
          fontface   = e$fontface,
          lineheight = e$lineheight
        )
      )
      
      items[[length(items) + 1]] <- g
    },error = function(e){ 
      warning(e) 
    })
    items
  }
  
  if(.theme.get.showlabels(theme)){
    outside     = .theme.get.outside(theme)
    showprimary = .theme.get.showprimary(theme)
    clockwise   = .theme.get.clockwise(theme)
    angle       = .get.angles(clockwise) + (!outside)*180
    angle.text  = .get.angles.ticklabels(clockwise) + .theme.get.rotation(self)
    #Generate the data
    df = ldply(X,function(x){
      ix =  which(seq.tlr == x)
      .get.grid.data(self,theme,data.extreme,X = x, major = TRUE, angle = angle[ix], angle.text = angle.text[ix])
    })
    for(name in unique(df$NameText))
      items = grobs(name=name,items=items,df=df[which(df$NameText  == name),,drop=F],outside,showprimary)
  }
  items
}

.render.grid <- function(self,data.extreme,panel_params,theme,items,X=self$required_scales){
  
  #Only Unique Entries
  X = unique(X)
  
  #Checks
  seq.tlr = self$required_scales
  if(any({!X %in% seq.tlr})) 
    stop('Invalid X')
  
  #Process the flags.
  clockwise     <- .theme.get.clockwise(theme)
  outside       <- .theme.get.outside(theme)
  showgrid.major<- .theme.get.showgrid.major(theme)
  showgrid.minor<- .theme.get.showgrid.minor(theme)
  
  #Get the Angles
  angle      <- .get.angles(clockwise) + (!outside)*180
  angle.text <- .get.angles.ticklabels(clockwise) + .theme.get.rotation(self)
  
  ##get the base data.
  df = expand.grid(ix=seq_along(seq.tlr),major=c(T,F))
  df = plyr::ddply(df,c('major','ix'),function(x){
    .get.grid.data(self,theme,data.extreme,X = seq.tlr[x$ix], major = x$major, angle = angle[x$ix], angle.text = angle.text[x$ix])
  })
  
  #If Nothing in 'd', return the curent list of items
  if(plyr::empty(df))
    return(items)
  
  grobs   <- function(name,items,df,showgrid.major=TRUE,showgrid.minor=TRUE){
    if((unique(df$Major) & showgrid.major) | (!unique(df$Major) & showgrid.minor)){
      tryCatch({
        e = calc_element(name,theme)
        g = lapply(1:nrow(df),function(ix){ 
          segmentsGrob(x0 = df$x[ix], 
                       x1 = df$xend.grid[ix], 
                       y0 = df$y[ix], 
                       y1 = df$yend.grid[ix],
                       default.units ="npc",
                       arrow         = NULL,
                       gp            = gpar(col     = e$colour %||% 'transparent', 
                                            lty     = e$linetype,
                                            lineend = 'butt',
                                            lwd     = (e$size %||% 0)*.pt)
          )
        })
        items = c(items,g)
      },error = function(e){ 
        warning(e) 
      })
    }
    items
  }
  
  #PROCESS TICKS AND LABELS
  if(showgrid.major | showgrid.minor){
    for(name in unique(df$NameGrid)){ 
      items <- grobs(name=name, items=items, df = df[ which(df$NameGrid  == name),,drop=FALSE], 
                     showgrid.major = showgrid.major, showgrid.minor = showgrid.minor)
    } 
  }
  
  items
}


.render.titles <- function(self,data.extreme,panel_params,theme,items){
  if(!.theme.get.showtitles(theme)) 
    return(items)
  
  #Determine the required scales
  seq.tlr = self$required_scales
  
  #Determine the length of the side
  sidelength = sqrt( diff(data.extreme$x[1:2])^2 + diff(data.extreme$y[1:2])^2)
  
  #Build the local function for building the grobs
  grobs = function(name,ix,items){
    tryCatch({
      
      e = calc_element(name,theme=theme,verbose=F)
      if(inherits(e, "element_blank"))
        return(items)
      
      ixc   <- c('x','y')
      latex <- calc_element('tern.plot.latex',theme)
      point <- as.numeric(data.extreme[ix,ixc])
      base  <- as.numeric(apply(data.extreme[-ix,ixc],2,mean))
      angle <- atan2((point[2]-base[2])*.ratio(),point[1]-base[1])
      n     <- regmatches(name,regexpr(".$",name)) 
      label <- c(self$scales[[n]]$name, self$labels_coord[[n]], self$labels_coord[[ self$mapping[[n]] ]], n)[1]
      grob  <- element_render(theme, name,
                              label = label_formatter(label,latex=latex),
                              x     = data.extreme$x[ix] + 0.05*sidelength*cos(angle),
                              y     = data.extreme$y[ix] + 0.05*sidelength*sin(angle),
                              hjust = e$hjust - 0.5*cos(angle),
                              vjust = e$vjust - 0.5*sin(angle))
      items[[length(items) + 1]] <- grob
    },error = function(e){
      #warning(e)
    })
    items
  }
  
  #process the axes
  for(ix in seq_along(seq.tlr))
    items = grobs(sprintf("tern.axis.title.%s",seq.tlr[ix]),ix,items)
  
  #Done, Return the list of items
  items
}

# The Arrows Parallel to the Axes
.render.arrows <- function(self,data.extreme,details,theme,items){
  
  if(!.theme.get.showarrows(theme))
    return(items)
  
  tryCatch({
    
    #clockwise or anticlockwise precession
    clockwise <- .theme.get.clockwise(theme)
    
    #The basic data.
    d.s <- data.extreme[ifthenelse(clockwise,c(2,3,1),c(3,1,2)),c('x','y')]
    d.f <- data.extreme[c(1,2,3),c('x','y')]
    rownames(d.s) <- rownames(d.f) #Correct rownames
    
    #Determine the length of the side
    sidelength = sqrt( diff(data.extreme$x[1:2])^2 + diff(data.extreme$y[1:2])^2)
    
    #arrow start and finish proportions
    arrowstart = calc_element('tern.axis.arrow.start', theme)
    arrowfinish= calc_element('tern.axis.arrow.finish',theme)
    
    #Ensure arrow start and finish length is 3.
    if(length(arrowstart) != 3 && length(arrowstart) >= 1)
      arrowstart <- rep(arrowstart[1],3)
    if(length(arrowfinish) != 3 && length(arrowfinish) >= 1)
      arrowfinish <- rep(arrowfinish[1],3)
    
    #Itterate over indexes 1:3
    for(i in c(1:3)){
      #Put in correct order.
      if(arrowfinish[i] < arrowstart[i]){
        warning(paste("Arrow size theme 'element tern.axis.arrow.finish[",i,"]' (",arrowfinish[i],") is < 'tern.axis.arrow.start[",i,"]' (",arrowstart[i],"), values will be swapped.",sep=""),call.=FALSE)
        #swapvalues
        tmp  = arrowstart[i]; arrowstart[i]  = arrowfinish[i]; arrowfinish[i] = tmp
      }
      #Check finish
      if(arrowfinish[i] > 1.0){
        warning(paste("Arrow size theme 'element tern.axis.arrow.finish[",i,"]' (",arrowfinish[i],") is > 1.0 and will be truncated",sep=""),call.=FALSE)
        arrowfinish[i] = 1.0
      }
      #Check start
      if(arrowstart[i] < 0.0){
        warning(paste("Arrow size theme 'element tern.axis.arrow.start[",i,"]' (",arrowstart[i],") is < 0.0 and will be truncated",sep=""),call.=FALSE)
        arrowstart[i] = 0.0
      }
    }
    
    #Cut down to relative proportion.
    dx   = (d.f - d.s)
    d.f  = d.f - (1 - arrowfinish)*dx
    d.s  = d.s +      arrowstart*dx
    d    = rbind(d.s,d.f)
    
    #Determine the start and end positions
    ixseq <- names(self$mapping)
    ixrow <- paste0("AT.",ixseq)
    ixcol <- c("x","y","xend","yend")
    ix    <- which(colnames(d) %in% ixcol[c(1:2)])
    d     <- cbind(d[1:3,ix],d[4:6,ix]);
    rownames(d) <- ixrow; colnames(d) <- ixcol
    
    #The arrow seperation in npc units.
    arrowsep      <- calc_element("tern.axis.arrow.sep",theme=theme,verbose=F)
    ticklength    <- max(calc_element("tern.axis.ticks.length.major",theme=theme,verbose=F),
                         calc_element("tern.axis.ticks.length.minor",theme=theme,verbose=F))
    
    #Ensure there are EXACTLY 3 values for each metric
    if(length(arrowsep)   != 3 && length(arrowsep)   >= 1){ arrowsep   = rep(arrowsep[1],3)   }
    if(length(ticklength) != 3 && length(ticklength) >= 1){ ticklength = rep(ticklength[1],3) }
    
    #Determine the Angles
    d[ixrow,"angle"]    <- .get.angles.arrows(clockwise)
    
    #get set of 3 arrowsep positions
    d[ixrow,"arrowsep"] <- arrowsep*sidelength
    
    #MOVE the Arrows Off the Axes.
    rotation = .theme.get.rotation(self)
    d[,ixcol[c(1,3)]]   <- d[,ixcol[c(1,3)]] + cos(pi*(d$angle + rotation)/180)*d$arrowsep #xcoordinates
    d[,ixcol[c(2,4)]]   <- d[,ixcol[c(2,4)]] + sin(pi*(d$angle + rotation)/180)*d$arrowsep #ycoorinates
    
    #Centerpoints, labels, arrowsuffix
    d$xmn = rowMeans(d[,ixcol[c(1,3)]])
    d$ymn = rowMeans(d[,ixcol[c(2,4)]])
    d$L   = unlist(lapply(ixseq,function(n){ .theme.get.label(self,n)  }))
    d$LA  = unlist(lapply(ixseq,function(n){ 
      #c(self$labels_coord[[ sprintf('%sarrow',n) ]],)[1]
      .theme.get.label(self,n,suffix='arrow')
    }))
    d$W   = unlist(lapply('W',function(n){  .theme.get.label(self,n,'')  }))
    d$A   = .get.angles.arrowmarker(clockwise)
    d$AL  = .valid.angle(d$A + .theme.get.rotation(self))
    
    ##Function to create new arrow grob
    .render.arrow <- function(name,ix,items){
      tryCatch({  
        e = calc_element(name,theme=theme,verbose=F)
        if(inherits(e, "element_blank"))
          return(items)
        
        g = segmentsGrob(x0 = d$x[ix], x1 = d$xend[ix], y0 = d$y[ix], y1 = d$yend[ix],
                         default.units ="npc",
                         arrow         = e$lineend,
                         gp            = gpar(col     = e$colour %||% 'transparent', 
                                              lty     = e$linetype,
                                              lineend = 'butt',
                                              lwd     = e$size)
        )
        items[[length(items) + 1]] <- g
      },error = function(e){
        warning(e)
      })
      items
    }
    
    #Function to greate new label grob
    .render.label <- function(name,ix,items){
      tryCatch({  
        e    = calc_element(name,theme=theme,verbose=F)
        if(inherits(e, "element_blank"))
          return(items)
        
        latex = calc_element('tern.plot.latex',theme)
        dA    = e$angle + d$AL - 180*(atan2((d$yend - d$y)*.ratio(), d$xend - d$x)/pi + as.numeric(clockwise))
        g     = textGrob( label = arrow_label_formatter(d$LA[ix],d$W[ix],latex=latex), 
                          x     = d$xmn[ix] + convertX(cos(pi*(d$angle[ix] + rotation)/180)*unit(2,'pt'),'npc',valueOnly = T), 
                          y     = d$ymn[ix] + convertY(sin(pi*(d$angle[ix] + rotation)/180)*unit(2,'pt'),'npc',valueOnly = T), 
                          hjust = e$hjust + 0.5*sin(dA[ix]*pi/180), 
                          vjust = e$vjust + 0.5*cos(dA[ix]*pi/180),
                          rot   = d$AL[ix] + e$angle, 
                          default.units="npc", 
                          gp   = gpar(col        = e$colour %||% 'transparent', 
                                      fontsize   = e$size,
                                      fontfamily = e$family, 
                                      fontface   = e$face, 
                                      lineheight = e$lineheight))
        items[[length(items) + 1]] <- g
      },error = function(e){ 
        warning(e) 
      })
      items
    }
    
    #process the axes
    for(i in 1:length(ixseq)){
      items <- .render.arrow(paste0("tern.axis.arrow.",     ixseq[i]),i,items) #Arrows
    }
    for(i in 1:length(ixseq)){
      items <- .render.label(paste0("tern.axis.arrow.text.",ixseq[i]),i,items) #Markers
    }
    
  },error=function(e){
    message(e)
  })
  items
}


.rotation = function (xy, angle, degrees=TRUE) {
  if(degrees) angle = pi*angle/180
  xy <- as.matrix(xy)
  ca <- cos(angle)
  sa <- sin(angle)
  xy.rot <- xy %*% t(matrix(c( ca, sa, 
                               -sa, ca), 2, 2))
  return(xy.rot)
}

.valid.angle = function(x){
  if(length(x) > 1)
    return(sapply(x,.valid.angle))
  if(x >  90) x = .valid.angle(x - 180)
  if(x <=-90) x = .valid.angle(x + 180)
  x
}

"%||%"  <- function(a, b) {if (!is.null(a)) a else b}

# Ternary / Cartesian Transformation
# 
# Functions to transform data from the ternary to cartesian spaces and vice-versa. 
# 
# @param data \code{data.frame} containing columns as required by the coordinate system. 
# Data will be scaled so that the rows sum to unity, in the event that the user has provided 
# data that does not.
# @param coord Coordinate system object, inheriting the \code{\link{CoordTern}} class, error will
# be thrown if a different coordinate system is sent to this method
# @param ... not used
# @param inverse logical if we are doing a forward (FALSE) or reverse (TRUE) transformation
# @param scale logical as to whether the transformed coordinates are scaled (or reverse scaled in the case of inverse 
# transformation) according to the training routine defined in the coordinate system.
# @param drop drop all non columns which are not involved in the transformation
# @author Nicholas Hamilton
# @examples 
# data(Feldspar)
# dfm = plyr::rename(Feldspar,c("Ab"="x","An"="y","Or"="z"))
# crd = coord_tern()
# fwd = tlr2xy(dfm,crd)
# rev = tlr2xy(fwd,crd,inverse = TRUE)
# @rdname ternary_transformation
# @name ternary_transformation
NULL

# @details \code{tlr2xy} transforms from the ternary to cartesian spaces, an inverse transformation 
# transforms between cartesian to ternary spaces
# @rdname ternary_transformation
# @export
tlr2xy <- function(data,coord,...,inverse=FALSE,scale=TRUE,drop=FALSE){
  
  #Run Some Checks
  if(!inherits(coord,"CoordTern")) 
    stop("argument 'coord' must be a CoordTern coordinate structure")
  if(class(data) != "data.frame")
    stop("argument 'data' must be of type 'data.frame'")
  if(!is.logical(inverse) | !is.logical(scale))
    stop("argument 'inverse' and 'scale' (both) must be logical")
  
  #Determine the Proposed Mapping and Reverse Mapping
  mapping = coord$mapping; ix.trl  = as.character(coord$mapping); ix.xy   = c("x","y")
  if(length(unique(ix.trl)) != 3 | length(unique(names(mapping))) != 3) stop('Mapping must have 3 unique named variable pairs',call.=FALSE)
  mappingRev = mapping; for( key in names(mappingRev) ) mappingRev[key] <- key; names(mappingRev) = as.character(mapping)
  
  #Determine the required and the destination position aesthetics
  ix.req    = if(!inverse){ ix.trl }else{ ix.xy  }
  ix.dst    = if(!inverse){ ix.xy  }else{ ix.trl }
  
  #Determine the aesthetic set groups
  setGroups = .get.sets(ix.req,names(data))
  
  #If there are more than one set group, process recursively, and re-combine
  if(length(setGroups) > 1){
    #Get the full set of combinations
    ix.full      = .combos(ix.req,setGroups)
    #Check the missing aesthetics
    .check.aes(coord,ix.full,names(data),inverse=inverse)
    #Split the data into transformation and non-transformation data
    data.notrans = data[,which(!names(data) %in% ix.full),drop=FALSE]
    data.totrans = data[,which( names(data) %in% ix.full),drop=FALSE]
    #Now for each group, conduct a transformation, when re-transformed, combine and return
    #Transformation is called recursively
    for(group in setGroups){
      ix  = .combos(ix.req,group)
      df  = data.totrans[,ix]
      re  = sprintf("^(%s)(%s)",paste(ix.req,collapse="|"),group)
      names(df) = gsub(re,"\\1",names(df))
      ##-------------------------------------------------------------
      df  = tlr2xy(df,coord,...,inverse=inverse,scale=scale)[,ix.dst]
      ##-------------------------------------------------------------
      names(df) = paste(ix.dst,group,sep="")
      if(!plyr::empty(df)) data.notrans = if(!plyr::empty(data.notrans)){ cbind(data.notrans,df) }else{ df }
    }
    return(data.notrans)
  }
  
  #Local function to adjust the range, depending on if inverse or not
  adjustRange <- function(input,lim,inv=inverse){
    if(is.null(lim)) lim=c(0,1)
    if( !diff(lim) ) lim=c(0,1)
    adl   = abs(diff(lim))
    ml    = min(lim)
    input = if(inv[1]){ input*adl + ml }else{ (input-ml)/adl }
    input
  }
  
  #Local function to scale ternary coordinates
  scaleCoordinatesToUnity <- function(input){
    s  <- rowSums(input[,ix.trl]);
    ix <- which(!as.logical(s))
    if(length(ix) > 0){ input[ix,ix.trl] <- 1/3; s[ix]  <- 1.0 } #Prevent Div By 0 error
    for(x in ix.trl){ input[,x] = input[,x]/s }
    input
  }
  
  #Forward transformation is ternary to cartesian
  if(!inverse[1]){
    #Check the missing aesthetics
    .check.aes(coord,ix.trl,names(data),inverse=FALSE)
    #If scale to composition sum of 1
    if(scale[1]){ data = scaleCoordinatesToUnity(data) }
    #Adjust for the Limits.
    for(x in mapping){ 
      data[,x] = adjustRange(data[,x],coord$limits[[ mappingRev[[x]] ]]) 
    }
    #Calculate
    data$y = data[,as.character(mapping['T'])]*tan(pi/3)*0.5 
    data$x = data[,as.character(mapping['R'])]+data$y*tan(pi/6)
    data = data[,-which(names(data) == 'z')]
    
    #Inverse transformation is cartesian to ternary
  }else{
    #Check the missing aesthetics
    .check.aes(coord,ix.xy,names(data),inverse=TRUE)
    #Calculate
    out.R = data[,ix.xy[1]] - data[,ix.xy[2]]*tan(pi/6)
    out.T = data[,ix.xy[2]]/(tan(pi/3)*0.5)
    out.L = 1.0 - out.R - out.T
    #Adjust for the Limits
    for(x in mappingRev){ 
      data[, mapping[[x]] ] = adjustRange( get(sprintf('out.%s',x)) ,coord$limits[[ x ]]) 
    }
  }
  
  #Done
  data
}

# @details \code{xy2tlr} transforms from the cartesian to ternary spaces, an inverse transformation 
# transforms between ternary to cartesian spaces, it is the reciprocal to \code{\link{tlr2xy}}, therefore
# an inverse transformation in \code{\link{xy2tlr}} function is the same as the forward 
# transformation in \code{\link{tlr2xy}} 
# @rdname ternary_transformation
# @export
xy2tlr <- function(data,coord,...,inverse=FALSE,scale=TRUE) 
  tlr2xy(data,coord,...,inverse=!inverse,scale=scale) 


#internal
#get the set groups
.get.sets <- function(vars,cols){
  re    = sprintf("^(%s)(.*)",paste(vars,collapse="|"))
  match = grepl(re,cols,perl=TRUE)
  ix    = which(match); if(length(ix) == 0) return(NULL)
  unique(gsub(re,"\\2",cols[ix]))
}

#Get combinations
.combos = function(a,b,collapse=""){ 
  apply(expand.grid(a,b), 1, paste, collapse=collapse) 
}

#Check the aesthetics
.check.aes = function(coord,ix,colNames,inverse=FALSE){
  missing = unique(setdiff(ix,colNames))
  if(length(missing) > 0){
    dir = c('tlr','xy'); if(inverse) dir = rev(dir)
    msg = sprintf("ggtern: %s requires the following missing aesthetics (%s) : %s", 
                  class(coord)[1], paste(dir,collapse="->"), paste(missing,collapse="', '"))
    stop(msg,call.=FALSE)
  }
}

#Expose some required functions from the parent ggplot2 namespace
.getFunctions <- function(){
  
  # OLD FUNCTIONS  
  #new_panel','train_layout','train_position','train_ranges','map_position','map_layout','reset_scales','facet_render',
  #xlabel','ylabel'
  
  .functions.ggplot2   = c('create_layout',
                           #expand_default', ## REMOVED
                           'plot_theme',
                           'element_render',
                           # 'message_wrap', ## Luis: removed. see function below
                           'set_last_plot','make_labels','build_guides','is.zero','add_ggplot','labelGrob',
                           'is.layer','is.facet','is.Coord','GeomSegment',
                           '.element_tree',
                           # 'el_def', ## NOW EXPORTED
                           'expand_limits_scale', ## NEW
                           'view_scale_primary', ## NEW
                           'view_scale_secondary', ## NEW
                           'combine_elements','aes_to_scale',
                           'is.Coord','is.facet','is.layer','make_labels','update_labels','update_guides',
                           # 'update_theme', ## REMOVED
                           'aes_to_scale',
                           'scales_add_missing','scales_list','scales_transform_df','scales_map_df','scales_train_df',
                           'predictdf',
                           # 'contour_lines', ## REMOVED
                           'check_required_aesthetics','snake_class',
                           'ggname','ggplot_gtable','camelize',
                           'element_grob.element_line','element_grob.element_rect','element_grob.element_text','element_grob.element_blank',
                           'plot_clone','compute_just','labelGrob',
                           'hexGrob',
                           # 'try_require', ## REMOVED
                           'hex_binwidth','hexBinSummarise',
                           'find_args','is.margin','justify_grobs')
  .functions.gridExtra  = c('latticeGrob')
  .functions          = rbind(data.frame(p='ggplot2',  f=unique(.functions.ggplot2)),
                              data.frame(p='gridExtra',f=unique(.functions.gridExtra)))
  
  structure(
    mapply(function(f,p){ getFromNamespace(f,p) },as.character(.functions$f), as.character(.functions$p)),
    class=c("internal")
  )
}

message_wrap <- function (...) {
  msg <- paste(..., collapse = "", sep = "")
  wrapped <- strwrap(msg, width = getOption("width") -
                       2)
  message(paste0(wrapped, collapse = "\n"))
}


ggint <- .getFunctions()

#Internal Functions
#
#@description INTERNAL FUNCTIONS: \code{ggtern} makes use of several non-exported internal functions, list are as follows:
#@keywords internal
#@rdname undocumented
#@name zzz-internal
NULL

# \code{ifthenelse} function takes input arguments \code{x}, \code{a} and \code{b} and returns \code{a} if \code{x} is \code{TRUE}, else, returns \code{b}
# @param x logical input to check
# @param a value to return if \code{x} is TRUE
# @param b value to return if \code{x} is FALSE
# @keywords internal
# @rdname undocumented
ifthenelse <- function(x,a,b){
  if(!is.logical(x))stop("x argument must be logical")
  if(x){a}else{b}
}

# \code{is.numericor} function takes input arguments \code{A} and \code{B} and returns \code{A} if \code{A} is numeric, else, returns \code{B}
# @param A value to return if numeric
# @param B numeric value to return if \code{A} is NOT numeric
# @keywords internal
# @rdname undocumented
is.numericor <- function(A,B){
  if(!is.numeric(B)){stop("b must be numeric")}
  if(is.numeric(A)){A}else{B}
}

"%||%"  <- function(a, b) {if (!is.null(a)) a else b}

# \code{find_global_tern} is a function that conducts a named search for the \code{name} object instance, within the \code{env} environment. 
# If an instance doesn't exist within the \code{env} environment, a search is then conducted within the \code{ggtern} and \code{ggplot2} 
# namespaces \emph{(in that order)}. This is a modified version of the original source as provided in \code{ggplot2}, which has the same functionality, however, the modification is such that the function
# now additionally searches within the \code{ggtern} namespace prior to the \code{ggplot2} namespace.
# @param name character name of object to search for
# @param env environment to search within as first priority
# @param mode the mode to search within
# @keywords internal
# @rdname undocumented
# @export
find_global_tern <- function (name, env=environment(),mode='any'){  
  if(!is.character(name)){stop("'name' must be provided as a character")}
  if(!inherits(environment(),"environment")){stop("'env' must inherit the environment class")}
  
  if (exists(name, envir = env, mode = mode)){ 
    return(get(name, envir = env, mode = mode))
  }
  
  nsenv <- asNamespace("ggtern")
  if(exists(name, envir = nsenv, mode=mode)){
    return(get(name, envir = nsenv, mode = mode))
  }
  
  nsenv <- asNamespace("ggplot2")
  if(exists(name, envir = nsenv, mode=mode)){
    return(get(name, envir = nsenv, mode = mode))
  }
  
  NULL
}

# Convert RGB to HEX Color
# 
# Function to convert rgb color to hex color
# @param r,g,b colors, numeric scalar between 0 and 255
# @keywords internal
# @author Nicholas Hamilton
# @examples 
# #Black
# rgb2hex(0,0,0)
# 
# #White
# rgb2hex(255,255,255)
# 
# #Red
# rgb2hex(255,0,0)
# 
# #Green
# rgb2hex(0,255,0) 
# 
# #Blue
# rgb2hex(0,0,255)
# 
# #Vectorised sequence of blue
# rgb2hex(0,0,seq(0,255,by=5))
# @export
rgb2hex = function(r = 0, g = 0, b = 0){
  df = data.frame(r, g, b)
  check = function(x, ix = NULL){
    nm = deparse(substitute(x))
    ix = as.character({ix %||% ''})
    if(!is.numeric(x))  stop(sprintf("'%s%s' must be numeric",             nm,ix),call. = FALSE)
    if(length(x) != 1)  stop(sprintf("'%s%s' must be scalar",              nm,ix),call. = FALSE)
    if(!is.finite(x))   stop(sprintf("'%s%s' must be finite",              nm,ix),call. = FALSE)
    if(x < 0 | x > 255) stop(sprintf("'%s%s' must be in the range [0,255]",nm,ix),call. = FALSE)
  }
  nr = nrow(df)
  sapply( c(1:nr), function(ix){
    n = if(nr > 1){ ix }else{ NULL }
    r = df$r[ix]; check(r,n)
    g = df$g[ix]; check(g,n)
    b = df$b[ix]; check(b,n)
    sprintf("#%.2x%.2x%.2x",r,g,b) 
  })
}

# Generate Axis Breaks
# 
# Calculates the Breaks for Major or Minor Gridlines based on the input limits.
# @param limits the scale limits
# @param isMajor major or minor grids
# @param n number of breaks
# @rdname breaks_tern
# @examples 
#  breaks_tern()
#  breaks_tern(limits = c(0,.5),FALSE,10)
# @export
breaks_tern <- function(limits = c(0,1), isMajor = TRUE, n = 5){
  if(is.null(limits) || !all(is.numeric(limits)))
    limits = c(0,1)
  
  if(diff(range(limits)) == 0){
    ret = if(isMajor) getOption("tern.breaks.default") else getOption("tern.breaks.default.minor")
    return(ret)
  }
  
  ret = pretty(limits,n = n)
  if(!isMajor){
    r = range(ret)
    d = diff(r)/(length(ret)-1)
    minor = seq(min(ret)-d/2,max(ret)+d/2,by = d)
    minor = minor[which(minor > min(limits) & minor < max(limits))]
    ret   = minor[which(!minor %in% ret)]
  }
  ret
}


# @rdname breaks_tern
# @name breaks_tern
# @usage NULL
# @format NULL
# @export
getBreaks = function(limits = c(0,1), isMajor = TRUE, n = 5){
  tern_dep("2.1.4","'getBreaks' has been superceded by the 'breaks' function")
  breaks_tern(limits,isMajor,n)
}

# Generate Axis Labels
# 
# Calculates the Labels for Major or Minor Gridlines based on the input limits.
# @param breaks numeric denoting the breaks to produce corresponding labels
# @inheritParams breaks_tern
# @param format the formatting string to be passed through to the \code{\link{sprintf}} function
# @param factor the multiplicative factor
# @examples 
#labels_tern()
# labels_tern(limits = c(0,.5))
# @author Nicholas Hamilton
# @rdname labels_tern
# @export
labels_tern = function(limits = c(0,1), breaks = breaks_tern(limits), format = "%g", factor = 100){
  if(!is.numeric(breaks)) 
    stop("'breaks' must be numeric",call.=FALSE)
  
  #Default Result
  result = factor[1]*breaks
  
  #Try and process...
  tryCatch({
    if(!is.numeric(factor)) 
      stop("'factor' must be numeric",call.=FALSE)
    result = sprintf(format,factor[1]*breaks)
    
    #Stop First Label interfering with the main label
    if(breaks[1] == min(limits))
      result[1] = ''
    
  },error=function(e){ })
  
  #Done
  result
}

# @rdname labels_tern
# @name labels_tern
# @usage NULL
# @format NULL
# @export
getLabels = function(limits = c(0,1), breaks = breaks_tern(limits), format = "%g", factor = 100){
  tern_dep("2.1.4","'getBreaks' has been superceded by the 'breaks' function")
  labels_tern(limits,breaks,format,factor)
}

# \code{tern_dep} is a function that gives a deprecation error, warning, or messsage, 
# depending on version number, it is based of the \code{\link[ggplot2]{gg_dep}} function which is
# used inside the \code{ggplot2} package
# @inheritParams ggplot2::gg_dep
# @keywords internal
# @rdname undocumented
tern_dep <- function(version, msg) {
  v <- as.package_version(version)
  cv <- packageVersion("ggtern")
  
  # If current major number is greater than last-good major number, or if
  #  current minor number is more than 1 greater than last-good minor number,
  #  give error.
  if (cv[[1,1]] > v[[1,1]]  ||  cv[[1,2]] > v[[1,2]] + 1) {
    stop(msg, " (Defunct; last used in version ", version, ")",
         call. = FALSE)
    
    # If minor number differs by one, give warning
  } else if (cv[[1,2]] > v[[1,2]]) {
    warning(msg, " (Deprecated; last used in version ", version, ")",
            call. = FALSE)
    
    # If only subminor number is greater, give message
  } else if (cv[[1,3]] > v[[1,3]]) {
    message(msg, " (Deprecated; last used in version ", version, ")")
  }
  
  invisible()
}

#internal
.makeValid <- function(x){
  x = x[[1]]
  if(class(x) == 'character'){
    x = gsub("%","'%'",x)
    x = gsub('([[:punct:]])\\1+', '\\1', x)
    x = gsub(" ","~",x)
  }
  x
}

# \code{arrow_label_formatter} is a function that formats the labels directly adjacent to the ternary arrows.
# @param label character label
# @param suffix chacater suffix behind each label
# @param sep the seperator between label and suffix 
# @param ... additional arguments
# @param latex logical as to whether latex formats should be parsed
# @keywords internal
# @rdname undocumented
arrow_label_formatter             = function(label,suffix=NULL,sep="/",...) UseMethod("arrow_label_formatter")
arrow_label_formatter.default     = function(label,suffix=NULL,sep="/",...) arrow_label_formatter.character( as.character(label), suffix, sep, ...)
arrow_label_formatter.call        = function(label,suffix=NULL,sep="/",...) arrow_label_formatter.expression(as.expression(label),suffix, sep, ...)    
arrow_label_formatter.expression  = function(label,suffix=NULL,sep="/",...){
  suffix = if(suffix  == "")   NULL else suffix
  sep    = if(is.null(suffix)) ""   else .trimAndPad(sep)
  parse(text=paste(as.character(label),suffix,sep))
}
arrow_label_formatter.character   = function(label,suffix=NULL,sep="/",latex = FALSE,...) {
  suffix = if(suffix  == "")   NULL else suffix
  sep    = if(is.null(suffix)) ""   else .trimAndPad(sep)
  result = paste(label,suffix,sep=sep)
  if(latex[1]) result = TeX(result)
  result
}
.trimAndPad <- function(x){
  x = gsub("^(\\s+)","",gsub("(\\s+)$","",x))
  if(nchar(x) == 1) x = sprintf(" %s ",x)
  x
}


# \code{label_formatter} is a function that formats / parses labels for use in the grid.
# @param label character label
# @param ... additional arguments
label_formatter = function(label,...){ arrow_label_formatter(label,suffix="",sep="",...) }


# \code{joinCharacterSeries} is a function will turn a character vector 
# from the format \code{c('a','b','c')} to a single string
# in the following format: \code{"'a','b' and 'c'"}
# @param x character vector
# @author Nicholas Hamilton
# @keywords internal
# @rdname undocumented
joinCharacterSeries <- function(x,lastWord='and'){
  if(!is.character(x) | !is.vector(x)) stop("'x' must be character vector",call.=FALSE)
  if(length(x) > 1){ x = paste(paste(x[-length(x)],collapse="', '"),x[length(x)],sep=sprintf("' %s '",lastWord)) }
  sprintf("'%s'",x)
}


# \code{identityInv} is a function which returns exactly the same as \code{\link{identity}} however
# it can be used within transformation logic via \code{do.call(...)} in the same way as for example
# \code{\link{ilrInv}} is to \code{\link{ilr}}.
# @param x input object
# @author Nicholas Hamilton
# @keywords internal
# @rdname undocumented
identityInv = function(z) identity(z)


# \code{getFormulaVars} is a function that returns a list of either dependent or independent variables used
# in an input formula
# @param x formula object
# @param dependent whether to return the dependent variables (TRUE) or the indpenedent variables (FALSE)
# @rdname undocumented
# @keywords internal
# @author Nicholas Hamilton
getFormulaVars = function(x,dependent=TRUE) {
  if(class(x) != 'formula') stop("x argument must be a formula",call.=FALSE)
  all.vars(x[[if(dependent) 3 else 2]])
}

# Function to add missing scales and other items to the plot and its coordinates sytem
# @param ggplot object
# @rdname undocumented
# @keywords internal
# @author Nicholas Hamilton
scales_add_missing_tern <- function(plot){
  
  #Run some checks
  stopifnot(inherits(plot,'ggplot'))
  stopifnot(inherits(plot$coordinates,'CoordTern'))
  
  #Ensure required scales have been added
  rs = plot$coordinates$required_scales
  
  aesthetics  = setdiff(rs, plot$scales$input())
  env = plot$plot_env
  for (aes in aesthetics) {
    scale_name <- paste("scale", aes, "continuous", sep = "_")
    scale_f <- find_global_tern(scale_name, env, mode = "function")
    plot$scales$add(scale_f())
  }
  
  #ggint$scales_add_missing(plot,rs,plot$plot_env) ##NH
  #plot$scales$scales = plot$scales$scales[!sapply(plot$scales$scales,is.null)] 
  #plot$scales$scales = compact(plot$scales$scales)
  
  #Push some details to the coordinates
  plot$coordinates$scales        = sapply(rs,plot$scales$get_scales) ##NH
  for(r in rs) 
    plot$coordinates$limits[[r]] = plot$scales$get_scales(r)$limits
  plot$coordinates$labels_coord  = plot$labels
  plot$coordinates$theme         = ggint$plot_theme(plot) #NH
  
  #done
  plot
}

# Function to add clipping mask if it isn't already present
# @param plot ggplot object
# @rdname undocumented
# @keywords internal
# @author Nicholas Hamilton
layers_add_or_remove_mask = function(plot){
  theme = ggint$plot_theme(plot) #NH
  mask  = calc_element('tern.panel.mask.show',theme)[1] %||% TRUE
  if(is.na(mask) || mask){
    if(!"GeomMask" %in% unlist(lapply(plot$layers,function(x){ class(x$geom) })))
      plot = plot + geom_mask()
  }else{
    plot$layers = plyr::compact(lapply(plot$layers,function(x){
      if(inherits(x$geom,'GeomMask')) return(NULL) else x
    }))
  }
  plot
} 

# @rdname ggtern_themes
# @export
theme_void <- function(base_size = 12, base_family = "") {
  theme_ggtern(base_size,base_family) %+replace%
    ggplot2::theme_void(base_size,base_family) %+replace%
    theme(
      #text                 = element_blank(),
      #line                 = element_blank(),
      #rect                 = element_blank(),
      tern.axis.text       = element_blank(),
      tern.axis.title      = element_blank()
    )
}

# @inheritParams gridExtra::arrangeGrob
# @inheritParams gridExtra::grid.arrange
# @rdname arrangeGrob 
# @aliases grid.arrange
# @export
grid.arrange = function (..., newpage = TRUE) {
  if (newpage) 
    grid.newpage()
  g <- arrangeGrob(...)
  grid.draw(g)
  invisible(g)
}

# Show or Hide Axis Ticklabels
# 
# Convenience functions to enable or disable the axis ticklabels
# 
# \code{theme_showlabels} is a function that apends to the current theme a flag to switch ON the axis ticklabels, whilst 
# \code{theme_hidelabels} or \code{theme_nolabels} (Alias) are functions that apends to the current theme a flag 
# to switch OFF the axis ticklabels
# @author Nicholas Hamilton
# @rdname theme_showlabels
# @name theme_showlabels
NULL

# @rdname theme_showlabels
# @export
theme_showlabels <- function(){.theme_showlabels(TRUE)}

# @rdname theme_showlabels
# @export
theme_hidelabels <- function(){.theme_showlabels(FALSE)}

# @rdname theme_showlabels
# @export
theme_nolabels   <- theme_hidelabels

#Internal Function
.theme_showlabels <- function(show){
  theme(tern.axis.text.show=show)
}

# @rdname ggtern_themes
# @export
theme_ggtern <- function(base_size = 11, base_family = ""){
  #Base ggplot2 theme
  baseTheme = get("theme_gray",asNamespace("ggplot2"))
  
  #Create Instance of the base theme
  base = baseTheme(base_size = base_size, base_family = base_family)
  
  #Start with the base theme 
  base %+replace%
    
    #Start with additional elements
    theme(
      
      ##TERNARY PANEL
      tern.panel.background          = element_rect(),  #Panel is the triangular region
      tern.plot.background           = element_rect(), #Plot  is the rectangular outer region
      tern.plot.latex                = getOption('tern.latex'), #Parse Labels as Latex
      
      ##AXIS ARROWS
      #tern.axis                     = element_line(),
      tern.axis.hshift               = getOption("tern.hshift"),
      tern.axis.vshift               = getOption("tern.vshift"),
      tern.axis.clockwise            = getOption("tern.clockwise"),
      
      tern.axis.line                 = element_line(),
      tern.axis.line.T               = element_line(),
      tern.axis.line.L               = element_line(),
      tern.axis.line.R               = element_line(),
      tern.axis.line.ontop           = getOption("tern.line.ontop"),
      
      #Axis Titles
      tern.axis.title                = element_text(),
      tern.axis.title.T              = element_text(),
      tern.axis.title.L              = element_text(),
      tern.axis.title.R              = element_text(),
      tern.axis.title.show           = getOption("tern.title.show"), 
      
      #Axis Text
      tern.axis.text                 = element_text(),
      tern.axis.text.T               = element_text(),
      tern.axis.text.L               = element_text(),
      tern.axis.text.R               = element_text(),
      tern.axis.text.show            = getOption("tern.text.show"),
      
      #Arrow
      tern.axis.arrow                = element_line(lineend = getOption('tern.arrow')),
      tern.axis.arrow.T              = element_line(),
      tern.axis.arrow.L              = element_line(),
      tern.axis.arrow.R              = element_line(),
      tern.axis.arrow.text           = element_text(),
      tern.axis.arrow.text.T         = element_text(),
      tern.axis.arrow.text.L         = element_text(),
      tern.axis.arrow.text.R         = element_text(),
      tern.axis.arrow.sep            = getOption("tern.arrow.sep"),
      tern.axis.arrow.show           = getOption("tern.arrow.show"),
      tern.axis.arrow.start          = getOption("tern.arrow.start"),
      tern.axis.arrow.finish         = getOption("tern.arrow.finish"),
      
      #Ticks
      tern.axis.ticks                = element_line(),
      tern.axis.ticks.major          = element_line(),
      tern.axis.ticks.major.T        = element_line(),
      tern.axis.ticks.major.L        = element_line(),
      tern.axis.ticks.major.R        = element_line(),
      tern.axis.ticks.length.major   = 1.0*base$axis.ticks.length,
      tern.axis.ticks.length.minor   = 0.5*base$axis.ticks.length,
      tern.axis.ticks.outside        = getOption("tern.ticks.outside"),
      tern.axis.ticks.primary.show   = getOption("tern.ticks.primary.show"),
      tern.axis.ticks.secondary.show = getOption("tern.ticks.secondary.show"),
      tern.axis.ticks.minor          = element_line(),
      tern.axis.ticks.minor.T        = element_line(),
      tern.axis.ticks.minor.L        = element_line(),
      tern.axis.ticks.minor.R        = element_line(),
      
      #Panel Grids
      #tern.panel.grid                = element_line(),
      tern.panel.grid.major          = element_line(),
      tern.panel.grid.major.T        = element_line(),
      tern.panel.grid.major.L        = element_line(),
      tern.panel.grid.major.R        = element_line(),
      tern.panel.grid.major.show     = getOption("tern.grid.major.show"),
      tern.panel.grid.minor          = element_line(),
      tern.panel.grid.minor.T        = element_line(),
      tern.panel.grid.minor.L        = element_line(),
      tern.panel.grid.minor.R        = element_line(),
      tern.panel.grid.minor.show     = getOption("tern.grid.minor.show"),
      tern.panel.grid.ontop          = getOption("tern.grid.ontop"),
      tern.panel.mask.show           = getOption("tern.mask.show"),
      tern.panel.expand              = getOption('tern.expand'),
      tern.panel.rotate              = getOption('tern.rotate')
    )
}

# Arrange multiple grobs on a page (ggtern version)
# 
# A very slight modification to the original function, removing the explicit direction to use the ggplotGrob function
# from the ggplot2 namespace
# @inheritParams gridExtra::arrangeGrob
# @author Nicholas Hamilton
# @rdname arrangeGrob
# @export
arrangeGrob = function (..., grobs = list(...), layout_matrix, vp = NULL, name = "arrange", 
                        as.table = TRUE, respect = FALSE, clip = "off", nrow = NULL, 
                        ncol = NULL, widths = NULL, heights = NULL, top = NULL, bottom = NULL, 
                        left = NULL, right = NULL, padding = unit(0.5, "line")) 
{
  n <- length(grobs)
  if (!is.null(ncol) && !is.null(widths)) {
    stopifnot(length(widths) == ncol)
  }
  if (!is.null(nrow) && !is.null(heights)) {
    stopifnot(length(heights) == nrow)
  }
  if (is.null(ncol) && !is.null(widths)) {
    ncol <- length(widths)
  }
  if (is.null(nrow) && !is.null(heights)) {
    nrow <- length(heights)
  }
  if (is.null(nrow) && !is.null(ncol)) {
    nrow <- ceiling(n/ncol)
  }
  if (is.null(ncol) && !is.null(nrow)) {
    ncol <- ceiling(n/nrow)
  }
  stopifnot(nrow * ncol >= n)
  if (is.null(nrow) && is.null(ncol) && is.null(widths) && 
      is.null(heights)) {
    nm <- grDevices::n2mfrow(n)
    nrow = nm[1]
    ncol = nm[2]
  }
  inherit.ggplot <- unlist(lapply(grobs, inherits, what = "ggplot"))
  inherit.trellis <- unlist(lapply(grobs, inherits, what = "trellis"))
  if (any(inherit.ggplot)) {
    stopifnot(requireNamespace("ggplot2", quietly = TRUE))
    toconv <- which(inherit.ggplot)
    #grobs[toconv] <- lapply(grobs[toconv], ggplot2::ggplotGrob)
    grobs[toconv] <- lapply(grobs[toconv], ggplotGrob) ##NH
  }
  if (any(inherit.trellis)) {
    stopifnot(requireNamespace("lattice", quietly = TRUE))
    toconv <- which(inherit.trellis)
    grobs[toconv] <- lapply(grobs[toconv], ggint$latticeGrob)
  }
  if (missing(layout_matrix)) {
    positions <- expand.grid(t = seq_len(nrow), l = seq_len(ncol))
    positions$b <- positions$t
    positions$r <- positions$l
    if (as.table) 
      positions <- positions[order(positions$t), ]
    positions <- positions[seq_along(grobs), ]
  }
  else {
    cells <- sort(unique(as.vector(layout_matrix)))
    range_cell <- function(ii) {
      ind <- which(layout_matrix == ii, arr.ind = TRUE)
      c(l = min(ind[, "col"]), r = max(ind[, "col"]), t = min(ind[, 
                                                                  "row"]), b = max(ind[, "row"]))
    }
    positions <- data.frame(do.call(rbind, lapply(cells, 
                                                  range_cell)))
    ncol <- max(positions$r)
    nrow <- max(positions$b)
  }
  if (is.null(widths)) 
    widths <- unit(rep(1, ncol), "null")
  if (is.null(heights)) 
    heights <- unit(rep(1, nrow), "null")
  if (!is.unit(widths)) 
    widths <- unit(widths, "null")
  if (!is.unit(heights)) 
    heights <- unit(heights, "null")
  gt <- gtable(name = name, respect = respect, heights = heights, ##NH
               widths = widths, vp = vp)
  gt <- gtable_add_grob(gt, grobs, t = positions$t, b = positions$b, 
                        l = positions$l, r = positions$r, z = seq_along(grobs), 
                        clip = clip)
  if (is.character(top)) {
    top <- textGrob(top)
  }
  if (is.grob(top)) {
    h <- grobHeight(top) + padding
    gt <- gtable_add_rows(gt, heights = h, 0)
    gt <- gtable_add_grob(gt, top, t = 1, l = 1, r = ncol(gt), 
                          z = Inf, clip = clip)
  }
  if (is.character(bottom)) {
    bottom <- textGrob(bottom)
  }
  if (is.grob(bottom)) {
    h <- grobHeight(bottom) + padding
    gt <- gtable_add_rows(gt, heights = h, -1)
    gt <- gtable_add_grob(gt, bottom, t = nrow(gt), l = 1, 
                          r = ncol(gt), z = Inf, clip = clip)
  }
  if (is.character(left)) {
    left <- textGrob(left, rot = 90)
  }
  if (is.grob(left)) {
    w <- grobWidth(left) + padding
    gt <- gtable_add_cols(gt, widths = w, 0)
    gt <- gtable_add_grob(gt, left, t = 1, b = nrow(gt), 
                          l = 1, r = 1, z = Inf, clip = clip)
  }
  if (is.character(right)) {
    right <- textGrob(right, rot = -90)
  }
  if (is.grob(right)) {
    w <- grobWidth(right) + padding
    gt <- gtable_add_cols(gt, widths = w, -1)
    gt <- gtable_add_grob(gt, right, t = 1, b = nrow(gt), 
                          l = ncol(gt), r = ncol(gt), z = Inf, clip = clip)
  }
  gt
}



