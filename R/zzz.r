# SET VERSION
build <- "Jacob"

# SET VERBOSITY

# dartR functions have a verbosity parameter that sets the level of reporting during the execution of the function. The verbosity level,
# set by parameter 'verbose' can be one of verbose 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary;
# 5, full report.

#options("dartR_verbose"=2)
# SET MESSAGES COLORS

# - For fatal errors use “error” which will print the message in red. Example usage: stop(error(“Fatal error”))
# - For warning messages use “warn” which will print the message in yellow. Example usage: cat(warn(“message”))
# - For reporting messages use “report” which will print the message in green. Example usage: cat(report(“message”))
# - For important messages use “important” which will print the message in blue. Example usage: cat(important(“message”))
# - For other messages as code use “code” which will print the message in cyan. Example usage: cat(code(“message”))
error <- crayon::red 
warn <- crayon::yellow
report <- crayon::green
important <- crayon::blue 
code <- crayon::cyan

# SET PLOTS COLORS

# function to replicate defaults colors of ggplot
discrete_palette <- function(n) {
  hues = seq(15, 375, length = n + 1)
  return(hcl(h = hues, l = 65, c = 100)[1:n])
}
# taken from wes_palette::Zissou1
diverging_palette <- colorRampPalette(c("#3B9AB2" ,"#78B7C5" ,"#EBCC2A" ,"#E1AF00" ,"#F21A00"))
# creating convergent to 0 palette
cool <- rainbow(50, start = rgb2hsv(col2rgb("cyan"))[1], end = rgb2hsv(col2rgb("blue"))[1])
warm <- rainbow(50, start = rgb2hsv(col2rgb("red"))[1], end = rgb2hsv(col2rgb("yellow"))[1])
cols <- c(rev(cool), rev(warm))
convergent_palette <- colorRampPalette(cols)
# taken from adegenet
viridis_palette <- colorRampPalette(c("#440154FF", "#482173FF", "#433E85FF", "#38598CFF",
                                      "#2D708EFF", "#25858EFF", "#1E9B8AFF", "#2BB07FFF",
                                      "#51C56AFF", "#85D54AFF", "#C2DF23FF", "#FDE725FF"))
two_colors <- c("#3B9AB2" ,"#78B7C5" )
three_colors <- c("#3B9AB2","deeppink","lemonchiffon")

# SET THEME FOR PLOTS

#' dartR theme
#'
#' This is the theme used as default for dartR plots.
#' This function controls all non-data display elements in the plots. 
#'
#' @param base_size base font size, given in pts.
#' @param base_family base font family
#' @param base_line_size base size for line elements
#' @param base_rect_size base size for rect elements
#' @examples
#' 
#' ggplot(data.frame(dummy=rnorm(1000)),aes(dummy)) +
#'  geom_histogram() +
#'  theme_dartR()
#'

# The half-line (base-fontsize / 2) sets up the basic vertical
# rhythm of the theme. Most margins will be set to this value.
# However, when we work with relative sizes, we may want to multiply
# `half_line` with the appropriate relative size. This applies in
# particular for axis tick sizes. And also, for axis ticks and
# axis titles, `half_size` is too large a distance, and we use `half_size/2`
# instead.
theme_dartR <- function(base_size = 11, base_family = "",
                        base_line_size = base_size / 22,
                        base_rect_size = base_size / 22) {
  half_line <- base_size / 2
  
  # Throughout the theme, we use three font sizes, `base_size` (`rel(1)`)
  # for normal, `rel(0.8)` for small, and `rel(1.2)` for large.
  
  # Elements in this first block aren't used directly, but are inherited by others
  t <- theme(
    line =               element_line(
      colour = "black", size = base_line_size,
      linetype = 1, lineend = "butt"
    ),
    rect =               element_rect(
      fill = "white", colour = "black",
      size = base_rect_size, linetype = 1
    ),
    text = element_text(family = base_family, face = "plain",colour = "black", size = base_size,lineheight = 0.9, hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), debug = FALSE),
    axis.line =          element_blank(),axis.line.x=NULL,axis.line.y = NULL,
    axis.text =          element_text(size = rel(1), colour = "black"),
    axis.text.x =        element_text(margin = margin(t = 0.8 * half_line / 2), vjust = 1),
    axis.text.x.top =    element_text(margin = margin(b = 0.8 * half_line / 2), vjust = 0),
    axis.text.y =        element_text(margin = margin(r = 0.8 * half_line / 2), hjust = 1),
    axis.text.y.right =  element_text(margin = margin(l = 0.8 * half_line / 2), hjust = 0),
    axis.ticks =         element_line(colour = "gray80"),
    axis.ticks.length =  unit(half_line / 2, "pt"),
    axis.ticks.length.x = NULL,
    axis.ticks.length.x.top = NULL,
    axis.ticks.length.x.bottom = NULL,
    axis.ticks.length.y = NULL,
    axis.ticks.length.y.left = NULL,
    axis.ticks.length.y.right = NULL,
    axis.title.x = element_text(margin = margin(t = half_line / 2),vjust = 1,face="bold"),
    axis.title.x.top = element_text(margin = margin(b = half_line / 2),vjust = 0),
    axis.title.y = element_text(angle = 90,margin = margin(r = half_line / 2),vjust = 1,face="bold"),
    axis.title.y.right = element_text(angle = -90,margin = margin(l = half_line / 2),vjust = 0),
    legend.background =  element_rect(colour = NA),
    legend.spacing =     unit(2 * half_line, "pt"),
    legend.spacing.x =    NULL,
    legend.spacing.y =    NULL,
    legend.margin =      margin(half_line, half_line, half_line, half_line),
    legend.key =         element_rect(fill = "transparent", colour = NA),
    legend.key.size =    unit(1.2, "lines"),
    legend.key.height =  NULL,
    legend.key.width =   NULL,
    legend.text =        element_text(size = rel(0.8)),
    legend.text.align =  NULL,
    legend.title =       element_text(hjust = 0),
    legend.title.align = NULL,
    legend.position =    "right",
    legend.direction =   NULL,
    legend.justification = "center",
    legend.box =         NULL,
    legend.box.margin =  margin(0, 0, 0, 0, "cm"),
    legend.box.background = element_blank(),
    legend.box.spacing = unit(2 * half_line, "pt"),
    panel.background =   element_rect(fill = "white", colour = NA),
    panel.border =       element_blank(),
    panel.grid =         element_line(colour = "gray80"),
    panel.grid.minor =   element_line(size = rel(0.5)),
    panel.spacing =      unit(half_line, "pt"),
    panel.spacing.x =    NULL,
    panel.spacing.y =    NULL,
    panel.ontop    =     FALSE,
    strip.background =   element_rect(fill = "white", colour = "black"),
    strip.text =         element_text(
      colour = "black",
      size = rel(1),
      face = "bold",
      margin = margin(0.8 * half_line, 0.8 * half_line, 0.8 * half_line, 0.8 * half_line)
    ),
    strip.text.x =      element_text(size = 14, face = "bold"),
    strip.text.y =       element_text(angle = -90),
    strip.text.y.left =  element_text(angle = 90),
    strip.placement =    "inside",
    strip.placement.x =  NULL,
    strip.placement.y =  NULL,
    strip.switch.pad.grid = unit(half_line / 2, "pt"),
    strip.switch.pad.wrap = unit(half_line / 2, "pt"),
    plot.background =    element_rect(colour = "white"),
    plot.title =         element_text( 
      face="bold",
      size = rel(1.2),
      hjust = 0.5, vjust = 1,
      margin = margin(b = half_line)
    ),
    plot.title.position = "panel",
    plot.subtitle =      element_text( 
      hjust = 0.5, vjust = 1,
      margin = margin(b = half_line)
    ),
    plot.caption =       element_text(
      size = rel(0.8),
      hjust = 1, vjust = 1,
      margin = margin(t = half_line)
    ),
    plot.caption.position = "panel",
    plot.tag =           element_text(
      size = rel(1.2),
      hjust = 0.5, vjust = 0.5
    ),
    plot.tag.position =  'topleft',
    plot.margin =        margin(half_line, half_line, half_line, half_line),
    complete = TRUE
  )
}

# WELCOME MESSAGE
.onAttach <- function(...) {
  
  packageStartupMessage(important("**** Welcome to dartR ****\n"))
  packageStartupMessage(report("Be aware that owing to CRAN requirements and compatibility reasons not all functions of the packages may run yet, as some dependencies could be missing. Hence for a most enjoyable experience we recommend to run the function "))
  packageStartupMessage(code("gl.install.vanilla.dartR()"))                 
  
  packageStartupMessage(report("This installs all missing and required packages for your version of dartR. \nFor citation information please use:"))
  packageStartupMessage(code("citation('dartR')"))
  
  options("dartR_verbose"=2)
  packageStartupMessage(report("Global verbosity is set to: 2\n"))
  
  packageStartupMessage(important("\n**** Have fun using dartR! ****"))
}

