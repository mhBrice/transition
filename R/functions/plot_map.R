plot_map <- function(ecoregion, col_reg = "transparent", box = c(1,2),
                     xy_pts = NULL, cex = .1, 
                     pch = 19, col = alpha("grey15",.3), bg = NULL,
                     axes = c(1,2)) {
  lim <- st_bbox(ecoregion)
  xax <-  seq(-80, -60, by = 5)
  yax <- seq(46, 52, by = 2)
  
  plot0(xlim = lim[c(1,3)]+c(.5,.2), ylim = lim[c(2,4)]+c(-.1,.1),
        grid.col = alpha("grey60", .3), fill = "white")
  
  plot(st_geometry(ecoregion), border = "grey50", 
       col = col_reg[ecoregion$SOUS_DOM6],
       axes = FALSE, add = TRUE)
  
  box2(box)

  
  axis(1, labels = FALSE, tcl = -.4)
  if(any(axes == 1))  axis(1, at = xax, labels = paste0(abs(xax), "°W"),
                           line = -.6, cex.axis = .75, tick = FALSE)
  axis(2, labels = FALSE, tcl = -.4)
  if(any(axes == 2))axis(2, at = yax, labels = paste0(yax, "°N"),
                         line = -.5, las = 1, cex.axis = .75, tick = FALSE)
  
  if(!is.null(xy_pts)){
    points(st_coordinates(xy_pts), pch = pch, col = col, bg = bg, cex = cex)
  }
}
