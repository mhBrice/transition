##################################
### Figure 1. REGION MAP ####
##################################

library(sf)
require(RColorBrewer)
require(scales)

reg_title <- c("Sugar maple-hickory\nSugar maple-basswood",
               "Sugar maple-yellow birch",
               "Balsam fir-yellow birch",
               "Balsam fir-white birch", "Spruce-moss")

col_reg <- brewer.pal(6,"Spectral")[c(1,1:3,5,6)]

n_reg <- table(BCDdf$ecoreg5)

lim <- st_bbox(ecoregion)

pdf("ms/figures/fig1_region.pdf",
    width = 5, height = 3)
# quartz(width = 5, height = 3)
par(mar=c(1.8,2.3,.3,0.3))

plot0(xlim = lim[c(1,3)]+c(.3,.1), ylim = lim[c(2,4)]+c(-.1,.1),
      grid.col = alpha("grey60", .3), fill = "white")

box2(1:2)
plot(st_geometry(ecoregion), border = "grey50",
     col = alpha(col_reg,.4)[ecoregion$SOUS_DOM6],
     axes=F, add=T)

axis(1, labels = F, tcl = -.4)
axis(1, at = seq(-80,-60,by=5), labels = paste0(abs(seq(-80,-60,by=5)), "°W"),
     line = -.3, cex.axis = .8, tick = F)
axis(2, labels = F, tcl = -.4)
axis(2, at = seq(46,52,by=2), labels = paste0(seq(46,52,by=2), "°N"),
     line = -.3, las = 1, cex.axis=.8, tick = F)

plot(st_geometry(xy), cex = .1, pch = 19, col = alpha("grey15",.3), add = T)

points(-61.8, 47.5, cex= 3, pch = 19, col = "white")
text(-59.7, 53.1, paste0(reg_title[5], " (", n_reg[5], ")"), cex = .75, xpd = NA)
text(-59.2, 49.6, paste0(reg_title[4], "\n(", n_reg[4], ")"), cex = .75, xpd = NA)
text(-63.8, 47.5, paste0(reg_title[3], " (", n_reg[3], ")"), cex = .75)
text(-64.9, 46.5, paste0(reg_title[2], " (", n_reg[2], ")"), cex = .75)
text(-66.8, 45.2, paste0(reg_title[1], " (", n_reg[1], ")"), cex = .75)

dev.off()
