##################################
### Figure 1. REGION MAP ####
##################################

library(sf)
library(RColorBrewer)
library(scales)

source('R/functions/plot_map.R')

### DATA ####

source('R/functions/prep_data.R')

reg_title <- c("Sugar maple-hickory &\nSugar maple-basswood",
               "Sugar maple-yellow birch",
               "Balsam fir-yellow birch",
               "Balsam fir-white birch", "Spruce-moss")

col_reg <- brewer.pal(6,"Spectral")[c(1,1:3,5,6)]

n_reg <- table(states_ba$ecoreg6)


pdf("res/fig1_region.pdf", width = 4.3, height = 2.4)
# quartz(width = 4.3, height = 2.4)
par(mar=c(1.5,2.1,.3,0.3))

plot_map(ecoregion, col_reg = alpha(col_reg,.4), xy_pts = xy)

points(-61.8, 47.4, cex = 3, pch = 15, col = "white")
text(-59.4, 53.15, paste0(reg_title[5], " (", n_reg[6], ")"), cex = .7, xpd = NA)
text(-59.1, 49.55, paste0(reg_title[4], "\n(", n_reg[5], ")"), cex = .7, xpd = NA)
text(-63.4, 47.5, paste0(reg_title[3], " (", n_reg[4], ")"), cex = .7)
text(-64.5, 46.5, paste0(reg_title[2], " (", n_reg[3], ")"), cex = .7)
text(-66.4, 45.2, paste0(reg_title[1], " (", n_reg[2]+n_reg[1], ")"), cex = .7)

dev.off()
