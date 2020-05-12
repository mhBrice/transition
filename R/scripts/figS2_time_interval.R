### DISTRIBUTION OF TIME INTERVALS ####

### PACKAGES & FUNCTIONS ####
source("R/functions/packages.R")
source("R/functions/cross_validation.R")

### DATA ####
source('R/functions/prep_data.R')


st_tr <- to_trans(states_ba, covar_names = "ecoreg6")

levels(st_tr$ecoreg6) <- c("Sugar maple-hickory \nSugar maple-basswood",
                           "Sugar maple-hickory \nSugar maple-basswood",
                           "Sugar maple-yellow birch",
                           "Balsam fir-yellow birch",
                           "Balsam fir-white birch",
                           "Spruce-moss")

pdf("res/figS2_hist_intervals.pdf", width = 4.5, height = 6.5)

par(mfrow = c(5, 1), mar = c(3.5, 4.5, .5, .6))
for(r in levels(st_tr$ecoreg6)) {
  hist(st_tr$t[st_tr$ecoreg6 == r],
       main = "", ylab = "", xlab = "",
       xlim = c(0, 40), ylim = c(0, 1200), breaks = seq(0, 40, 1),
       xpd = NA, las = 1, xaxs = "i", yaxs = "i",
       cex.axis = .95, tcl = -.5)
  mtext(r, 3, adj = .9, font = 2, cex = 0.8, line = -1, padj = .5)

}
mtext("Time intervals between survey (years)", 1, outer = TRUE,
      line = -1.2, cex = 0.85, font = 2)
mtext("Frequency of forest plots", 2, outer = TRUE,
      line = -1.3, cex = 0.85, font = 2)

dev.off()
