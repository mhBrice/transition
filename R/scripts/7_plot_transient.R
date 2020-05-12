### PLOT TRANSIENT MEASURES ####

# Figure 6
# Figure S8
# Figure S9

### PACKAGES & FUNCTIONS ####

source("R/functions/packages.R")

source('R/functions/markov.R')

### DATA ####

source('R/functions/prep_data.R')

# Load msm results

load("res/msm_all75.rda")

msm_glb <- msm_all75[["msm_glb"]]

### DATA PREP ####

# Mean environmental conditions at the ecotone

ll_ecotone <- which(states_ba$ecoreg3=="Mixed")

mu_ecotone <- as.list(apply(states_ba[ll_ecotone, c("sTP", "sCMI", "DRAIN", "PH_HUMUS")],
                            2, mean))

# Temperature boundaries of the ecotone

tp_ecotone <- quantile(states_ba$sTP[ll_ecotone], c(.2, .8))

# Temperature gradient

tp_grad <- seq(-1.9, 1.6, len = 50)

# Expand data frame with temperature gradient and disturbances levels

dfl <- expand.grid(sTP = tp_grad, logging = c(0, 1, 2),
                   sCMI = mu_ecotone[[2]],
                   DRAIN = mu_ecotone[[3]],
                   PH_HUMUS = mu_ecotone[[4]])

dfn <- expand.grid(sTP = tp_grad, natural = c(0, 1, 2),
                   sCMI = mu_ecotone[[2]],
                   DRAIN = mu_ecotone[[3]],
                   PH_HUMUS = mu_ecotone[[4]])


### FIGURE 6. TRANSIENT METRICS #####

transients_log <- transient_index(mod = msm_glb, df = dfl)
transients_nat <- transient_index(mod = msm_glb, df = dfn)


# Layout matrix

mat <- matrix(c(1:8), nrow = 4)
mat <- rbind(mat, c(9, 9))

pdf("res/fig6_transients.pdf", width = 5.5, height = 6.8)

layout(mat, heights = c(.12, 1, 1, 1, .21), widths = c(1, 1))

par(oma = c(0, 2.2, 0, 0))

# Metrics for natural disturbances

par(mar = c(0, .5, 0, .3))
plot0(text = "Natural", cex = 1.2, font = 2, xpd = NA)

par(mar = c(1,1,.1,.3))
plot_transient(index = matrix(transients_nat$soj_comm, ncol = 3),
               tp_grad = tp_grad, tp_ecotone = tp_ecotone,
               ylab = "Turnover time (y)",
               unscale = sc_sTP, axes = 2, ylim = c(0,350))
mtext("(a)", 3, adj = .97, line = -1.2, cex = 0.8, font = 2)

plot_transient(index = matrix(transients_nat$entropy_comm, ncol = 3),
               tp_grad = tp_grad, tp_ecotone = tp_ecotone,
               ylab = "Entropy",
               unscale = sc_sTP, axes = 2, ylim = c(0,.9))
mtext("(c)", 3, adj = .97, line = -1.2, cex = 0.8, font = 2)

plot_transient(index = matrix(transients_nat$halflife, ncol = 3),
               tp_grad = tp_grad, tp_ecotone = tp_ecotone,
               ylab = "Half-life to equilibrium (y)",
               unscale = sc_sTP, ylim = c(0,180))
mtext("(e)", 3, adj = .97, line = -1.2, cex = 0.8, font = 2)

# Metrics for logging

par(mar = c(0,.5,0,.3))
plot0(text = "Logging", cex = 1.2, font = 2, xpd = NA)

par(mar = c(1,1,.1,.3))
plot_transient(index = matrix(transients_log$soj_comm, ncol = 3),
               tp_grad = tp_grad, tp_ecotone = tp_ecotone,
               axes = NULL,
               unscale = sc_sTP, ylim = c(0, 350))
mtext("(b)", 3, adj = .97, line = -1.2, cex = 0.8, font = 2)

plot_transient(index = matrix(transients_log$entropy_comm, ncol = 3),
               tp_grad = tp_grad, tp_ecotone = tp_ecotone,
               unscale = sc_sTP, axes = NULL, ylim = c(0, .9))
mtext("(d)", 3, adj = .97, line = -1.2, cex = 0.8, font = 2)

plot_transient(index = matrix(transients_log$halflife, ncol = 3),
               tp_grad = tp_grad, tp_ecotone = tp_ecotone,
               axes = 1,
               unscale = sc_sTP,  ylim = c(0, 180))
mtext("(f)", 3, adj = .97, line = -1.2, cex = 0.8, font = 2)

# Axis & Legend

par(mar = c(.1, .5, .7, .5))
plot0()
mtext("Mean temperature of the growing season (°C)", 3, line = -1, cex = 0.9)

legend("bottom", legend = c("No or minor", "Moderate", "Major"), cex = 1.15,
       col = "black", horiz = TRUE, inset = c(0, -.2),
       lty = 1:3, lwd = 1.4,
       xpd = NA, bty = "n", yjust = .5, seg.len = 2.2)

dev.off()




### FIGURE S8. State contribution to turnover time ####

# Layout matrix

mat <- matrix(c(1:10), 5)
mat <- rbind(mat, c(11, 11))


pdf("res/figS8_contrib2turnover.pdf", width = 5.5, height = 6.8)
#quartz(width = 5.5, height = 6.8)

par(oma = c(0, 2, 0, 0))
layout(mat, heights = c(.12, 1, 1, 1, 1, .25))

# Metrics for natural disturbances

par(mar = c(0, .5, 0, .3))
plot0(text = "Natural", cex = 1.2, font = 2, xpd = NA)

par(mar = c(1,1,.1,.3))

for(s in 1:4) {
  axes = 2
  if(s==1) ylim = c(0, 350)
  if(s==2) ylim = c(0, 100)
  if(s==3) ylim = c(0, 200)
  if(s==4) {
    ylim = c(0, 200)
    axes = c(1, 2)
  }
  plot_transient(index = matrix(transients_nat$soj_st_contrib[s,], ncol = 3),
                 tp_grad = tp_grad, tp_ecotone = tp_ecotone,
                 main = states[s],
                 unscale = sc_sTP, axes = axes, ylim = ylim)
}
mtext("State contribution to turnover time", 2, outer = TRUE, line = .8, cex = .87)

# Metrics for logging

par(mar = c(0, .5, 0, .3))
plot0(text = "Logging", cex = 1.2, font = 2, xpd = NA)

par(mar = c(1, 1, .1, .3))
for(s in 1:4) {
  axes = 0
  if(s==1) ylim = c(0, 350)
  if(s==2) ylim = c(0, 100)
  if(s==3) ylim = c(0, 200)
  if(s==4) {
    ylim = c(0, 200)
    axes = 1
  }
  plot_transient(index = matrix(transients_log$soj_st_contrib[s,], ncol = 3),
                 tp_grad = tp_grad, tp_ecotone = tp_ecotone,
                 unscale = sc_sTP, axes = axes, ylim = ylim)
}

# Axis & Legend

par(mar = c(0, .5, 0, .5))
plot0()
text(0, .25,"Mean temperature of the growing season (°C)", cex = 1.3, xpd = NA)
legend(0, .1, legend = c("No or minor", "Moderate", "Major"), cex = 1.1,
       col = "black", horiz = TRUE,
       lty = 1:3, lwd = 1.4,
       xpd = NA, bty = "n", xjust = 0.5, seg.len = 2.2)

dev.off()



### FIGURE S9. State ontribution to entropy ####

pdf("res/figS9_contrib2entropy.pdf", width = 5.5, height = 6.8)

par(oma = c(0, 2, 0, 0))
layout(mat, heights = c(.12, 1, 1, 1, 1, .25))

# Metrics for natural disturbances

par(mar = c(0, .5, 0, .3))
plot0(text = "Natural", cex = 1.2, font = 2, xpd = NA)

par(mar = c(1, 1, .1, .3))

for(s in 1:4) {
  axes = 2
  if(s==4) axes = c(1, 2)

  plot_transient(index = matrix(transients_nat$entropy_st_contrib[s,], ncol = 3),
                 tp_grad = tp_grad, tp_ecotone = tp_ecotone,
                 main = states[s],
                 unscale = sc_sTP, axes = axes, ylim = c(0, 1))
}
mtext("State contribution to entropy", 2, outer = TRUE, line = .8, cex = .87)

# Metrics for logging

par(mar = c(0, .5, 0, .3))
plot0(text = "Logging", cex = 1.2, font = 2, xpd = NA)

par(mar = c(1,1,.1,.3))
for(s in 1:4) {
  axes = 0
  if(s==4) axes = 1
  plot_transient(index = matrix(transients_log$entropy_st_contrib[s,], ncol = 3),
                 tp_grad = tp_grad, tp_ecotone = tp_ecotone,
                 unscale = sc_sTP, axes = axes, ylim = c(0, 1))
}

# Axis & Legend

par(mar = c(0, .5, 0, .5))
plot0()
text(0, .25, "Mean temperature of the growing season (°C)", cex = 1.3, xpd = NA)
legend(0, .1, legend = c("No or minor", "Moderate", "Major"), cex = 1.1,
       col = "black", horiz = TRUE,
       lty = 1:3, lwd = 1.4,
       xpd = NA, bty = "n", xjust = 0.5, seg.len = 2.2)


dev.off()
