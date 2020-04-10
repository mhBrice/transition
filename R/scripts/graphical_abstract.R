### Graphical abstract ####

### PACKAGES & FUNCTIONS ####

source("R/functions/packages.R")

source('R/functions/markov.R')

### DATA ####

source('R/functions/prep_data.R')

# Load msm results

load("res/msm_all75.rda")

msm_glb <- msm_all75[["msm_glb"]]

### PREP DATA ####

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

# Plot 

png("res/graphical_abstract.png", width = 5.5, height = 3.4, unit = "in", res = 300)
#quartz(width = 5.5, height = 3.4)
layout(matrix(c(1,2),1), widths = c(1,.3))
par(mar=c(3.3,2.9,1.2,.5))
plot_ss(mod = msm_glb, df = dfl, tp_ecotone = tp_ecotone, 
        dist = "logging", cex.axis = .8, cex.st = .9,
        unscale = sc_sTP,
        ylab = "", xlab = "")
mtext("State proportion at equilibrium", 2, line = 2, font = 2, cex = .95)
mtext("Mean temperature of the growing season\n Latitunal gradient", 
      1, line = 2.3, font = 2, cex = .95)


par(mar=c(3,0,5,0))
plot0()
mtext("Logging intensity", line = -.2, font = 2, cex = .87, at = -1.1, adj = 0, xpd = NA)
legend("topleft", legend = c("No or minor", "Moderate", "Major"), 
       cex = .85, pt.cex = 1,
       col = c("grey75","grey45","grey15"), pch = 19, lty = 1:3, lwd = 1.6, 
       seg.len = 2.8, x.intersp = 0.9, inset = c(-.03,0), bty = "n")

mtext("Current ecotone", 1, line = -5, font = 2, cex = .87, at = -1.1, adj = 0, xpd = NA)
polygon(x = c(-1,.8,.8,-1), y = c(-.2,-.2,-.6,-.6),
        col = alpha("grey", .2), border = NA)
dev.off()
