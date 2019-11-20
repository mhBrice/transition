### Plot steady state ####

### PACKAGES ####
library(graphicsutils)
library(dplyr)
library(msm)
library(sf)
library(scales)

### FUNCTIONS ####

source('R/functions/markov.R')
source('R/functions/plot_SS_change.R')

### DATA ####

source('R/functions/prep_data.R')

# Load msm results

load("res/msm_all75.rda")

msm_glb <- msm_all75[["msm_glb"]]

### 
ll_mixed <- which(states_ba$ecoreg3=="Mixed")
mixed_mean <- as.list(apply(states_ba[ll_mixed, 
                                      c("sTP", "sCMI", "DRAIN", "PH_HUMUS")],
                            2, mean))

tp_mixed <- quantile(states_ba$sTP[ll_mixed], c(.2,.8))


tp_grad <- seq(-1.9, 1.6, len = 50)

dfl <- expand.grid(sTP = tp_grad, logging = c(0, 1, 2), 
                  sCMI = mixed_mean[[2]],
                  DRAIN = mixed_mean[[3]],
                  PH_HUMUS = mixed_mean[[4]])

dfn <- expand.grid(sTP = tp_grad, natural = c(0, 1, 2), 
                   sCMI = mixed_mean[[2]],
                   DRAIN = mixed_mean[[3]],
                   PH_HUMUS = mixed_mean[[4]])

init <- states_ba %>% 
  group_by(ID_PE) %>% 
  arrange(year_measured) %>% 
  slice(1) %>% filter(year_measured<15)
init <- table(states_ba$states_ba[states_ba$year_measured<10])
init <- init/sum(init)

dn <- list(c(mixed_mean, natural1 = 0), 
           c(mixed_mean, natural1 = 1), 
           c(mixed_mean, logging1 = 1),
           c(mixed_mean, natural2 = 1), 
           c(mixed_mean, logging2 = 1))
qmats <- lapply(dn, function(x)
  qmatrix.msm(msm_glb, covariates = x, ci = "none"))

SSn <- lapply(qmats, function(x) steady_state(qmat = x))
SSn <- rbind(init, do.call(rbind, SSn))

dl <- list(c(mixed_mean, logging1 = 0), 
           c(mixed_mean, logging1 = 1), 
           c(mixed_mean, logging2 = 1))
qmats <- lapply(dl, function(x)
  qmatrix.msm(msm_glb, covariates = x, ci = "none"))

SSl <- lapply(qmats, function(x) steady_state(qmat = x))
SSl <- rbind(init, do.call(rbind, SSl))


### STEADY STATE ALONG TP GRADIENT ####


mat <- matrix(c(4,4,4,5,1:3,5), 4)

pdf("res/fig6_steady.pdf", width = 3.6, height = 7)
#quartz(width = 3.6, height = 7.2)
layout(mat, widths = c(0.08, 1), heights = c(1,1,1,0.2))
par(mar = c(2,2.2,1,0.2), oma = c(0,0,0.1,0))

barplot_index(index = SSn, ylim = c(0,1))
mtext(letters[1], 3, adj = -.17)
# Steady state for natural disturbances
x=plot_SS_change(mod = msm_glb, df = dfn, tp_mixed = tp_mixed, 
               dist = "natural",
               unscale = sc_sTP, axes = 2, 
               xlab = NULL, ylab = NULL, main = "Natural")
mtext(letters[2], 3, adj = -.17)

# Steady state for logging
x=plot_SS_change(mod = msm_glb, df = dfl, tp_mixed = tp_mixed, 
               dist = "logging",
               unscale = sc_sTP,
               ylab = NULL, xlab = "Mean temperature of the growing season",
               main = "Logging")
mtext(letters[3], 3, adj = -.17)


# Axis
par(mar = c(1,0,1.5,0))
plot0()
mtext("State proportion at equilibrium", 2, line = -2, cex =.85)

# Legend
par(mar = c(0,2.2,0,0.2))
plot0()
legend("bottom", legend = c("Minor", "Moderate", "Major"), cex = 1.1,
       col = c("grey75","grey45","grey15"), lty = 1:3, lwd = 1.5, 
       seg.len = 2, x.intersp = 0.9,
       horiz = TRUE, bty = "n")

dev.off()


### SUPP - STEADY STATE ALONG TP GRADIENT + DISTURBANCE GRADIENT ####

d_grad <- seq(0, 1, len = 11)
nat1_grad <- list(natural1 = d_grad)
nat2_grad <- list(natural2 = d_grad)
log1_grad <- list(logging1 = d_grad)
log2_grad <- list(logging2 = d_grad)

mat <- matrix(c(5,5,0,0,1:2,6,7,3,4,6,7), 4)

pdf("res/figSupp_SS_dfreq.pdf", width = 7, height = 8)
#quartz(width = 7, height = 8)
layout(mat, widths = c(0.06, 1,1), heights = c(1,1,0.3,.8))
par(mar = c(1,2,1,0.3), oma = c(2,0,0,0))

# Steady state for natural disturbances
int_nat1 <- plot_SS_change(mod = msm_glb, 
               other_covar = mixed_mean, 
               d_grad = nat1_grad, 
               tp_grad = tp_grad, 
               xlab = NULL, ylab = NULL, main = "Moderate natural", axes = 2,
               unscale = sc_sTP, tp_ecoreg = tp_ecoreg)
mtext("a", 3, adj = .05, line = -1.5, cex = 0.9)

int_nat2 <- plot_SS_change(mod = msm_glb, 
               other_covar = mixed_mean, 
               d_grad = nat2_grad, 
               tp_grad = tp_grad, 
               xlab = NULL, ylab = NULL, main = "Major natural", 
               unscale = sc_sTP, tp_ecoreg = tp_ecoreg)
mtext("c", 3, adj = .05, line = -1.5, cex = 0.9)

# Steady state for logging
int_log1 <- plot_SS_change(mod = msm_glb, 
                           other_covar = mixed_mean, 
                           d_grad = log1_grad, 
                           tp_grad = tp_grad, 
                           xlab = NULL, ylab = NULL, axes = 0,
                           main = "Moderate logging", 
                           unscale = sc_sTP, tp_ecoreg = tp_ecoreg)
mtext("b", 3, adj = .05, line = -1.5, cex = 0.9)

int_log2 <- plot_SS_change(mod = msm_glb, 
                           other_covar = mixed_mean, 
                           d_grad = log2_grad, 
                           tp_grad = tp_grad, 
                           xlab = NULL, ylab = NULL, axes = 1,
                           main = "Major logging", 
                           unscale = sc_sTP, tp_ecoreg = tp_ecoreg)
mtext("d", 3, adj = .05, line = -1.5, cex = 0.9)

# Axis
par(mar = c(1,0,1.5,0))
plot0()
text(0,0, labels = "State proportion at equilibrium", srt = 90, cex= 1.4, xpd = NA)

# Legend
par(mar = c(.5,6,2.5,5))
plot0()
mtext("Mean temperature of the growing season", 3, line = .7, cex = .9)

legend(-1, .5, xjust=0, yjust = .6, 
       legend = c("Boreal", "Temperate + Mixed"), cex = 1.2,
       col = st_col[c(1,4)], lwd = 1.5,
       xpd = NA, bty = "n", seg.len = 2.5)
for(i in d_grad) {
  points(i/4, .1, pch = 21, col = "grey", bg = alpha("black", i), cex = 1, lwd = .5)
}
arrows(x0 = 0, x1 = 1/4, y0 = -.2, y1 = -.2, length = .05)
text(-.03, -.2, "0", cex = 1.1)
text(.28, -.2, "1", cex = 1.1)
text(.6, 0, "Disturbance frequency", cex = 1.2)


par(lwd = 1.5, mar =c(3,10,0,11), las =1)
plot(int_nat1 ~ d_grad, type = "l", col = "grey65", ylim = c(12,13.1),
     ylab = "", xlab = "", xaxs = "i")
mtext("Frequency of disturbances", 1, line = 2.5, cex = 0.85)
mtext("Position of the B-T transition\nalong the temperature gradient", 2, 
      line = 2.7, las = 0, cex = 0.85)
lines(int_nat2 ~ d_grad, col = "grey65", lty = 2)
lines(int_log1 ~ d_grad)
lines(int_log2 ~ d_grad, lty = 2)
legend("bottomleft", 
       legend = c("Moderate natural", "Major natural", 
                  "Moderate logging", 'Major logging'),
       col = c("grey65", "grey65", "black", "black"), 
       lty = c(1,2,1,2), bty = "n")
mtext("e", 3, adj = .05, line = -1.5, cex = 0.9)

dev.off()




### MARKOV METRICS #####
transients_log <- index_gradient(mod = msm_glb, df = dfl)
transients_nat <- index_gradient(mod = msm_glb, df = dfn)


# Layout matrix
mat <- matrix(c(1:8),4)
mat <- rbind(mat, c(9,9))

pdf("res/fig7_transients.pdf", width = 5.5, height = 6.8)
#quartz(width = 5.5, height = 6.8)
layout(mat, heights = c(.12,1,1,1,.21), widths = c(1,1))
par(oma = c(0,2.2,0,0))

# Metrics for natural disturbances
par(mar = c(0,.5,0,.3))
plot0(text = "Natural", cex = 1.2, font = 2, xpd = NA)

par(mar = c(1,1,.1,.3))
plot_transients(index = matrix(transients_nat$soj_comm, ncol = 3), 
                tp_grad = tp_grad, tp_mixed = tp_mixed, 
                ylab = "Turnover time (y)", 
                unscale = sc_sTP, axes = 2, ylim = c(0,350))
mtext("a", 3, adj = .97, line = -1.2, cex = 0.9)

plot_transients(index = matrix(transients_nat$entropy_comm, ncol = 3), 
                tp_grad = tp_grad, tp_mixed = tp_mixed,
                ylab = "Entropy", 
                unscale = sc_sTP, axes = 2, ylim = c(0,.9))
mtext("c", 3, adj = .97, line = -1.2, cex = 0.9)

plot_transients(index = matrix(transients_nat$halflife, ncol = 3), 
                tp_grad = tp_grad, tp_mixed = tp_mixed,
                ylab = "Half-life to equilibrium (y)",
                unscale = sc_sTP, ylim = c(0,180))
mtext("e", 3, adj = .97, line = -1.2, cex = 0.9)

# Metrics for logging
par(mar = c(0,.5,0,.3))
plot0(text = "Logging", cex = 1.2, font = 2, xpd = NA)

par(mar = c(1,1,.1,.3))
plot_transients(index = matrix(transients_log$soj_comm, ncol = 3), 
                tp_grad = tp_grad, tp_mixed = tp_mixed, 
                axes = NULL,
                unscale = sc_sTP, ylim = c(0,350))
mtext("b", 3, adj = .97, line = -1.2, cex = 0.9)

plot_transients(index = matrix(transients_log$entropy_comm, ncol = 3), 
                tp_grad = tp_grad, tp_mixed = tp_mixed,
                unscale = sc_sTP, axes = NULL, ylim = c(0,.9))
mtext("d", 3, adj = .97, line = -1.2, cex = 0.9)

plot_transients(index = matrix(transients_log$halflife, ncol = 3), 
                tp_grad = tp_grad, tp_mixed = tp_mixed,
                axes = 1,
                unscale = sc_sTP,  ylim = c(0,180))
mtext("f", 3, adj = .97, line = -1.2, cex = 0.9)

par(mar = c(.1,.5,.7,.5))
plot0()
mtext("Mean temperature of the growing season", 3, line = -1, cex = 0.9)

legend("bottom", legend = c("Minor", "Moderate", "Major"), cex = 1.15,
       col = "black", horiz = TRUE, inset = c(0,-.2),
       lty = 1:3, lwd = 1.4,
       xpd = NA, bty = "n", yjust=.5, seg.len = 2.2)

dev.off()




### SUPPMAT - State contribution to turnover time ####

# Layout matrix
mat <- matrix(c(1:10), 5)
mat <- rbind(mat, c(11, 11))


pdf("res/figSupp_contrib2turnover.pdf", width = 5.5, height = 6.8)
#quartz(width = 5.5, height = 6.8)
par(oma = c(0,2,0,0))
layout(mat, heights = c(.12,1,1,1,1,.25))

# Metrics for natural disturbances
par(mar = c(0,.5,0,.3))
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
  plot_transients(index = matrix(transients_nat$soj_st_contrib[s,], ncol = 3), 
                  tp_grad = tp_grad, tp_mixed = tp_mixed,
                      main = states[s],
                      unscale = sc_sTP, axes = axes, ylim = ylim)
}
mtext("State contribution to turnover time", 2, outer = T, line = .8, cex = .87)

# Metrics for logging
par(mar = c(0,.5,0,.3))
plot0(text = "Logging", cex = 1.2, font = 2, xpd = NA)

par(mar = c(1,1,.1,.3))
for(s in 1:4) {
  axes = 0
  if(s==1) ylim = c(0, 350)
  if(s==2) ylim = c(0, 100)
  if(s==3) ylim = c(0, 200)
  if(s==4) { 
    ylim = c(0, 200) 
    axes = 1
  }
  plot_transients(index = matrix(transients_log$soj_st_contrib[s,], ncol = 3),
                      tp_grad = tp_grad, tp_mixed = tp_mixed,
                      unscale = sc_sTP, axes = axes, ylim = ylim)
}


par(mar = c(0,.5,0,.5))
plot0()
text(0, .25,"Mean temperature of the growing season", cex = 1.3, xpd = NA)
legend(0, .1, legend = c("Minor", "Moderate", "Major"), cex = 1.1,
       col = "black", horiz = TRUE,
       lty = 1:3, lwd = 1.4,
       xpd = NA, bty = "n", xjust = 0.5, seg.len = 2.2)

dev.off()



### SUPPMAT - State ontribution to entropy ####

pdf("res/figSupp_contrib2entropy.pdf", width = 5.5, height = 6.8)
#quartz(width = 5.5, height = 6.8)
par(oma = c(0,2,0,0))
layout(mat, heights = c(.12,1,1,1,1,.25))

# Metrics for natural disturbances
par(mar = c(0,.5,0,.3))
plot0(text = "Natural", cex = 1.2, font = 2, xpd = NA)

par(mar = c(1,1,.1,.3))

for(s in 1:4) {
  axes = 2
  if(s==4) axes = c(1, 2) 

  plot_transients(index = matrix(transients_nat$entropy_st_contrib[s,], ncol = 3), 
                  tp_grad = tp_grad, tp_mixed = tp_mixed,
                  main = states[s],
                  unscale = sc_sTP, axes = axes, ylim = c(0,1))
}
mtext("State contribution to entropy", 2, outer = T, line = .8, cex = .87)

# Metrics for logging
par(mar = c(0,.5,0,.3))
plot0(text = "Logging", cex = 1.2, font = 2, xpd = NA)

par(mar = c(1,1,.1,.3))
for(s in 1:4) {
  axes = 0
  if(s==4) axes = 1
  plot_transients(index = matrix(transients_log$entropy_st_contrib[s,], ncol = 3),
                  tp_grad = tp_grad, tp_mixed = tp_mixed,
                  unscale = sc_sTP, axes = axes, ylim = c(0,1))
}


par(mar = c(0,.5,0,.5))
plot0()
text(0, .25,"Mean temperature of the growing season", cex = 1.3, xpd = NA)
legend(0, .1, legend = c("Minor", "Moderate", "Major"), cex = 1.1,
       col = "black", horiz = TRUE,
       lty = 1:3, lwd = 1.4,
       xpd = NA, bty = "n", xjust = 0.5, seg.len = 2.2)


dev.off()

