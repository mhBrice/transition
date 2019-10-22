### Plot steady state ####

### PACKAGES ####
library(graphicsutils)
library(dplyr)
library(msm)
library(sf)

### FUNCTIONS ####

source('R/functions/markov.R')
source('R/functions/plot_SS_change.R')

### DATA ####

source('R/functions/prep_data.R')

# Load msm results

load("res/msm_all75.rda")

msm_glb <- msm_all75[["msm_glb"]]

### 
envmean <- aggregate(states_ba[,c("CMI", "DRAIN", "PH_HUMUS")], 
                     by = list(states_ba$ecoreg3), mean)
mixed_mean <- envmean[3,-1]

states_ba$ecoreg5 <- states_ba$ecoreg6
levels(states_ba$ecoreg5) <- c("Sugar maple-basswood", 
                               "Sugar maple-basswood",
                               "Sugar maple-yellow birch",
                               "Balsam fir-yellow birch",
                               "Balsam fir-white birch",
                               "Spruce-moss")

tp_ecoreg_min <- aggregate(states_ba$sTP, 
                       by = list(states_ba$ecoreg5), function(x) quantile(x, .25))
tp_ecoreg_max <- aggregate(states_ba$sTP, 
                           by = list(states_ba$ecoreg5), function(x) quantile(x, .75))
tp_ecoreg <- cbind(tp_ecoreg_min[,2], tp_ecoreg_max[,2])

reg_title <- c("Sugar maple-hickory & \nSugar maple-basswood",
               "Sugar maple-yellow birch",
               "Balsam fir-yellow birch",
               "Balsam fir-white birch", 
               "Spruce-moss")

col_reg = c("#D53E4F", "#FC8D59", "#FEE08B", "#99D594", "#3288BD")

### Covariates ####
quantile(states_ba$sTP, c(.03,.97))
tp_grad <- seq(-1.9, 1.6, len = 50)

nat_level <- list(list(natural1 = 0, natural2 = 0), 
                  list(natural1 = 1, natural2 = 0), 
                  list(natural1 = 0, natural2 = 1))
log_level <- list(list(logging1 = 0, logging2 = 0), 
                  list(logging1 = 1, logging2 = 0), 
                  list(logging1 = 0, logging2 = 1))



covar_d <- list(c(mixed_mean), 
                c(natural1 = 1, mixed_mean),
                c(natural2 = 1, mixed_mean), 
                c(logging1 = 1, mixed_mean),
                c(logging2 = 1, mixed_mean))


### STEADY STATE ALONG TP GRADIENT ####

mat <- matrix(c(4,4,1:2,3,3), 2)

pdf("res/fig4_SS_gradient.pdf", width = 6, height = 5.5)
#quartz(width = 6, height = 5.5)
layout(mat, widths = c(0.06, 1, .5))
par(mar = c(1,2,1.5,0), oma = c(2,0,0,0))

# Steady state for natural disturbances
plot_SS_change(mod = msm_glb, 
               other_covar = mixed_mean, 
               d_grad = nat_level, 
               tp_grad = tp_grad, 
               xlab = NULL, ylab = NULL, main = "Natural", axes = 2,
               unscale = sc_sTP, tp_ecoreg = tp_ecoreg)



# Steady state for logging
plot_SS_change(mod = msm_glb, 
               other_covar = mixed_mean, 
               d_grad = log_level, 
               tp_grad = tp_grad, 
               ylab = NULL, xlab = "Mean temperature of the growing season",
               main = "Logging",
               unscale = sc_sTP, tp_ecoreg = tp_ecoreg)


# Legend
par(mar = c(7,0,7,0))
plot0(xaxs="i")
legend(-1,.5, legend = c("Minor", "Moderate", "Major"), cex = 1.1,
       pch = 21, col = "black", pt.bg = c("white", "grey55", "black"), 
       pt.cex = 2, pt.lwd = 1,
       lty = 1:3, lwd = 1.5,
       xpd = NA, bty = "n", xjust = -.1, yjust = .5, seg.len = 2.5)
legend(-1, 0, legend = c("Boreal", "Temperate + Mixed"), cex = 1.1,
       col = st_col[c(1,4)], lwd = 1.5,
       xpd = NA, bty = "n", yjust = .5, seg.len = 2.5)
legend(-1,-.7, legend = rev(reg_title), cex = 1.1,
       col = rev(col_reg), lwd = 2.5,
       xpd = NA, bty = "n", yjust = .5, seg.len = 2.5)

# Axis
par(mar = c(1,0,1.5,0))
plot0()
text(0,0, labels = "State proportion at equilibrium", srt = 90, cex= 1.4, xpd = NA)

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

index_res <- index_gradient(mod = msm_glb, covar_d = covar_d, tp_grad = tp_grad)

# Find min and max
# tp_grad[which.max(index_res$halflife_df[,1])]  * sc_sTP[2] + sc_sTP[1]


# Layout matrix
mat <- matrix(c(1:8),4)
mat <- rbind(mat, c(9,9))
mat <- cbind(mat, c(0, 10,11,0, 0))

pdf("res/fig7_index_gradient.pdf", width = 6.7, height = 6.8)
#quartz(width = 6.7, height = 6.8)
layout(mat, heights = c(.12,1,1,1,.14), widths = c(1,1,0.8))
par(oma = c(0,2,0,0))

# Metrics for natural disturbances
par(mar = c(0,.5,0,.3))
plot0(text = "Natural", cex = 1.2, font = 2, xpd = NA)

par(mar = c(1,1,.1,.3))
plot_index_gradient(index = index_res$soj_df[,1:3], tp_grad = tp_grad, 
                    ylab = "Turnover time (y)", 
                    unscale = sc_sTP,axes = 2, ylim = c(0,370),
                    tp_ecoreg = tp_ecoreg)
mtext("a", 3, adj = .97, line = -1.2, cex = 0.9)

plot_index_gradient(index = index_res$entropy_df[,1:3], tp_grad = tp_grad, 
                    ylab = "Entropy", 
                    unscale = sc_sTP, axes = 2, ylim = c(0,.9), 
                    tp_ecoreg = tp_ecoreg)
mtext("c", 3, adj = .97, line = -1.2, cex = 0.9)

plot_index_gradient(index = index_res$halflife_df[,1:3], tp_grad = tp_grad, 
                    ylab = "Half-life to equilibrium (y)",
                    unscale = sc_sTP, ylim = c(0,160), 
                    tp_ecoreg = tp_ecoreg)
mtext("e", 3, adj = .97, line = -1.2, cex = 0.9)

# Metrics for logging
par(mar = c(0,.5,0,.3))
plot0(text = "Logging", cex = 1.2, font = 2, xpd = NA)

par(mar = c(1,1,.1,.3))
plot_index_gradient(index = index_res$soj_df[,c(1,4,5)], tp_grad = tp_grad, 
                    axes = NULL,
                    unscale = sc_sTP, ylim = c(0,370), 
                    tp_ecoreg = tp_ecoreg)
mtext("b", 3, adj = .97, line = -1.2, cex = 0.9)

plot_index_gradient(index = index_res$entropy_df[,c(1,4,5)], tp_grad = tp_grad, 
                    unscale = sc_sTP, axes = NULL, ylim = c(0,.9), 
                    tp_ecoreg = tp_ecoreg)
mtext("d", 3, adj = .97, line = -1.2, cex = 0.9)

plot_index_gradient(index = index_res$halflife_df[,c(1,4,5)], tp_grad = tp_grad, 
                    axes = 1,
                    unscale = sc_sTP,  ylim = c(0,160), 
                    tp_ecoreg = tp_ecoreg)
mtext("f", 3, adj = .97, line = -1.2, cex = 0.9)

par(mar = c(.1,.5,.7,.5))
plot0(text = "Mean temperature of the growing season", cex = 1.2)

par(mar = c(1,.1,.1,0.1))
plot0(xaxs = "i")
legend(-1.1, 0, legend = c("Minor", "Moderate", "Major"), cex = 1.1,
       col = "black", 
       lty = 1:3, lwd = 1.4,
       xpd = NA, bty = "n", yjust=.5, seg.len = 2.2)
plot0(xaxs = "i")
legend(-1.1, 0, legend = rev(reg_title), cex = 1.1,
       col = rev(col_reg), lwd = 2.2,
       xpd = NA, bty = "n", yjust=.5, seg.len = 2.2)

dev.off()




### SUPPMAT - State contribution to turnover time ####

# Layout matrix
mat <- matrix(c(1:10), 5)
mat <- rbind(mat, c(11, 11))
mat <- cbind(mat, c(0,0, 12,13,0, 0))


pdf("res/figSupp_contrib2turnover.pdf", width = 6.7, height = 6.8)
#quartz(width = 6.7, height = 6.8)
par(oma = c(0,2,0,0))
layout(mat, heights = c(.12,1,1,1,1,.15), widths = c(1,1,0.8))


# Metrics for natural disturbances
par(mar = c(0,.5,0,.3))
plot0(text = "Natural", cex = 1.2, font = 2, xpd = NA)

par(mar = c(1,1,.1,.3))

for(s in 1:4) {
  axes = 2
  if(s==1) ylim = c(0, 350)
  if(s %in% c(2,3)) ylim = c(0, 100) 
  if(s==4) { 
    ylim = c(0, 200) 
    axes = c(1, 2) 
    }
  plot_index_gradient(index = do.call(cbind, lapply(index_res$soj_st_ls[c(1:3)], function(x) as.numeric(x[,s]))), tp_grad = tp_grad, 
                      main = states[s],
                      unscale = sc_sTP, axes = axes, ylim = ylim,
                      tp_ecoreg = tp_ecoreg)
}
mtext("State contribution to turnover time", 2, outer = T, line = .8, cex = .87)

# Metrics for logging
par(mar = c(0,.5,0,.3))
plot0(text = "Logging", cex = 1.2, font = 2, xpd = NA)

par(mar = c(1,1,.1,.3))
for(s in 1:4) {
  axes = 0
  if(s==1) ylim = c(0, 350)
  if(s %in% c(2,3)) ylim = c(0, 100) 
  if(s==4) { 
    ylim = c(0, 200) 
    axes = 1
  }
  plot_index_gradient(index = do.call(cbind, lapply(index_res$soj_st_ls[c(1,4,5)], function(x) as.numeric(x[,s]))), tp_grad = tp_grad,
                      unscale = sc_sTP, axes = axes, ylim = ylim,
                      tp_ecoreg = tp_ecoreg)
}


par(mar = c(.1,.5,.7,.5))
plot0(text = "Mean temperature of the growing season", cex = 1.3, xpd = NA)

par(mar = c(1,.1,.1,0.1))
plot0(xaxs = "i")
legend(-1.1, 0, legend = c("Minor", "Moderate", "Major"), cex = 1.1,
       col = "black", 
       lty = 1:3, lwd = 1.4,
       xpd = NA, bty = "n", yjust=.5, seg.len = 2.2)
plot0(xaxs = "i")
legend(-1.1, 0, legend = rev(reg_title), cex = 1.1,
       col = rev(col_reg), lwd = 2.2,
       xpd = NA, bty = "n", yjust=.5, seg.len = 2.2)

dev.off()



### SUPPMAT - State ontribution to entropy ####

pdf("res/figSupp_contrib2entropy.pdf", width = 6.7, height = 6.8)
#quartz(width = 6.7, height = 6.8)
par(oma = c(0,2,0,0))
layout(mat, heights = c(.12,1,1,1,1,.15), widths = c(1,1,0.8))


# Metrics for natural disturbances
par(mar = c(0,.5,0,.3))
plot0(text = "Natural", cex = 1.2, font = 2, xpd = NA)

par(mar = c(1,1,.1,.3))

for(s in 1:4) {
  axes = 2
  if(s==4) axes = c(1, 2) 

  plot_index_gradient(index = do.call(cbind, lapply(index_res$entropy_st_ls[c(1:3)], function(x) as.numeric(x[,s]))), tp_grad = tp_grad, 
                      ylab = NULL, xlab = NULL, main = states[s],
                      unscale = sc_sTP, axes = axes, ylim = c(0,1),
                      tp_ecoreg = tp_ecoreg)
}
mtext("State contribution to entropy", 2, outer = T, line = .8, cex = .87)

# Metrics for logging
par(mar = c(0,.5,0,.3))
plot0(text = "Logging", cex = 1.2, font = 2, xpd = NA)

par(mar = c(1,1,.1,.3))
for(s in 1:4) {
  axes = 0
  if(s==4) axes = 1
  plot_index_gradient(index = do.call(cbind, lapply(index_res$entropy_st_ls[c(1,4,5)], function(x) as.numeric(x[,s]))), tp_grad = tp_grad, 
                      ylab = NULL, xlab = NULL,
                      unscale = sc_sTP, axes = axes, ylim = c(0,1),
                      tp_ecoreg = tp_ecoreg)
}


par(mar = c(.1,.5,.7,.5))
plot0(text = "Mean temperature of the growing season", cex = 1.3, xpd = NA)

par(mar = c(1,.1,.1,0.1))
plot0(xaxs = "i")
legend(-1.1, 0, legend = c("Minor", "Moderate", "Major"), cex = 1.1,
       col = "black", 
       lty = 1:3, lwd = 1.4,
       xpd = NA, bty = "n", yjust=.5, seg.len = 2.2)
plot0(xaxs = "i")
legend(-1.1, 0, legend = rev(reg_title), cex = 1.1,
       col = rev(col_reg), lwd = 2.2,
       xpd = NA, bty = "n", yjust=.5, seg.len = 2.2)

dev.off()

library(popdemo)

p =pmatrix.msm(msm_glb, t = 1, covariates = c(covar_d[[1]], sTP = tp_grad[20]), ci = "none")
class(p)="matrix"
convt(p)

dr(p, return.time=T)

KeyfitzD(p, vector =c(.25,.25,.25,.25))

         