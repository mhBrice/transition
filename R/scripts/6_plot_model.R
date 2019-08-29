### Plot best model ####

### PACKAGES ####
library(graphicsutils)
library(dplyr)
library(msm)

### FUNCTIONS ####

source('R/functions/plot_msm.R')

### DATA ####

source('R/scripts/3_prep_trans.R')

# Load msm results

load("res/msm_all75_drainph.rda")

msm_glb <- msm_all75[["msm_glb"]]

### Estimated ratio of transition intensities
x = covar_nat[[2]]
unique(states_ba$DRAIN)
x$DRAIN = -1.1915525
qratio.msm(msm_glb, ind1 = c(2,4), ind2 = c(4,2), covariates =x)
qratio.msm(msm_glb, ind1 = c(3,1), ind2 = c(1,3))
qratio.msm(msm_glb, ind1 = c(1,2), ind2 = c(2,1), covariates = x)
0.0237933/0.0139635

### PLOT COEFFICIENTS BEST MODEL ####

varnames <- c("Temperature", "CMI", 
              "Drainage", "pH", 
              "Natural 1", "Natural 2", 
              "Logging 1", "Logging 2")

pdf("res/fig2_HR.pdf",
    width = 7, height = 6)
#quartz(width = 7, height = 6)
plot_risk(msm_glb, varnames = varnames)
dev.off()

# quartz(width = 4, height = 7)
# plot_risk2(msm5, varnames = varnames)


### PLOT INFLUENCE OF COVARIATE ON STEADY STATE ####



### PREVALENCE OF STATES THROUGH TIME ####

prevalence.msm(msm_glb, times = c(10,25,40), covariates = 'population')


### PLOT TRANSITION PROBABILITY ####

aggregate(states_ba[,c("natural", "logging")], 
          by = list(states_ba$ecoreg3), table)

envmean <- aggregate(states_ba[,c("sTP", "CMI", "DRAIN", "PH_HUMUS")], 
                     by = list(states_ba$ecoreg3), mean)
mixed_mean <- envmean[3,-1]

covar_nat <- list(c(mixed_mean), 
                  c(natural1 = 1, mixed_mean),
                  c(natural2 = 1, mixed_mean))

covar_log <- list(c(mixed_mean), 
                  c(logging1 = 1, mixed_mean),
                  c(logging2 = 1, mixed_mean))



mat <- matrix(c(0:14,0,0,15,15,0), 5, 4)
mat <- rbind(mat, c(0,16,16,0))


pdf("res/figSupp_proba.pdf",
    width = 7, height = 7)
#quartz(width = 7, height = 7)
layout(mat, widths = c(.6,1,1,.45), heights = c(.17,1,1,1,1,.2))

par(mar=c(9,0,1,2))

for(st in states) {
  plot0(fill = "grey95")
  text(.9, 0, labels = paste("From", st), cex = 1.1, font = 2, adj = 1)
}


par(mar=c(0,2,0,0))
plot0(text = "Natural", fill = "grey95", font = 2, cex = 1.2)
par(mar=c(1,2,1,.5))
plot_pmatrix(msm_glb, covar = covar_nat)

par(mar=c(0,2,0,0))
plot0(text = "Logging", fill = "grey95", font = 2, cex = 1.2)
par(mar=c(1,2,1,.5))
plot_pmatrix(msm_glb, covar = covar_log, yaxis = F)

mtext("Probability of transition", 2, line = -9.8, outer = T, cex= .75, font =2)

# Legend
par(mar=c(5,0,4,0))
plot0()
legend("top", legend = states, title = "Transition to", 
       col = st_col, cex = 1.1, lwd = 1.3, bty = "n")
legend("bottom", legend = c("Minor", "Moderate", "Major"), title = "Disturbance", 
       col = "grey15", cex = 1.1, lwd = 1.3, lty = 1:3, bty = "n")

par(mar=c(0,0,1,0))
plot0(text = "Time (Years)", font = 2, cex = 1.1)

dev.off()


prevalence.msm(msm_glb, times = c(0,10,20,30,40), ci = "none", plot = T)
