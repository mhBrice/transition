### Plot best model ####

### PACKAGES ####
library(graphicsutils)
library(dplyr)
library(msm)
library(sf)
library(knitr)
library(kableExtra)

### FUNCTIONS ####

source('R/functions/plot_msm.R')

### DATA ####

source('R/functions/prep_data.R')

# Load msm results

load("res/msm_all75.rda")
msm_all <- msm_all75
msm_glb <- msm_all[["msm_glb"]]


### Estimated ratio of transition intensities
envmean <- aggregate(states_ba[,c("sTP", "sCMI", "DRAIN", "PH_HUMUS")], 
                     by = list(states_ba$ecoreg3), mean)
mixed_mean <- envmean[3,-1]

covar_nat <- list(c(mixed_mean), 
                  c(natural1 = 1, mixed_mean),
                  c(natural2 = 1, mixed_mean))

covar_log <- list(c(mixed_mean), 
                  c(logging1 = 1, mixed_mean),
                  c(logging2 = 1, mixed_mean))

# Probability matrix
covar <- list(c(mixed_mean), 
              c(natural1 = 1, mixed_mean),
              c(natural2 = 1, mixed_mean), 
              c(logging1 = 1, mixed_mean),
              c(logging2 = 1, mixed_mean))

p_list <- lapply(covar, 
                 function(x) pmatrix.msm(msm_glb, t=10, 
                                         covariates = as.list(x), 
                                         ci = "normal"))


### TABLE OF MODEL RESULTS ####


stats_msm <- lapply(msm_all, 
                    function(x) c("Covariates" = NA,
                                  "Number of parameters" = x$paramdata$npars-length(x$fixedpars),
                                  "-2 Log-likelihood" = round(x$minus2loglik,1), 
                                  "Delta AIC" = round(AIC(x)-AIC(msm_all$msm_glb),1), 
                                  "LR test" = NA)) %>%
  do.call(rbind,.) %>% as.data.frame(.)

rownames(stats_msm) <- c("Baseline", "Climate", "Soil", "Disturbances", "Full")

stats_msm[,"Covariates"] <- c("Intercept",
                              paste("Temperature,", "CMI"),
                              paste("Drainage,", "pH"),
                              paste("Natural,", "Logging"),
                              "All")
stats_msm <- stats_msm[order(stats_msm$`Delta AIC`, decreasing = TRUE),]
stats_msm$`LR test` <- c("---", rep("< 0.001", 4))

kable(stats_msm, format = "latex",
      booktabs = T, linesep = "") %>% 
  column_spec(1, bold = TRUE) %>% 
  row_spec(c(0,5), bold = TRUE) 


### PLOT COEFFICIENTS BEST MODEL ####

varnames <- c("Temperature", "CMI", 
              "Drainage", "pH", 
              "Natural 1", "Natural 2", 
              "Logging 1", "Logging 2")

pdf("res/fig4_HR.pdf",
    width = 7, height = 6)
#quartz(width = 7, height = 6)
plot_risk(msm_glb, varnames = varnames)
dev.off()



### Plot transition matrix #####

dist_title <- c("Minor", "Moderate natural", "Major natural", 
                "Moderate logging", "Major logging")

mat <- matrix(c(0,1,1,0,
                2,2,4,4,
                3,3,5,5),3,byrow = T)

pdf("res/fig5_pmatrix.pdf", width = 4.2, height = 6.2)
#quartz(width = 4.2, height = 6.2)
layout(mat)
par(mar=c(.5,2,3.5,1))
for(i in 1:length(p_list)) {
  tr <- p_list[[i]]
  tr <- tr[["estimates"]]
  #dimnames(tr) <- list(states,states)
  if(i == 1) labels = TRUE else labels = FALSE
  plot_trans(pmat=tr, states_lab = c("B", "M", "P", "T"), 
             labels = labels, main = dist_title[i])
}
dev.off()

### PLOT TRANSITION PROBABILITY ####

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


