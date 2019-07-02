### Plot best model ####

### DATA ####

source('R/scripts/3_prep_trans.R')

# Load msm results

load("res/msm_all.rda")

msm_glb <- msm_all[["msm_glb"]]


### PLOT COEFFICIENTS BEST MODEL ####

varnames <- c("Temperature", "CMI", 
              "Natural 1", "Natural 2", 
              "Logging 1", "Logging 2", 
              "Moder", "Mor")

pdf("res/res_state/figx_HR.pdf",
    width = 7, height = 6)
#quartz(width = 7, height = 6)
plot_risk(msm5, varnames = varnames)
dev.off()

quartz(width = 4, height = 7)
plot_risk2(msm5, varnames = varnames)





mean(states_ba$sTP)
mf <- model.matrix(~sTP + CMI + natural + logging + TYPEHUMUS, states_ba)

x=qmatrix.msm(msm5, covariates = list(sTP=mean(mm.msm5[,"sTP"]), CMI=mean(mm.msm5[,"CMI"]),
                                      natural1 = mean(mm.msm5[,"natural1"]),
                                      natural2 = mean(mm.msm5[,"natural2"]),
                                      logging1 = 0,
                                      logging2 = 1,
                                      TYPEHUMUSMD = mean(mm.msm5[,"TYPEHUMUSMD"]),
                                      TYPEHUMUSMR = mean(mm.msm5[,"TYPEHUMUSMR"])))
eig=eigen(t(x$estimates))
eig$vectors[,4]/sum(eig$vectors[,4])


### PREVALENCE OF STATES THROUGH TIME ####

#prevalence.msm(msm5, times = c(10,25,40), covariates = 'mean')



# Standard errors
# par_se<-sqrt(diag(solve(msm5$opt$hessian)))

### PLOT TRANSITION PROBABILITY ####

pmatrix.msm(msm5, t=10, covariates = list(natural="0", logging = "0"))


covar_nat <- list(list(natural1 = 0, natural2 = 0, logging1 = 0, logging2 = 0), 
                  list(natural1 = 1, natural2 = 0, logging1 = 0, logging2 = 0),
                  list(natural1 = 0, natural2 = 1, logging1 = 0, logging2 = 0))

covar_log <- list(list(natural1 = 0, natural2 = 0, logging1 = 0, logging2 = 0), 
                  list(natural1 = 0, natural2 = 0, logging1 = 1, logging2 = 0),
                  list(natural1 = 0, natural2 = 0, logging1 = 0, logging2 = 1))

mm.msm5 <- model.matrix(form_all, states_ba)

mat <- matrix(c(0:14,0,0,15,15,0), 5, 4)
mat <- rbind(mat, c(0,16,16,0))


pdf("res/res_state90/figx_proba.pdf",
    width = 7, height = 6.5)
#quartz(width = 7, height = 6.5)
layout(mat, widths = c(.6,1,1,.5), heights = c(.2,1,1,1,1,.2))

par(mar=c(8,0,1,2))

for(st in states) {
  plot0(fill = "grey95")
  text(.9, 0, labels = paste("From", st), cex = 1.1, font = 2, adj = 1)
}


par(mar=c(0,2,0,.5))
plot0(text = "Natural disturbances", fill = "grey95", font = 2, cex = 1.2)
par(mar=c(1,2,1,.5))
plot_pmatrix(msm5, covar = covar_nat, mm = mm.msm5)

par(mar=c(0,2,0,.5))
plot0(text = "Harvesting", fill = "grey95", font = 2, cex = 1.2)
par(mar=c(1,2,1,.5))
plot_pmatrix(msm5, covar = covar_log, yaxis = F, mm = mm.msm5)

mtext("Probability of transition", 2, line = -9.8, outer = T, cex= .75, font =2)

# Legend
par(mar=c(5,0,4,0))
plot0()
legend("top", legend = states, title = "Transition to", 
       col = st_col, cex = 1.1, lwd = 1.2, bty = "n")
legend("bottom", legend = c("Minor", "Moderate", "Major"), title = "Disturbance", 
       col = "grey15", cex = 1.1, lwd = 1.2, lty = 3:1, bty = "n")

par(mar=c(0,0,1,0))
plot0(text = "Time (Years)", font = 2, cex = 1.1)

dev.off()


