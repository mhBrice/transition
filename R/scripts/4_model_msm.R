### CONTINUOUS TIME MARKOV MODEL ####

### PACKAGES ####

library(msm)

library(graphicsutils)
library(dplyr)

### FUNCTIONS ####

source('R/functions/plot_msm.R')

### DATA ####
source('R/scripts/1_dataFormatting_env.R')
source('R/scripts/2_dataFormatting_transition.R')


source('R/scripts/3_prep_trans.R')


# voisinage -> calculer % de chaque etat par région
# factor is better than continuous variable!
# try opt.method="SANN" to find true maximum (minimum): simulated annealing is slower as it requires many evaluations of the likelihood surface to arrive at final estimates, can be used. It handles the problem by periodically making random jumps to new parameter values as part of the routine. 


# changer définition des états = tester boréal pure, diviser mixte et augmenter le seuil de B et T...
# comparer matrix de P observé vs modélisé

(stable <- statetable.msm(states_num, plot_id, data = states_ba))

# Final proportion
colSums(stable)/sum(stable)

# Initial proportion
rowSums(stable)/sum(stable)


### Q matrix with allowed transition ####

Q  <-  rbind(c(0.7, 0.1, 0.2, 0),
             c(0.1, 0.6, 0.1, 0.2),
             c(0.3, 0.1, 0.5, 0.1),
             c(0, .1, .1, .8))

rownames(Q) <- colnames(Q) <- states

Q.crude  <- crudeinits.msm(states_num ~ year_measured, plot_id, data=states_ba, qmatrix=Q)


### CANDIDATE MODELS ####

# ---------------------------------------------------------------------#
# msm0 - Null model
# msm1 - Climate model
# msm2 - Climate and climate change model
# msm3 - Disturbance model (all, natural + harvest, histo vs recent)
# msm4 - Soil model
# msm5 - Complete model
# ---------------------------------------------------------------------#

### msm0 - Null model ####

msm0 <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
                 qmatrix = Q, 
                 gen.inits = TRUE,
                 obstype = 1, 
                 control = list(trace = 1, maxit = 5000, fnscale = 30800),
                 opt.method = "optim")

### msm1 - Climate model ####

covar <- c("sTP", "CMI")
covar_p <- covar

form_all <- as.formula(paste0("~ ", paste(covar, collapse = "+")))
form_p <- as.formula(paste0("~ ", paste(covar_p, collapse = "+")))

covariates =  list(
  # From boreal
  "1-2" =  form_all,
  "1-3" =  form_p,
  # From mixed
  "2-1" = form_all,
  "2-3" = form_p,
  "2-4" = form_all,
  # From pioneer
  "3-1" = form_all,
  "3-2" = form_all,
  "3-4" = form_all,
  # From temperate
  "4-2" = form_all,
  "4-3" = form_p
)

msm1.0 <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
            qmatrix = Q, 
            gen.inits=TRUE,
            obstype = 1, 
            control = list(trace=1, maxit=5000, fnscale=28990),
            opt.method = "optim", 
            covariates = covariates)

1 - msm1.0$minus2loglik/msm0$minus2loglik
lrtest.msm(msm0, msm1.0)


### msm2 - Climate and climate change model ####

covar <- c("sTP", "CMI", "TP_xmin")
covar_p <- covar

form_all <- as.formula(paste0("~ ", paste(covar, collapse = "+")))
form_p <- as.formula(paste0("~ ", paste(covar_p, collapse = "+")))

covariates =  list(
  # From boreal
  "1-2" =  form_all,
  "1-3" =  form_p,
  # From mixed
  "2-1" = form_all,
  "2-3" = form_p,
  "2-4" = form_all,
  # From pioneer
  "3-1" = form_all,
  "3-2" = form_all,
  "3-4" = form_all,
  # From temperate
  "4-2" = form_all,
  "4-3" = form_p
)

msm2.3 <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
            qmatrix = Q, 
            gen.inits=TRUE,
            obstype = 1, 
            control = list(trace=1, maxit=5000, fnscale=28990),
            opt.method = "optim", 
            covariates = covariates)

1 - msm2.1$minus2loglik/msm0$minus2loglik
lrtest.msm(msm2, msm2.2)


### msm3 - Disturbance model (all, natural + harvest, histo vs recent) ####

covar <- c("natural_lag", "logging_lag")
covar_p <- covar

form_all <- as.formula(paste0("~ ", paste(covar, collapse = "+")))
form_p <- as.formula(paste0("~ ", paste(covar_p, collapse = "+")))

covariates =  list(
  # From boreal
  "1-2" =  form_all,
  "1-3" =  form_p,
  # From mixed
  "2-1" = form_all,
  "2-3" = form_p,
  "2-4" = form_all,
  # From pioneer
  "3-1" = form_all,
  "3-2" = form_all,
  "3-4" = form_all,
  # From temperate
  "4-2" = form_all,
  "4-3" = form_p
)

msm3.1 <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
            qmatrix = Q, 
            gen.inits=TRUE,
            obstype = 1, 
            control = list(trace=1, maxit=5000, fnscale=27990),
            opt.method = "optim", 
            covariates = covariates)

1 - msm3$minus2loglik/msm0$minus2loglik
1 - msm3.1$minus2loglik/msm0$minus2loglik

lrtest.msm(msm0, msm2, msm3, msm3.1)
AIC(msm0, msm2,msm3)

### msm4 - Soil model ####

covar <- c("TYPEHUMUS")
covar_p <- covar

form_all <- as.formula(paste0("~ ", paste(covar, collapse = "+")))
form_p <- as.formula(paste0("~ ", paste(covar_p, collapse = "+")))

covariates =  list(
  # From boreal
  "1-2" =  form_all,
  "1-3" =  form_p,
  # From mixed
  "2-1" = form_all,
  "2-3" = form_p,
  "2-4" = form_all,
  # From pioneer
  "3-1" = form_all,
  "3-2" = form_all,
  "3-4" = form_all,
  # From temperate
  "4-2" = form_all,
  "4-3" = form_p
)

msm4 <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
            qmatrix = Q, 
            gen.inits=TRUE,
            obstype = 1, 
            control = list(trace=1, maxit=5000, fnscale=28990),
            opt.method = "optim", 
            covariates = covariates)

1 - msm4$minus2loglik/msm0$minus2loglik
lrtest.msm(msm0, msm4, msm4.1)

### msm5 - Complete model ####

covar <- c("sTP", "CMI", "natural", "logging", "TYPEHUMUS")
covar_p <- c("natural", "logging")

form_all <- as.formula(paste0("~ ", paste(covar, collapse = "+")))
form_p <- as.formula(paste0("~ ", paste(covar_p, collapse = "+")))

covariates =  list(
  # From boreal
  "1-2" =  form_all,
  "1-3" =  form_p,
  # From mixed
  "2-1" = form_all,
  "2-3" =  form_p,
  "2-4" = form_all,
  # From pioneer
  "3-1" = form_all,
  "3-2" = form_all,
  "3-4" = form_all,
  # From temperate
  "4-2" = form_all,
  "4-3" = form_p
)

msm5 <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
                 qmatrix = Q, 
                 gen.inits = TRUE,
                 obstype = 1, 
                 control = list(trace=1, maxit=5000, fnscale=28800),
                 opt.method = "optim", 
                 covariates = covariates)


msm_all <- list(msm0=msm0, msm1=msm1, msm3=msm3, msm4=msm4, msm5=msm5)
save(msm_all, file = "res/msm_all.rda")

### MODEL COMPARISON ####
rsq <- round((1 - (msm5$minus2loglik-msm5$paramdata$nopt)/msm0$minus2loglik)*100, 2)

lrtest.msm(msm0, msm1, msm2, msm3, msm4, msm5)
AIC(msm0, msm1, msm2, msm3, msm4, msm5)


### BASELINE HAZARD ####

qmatrix.msm(msm5)

### PLOT COEFFICIENTS BEST MODEL ####

varnames <- c("Temperature", "CMI", 
              "Natural 1", "Natural 2", 
              "Logging 1", "Logging 2", 
              "Moder", "Mor")

pdf("res/res_state90/figx_HR.pdf",
    width = 7, height = 6)
#quartz(width = 7, height = 6)
plot_risk(msm5, varnames = varnames)
dev.off()

quartz(width = 4, height = 7)
plot_risk2(msm5, varnames = varnames)




x=qmatrix.msm(msm5, covariates = "mean")
eig=eigen(t(x$estimates))
eig$vectors[,4]/sum(eig$vectors[,4])

mean(states_ba$sTP)
mm.msm5=model.matrix(~sTP + CMI + natural + logging + TYPEHUMUS, states_ba)

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

prevalence.msm(msm5, times = c(10,25,40), covariates = 'mean')

### Goodness-of-fit test ####

# p0 <- pearson.msm(state.msm0)
# p1 <- pearson.msm(state.msm1)
# p2 <- pearson.msm(state.msm2)
# p3 <- pearson.msm(state.msm)

# Optimal parameters
# par_opt <- msm5$opt$par


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





