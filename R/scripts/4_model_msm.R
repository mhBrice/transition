### CONTINUOUS TIME MARKOV MODEL ####

### PACKAGES ####

library(msm)

library(graphicsutils)
library(dplyr)

### FUNCTIONS ####

source('R/functions/plot_msm.R')

### DATA ####
# source('R/scripts/1_dataFormatting_env.R')
# source('R/scripts/2_dataFormatting_transition.R')


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

### msm_c - Climate model ####

covar <- c("sTP", "CMI")
covar_p <- ~1

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

msm_c <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
            qmatrix = Q, 
            gen.inits=TRUE,
            obstype = 1, 
            control = list(trace=1, maxit=5000, fnscale=28990),
            opt.method = "optim", 
            covariates = covariates)

1 - msm_c$minus2loglik/msm0$minus2loglik
lrtest.msm(msm0, msm_c)


### msm_r - Ecoregion model ####

covar <- c("ecoreg3")
covar_p <- ~1

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

msm_r <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
            qmatrix = Q, 
            gen.inits=TRUE,
            obstype = 1, 
            control = list(trace=1, maxit=5000, fnscale=28990),
            opt.method = "optim", 
            covariates = covariates)

1 - msm_r$minus2loglik/msm0$minus2loglik
lrtest.msm(msm0, msm_r)

### msm_s - Soil model ####

covar <- "TYPEHUMUS"
covar_p <- ~1

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

msm_s <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
             qmatrix = Q, 
             gen.inits=TRUE,
             obstype = 1, 
             control = list(trace=1, maxit=5000, fnscale=28990),
             opt.method = "optim", 
             covariates = covariates)

1 - msm_s$minus2loglik/msm0$minus2loglik
lrtest.msm(msm0, msm_s)



### msm_d - Disturbance model ####

covar <- c("natural", "logging")
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

msm_d <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
              qmatrix = Q, 
              gen.inits=TRUE,
              obstype = 1, 
              control = list(trace=1, maxit=5000, fnscale=27990),
              opt.method = "optim", 
              covariates = covariates)

1 - msm_d$minus2loglik/msm0$minus2loglik

lrtest.msm(msm0, msm_c, msm_d)


### msm_dlag - Disturbance model with lag for other transition than pioneer ####

covar <- c("natural_lag", "logging_lag") 
covar_p <- c("natural", "logging") 

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

msm_dlag <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
            qmatrix = Q, 
            gen.inits=TRUE,
            obstype = 1, 
            control = list(trace=1, maxit=5000, fnscale=27990),
            opt.method = "optim", 
            covariates = covariates)

1 - msm_dlag$minus2loglik/msm0$minus2loglik

lrtest.msm(msm0, msm_c, msm_cc, msm_d, msm_dlag)

### msm_glb - Complete model ####

covar <- c("sTP", "CMI", "TYPEHUMUS", "ecoreg3", "natural", "logging")
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



msm_glb <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
                 qmatrix = Q, 
                 gen.inits = TRUE,
                 obstype = 1, 
                 control = list(trace=1, maxit=5000, fnscale=28800),
                 opt.method = "optim", 
                 covariates = covariates)

1 - msm_glb$minus2loglik/msm0$minus2loglik

msm_all <- list(msm0 = msm0, 
                msm_c = msm_c, 
                msm_cc = msm_cc, 
                msm_s = msm_s, 
                msm_d = msm_d, 
                msm_dlag = msm_dlag, 
                msm_glb = msm_glb)
save(msm_all, file = "res/msm_all.rda")



### MODEL COMPARISON ####
rsq <- round((1 - (msm_glb$minus2loglik-msm_glb$paramdata$nopt)/msm0$minus2loglik)*100, 2)

lrtest.msm(msm0, msm_c, msm_cc, msm_s, msm_d, msm_dlag, msm_glb)
AIC(msm0, msm_c, msm_cc, msm_s, msm_d, msm_dlag, msm_glb)





