### CONTINUOUS TIME MARKOV MODEL ####

### PACKAGES ####

library(msm)
library(dplyr)
library(data.table)

### FUNCTIONS ####

source("R/functions/cross_validation.R")

### DATA ####
# source('R/scripts/1_dataFormatting_env.R')
# source('R/scripts/2_dataFormatting_transition.R')

source('R/scripts/3_prep_trans.R')

states_ba$states_num = states_ba$states_num75

states_ba$states = states_ba$states_ba75

#levels(states_ba$TYPEHUMUS) = c("MU", "MD", "MR", "SO", "SO")
states_ba <- states_ba %>%
  filter(TYPEHUMUS %in% c("MU", "MD", "MR", "TO")) %>% # Subset plots for humus type
  mutate(TYPEHUMUS = droplevels(TYPEHUMUS))

# voisinage -> calculer % de chaque etat par région
# factor is better than continuous variable!
# try opt.method="SANN" to find true maximum (minimum): simulated annealing is slower as it requires many evaluations of the likelihood surface to arrive at final estimates, can be used. It handles the problem by periodically making random jumps to new parameter values as part of the routine. 


# changer définition des états = tester boréal pure, diviser mixte et augmenter le seuil de B et T...
# comparer matrix de P observé vs modélisé

(st_table <- statetable.msm(states_num, plot_id, data = states_ba))

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
# msm_c - Climate model
# msm_d - Disturbance model (natural + harvest)
# msm_dlag - Disturbance model (natural + harvest)
# msm_s - Soil model
# msm_glb - Complete model
# ---------------------------------------------------------------------#

# strata <- states_ba %>% group_by(plot_id) %>% slice(1) %>% ungroup() %>% select(ecoreg3)
# id <- unique(states_ba$plot_id)
# 
# fold <- kfold(strata = strata$ecoreg3, id = id)
# 
# states_trans <- to_trans(states_ba)
# 
# x=lapply(fold, function(x) check_fold(x, data = states_ba, data_trans=states_trans))

# saveRDS(fold, "res/fold.rds")
fold <- readRDS("res/fold.rds")

cv_msm0 <- cv_msm(data = states_ba, fold = fold, Q = Q, 
                  covar_form = NULL, covar_names = NULL)
### msm0 - Null model ####

msm0 <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
                 qmatrix = Q, 
                 gen.inits = TRUE,
                 obstype = 1, 
                 control = list(trace = 1, maxit = 5000, fnscale = 30800),
                 opt.method = "optim")



### msm_c - Climate model ####


covariates_c <- make_forms(covar = c("sTP", "CMI"), covar_p = ~1)

cv_msm_c <- cv_msm(data = states_ba, fold = fold, Q = Q, 
                  covar_form = covariates_c, covar_names = c("sTP", "CMI"))

msm_c <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
            qmatrix = Q, 
            gen.inits=TRUE,
            obstype = 1, 
            control = list(trace=1, maxit=5000, fnscale=28990),
            opt.method = "optim", 
            covariates = covariates_c)

1 - msm_c$minus2loglik/msm0$minus2loglik
lrtest.msm(msm0, msm_c)


### msm_s - Soil model ####

covariates_s <- make_forms(covar = "TYPEHUMUS", covar_p = ~1)

cv_msm_s <- cv_msm(data = states_ba, fold = fold, Q = Q, 
                   covar_form = covariates_s, covar_names = "TYPEHUMUS")

msm_s <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
             qmatrix = Q, 
             gen.inits=TRUE,
             obstype = 1, 
             control = list(trace=1, maxit=5000, fnscale=28990),
             opt.method = "optim", 
             covariates = covariates_s)

1 - msm_s$minus2loglik/msm0$minus2loglik
lrtest.msm(msm0, msm_s)



### msm_d - Disturbance model ####

covariates_d <- make_forms(covar = c("natural", "logging"))

cv_msm_d <- cv_msm(data = states_ba, fold = fold, Q = Q, 
                   covar_form = covariates_d, covar_names = c("natural", "logging"))


msm_d <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
              qmatrix = Q, 
              gen.inits=TRUE,
              obstype = 1, 
              control = list(trace=1, maxit=5000, fnscale=27990),
              opt.method = "optim", 
              covariates = covariates_d)

1 - msm_d$minus2loglik/msm0$minus2loglik

lrtest.msm(msm0, msm_c, msm_d)


### msm_glb - Complete model ####

covariates_glb <- make_forms(covar = c("sTP", "CMI", "TYPEHUMUS",  "natural", "logging"),
           covar_p = c("natural", "logging"))


cv_msm_glb <- cv_msm(data = states_ba, fold = fold, Q = Q, 
                   covar_form = covariates_glb, 
                   covar_names = c("sTP", "CMI", "TYPEHUMUS",  "natural", "logging"))


msm_glb <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
                 qmatrix = Q, 
                 gen.inits = TRUE,
                 obstype = 1, 
                 control = list(trace=1, maxit=5000, fnscale=28800),
                 opt.method = "optim", 
                 covariates = covariates_glb)

1 - msm_glb$minus2loglik/msm0$minus2loglik

### Save all models ####
msm_all75 <- list(msm0 = msm0, 
                msm_c = msm_c, 
                msm_s = msm_s, 
                msm_d = msm_d, 
                msm_glb = msm_glb)
save(msm_all75, file = "res/msm_all75.rda")

### Save all cv models ####
cv_msm_all75 <- list(cv_msm0 = cv_msm0, 
                cv_msm_c = cv_msm_c, 
                cv_msm_s = cv_msm_s, 
                cv_msm_d = cv_msm_d, 
                cv_msm_glb = cv_msm_glb)
save(cv_msm_all75, file = "res/cv_msm_all75.rda")



### MODEL COMPARISON ####
# rsq <- round((1 - (msm_glb$minus2loglik-msm_glb$paramdata$nopt)/msm0$minus2loglik)*100, 2)
# 
# lrtest.msm(msm0, msm_c, msm_s, msm_d, msm_glb)
# AIC(msm0, msm_c, msm_s, msm_d, msm_dlag, msm_glb)





