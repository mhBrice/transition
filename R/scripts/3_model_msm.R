### CONTINUOUS TIME MARKOV MODEL ####

### PACKAGES ####

library(msm)
library(dplyr)
library(data.table)
library(DescTools)
library(sf)
library(caret)

### FUNCTIONS ####

source("R/functions/cross_validation.R")

### DATA ####

source('R/functions/prep_data.R')


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



### msm0 - Null model ####

msm0 <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
                 qmatrix = Q, 
                 gen.inits = TRUE,
                 obstype = 1, 
                 control = list(trace = 1, maxit = 5000, fnscale = 30800),
                 opt.method = "optim")


### msm_c - Climate model ####

covariates_c <- make_forms(covar = c("sTP", "sCMI"), covar_p = ~1)

msm_c <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
             qmatrix = Q, 
             gen.inits=TRUE,
             obstype = 1, 
             control = list(trace=1, maxit=5000, fnscale=37000),
             opt.method = "optim", 
             covariates = covariates_c)


1 - msm_c$minus2loglik/msm0$minus2loglik
lrtest.msm(msm0, msm_c)


### msm_s - Soil model ####
covariates_s <- make_forms(covar = c("DRAIN","PH_HUMUS"), covar_p = ~1)


msm_s <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
             qmatrix = Q, 
             gen.inits=TRUE,
             obstype = 1, 
             control = list(trace=1, maxit=5000, fnscale=38000),
             opt.method = "optim", 
             covariates = covariates_s)


1 - msm_s$minus2loglik/msm0$minus2loglik
lrtest.msm(msm0, msm_s)



### msm_d - Disturbance model ####

covariates_d <- make_forms(covar = c("natural", "logging"))

msm_d <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
              qmatrix = Q, 
              gen.inits=TRUE,
              obstype = 1, 
              control = list(trace=1, maxit=5000, fnscale=38000),
              opt.method = "optim", 
              covariates = covariates_d)

# hr <- hazard.msm(msm_d)
# lhr <- lapply(hr, function(x) log(x[,"HR"]))
# 

1 - msm_d$minus2loglik/msm0$minus2loglik

lrtest.msm(msm0, msm_c, msm_d)



### msm_glb - Complete model ####

covariates_glb <- make_forms(covar = c("sTP", "sCMI", "DRAIN", "PH_HUMUS", "natural", "logging"),
           covar_p = c("natural", "logging"))


msm_glb <- msm(states_num ~ year_measured, subject = plot_id, 
               data = states_ba,
               qmatrix = Q, 
               gen.inits = TRUE,
               obstype = 1, 
               control = list(trace=1, maxit=5000, fnscale=36000),
               opt.method = "optim", 
               covariates = covariates_glb)
hr <- hazard.msm(msm_glb)
lhr <- lapply(hr, function(x) log(x[,"HR"]))
#
1 - msm_glb$minus2loglik/msm0$minus2loglik

### Save all models ####
msm_all75 <- list(msm0 = msm0, 
                msm_c = msm_c, 
                msm_s = msm_s, 
                msm_d = msm_d, 
                msm_glb = msm_glb)
save(msm_all75, file = "res/msm_all75.rda")

#################################
### 10-FOLD CROSS-VALIDATION ####
#################################
strata <- data.table(states_ba) 
strata <- strata[ , list(ecoreg = first(ecoreg6), 
                         sTP = first(sTP),
                         sCMI = first(sCMI),
                         PH_HUMUS = first(PH_HUMUS),
                         DRAIN = first(DRAIN),  
                         natural = max(as.numeric(natural)),  
                         logging = max(as.numeric(logging))), 
                  by = plot_id]

id <- strata$plot_id
strata$strata <- paste0(strata$ecoreg, strata$sTP, strata$sCMI, strata$PH_HUMUS,
                          strata$DRAIN, strata$natural, strata$logging)
fold <- kfold(strata = strata$strata, id = id, k = 10)

states_trans <- to_trans(states_ba, covar_names = c("sTP", "sCMI", "PH_HUMUS",  "natural", "logging"))
lhr[1:4] = lapply(lhr[1:4], function(x) x[-which(x==0)])

x=lapply(fold, function(x) check_fold(x, data = states_ba, data_trans=states_trans))

saveRDS(fold, "res/fold.rds")
fold <- readRDS("res/fold.rds")
### msm0 - Null model ####
cv_msm0 <- cv_msm(data = states_ba, fold = fold, Q = Q, 
                  covar_form = NULL, covar_names = NULL)

### msm_c - Climate model ####
cv_msm_c <- cv_msm(data = states_ba, fold = fold, Q = Q, 
                   covar_form = covariates_c, covar_names = c("sTP", "sCMI"))


### msm_s - Soil model ####
cv_msm_s <- cv_msm(data = states_ba, fold = fold, Q = Q, 
                   covar_form = covariates_s, covar_names = c("DRAIN","PH_HUMUS"))



### msm_d - Disturbance model ####
cv_msm_d <- cv_msm(data = states_ba, fold = fold, Q = Q, 
                   covar_form = covariates_d, covar_names = c("natural", "logging"))



### msm_glb - Complete model ####

cv_msm_glb <- cv_msm(data = states_ba, fold = fold, Q = Q, 
                     covar_form = covariates_glb, 
                     covar_names = c("sTP", "sCMI", "DRAIN","PH_HUMUS","natural", "logging"),
                     covinits = lhr[c(4:8)])

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





