### CONTINUOUS TIME MARKOV MODEL ####

### PACKAGES & FUNCTIONS ####

source("R/functions/packages.R")

source("R/functions/cross_validation.R")

### DATA ####

source('R/functions/prep_data.R')

# states_ba$states_num=states_ba$states_num85
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
# msm_c - Climate model (sTP + sCMI)
# msm_d - Disturbance model (natural + logging)
# msm_s - Soil model (PH_HUMUS + DRAIN)
# msm_glb - Complete model (sTP + sCMI + PH_HUMUS + DRAIN + natural + logging)
# ---------------------------------------------------------------------#



### Null model ####

msm0 <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
                 qmatrix = Q, 
                 gen.inits = TRUE,
                 obstype = 1, 
                 control = list(trace = 1, maxit = 5000, fnscale = 30800),
                 opt.method = "optim")


### Climate model ####

covariates_c <- make_forms(covar = c("sTP", "sCMI"), covar_p = ~1)

msm_c <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
             qmatrix = Q, 
             gen.inits=TRUE,
             obstype = 1, 
             control = list(trace=1, maxit=5000, fnscale=37000),
             opt.method = "optim", 
             covariates = covariates_c)




### Soil model ####
covariates_s <- make_forms(covar = c("DRAIN","PH_HUMUS"), covar_p = ~1)


msm_s <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
             qmatrix = Q, 
             gen.inits=TRUE,
             obstype = 1, 
             control = list(trace=1, maxit=5000, fnscale=38000),
             opt.method = "optim", 
             covariates = covariates_s)



### Disturbance model ####

covariates_d <- make_forms(covar = c("natural", "logging"))

msm_d <- msm(states_num ~ year_measured, subject = plot_id, data = states_ba,
              qmatrix = Q, 
              gen.inits=TRUE,
              obstype = 1, 
              control = list(trace=1, maxit=5000, fnscale=38000),
              opt.method = "optim", 
              covariates = covariates_d)



### Complete model ####

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


### Save all models ####

msm_all75 <- list(msm0 = msm0, 
                msm_c = msm_c, 
                msm_s = msm_s, 
                msm_d = msm_d, 
                msm_glb = msm_glb)
save(msm_all75, file = "res/msm_all75.rda")
#save(msm_glb, file = "res/msm_glb85.rda")


### MODEL COMPARISON ####

# AIC
lapply(msm_all75, AIC)

# Likelihood ratio tests
lapply(msm_all75, function(x) lrtest.msm(msm_all75$msm0, x))

# Pseudo R2
lapply(msm_all75, function(x) 1 - x$minus2loglik/msm_all75$msm0$minus2loglik) 


