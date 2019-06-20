### PARAMETER ESTIMATION ####


setwd("~/Documents/GitHub/Doctorat/community_transition")

### DATA ####

#source('R/functions/multistate_LL.R')
load("data/data_gensa.rda")

dat <- data_gensa$states_trans
mm <- data_gensa$mm
mm_p <- data_gensa$mm_p
par_init <- data_gensa$par_init
par_lo <- data_gensa$par_lo
par_hi <- data_gensa$par_hi
names_params <- names(par_init)


### PACKAGES & FUNCTIONS ####
library(expm)
library(GenSA)
library(parallel)

source('R/functions/multistate_LL.R')

# Maximum likelihood estimation using GenSA

np <- detectCores()
cl <- makeForkCluster(np)

system.time({
  estim.pars <- GenSA(par = par_init, lower = par_lo, upper = par_hi,
                      fn = multistate_LL, 
                      control = list(verbose = TRUE, smooth = FALSE, 
                                     nb.stop.improvement = 1000, 
                                     maxit = 10000, temperature = 7000), 
                      dat = dat, 
                      mm = mm, mm_p = mm_p, 
                      names_params = names_params) 
})

stopCluster(cl)

names(estim.pars$par) <- names(par_init)
write.table(estim.pars$par,
            file = paste("../estimated_params/GenSA_", Sys.Date(), ".txt", sep=""), 
            quote=FALSE, col.names=FALSE)
save(estim.pars, par_lo, par_hi, 
     file = paste("../estimated_params/GenSA_", Sys.Date(), ".RData", sep=""))


# states_test <- states_trans[1:100,]
# # Maximum likelihood estimation using optim:
# max <- optim(par = par_init, fn = multistate_LL, method = c("Nelder-Mead"),
#              control = list(maxit=50000), hessian=TRUE,
#              dat = dat,
#              mm = mm, mm_p = mm_p, names_params=names(par_init))
# 
# 
# p<-max$par
