### MULTISTATE MODEL #####

setwd("~/Documents/GitHub/Doctorat/community_transition")

### PACKAGES & FUNCTIONS ####
require(dplyr)
require(nnet)

source('R/functions/multistate_prep.R')

## DATA ####

states_trans <- readRDS("data/states_ba_eachtrans.RDS")


### COVARIATES ####

# Scale covariates
var_to_scale <- c("sTP", "sPP", "slope_sTP", "slope_sPP", 
                  "CMI", "slope_CMI",
                  "all_mort_100", "nat_mort_100", "harvest_100",
                  "neigh_T", "neigh_B",
                  "age_mean")

states_trans[, var_to_scale] <- scale(states_trans[, var_to_scale])


### Covariates

covar <- c("sTP", "sPP", 
           #"slope_sTP", "slope_sPP", 
           "nat_mort_100", "harvest_100", 
           "neigh_T", "neigh_B" 
           #"TYPEHUMUS"
           )
covar_p <- c("nat_mort_100", "harvest_100")

if("TYPEHUMUS" %in% covar) {
  states_trans <- states_trans %>% 
    filter(!(TYPEHUMUS %in% c("AN", "NA"))) %>%
    filter(!(is.na(TYPEHUMUS))) %>%
    droplevels()
}

# covariates to try
# disturbance effect: disturb100 vs nat_mort100 + harvest100
# neighborhood effect: proportion of B and T in a region or in a buffer
# soil effect: EPMATORG+PH_HUMUS vs TYPEHUMUS
# time: none vs linear vs gompezt or weibull

# Create model marix

mm_list <- my_mm(covar = list(covar, covar_p), data = states_trans)


### INITIALIZE PARAMETERS ####

# Set Initial parameters using multinomial models
par_init <- init_param(covar = covar, covar_p = covar_p, dat_trans = states_trans)
length(par_init)
range(par_init)

par_lo <- rep(-15, length(par_init))
par_hi <- rep(15, length(par_init))

data_gensa <- list(states_trans = select(states_trans, plot_id:To),
                   mm = mm_list[[1]],
                   mm_p = mm_list[[2]],
                   par_init = par_init,
                   par_lo = par_lo,
                   par_hi = par_hi)

save(data_gensa, file = "data/data_gensa.rda")
