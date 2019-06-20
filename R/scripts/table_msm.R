### TABLE OF MODEL RESULTS ####


### PACKAGES ####

library(msm)

library(graphicsutils)
library(dplyr)
library(knitr)

### DATA ####

load("res/msm_all.rda")

stats_msm <- lapply(msm_all, 
                    function(x) c("Covariates" = NA,
                                  "Nb of parameters" = x$paramdata$npars,
                                  "Delta AIC" = AIC(x)-AIC(msm_all$msm5), 
                                  "-2 Log-likelihood" = x$minus2loglik, 
                                  "McFadden R2" = (1 - x$minus2loglik/msm_all$msm0$minus2loglik)*100)) %>%
  do.call(rbind,.) %>% as.data.frame(.)

rownames(stats_msm) <- c("Null", "Climate", "Disturbances", "Soil", "Full")

stats_msm[,"Covariates"] <- c("Intercept",
                              paste("Temperature,", "CMI"),
                              paste("Natural,", "Harvesting"),
                              paste("Humus type"),
                              "All")

kable(stats_msm, digits = 2)
