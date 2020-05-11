#### MODEL VALIDATION #####

# Figure S3

### PACKAGES & FUNCTIONS ####

source("R/functions/packages.R")

source('R/functions/multistate_performance.R')

### DATA ####

source('R/functions/prep_data.R')

# Load complete msm results

load("res/msm_all75.rda")
msm_all <- msm_all75
msm_glb <- msm_all$msm_glb

### Q matrix with allowed transition ####

Q  <-  rbind(c(0.7, 0.1, 0.2, 0),
             c(0.1, 0.6, 0.1, 0.2),
             c(0.3, 0.1, 0.5, 0.1),
             c(0, .1, .1, .8))

rownames(Q) <- colnames(Q) <- states

Q.crude  <- crudeinits.msm(states_num ~ year_measured, plot_id, data=states_ba, qmatrix=Q)


### CREATE STRATIFIED FOLDS ####

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

states_trans <- to_trans(states_ba,
                         covar_names = c("sTP", "sCMI", "PH_HUMUS", "DRAIN", "natural", "logging"))

x=lapply(fold, function(x) check_fold(x, data = states_ba, data_trans=states_trans))

saveRDS(fold, "res/fold.rds")

fold <- readRDS("res/fold.rds")

### CROSS-VALIDATION ####

### Null model ####

cv_msm0 <- cv_msm(data = states_ba, fold = fold, Q = Q,
                  covar_form = NULL, covar_names = NULL)

### Climate model ####

covariates_c <- make_forms(covar = c("sTP", "sCMI"), covar_p = ~1)

cv_msm_c <- cv_msm(data = states_ba, fold = fold, Q = Q,
                   covar_form = covariates_c, covar_names = c("sTP", "sCMI"))


### Soil model ####

covariates_s <- make_forms(covar = c("DRAIN","PH_HUMUS"), covar_p = ~1)

cv_msm_s <- cv_msm(data = states_ba, fold = fold, Q = Q,
                   covar_form = covariates_s, covar_names = c("DRAIN","PH_HUMUS"))



### Disturbance model ####

covariates_d <- make_forms(covar = c("natural", "logging"))

cv_msm_d <- cv_msm(data = states_ba, fold = fold, Q = Q,
                   covar_form = covariates_d, covar_names = c("natural", "logging"))



### Complete model ####

# use covinits to help and accelerate model convergence
hr <- hazard.msm(msm_glb)
lhr <- lapply(hr, function(x) log(x[,"HR"]))

covariates_glb <- make_forms(covar = c("sTP", "sCMI", "DRAIN", "PH_HUMUS", "natural", "logging"),
                             covar_p = c("natural", "logging"))

cv_msm_glb <- cv_msm(data = states_ba, fold = fold, Q = Q,
                     covar_form = covariates_glb,
                     covar_names = c("sTP", "sCMI", "DRAIN","PH_HUMUS","natural", "logging"),
                     covinits = lhr[c(4:8)])

### SAVE ####


cv_msm_all75 <- list(cv_msm0 = cv_msm0,
                     cv_msm_c = cv_msm_c,
                     cv_msm_s = cv_msm_s,
                     cv_msm_d = cv_msm_d,
                     cv_msm_glb = cv_msm_glb)

save(cv_msm_all75, file = "res/cv_msm_all75.rda")

### COMPUTE MODEL PERFORMANCE #####

# Load cross validated msm results

load("res/cv_msm_all75.rda")

mod_names <- c("Baseline", "Climate", "Soil", "Disturbances", "Full")

names_cvmod <- names(cv_msm_all75)
n_cvmod <- length(cv_msm_all75)
n_fold <- length(cv_msm_all75$cv_msm0)


cv_msm_train <- list()
cv_pred_valid <- list()
for(cv in names_cvmod) {
  tmp <- cv_msm_all75[[cv]]
  cv_msm_train[[cv]] <- lapply(tmp, function(x) x$msm_train)
  cv_pred_valid[[cv]] <- lapply(tmp, function(x) x$pred_valid)
}

col_mod <- c("grey55", "#FDAE61", "#68340e", "#D7191C",  "#780a56")



### 1. MULTI-CLASS AUC  ####
# generalization from Hand and Till

gener_auc <- pair_auc <- list()
for(cv in names_cvmod) {
  multi_auc <- lapply(cv_pred_valid[[cv]], function(x) multiclass.auc(x[,4:7], x$to))
  gener_auc[[cv]] <- unlist(lapply(multi_auc, function(x) x[1]))
  pair_auc[[cv]] <- lapply(multi_auc, function(x) attr(x, "pair_AUCs"))
}

gener_auc_summ <- signif(apply(do.call(cbind.data.frame, gener_auc), 2, mean),3)

pair_auc_summ <- lapply(pair_auc, function(x) do.call(rbind, x))


### 2. LOGARITHMIC SCORING RULES  ####

logscore_cv <-logscore_glb_cv <- list()

for(cv in names_cvmod) {
  tmp <- cv_pred_valid[[cv]]

  logscore_f <-logscore_glb_f <-  c()

  for(f in 1:n_fold) {
    to <- model.matrix(~tmp[[f]]$to + 0)
    colnames(to) <- states

    pred <- tmp[[f]][,states]

    logscore_f <- rbind(logscore_f, colMeans(logscore_state(to, pred)))

    q <- to*pred
    q <- q[q>0]
    logscore_glb_f <- c(logscore_glb_f, -log(q))

  }

  logscore_cv[[cv]] <- logscore_f
  logscore_glb_cv[[cv]] <- logscore_glb_f

}


lapply(logscore_cv, function(x) apply(x, 2, mean))

logscore_summ <- unlist(lapply(logscore_glb_cv, mean))




### PLOT CV SCORES ####

st_auc <- gsub("[^A-Z/]", "", colnames(pair_auc_summ$cv_msm0))

pdf("res/figS3_cv.pdf", width = 7.7, height = 3.2)
#quartz(width = 7.7, height = 3.2)
par(mfrow = c(1,2),mar = c(1.5, 2.5, .5, .3))
plot_score(pair_auc_summ, states = st_auc,
           ylab = "Pairwise AUCs", stats = gener_auc_summ, ylim = c(.8,1))
mtext("a",3, adj = -.11, line = -.5)
plot_score(logscore_cv, states = c("B", "M", "P", "T"),
           ylab = expression("Logarithmic Score [0,"*infinity*"["),
           ylim = c(0.14,.38),
           stats = signif(logscore_summ, 3), text.width = 1.2)
mtext("b",3, adj = -.11, line = -.5)
dev.off()

