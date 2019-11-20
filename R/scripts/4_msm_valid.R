#### MODEL VALIDATION #####

### PACKAGES ####

library(dplyr)
library(msm)

# library(caret)
library(DescTools)
library(sf)
library(scoring)

library(RColorBrewer)
library(graphicsutils)
library(knitr)

source('R/functions/multistate_performance.R')

### DATA ####

source('R/functions/prep_data.R')

# Load complete msm results

load("res/msm_all75.rda")
msm_all <- msm_all75
msm_glb <- msm_all$msm_glb

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

### 1. LIKELIHOOD RATIO TEST + PSEUDO R2 ####

# Complete model
lapply(msm_all75, function(x) lrtest.msm(msm_all75$msm0, x))

lapply(msm_all75, function(x) 1 - x$minus2loglik/msm_all75$msm0$minus2loglik) 


### 2. LOGARITHMIC SCORING RULES  ####

logscore_cv <-logscore_glb_cv <- list()

for(cv in names_cvmod) {
  tmp <- cv_pred_valid[[cv]]
  
  logscore_f <-logscore_glb_f <-  c()

  for(f in 1:n_fold) {
    to <- Dummy(tmp[[f]]$to, "full")
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


ttest_MSE <- list()
for(cv in names_cvmod[-1]) {
  ttest_state <- c()
  for(s in states) {
    tmp <- t.test(logscore_cv[[cv]][,s], logscore_cv$cv_msm0[,s], 
           alternative = "less", paired = T)$p.value
    ttest_state <- c(ttest_state, tmp)

  }
  names(ttest_state) <- states
  ttest_MSE[[cv]] <- ttest_state
}
  




# LSS_cv <- lapply(logscore_cv, function(x) 1-x/logscore_cv$cv_msm0)
# LSS_cv$cv_msm0=NULL

# pdf("res/fig3_cv_LSS.pdf", width = 4, height = 3.3)
# #quartz(width = 4, height = 3.3)
# par(mar = c(1.5, 2.8, .5, .5))
# plot_score(LSS_cv, ylab = expression("Logarithmic Skill Score (1-LS/LS"["baseline"]*")"), 
#            mod_names = mod_names[-1], col_mod = col_mod[-1], ylim = c(0,.25),
#            stats = signif(logscore_summ[-1], 3), text.width = 1.2)
# dev.off()


### 3. AUC  ####
# generalization from Hand and Till

gener_auc <- pair_auc <- list()
for(cv in names_cvmod) {
  multi_auc <- lapply(cv_pred_valid[[cv]], function(x) multiclass.auc(x[,4:7], x$to))
  gener_auc[[cv]] <- unlist(lapply(multi_auc, function(x) x[1]))
  pair_auc[[cv]] <- lapply(multi_auc, function(x) attr(x, "pair_AUCs"))
}

gener_auc_summ <- signif(apply(do.call(cbind.data.frame, gener_auc), 2, mean),3)

pair_auc_summ=lapply(pair_auc, function(x) do.call(rbind, x))

### Plot of pairs AUC ####
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

