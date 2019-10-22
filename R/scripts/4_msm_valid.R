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
msm_glb <- msm_all75$msm_glb

# Load cross validated msm results

load("res/cv_msm_all75_drainph.rda")

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
  


pdf("res/fig3_cv_LS.pdf", width = 4.2, height = 3.3)
#quartz(width = 4.2, height = 3.3)
par(mar = c(1.5, 2.8, .5, 1.2))
plot_score(logscore_cv, ylab = expression("Logarithmic Score [0,"*infinity*"["), ylim = c(0.15,.35),
           stats = signif(logscore_summ, 3), text.width = 1.2)
dev.off()

LSS_cv <- lapply(logscore_cv, function(x) 1-x/logscore_cv$cv_msm0)
LSS_cv$cv_msm0=NULL

pdf("res/fig3_cv_LSS.pdf", width = 4, height = 3.3)
#quartz(width = 4, height = 3.3)
par(mar = c(1.5, 2.8, .5, .5))
plot_score(LSS_cv, ylab = expression("Logarithmic Skill Score (1-LS/LS"["baseline"]*")"), 
           mod_names = mod_names[-1], col_mod = col_mod[-1], ylim = c(0,.25),
           stats = signif(logscore_summ[-1], 3), text.width = 1.2)
dev.off()

########################
xy <- st_read("data/plot_xy32198_may2018.gpkg") %>% filter(ID_PE %in% states_ba$ID_PE)
st_crs(xy) <- 32198
ecoreg_df <- readRDS("data/ecoreg_df.RDS") %>% filter(ID_PE %in% states_ba$ID_PE)


brier_plot <- as.data.frame(aggregate(brier_cv$cv_msm_glb$brier, by = list(brier_cv$cv_msm_glb$plot_id),mean))

boxplot(xz2$brier~xz2$ecoreg, outline=F)

xz2 <- xy %>% mutate(brier = (brier_plot$x)) %>% mutate(ecoreg = (ecoreg_df$ecoreg11)) #%>% filter(brier>1)

myPal <- colorRampPalette(c(alpha("grey",.2), alpha("blue",.4),alpha("orange",.4), alpha("red",.4)), alpha=T)
plot(xz2["brier"], pal = myPal, nbreaks = 3, pch = 19, cex = .5)

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

pdf("res/fig3_cv_auc.pdf", width = 4, height = 3.3)
#quartz(width = 4, height = 3.3)
par(mar = c(1.5, 2.5, .5, .3))
plot_score(pair_auc_summ, states = st_auc, 
           ylab = "Pairwise AUCs", stats = gener_auc_summ, ylim = c(.8,1))
dev.off()


### TABLE OF MODEL RESULTS ####


stats_msm <- lapply(msm_all75, 
                    function(x) c("Covariates" = NA,
                                  "Nb of parameters" = x$paramdata$npars-length(x$fixedpars),
                                  "-2 Log-likelihood" = round(x$minus2loglik,1), 
                                  "Delta AIC" = round(AIC(x)-AIC(msm_all75$msm_glb),1), 
                                  "LR test" = "< 0.001",
                                  "mAUC" = NA,
                                  "Log Score" = NA)) %>%
  do.call(rbind,.) %>% as.data.frame(.)

rownames(stats_msm) <- c("Baseline", "Climate", "Soil", "Disturbances", "Full")

stats_msm[,"Covariates"] <- c("Intercept",
                              paste("Temperature,", "CMI"),
                              paste("Drainage,", "pH"),
                              paste("Natural,", "Harvesting"),
                              "All")
stats_msm[,"mAUC"] <- gener_auc_summ
stats_msm[,"Log Score"] <- round(logscore_summ,3)

kable(stats_msm[,5], format = "markdown")





# 
# ### CHI2 ON TRANSITION MATRIX ####
# 
# G_stats <- function(observ, expect) {
#   g <- 0
#   for (i in 1:nrow(observ)) {
#     for (j in 1:ncol(observ)) {
#       if (observ[i, j] != 0 & expect[i, j] != 0) 
#         g <- g + observ[i, j] * log(observ[i, j]/expect[i, j])
#     }
#   }
#   g <- 2*g
#   g
# }
# 
# Chi_stats <- function(observ, expect) {
#   g <- 0
#   for (i in 1:nrow(observ)) {
#     for (j in 1:ncol(observ)) {
#       if (expect[i, j] != 0) 
#         g <- g +((observ[i, j] - expect[i, j])^2)/expect[i, j]
#     } 
#   }
#   g
# }
# perm <- 1000
# chi_perm <- list()
# g_perm <- list()
# for(cv in names_cvmod) {
#   
#   g_f <- chi_f <- c()
#   for(f in 1:n_fold) {
#     mod_f <- cv_pred_valid[[cv]][[f]]
#     
#     # Observed transition matrix
#     obs_trmat <- table(mod_f$from, mod_f$to)
#     
#     p_tmp <- apply(mod_f[,states], 1, function(x) sample(states, size = perm, replace = T, prob = x))
#     p_hard <- as.data.frame(t(p_tmp))
#     
#     pred_trmat <- apply(p_hard, 2, function(x) table(mod_f$from, x))
#     
#     chi_tmp <- apply(pred_trmat, 2, function(x) Chi_stats(obs_trmat, matrix(x,4)))
#     chi_f <- c(chi_f, chi_tmp)
#     
#     g_tmp <- apply(pred_trmat, 2, function(x) G_stats(obs_trmat, matrix(x,4)))
#     g_f <- c(g_f, g_tmp)
#   }
#   chi_perm[[cv]] <- chi_f
#   g_perm[[cv]] <- g_f
# 
# }
# 
# colMeans(chi_perm2)
# colMeans(apply(chi_perm2, 2, function(x) 1-(x/chi_perm2[,1])))
# 
# chi_perm2 = do.call(cbind, chi_perm)
# quartz()
# boxplot(abs(chi_perm2))
# 
# 
# 
# ### Hard classifier ####
# 
# # From probability to state predictions = translate P into number of times a state was picked over 100 trials
# rep <- 100
# 
# # Observed
# obs_rep <- cbind.data.frame(from = rep(from, each = rep), 
#                             to = rep(to, each = rep))
# 
# 
# # Expected
# p_hard <- list()
# for(i in names_mod) {
#   p_tmp <- apply(p_df[[i]], 1, function(x) sample(states, size = rep, replace = T, prob = x))
#   p_hard[[i]] <- as.factor(as.vector(p_tmp))
# }
# 
# 
# 
# ### 1. Accuracy ####
# 
# (accu <- lapply(p_hard, function(x) calculate.accuracy(predictions = x, ref.labels = obs_rep$to)))
# 
# ### 2. Micro and macro averages of the F1-score from confusion matrix ####
# 
# ### One-vs-all confusion matrices
# cm_all <- list()
# for(i in 1:n_mod) {
#   p_tmp <- p_hard[[i]]
#   cm <- list()
#   
#   for (s in states) {
# 
#     # in the i-th iteration, use the i-th class as the positive class
#     cm_tmp <- apply(p_tmp, 2, function(x) confusionMatrix(as.factor(x), to, 
#                               positive = s))
#     
#     cm_accurary <- lapply(cm_tmp, function(x) x$overall[1])
#     cm_metrics <- lapply(cm_tmp, function(x) x$byClass[, metrics])
#     cm[[s]] <- list(cm_accurary = cm_accurary, cm_metrics = cm_metrics)
#     # cm[[s]] <- confusionMatrix(as.factor(p_tmp), obs_rep$to, 
#     #                            positive = positive.class)
#   }
#   cm_all[[names_mod[i]]] <- cm
# }
# 
# apply(simplify2array(cm_all$msm_glb$Temperate$cm_metrics), 1:2, mean)
# 
# mean(unlist(cm_all$msm_glb$Temperate$cm_metrics))
# 
# metrics <- c("Precision", "Recall")
# lapply(cm_all, function(x) x[[1]]$byClass[, metrics])
# 
# # Micro F1
# (micro.f1 <- lapply(cm_all, get.micro.f1))
# 
# 
# # Macro F1
# (macro.f1 <- lapply(cm_all, get.macro.f1))
# 
# 
# ### Soft classifier ####
# 
# 
# 
# # Radar chart
# # quartz()
# # par(mfrow=c(2,2), mar = c(1,2,1,5), oma = c(5,0,0,0))
# # for(s in states) {
# #   gener_auc <- unlist(lapply(multi_auc2[[s]], function(x) x[1]))
# #   pair_auc <- lapply(multi_auc2[[s]], function(x) attr(x, "pair_AUCs"))
# #   
# #   plot0(xlim = c(1, 6), ylim = c(.4, 1), grid.col = "grey", yaxs = "i")
# #   axis(2, tick = F, line = -1, las = 1, cex.axis = .8)
# #   mtext("Pairwise AUCs", 2, line = 1.7)
# #   if(s %in% states[3:4]) text(1:6, .4, names(pair_auc$msm0)[ord], cex = 0.9,
# #                               srt = 90, xpd = NA, adj = 1)
# #   
# #   
# #   for(i in 1:n_mod) {
# #     x <- pair_auc[[i]][ord]
# #     points(x, type = "b", col = colors_border[i], pch = 19, cex = .8)
# #     
# #   }
# #   
# #   legend(x=5.7, y = 1.05, legend = rev(paste(names(pair_auc), round(gener_auc,3))), bty = "n", pch = 20 , 
# #          col = rev(colors_border), 
# #          cex = 1, pt.cex=1, xpd = NA)
# # }
# 
# 
# 
# ### 4. Chi2 ####
# 
# xy <- st_read("data/plot_xy32198_may2018.gpkg") %>% filter(ID_PE %in% states_ba$ID_PE)
# st_crs(xy) <- 32198
# 
# o_df <- Dummy(to, "full")
# 
# chi_ind <- list()
# chi_sum <- list()
# for(i in names_mod) {
#   chi_ind[[i]] <- chi_contrib(observ = o_df, expect = p_df[[i]])
#   chi_sum[[i]] <- sum(chi_ind[[i]])
# }
# xz2 <- xy %>% mutate(chi = (chi_ind$msm_glb)) %>% filter(to!=from) %>% filter(chi<10)
# 
# plot(xz2["chi"])
# 
# # Permutation within rows of observed df
# 
# o_perm <- vegan::permatswap(o_df, fixedmar = "both", shuffle = "samp", times = 999)
# 
# chi_ind_perm <- list()
# chi_sum_perm <- list()
# for(i in names_mod) {
#   chi_ind_perm[[i]] <- lapply(o_perm$perm, function(x) chi_contrib(observ = x, expect = p_df[[i]]))
#   chi_sum_perm[[i]] <- do.call(cbind, lapply(chi_ind_perm[[i]], sum))
# }
# 
# 
# hist(chi_sum_perm$msm0, breaks = 30, xlim = range(chi_sum_perm$msm0, chi_sum$msm0), main = "True and permutation Chi2")
# abline(v = chi_sum$msm0, col = "red3", lwd = 1.2)
# 
# ## compute permutation p-value
# pval.row <- (sum(chi_sum_perm <= chi_sum) + 1) / (length(chi_sum_perm) + 1)
# 
# x=as.data.frame(p_df-p_df0)
# x$pred <- states[apply(x, 1, which.max)]
# 
# x=cbind.data.frame(from,to, x)
# 
# x2 = x[which(x$from != x$to),]
# 
# cbind(table(x$pred), table(x$to))
# 
# (Xsq = chisq.test(table(x$pred), p = table(x$to)/length(x$to)))
# 
# 
# 
