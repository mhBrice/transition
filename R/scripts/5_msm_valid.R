#### MODEL VALIDATION #####

### PACKAGES ####

library(dplyr)
library(msm)

library(caret)
library(DescTools)
library(sf)

library(RColorBrewer)
library(graphicsutils)
library(fmsb) # radarchart

source('R/functions/multistate_performance.R')

### DATA ####

source('R/scripts/3_prep_trans.R')

# Load msm results

load("res/msm_all85.rda")
load("res/msm_all75.rda")

lapply(msm_all, function(x) 1 - x$minus2loglik/msm_all$msm0$minus2loglik) 
msm_all$msm_glbreg=NULL
names_mod <- names(msm_all)
n_mod <- length(msm_all)





# Test data
covar_names <- attr(terms(msm_all$msm_glb$covariates), "term.labels")
# 
states_ba <- states_ba %>%
  filter(TYPEHUMUS %in% c("MU", "MD", "MR", "TO")) %>% # Subset plots for humus type
  mutate(TYPEHUMUS = droplevels(TYPEHUMUS))

sample_plot <- sample(unique(states_ba$plot_id), 3000)

st_from <- states_ba %>% group_by(plot_id) %>% 
  mutate(t = lead(year_measured) - year_measured) %>% 
  slice(n()-1) %>% dplyr::select(plot_id, from = states_ba75, t, covar_names) %>% 
  filter(plot_id %in% sample_plot)

st_to <- states_ba %>% group_by(plot_id) %>% 
  slice(n()) %>% dplyr::select(plot_id, to = states_ba75) %>% 
  filter(plot_id %in% sample_plot)

to <- st_to$to
from <- st_from$from
t <- st_from$t
covariates <- st_from[,covar_names]

### Make predictions for the testing dataset ####
p_df <- list()
for(i in names_mod) {
  mod <- msm_all[[i]]
  
  covar_names <- attr(terms(mod$covariates), "term.labels")
  
  p_tmp <- msm_pred(mod = mod, covariates = covariates[covar_names], t = t, from = from)
  
  p_df[[i]] <- p_tmp
}


### Hard classifier ####

# From probability to state predictions = translate P into number of times a state was picked over 100 trials
rep <- 100

# Observed
obs_rep <- cbind.data.frame(from = rep(from, each = rep), 
                            to = rep(to, each = rep))


# Expected
p_hard <- list()
for(i in names_mod) {
  p_tmp <- apply(p_df[[i]], 1, function(x) sample(states, size = rep, replace = T, prob = x))
  p_hard[[i]] <- as.factor(as.vector(p_tmp))
}



### 1. Accuracy ####

(accu <- lapply(p_hard, function(x) calculate.accuracy(predictions = x, ref.labels = obs_rep$to)))

### 2. Micro and macro averages of the F1-score from confusion matrix ####

### One-vs-all confusion matrices
cm_all <- list()
for(i in 1:n_mod) {
  p_tmp <- p_hard[[i]]
  cm <- list()
  
  for (s in states) {

    # in the i-th iteration, use the i-th class as the positive class
    cm_tmp <- apply(p_tmp, 2, function(x) confusionMatrix(as.factor(x), to, 
                              positive = s))
    
    cm_accurary <- lapply(cm_tmp, function(x) x$overall[1])
    cm_metrics <- lapply(cm_tmp, function(x) x$byClass[, metrics])
    cm[[s]] <- list(cm_accurary = cm_accurary, cm_metrics = cm_metrics)
    # cm[[s]] <- confusionMatrix(as.factor(p_tmp), obs_rep$to, 
    #                            positive = positive.class)
  }
  cm_all[[names_mod[i]]] <- cm
}

apply(simplify2array(cm_all$msm_glb$Temperate$cm_metrics), 1:2, mean)

mean(unlist(cm_all$msm_glb$Temperate$cm_metrics))

metrics <- c("Precision", "Recall")
lapply(cm_all, function(x) x[[1]]$byClass[, metrics])

# Micro F1
(micro.f1 <- lapply(cm_all, get.micro.f1))


# Macro F1
(macro.f1 <- lapply(cm_all, get.macro.f1))


### Soft classifier ####

### 3. AUC generalization from Hand and Till ####

(multi_auc <- lapply(cv_msm_d, function(x) multiclass.auc(x[,-c(1,2)], x$to)))

multi_auc2 <- list()
for(s in states) {
  l_from <- which(from == s)
  (M <- lapply(p_df, function(x) multiclass.auc(x[l_from,], to[l_from])))
  multi_auc2[[s]] <- M
}

### Plot of pairs AUC ####

gener_auc <- unlist(lapply(multi_auc, function(x) x[1]))
pair_auc <- lapply(multi_auc, function(x) attr(x, "pair_AUCs"))

# Radar chart
dat = do.call(rbind, pair_auc)
dat = rbind(rep(1,6) , rep(0.8,6) , dat)
dat = as.data.frame(dat)


colors_border <- brewer.pal(n_mod,"Spectral")

# quartz()
# radarchart(dat, axistype = 1, 
#            #custom polygon
#            pcol = colors_border, plwd = 1.5, plty = 2,
#            #custom the grid
#            cglcol = "grey", cglty = 1, axislabcol = "grey", 
#            caxislabels = seq(0.8, 1, length.out = 5), cglwd = 0.8,
#            #custom labels
#            vlcex = 0.8)
# legend(x=0.7, y = 1, legend = names(pair_auc), bty = "n", pch = 20 , col = colors_border, 
#       cex = 1, pt.cex=1, title = "Models")

# Parallel plot
ord <- order(pair_auc$msm_glb, decreasing = T)

quartz(width = 7, height = 5)
par(mar = c(7, 3.2, 1, 8))
plot0(xlim = c(1, 6), ylim = c(.7, 1), grid.col = "grey", yaxs = "i")
axis(2, tick = F, line = -1, las = 1, cex.axis = .8)
mtext("Pairwise AUCs", 2, line = 1.7)

for(i in 1:n_mod) {
  x <- pair_auc[[i]][ord]
  points(x, type = "b", col = colors_border[i], pch = 19, cex = .8)
}
text(1:6, 0.7, names(pair_auc$msm0)[ord], cex = 0.9,
     srt = 90, xpd = NA, adj = 1)
legend(x = 6.1, y = 0.9, legend = rev(names(pair_auc)), 
       col = rev(colors_border), lwd =1, cex = .8, bty = "n", xpd = NA)
legend(x = 7.1, y = 0.9, legend = rev(round(gener_auc,3)), cex = .8, bty = "n", xpd = NA)



# Radar chart
# quartz()
# par(mfrow=c(2,2), mar = c(1,2,1,5), oma = c(5,0,0,0))
# for(s in states) {
#   gener_auc <- unlist(lapply(multi_auc2[[s]], function(x) x[1]))
#   pair_auc <- lapply(multi_auc2[[s]], function(x) attr(x, "pair_AUCs"))
#   
#   plot0(xlim = c(1, 6), ylim = c(.4, 1), grid.col = "grey", yaxs = "i")
#   axis(2, tick = F, line = -1, las = 1, cex.axis = .8)
#   mtext("Pairwise AUCs", 2, line = 1.7)
#   if(s %in% states[3:4]) text(1:6, .4, names(pair_auc$msm0)[ord], cex = 0.9,
#                               srt = 90, xpd = NA, adj = 1)
#   
#   
#   for(i in 1:n_mod) {
#     x <- pair_auc[[i]][ord]
#     points(x, type = "b", col = colors_border[i], pch = 19, cex = .8)
#     
#   }
#   
#   legend(x=5.7, y = 1.05, legend = rev(paste(names(pair_auc), round(gener_auc,3))), bty = "n", pch = 20 , 
#          col = rev(colors_border), 
#          cex = 1, pt.cex=1, xpd = NA)
# }



### 4. Chi2 ####

xy <- st_read("data/plot_xy32198_may2018.gpkg") %>% filter(ID_PE %in% states_ba$ID_PE)
st_crs(xy) <- 32198

o_df <- Dummy(to, "full")

chi_ind <- list()
chi_sum <- list()
for(i in names_mod) {
  chi_ind[[i]] <- chi_contrib(observ = o_df, expect = p_df[[i]])
  chi_sum[[i]] <- sum(chi_ind[[i]])
}
xz2 <- xy %>% mutate(chi = (chi_ind$msm_glb)) %>% filter(to!=from) %>% filter(chi<10)

plot(xz2["chi"])

# Permutation within rows of observed df

o_perm <- vegan::permatswap(o_df, fixedmar = "both", shuffle = "samp", times = 999)

chi_ind_perm <- list()
chi_sum_perm <- list()
for(i in names_mod) {
  chi_ind_perm[[i]] <- lapply(o_perm$perm, function(x) chi_contrib(observ = x, expect = p_df[[i]]))
  chi_sum_perm[[i]] <- do.call(cbind, lapply(chi_ind_perm[[i]], sum))
}


hist(chi_sum_perm$msm0, breaks = 30, xlim = range(chi_sum_perm$msm0, chi_sum$msm0), main = "True and permutation Chi2")
abline(v = chi_sum$msm0, col = "red3", lwd = 1.2)

## compute permutation p-value
pval.row <- (sum(chi_sum_perm <= chi_sum) + 1) / (length(chi_sum_perm) + 1)

x=as.data.frame(p_df-p_df0)
x$pred <- states[apply(x, 1, which.max)]

x=cbind.data.frame(from,to, x)

x2 = x[which(x$from != x$to),]

cbind(table(x$pred), table(x$to))

(Xsq = chisq.test(table(x$pred), p = table(x$to)/length(x$to)))

# CHi2 on transition matrix

# Observed transition matrix
obs_trmat <- table(from, to)
obs_trmat2 <- car::logit(obs_trmat/rowSums(obs_trmat),adjust = 0.001)

# predicted transition matrix
# Expected
G_stats <- function(observ, expect) {
  g <- 0
  for (i in 1:nrow(observ)) {
    for (j in 1:ncol(observ)) {
      if (observ[i, j] != 0 & expect[i, j] != 0) 
        g <- g + observ[i, j] * log(observ[i, j]/expect[i, j])
    }
  }
  g <- 2*g
  g
}

Chi_stats <- function(observ, expect) {
  g <- 0
  for (i in 1:nrow(observ)) {
    for (j in 1:ncol(observ)) {
      if (expect[i, j] != 0) 
        g <- g +((observ[i, j] - expect[i, j])^2)/expect[i, j]
    } 
  }
  g
}
perm <- 1000
p_hard <- list()
pred_trmat <- list()
chi_perm <- list()
g_perm <- list()
for(i in names_mod) {
  p_tmp <- apply(p_df[[i]], 1, function(x) sample(states, size = perm, replace = T, prob = x))
  p_hard[[i]] <- as.data.frame(t(p_tmp))
  
  pred_trmat[[i]] <- apply(p_hard[[i]], 2, function(x) table(from, x))

  chi_perm[[i]]  <- apply(pred_trmat[[i]], 2, function(x) Chi_stats(obs_trmat, matrix(x,4)))

  g_perm[[i]]  <- apply(pred_trmat[[i]], 2, function(x) G_stats(obs_trmat, matrix(x,4)))

}

chi_perm2 = do.call(cbind, g_perm)
quartz()
boxplot(abs(chi_perm2))

chi_perm3 = reshape2::melt(chi_perm2)
chi_perm3=chi_perm3[which(chi_perm3$value<1000),]
quartz()
beeswarm::beeswarm(chi_perm3$value~chi_perm3$Var2, spacing = .5, cex = .5, col = alpha("black",0.1), pch = 19)



