#### MODEL VALIDATION #####

### PACKAGES ####

library(dplyr)
library(msm)

library(caret)
library(DescTools)
library(sf)

source('R/functions/multistate_performance.R')

### DATA ####

source('R/scripts/3_prep_trans.R')

# Load msm results

load("res/msm_all.rda")

msm_all[["msm_glb2"]]=msm_glb2

names_mod <- names(msm_all)
n_mod <- length(msm_all)


# Test data
covar_names <- attr(terms(msm_all$msm_glb2$covariates), "term.labels")

covar_names <- c(covar_names, "TP_xmin")
st_from <- states_ba %>% group_by(plot_id) %>% 
  mutate(t = lead(year_measured) - year_measured) %>% 
  slice(n()-1) %>% select(plot_id, from = states_ba, t, covar_names)

st_to <- states_ba %>% group_by(plot_id) %>% 
  slice(n()) %>% select(plot_id, to = states_ba)

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
accu <- lapply(p_hard, calculate.accuracy(predictions = x, ref.labels = obs_rep$to))


accu <- list()
for(i in 1:n_mod) {
  accu[[names_mod[i]]] <- calculate.accuracy(predictions = p_hard[[i]], ref.labels = obs_rep$to)
}
accu

### 2. Micro and macro averages of the F1-score from confusion matrix ####

### One-vs-all confusion matrices
cm_all <- list()
for(i in 1:n_mod) {
  p_tmp <- p_hard[[i]]
  cm <- vector("list", length(states))
  
  for (s in seq_along(cm)) {
    positive.class <- states[s]
    # in the i-th iteration, use the i-th class as the positive class
    cm[[s]] <- confusionMatrix(as.factor(p_tmp), obs_rep$to, 
                               positive = positive.class)
  }
  cm_all[[names_mod[i]]] <- cm
}

metrics <- c("Precision", "Recall")
lapply(cm_all, function(x) x[[1]]$byClass[, metrics])

# Micro F1
(micro.f1 <- lapply(cm_all, get.micro.f1))


# Macro F1
(macro.f1 <- lapply(cm_all, get.macro.f1))


### Soft classifier ####

### 3. AUC generalization from Hand and Till ####

(multi_auc <- lapply(p_df, function(x) multiclass.auc(x, to)))

multi_auc2 <- list()
for(s in states) {
  l_from <- which(from == s)
  (M <- lapply(p_df, function(x) multiclass.auc(x[l_from,], to[l_from])))
  multi_auc2[[s]] <- M
}

### Plot of pairs AUC

gener_auc <- unlist(lapply(multi_auc, function(x) x[1]))
pair_auc <- lapply(multi_auc, function(x) attr(x, "pair_AUCs"))

# Radar chart
dat = do.call(rbind, pair_auc)
dat=rbind(rep(1,6) , rep(0.8,6) , dat)
dat = as.data.frame(dat)


colors_border <- brewer.pal(8,"Spectral")

radarchart(dat, axistype = 1, 
           #custom polygon
           pcol = colors_border, plwd = 1.5, plty = 2,
           #custom the grid
           cglcol = "grey", cglty = 1, axislabcol = "grey", 
           caxislabels = seq(0.8, 1, length.out = 5), cglwd = 0.8,
           #custom labels
           vlcex = 0.8)
legend(x=0.7, y = 1, legend = names(pair_auc), bty = "n", pch = 20 , col = colors_border, 
      cex = 1, pt.cex=1, title = "Models")

# Parallel plot
ord <- order(pair_auc$msm_glb2, decreasing = T)

quartz(width = 7, height = 5)
par(mar = c(7, 3.2, 1, 8))
plot0(xlim = c(1, 6), ylim = c(.8, 1), grid.col = "grey", yaxs = "i")
axis(2, tick = F, line = -1, las = 1, cex.axis = .8)
mtext("Pairwise AUCs", 2, line = 1.7)

for(i in 1:n_mod) {
  x <- pair_auc[[i]][ord]
  points(x, type = "b", col = colors_border[i], pch = 19, cex = .8)
}
text(1:6, 0.795, names(pair_auc$msm0)[ord], cex = 0.9,
     srt = 90, xpd = NA, adj = 1)
lgd <- legend(x = 6.5, y = 0.9, legend = rev(names(pair_auc)), 
              col = rev(colors_border), lwd =1, cex = .8, bty = "n", xpd = NA)


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

o_perm <- vegan::permatfull(o_df, fixedmar = "rows", shuffle = "samp", times = nperm)

chi_ind_perm <- lapply(o_perm$perm, function(x) chi_contrib(observ = x, expect = p_df))
chi_sum_perm <- do.call(cbind, lapply(chi_ind_perm, sum))

hist(chi_sum_perm, breaks = 30, xlim = range(chi_sum_perm, chi_sum), main = "True and permutation Chi2")
abline(v = chi_sum, col = "red3", lwd = 1.2)

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
obs_trmat <- table(x$from, x$to)

# predicted transition matrix

pred_trmat <- table(x$from, x$pred)

chisq.test(pred_trmat, p=obs_trmat/rowSums(obs_trmat))




