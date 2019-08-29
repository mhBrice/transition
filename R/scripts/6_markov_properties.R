### 


library(msm)
library(graphicsutils)
library(scales)

source('R/functions/markov.R')
source('R/functions/plot_msm.R')

### DATA ####

source('R/scripts/prep_trans.R')

# Load msm results

load("res/msm_all75_drainph.rda")

msm_glb <- msm_all75[["msm_glb"]]

### PREDICTED TRANSITION PROBABILITY MATRIX AND JUMP MATRIX ####

envmean <- aggregate(states_ba[,c("sTP", "CMI", "DRAIN", "PH_HUMUS")], 
                     by = list(states_ba$ecoreg3), mean)
mixed_mean <- envmean[3,-1]

# Probability matrix

p0 <- pmatrix.msm(msm_glb, t=10, covariates = as.list(mixed_mean), ci = "normal")

pN1 <- pmatrix.msm(msm_glb, t=10, covariates = as.list(c(mixed_mean, natural1=1)), ci = "normal")

pN2 <- pmatrix.msm(msm_glb, t=10, covariates = as.list(c(mixed_mean, natural2=1)), ci = "normal")

pL1 <- pmatrix.msm(msm_glb, t=10, covariates = as.list(c(mixed_mean, logging1=1)), ci = "normal")

pL2 <- pmatrix.msm(msm_glb, t=10, covariates = as.list(c(mixed_mean, logging2=1)), ci = "normal")

# Jump matrix 

jump0 <- pnext.msm(msm_glb, covariates = as.list(mixed_mean))

jumpN1 <- pnext.msm(msm_glb, covariates = as.list(c(mixed_mean, natural1=1)))

jumpN2 <- pnext.msm(msm_glb, covariates = as.list(c(mixed_mean, natural2=1)))

jumpL1 <- pnext.msm(msm_glb, covariates = as.list(c(mixed_mean, logging1=1)))

jumpL2 <- pnext.msm(msm_glb, covariates = as.list(c(mixed_mean, logging2=1)))

p_list <- list(p0, pN1, pN2, pL1, pL2)
j_list <- list(jump0, jumpN1, jumpN2, jumpL1, jumpL2)

### COMPUTE MARKOV PROPERTIES ####

markov_res <- list()
for(i in 1:length(p_list)) {
  p <- p_list[[i]]
  j <- j_list[[i]]
  succ <- succession(pmat = p, jump = j, ci = T)
  
  markov_res[[i]] <- succ 
}


### Plot transition matrix #####

dist_title <- c("Minor", "Moderate natural", "Major natural", 
                "Moderate logging", "Major logging")

mat <- matrix(c(0,1,1,0,
                2,2,4,4,
                3,3,5,5),3,byrow = T)

pdf("res/fig4_pmatrix.pdf", width = 4.2, height = 6.2)
#quartz(width = 4.2, height = 6.2)
layout(mat)
par(mar=c(.5,2,3.5,1))
for(i in 1:length(p_list)) {
  tr <- p_list[[i]]
  tr <- tr[["estimates"]]
  #dimnames(tr) <- list(states,states)
  if(i == 1) labels = TRUE else labels = FALSE
  plot_trans(pmat=tr, states_lab = c("B", "M", "P", "T"), 
             labels = labels, main = dist_title[i])
}
dev.off()

### Steady states ####

pdf("res/figx_steady_bp.pdf",
    width = 5, height = 4)
#quartz(width = 5, height = 4)
par(mar=c(3,3.3,1,.1))
tr <- do.call(rbind, lapply(markov_res, function(x) x$steady))
bp <- barplot(t(tr), ylim=c(0,1), beside = T, 
              las = 1, cex.axis = 0.8,
              legend.text = states, col = alpha(st_col,.8), border = "white",
              args.legend=list(x="topleft", bty = "n", cex = .8, border="white"))


# Axis labels
mtext("State relative frequency at equilibrium", 2, line = 2.2, cex = .9)

mtext("Minor", 1, at = mean(bp[,1]), line = .3, cex = .8)

lines(x=range(bp[,2:3]), y = c(-.1,-.1), lwd = 2, xpd = NA)
mtext(c("Moderate", "Major"), 1, at = apply(bp[,2:3], 2, mean), line = .3, cex = .8)
mtext("Natural disturbances", 1, at = mean(bp[,2:3]), line = 1.8, font = 2, cex = .8)

lines(x=range(bp[,4:5]), y = c(-.1,-.1), lwd = 2, xpd = NA)
mtext(c("Moderate", "Major"), 1, at = apply(bp[,4:5], 2, mean), line = .3, cex = .8)
mtext("Harvesting", 1, at = mean(bp[,4:5]), line = 1.8, font = 2, cex = .8)

dev.off()

### Steady states over time ####

init_prop <- states_ba %>% group_by(plot_id) %>% slice(1) 

init_prop <- table(init_prop$states)

init_prop <- init_prop/sum(init_prop)


dist_list <- list(as.list(as.list(mixed_mean)),
                  as.list(c(mixed_mean, natural1=1)),
                  as.list(c(mixed_mean, natural2=1)),
                  as.list(c(mixed_mean, logging1=1)),
                  as.list(c(mixed_mean, logging2=1)))

mat <- matrix(c(1,2,3,0,4,5),2,byrow = T)

time <- 1:500
ci <- "none"

pdf("res/res_state90/figx_steady.pdf",
    width = 9, height = 5.5)
#quartz(width = 9, height = 5.5)
layout(mat)
par(mar=c(3.5,3.5,3.5,1))

for(d in 1:5) {
  
  plot0(xlim=range(time), ylim=c(0,1))
  axis(1)
  axis(2, las = 0)
  mtext("Time (years)", 1, line = 2.2, cex = .9)
  mtext("State proportion", 2, line = 2.2, cex = .9)
  mtext(dist_title[d], 3, line = 2, font = 2)
  
  dd <- dist_list[[d]]
  
  # Probability matrix
  p_list <- list()
  for(t in time) {
    p <- pmatrix.msm(msm_glb, t=t, covariates = dd, ci = ci)
    p <- p[1:4,,]
    p <- provideDimnames(p, base = list(states,states))
    p_list[[t]] <- p
  }

    # half life
  half <- markov_res[[d]]$damping
  mtext(paste("Half-life:", round(half,2), "years"), 3, line = .5, cex = .8)
  
  # Proportion
  if(ci == "none") {
    prop <- lapply(p_list, FUN = function(x) init_prop %*% x)
    prop_estimate <- do.call(rbind, prop)
  } else {
    prop <- lapply(p_list, FUN = function(x) apply(x, 3, function(x2) init_prop %*% x2))
    prop <- lapply(prop, function(x) provideDimnames(x, base=list(states)))
    prop_estimate <- do.call(rbind, lapply(prop, function(x) x[,"estimate"]))
    prop_lower <- do.call(rbind, lapply(prop, function(x) x[,"lower"]))
    prop_upper <- do.call(rbind, lapply(prop, function(x) x[,"upper"]))
  }

  
  for(st in 1:4) {
    lines(time, prop_estimate[,st], col= alpha(st_col[st],.7), pch = 19, lwd = 2)
    # confidence interval
    if(ci != "none") {
      polygon(c(time, rev(time)), 
              c(prop_lower[,st], rev(prop_upper[,st])), 
              col = alpha(st_col[st],.3), border = NA) 
    }

  }
  
}
dev.off()

### TOTAL LENGTH OF STAY ####
totlos.msm(msm_glb, start=init_prop, fromt=0, tot=100, covariates = covar_log[[2]])

### Expected first passage time ####

efpt.msm(msm_glb, tostate = 4)

### Estimated ratio of transition intensities ####
# M-T is xx times as likely as T-M under the no disturbance scenario
qratio.msm(msm_glb, ind1 = c(2,4), ind2 = c(4,2), covariates = dist_list[[4]])
qratio.msm(msm_glb, ind1 = c(1,2), ind2 = c(2,1), covariates = dist_list[[1]])

# quartz(width = 11, height = 3.5)
# par(mfrow=c(1,3), mar=c(3,5.5,5,1), oma = c(0,1,0,0))
# for(i in 1:3) {
#   tr = rbind(markov_res[[i]]$steady, 
#              rowSums(markov_res[[i]]$freq_mat)/sum(markov_res[[i]]$freq_mat))
#   barplot(tr, ylab = "Steady vs Observed", beside = T, las = 1, ylim = c(0,.5),
#           legend.text = c("Steady", "Observed"), args.legend=list(x="topleft", bty = "n"))
#   mtext(paste(dist_title[i], "disturbance"), 3, line = 3)
# }


### Persistence -probability to stay- and recurrence time ####

mat <- matrix(1:12, 4, byrow = T)
mat <- cbind(rep(13,4), mat[,1], rep(14,4), mat[,2], rep(15,4), mat[,3])
mat <- rbind(mat, c(0,16,0,17,0,18))
  
ps <- lapply(markov_res, function(x) x$p_persist)
ps_estimate <- do.call(rbind, lapply(ps, function(x) x[,"estimate"]))
ps_lower <- do.call(rbind, lapply(ps, function(x) x[,"lower"]))
ps_upper <- do.call(rbind, lapply(ps, function(x) x[,"upper"]))

recur <- do.call(rbind, lapply(markov_res, function(x) x$recurrence_time))

H_st <- lapply(markov_res, function(x) x$entropy_st)
H_estimate <- do.call(rbind, lapply(H_st, function(x) x[,"estimate"]))
# H_lower <- do.call(rbind, lapply(H_st, function(x) x[,"lower"]))
# H_upper <- do.call(rbind, lapply(H_st, function(x) x[,"upper"])) # weird

H_rel <- do.call(rbind, lapply(markov_res, function(x) x$entropy_rel))

bp <- barplot(ps_estimate[,1], beside = T, plot = F)

pdf("res/res_state85/figx_persist-recur.pdf", width = 7, height = 6)
#quartz(width = 7, height = 6)
layout(mat, heights = c(1,1,1,1,.3), widths = c(.08,1,.08,1, .08, 1))

par(mar=c(.5, 2.1, .5, .6))
for(st in 1:4) {
  # Persistence plot
  plot0(bp, ylim = c(0,1.01), 
        bty = "l", frame.plot = T, yaxs = "i")
  grid(nx = NA, ny = NULL, col = alpha("lightgray", .5), lty = 1, lwd = .8)
  mtext(states[st], 1, adj = 0.05, line = -1.1, cex = .7)
  axis(2, cex.axis = .9, las = 1)
  
  arrows(x0 = bp, 
         y0 = ps_lower[,st],
         y1 = ps_upper[,st], 
         angle = 90, code = 1, 
         length = 0, col = st_col[st], lwd = 1.2, xpd = NA)

  points(bp, ps_estimate[,st], pch = 19, col = st_col[st], cex = 1.2)
  points(bp, ps_estimate[,st], pch = 19, col = alpha(st_col[st], .8), type = "l")
  
  # Recurrence plot
  plot0(xlim = range(bp), ylim = c(0,51), 
        bty = "l", frame.plot = T, yaxs = "i")
  
  grid(nx = NA, ny = NULL, col = alpha("lightgray", .5), lty = 1, lwd = .8)
  mtext(states[st], 1, adj = 0.05, line = -1.1, cex = .7)
  axis(2, cex.axis = .9, las = 1)
  
  points(bp, recur[,st], pch = 19, col = st_col[st], cex = 1.2, xpd = NA)
  points(bp, recur[,st], pch = 19, col = alpha(st_col[st], .8), type = "l")
  
  # Entropy plot
  plot0(xlim = range(bp), ylim = c(0,1), 
        bty = "l", frame.plot = T, yaxs = "i")
  
  grid(nx = NA, ny = NULL, col = alpha("lightgray", .5), lty = 1, lwd = .8)
  mtext(states[st], 1, adj = 0.05, line = -1.1, cex = .7)
  axis(2, cex.axis = .9, las = 1)
  
  points(bp, H_estimate[,st], pch = 19, col = st_col[st], cex = 1.2)
  points(bp, H_estimate[,st], pch = 19, col = alpha(st_col[st], .8), type = "l")
  
}

# Y axis
par(mar=c(0,0,0,0))
plot0()
text(0,0, expression(bold(paste("Probability of persistence ", (p[rr])))), cex = 1.1, srt = 90)

plot0()
text(0,0, "Recurrence time (decades)", cex = 1.1, srt = 90, font = 2)

plot0()
text(0,0, "Entropy", cex = 1.1, srt = 90, font = 2)



# X Axis
par(mar=c(0.5,2.1,0.1,.6))
for(i in 1:3){
  plot0(bp, ylim = c(-1,1), xaxs = "r")
  text(bp, .8, c("0", "1", "2", "1", "2"), cex = 1.1, xpd = NA, adj = .5)
  
  lines(x=range(bp[2:3])+c(-.2,.2), y = c(0,0), lwd = 2, xpd = NA)
  lines(x=range(bp[4:5])+c(-.2,.2), y = c(0,0), lwd = 2, xpd = NA)
  
  text(c(mean(bp[2:3]), mean(bp[4:5])), -.8, c("Natural", "Harvesting"), 
       font = 2, cex = 1.1, xpd = NA, adj = .5)
}

dev.off()

### ENTROPY OF CONTINUOUS-TIME MARKOV CHAINS ####



# without ci
# tr <- do.call(rbind, lapply(markov_res, function(x) x$recurrence_time))
# tr <- do.call(bind, lapply(tr, function(x) x[,"estimate"]))
# quartz(width = 6, height = 4)
# par(mar=c(3,4.5,1,7))
# barplot(t(tr), ylab = expression(paste("Probability of persistence ", (p[ii]))), 
#         beside = T, las = 1, ylim = c(0,1),
#         legend.text = states, col = alpha(st_col,.8), border = "white", 
#         args.legend=list(x="right", bty = "n", inset=c(-.35,0), title = "Disturbance", cex = .8))

# Recurrence time - measures resilience time to go back to a state after a disturbance
# quartz(width = 6, height = 4)
# par(mar=c(3,4.5,1,7))
# tr <- do.call(rbind, lapply(markov_res, function(x) x$recurrence_time))
# 
# barplot(tr, ylab = "Recurrence time (decades)", beside = T, las = 1,
#         legend.text = dist_title, 
#         args.legend=list(x="right", bty = "n", inset=c(-.35,0), title = "Disturbance", cex = .8))

# quartz(width = 6, height = 4)
# par(mar=c(3,4.5,1,7))
# tr <- do.call(rbind, lapply(markov_res, function(x) x$p_replaceof))
# 
# barplot(tr, ylab = "Replacement of", beside = T, las = 1,
#         legend.text = dist_title, 
#         args.legend=list(x="right", bty = "n", inset=c(-.35,0), title = "Disturbance", cex = .8))



################################
### RAW STATE TRANSITION ####
################################

states_d0 <- states_trans %>%
  mutate(disturb=as.numeric(as.character(disturb))) %>%
  group_by(plot_id) %>%
  filter(max(disturb)==0) 

states_d1 <- states_trans %>%
  mutate(disturb=as.numeric(as.character(disturb))) %>%
  group_by(plot_id) %>%
  filter(max(disturb)==1)

states_d2 <- states_trans %>%
  mutate(disturb=as.numeric(as.character(disturb))) %>%
  group_by(plot_id) %>%
  filter(max(disturb)==2)







states_list <- list(states_d0, states_d1, states_d2)
markov_res <- list()
for(i in 1:3) {
  st <- states_list[[i]]
  tr_ls <- trans_mat(st$From, st$To)
  
  tr <- tr_ls$trans
  freq_mat <- tr_ls$raw
  
  steady <- steady_state(tr, plot = F, freq = freq_mat)
  
  succ <- succession(tr, disturb = "Pioneer")
  
  markov_res[[i]] <- list(trans_mat = tr, freq_mat = freq_mat,
                          steady = steady, succ = succ)
}


dist_title <- c("Minor", "Moderate", "Major")

quartz(width = 11, height = 3.8)
par(mfrow=c(1,3), mar=c(1,5.5,5,1), oma = c(0,1,0,0))
for(i in 1:3) {
  tr = markov_res[[i]]$trans_mat
  plot_trans(tr)
  mtext(paste(dist_title[i], "disturbance"), 3, line = 3)
}

quartz(width = 6, height = 4)
par(mar=c(3,5,1,5))
tr <- rbind(markov_res[[1]]$steady, 
            markov_res[[2]]$steady,
            markov_res[[3]]$steady)
barplot(tr, ylim=c(0,.6), ylab = "Steady state", beside = T, las = 1,
        legend.text = dist_title, 
        args.legend=list(x="right", bty = "n", inset=c(-.25,0), title = "Disturbance"))


quartz(width = 11, height = 3.5)
par(mfrow=c(1,3), mar=c(3,5.5,5,1), oma = c(0,1,0,0))
for(i in 1:3) {
  tr = rbind(markov_res[[i]]$steady, 
             rowSums(markov_res[[i]]$freq_mat)/sum(markov_res[[i]]$freq_mat))
  barplot(tr, ylab = "Steady vs Observed", beside = T, las = 1, ylim = c(0,.5),
          legend.text = c("Steady", "Observed"), args.legend=list(x="topleft", bty = "n"))
  mtext(paste(dist_title[i], "disturbance"), 3, line = 3)
}



quartz(width = 6, height = 4)
par(mar=c(3,5,1,5))
tr <- rbind(markov_res[[1]]$succ$p_replaceby, 
            markov_res[[2]]$succ$p_replaceby,
            markov_res[[3]]$succ$p_replaceby)
barplot(tr, ylab = "Replacement by", beside = T, las = 1,
        legend.text = dist_title, 
        args.legend=list(x="right", bty = "n", inset=c(-.25,0), title = "Disturbance"))


quartz(width = 6, height = 4)
par(mar=c(3,5,1,5))
tr <- rbind(markov_res[[1]]$succ$p_replaceof, 
            markov_res[[2]]$succ$p_replaceof,
            markov_res[[3]]$succ$p_replaceof)
barplot(tr, ylab = "Replacement of", beside = T, las = 1,
        legend.text = dist_title, 
        args.legend=list(x="right", bty = "n", inset=c(-.25,0), title = "Disturbance"))


quartz(width = 6, height = 4)
par(mar=c(3,5,1,5))
tr <- rbind(markov_res[[1]]$succ$p_persist, 
            markov_res[[2]]$succ$p_persist,
            markov_res[[3]]$succ$p_persist)
barplot(tr, ylab = "Persistence", beside = T, las = 1,
        legend.text = dist_title, 
        args.legend=list(x="right", bty = "n", inset=c(-.25,0), title = "Disturbance"))

# turnover rate = replacement of
quartz(width = 6, height = 4)
par(mar=c(3,5,1,5))
tr <- rbind(markov_res[[1]]$succ$turn_rate, 
            markov_res[[2]]$succ$turn_rate,
            markov_res[[3]]$succ$turn_rate)
barplot(tr, ylab = "Turnover rate", beside = T, las = 1,
        legend.text = dist_title, 
        args.legend=list(x="right", bty = "n", inset=c(-.25,0), title = "Disturbance"))

quartz(width = 6, height = 4)
par(mar=c(3,5,1,5))
tr <- rbind(markov_res[[1]]$succ$recurrence_time, 
            markov_res[[2]]$succ$recurrence_time,
            markov_res[[3]]$succ$recurrence_time)
barplot(tr, ylab = "Recurrence time", beside = T, las = 1,
        legend.text = dist_title, 
        args.legend=list(x="right", bty = "n", inset=c(-.25,0), title = "Disturbance"))



# visualize how the probabilities change as the number of steps increases

timestep <- 1:10
quartz(width = 11, height = 4)
par(mfrow=c(1,3), mar=c(4,4,3,1))
for(i in 1:3) {
  initState <- rowSums(markov_res[[i]]$freq_mat)/sum(markov_res[[i]]$freq_mat)
  pij <- c()
  for(k in timestep) {
    pij <- rbind.data.frame(pij, t(initState) %*% (markov_res[[i]]$trans_mat %^% k))
  }
  colnames(pij) <- st_name
  
  plot(pij$Boreal ~ timestep, type = "l", col = st_color[1], lwd = 2, ylim = c(0,1),
       ylab = "Proportion")
  lines(pij$Mixed ~ timestep, col = st_color[2], lwd = 2)
  lines(pij$Pioneer ~ timestep, col = st_color[3], lwd = 2)
  lines(pij$Temperate ~ timestep, col = st_color[4], lwd = 2)
  mtext(paste("Disturbance", dist_title[i]), 3, line = 1)
}
legend("topright", legend = st_name, col = st_color, lwd = 2, bty = "n")
