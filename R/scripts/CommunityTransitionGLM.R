#### Model of community transition ####

setwd("~/Documents/GitHub/Doctorat/community_transition")

### PACKAGES ####
require(dplyr)
require(reshape2)

#require(labdsv)
#require(factoextra)
require(vegan)
require(nnet)

require(sf)
#require(nngeo)
require(effects)

require(ggplot2)
require(RColorBrewer)
require(scales)
#require(colorRamps)
require(graphicsutils)

source('R/functions/fig_prob_ima.R')


## DATA ####

# Species
states_ba <- readRDS("data/states_envba.RDS")
states_trans <- readRDS("data/states_ba_trans.RDS")

ecoregion <- st_read("../data/map/ecoregion_simple.gpkg", quiet = T)
ecoregion <- st_transform(ecoregion, 32198)

ecoregion$SOUS_DOM6 <- factor(ecoregion$SOUS_DOM6, rev(c("Spruce-moss",
                                                         "Balsam fir-white birch",
                                                         "Balsam fir-yellow birch",
                                                         "Sugar maple-yellow birch",
                                                         "Sugar maple-basswood",
                                                         "Sugar maple-bitternut hickory")))

# Plotting parameters
st_name <- c("Boreal", "Mixed", "Pioneer", "Temperate")
st_color <- c("#158282", "#A1BD93","#FEAC19", "#D43650")

## TRANSITION TABLE

# All transitions
states_trans2 <- states_trans %>%
  #filter(disturb==0) %>%
  group_by(plot_id) %>%
  slice(n())

tr <- table(states_trans2$From, states_trans2$To)

x=decostand(tr, method = "total", 1)

eigen(x)

# The closer the second eigenvalue is in magnitude to the first, the slower the rate of convergence

markov2 <- new('markovchain',
               transitionMatrix = matrix(x,4), 
               states = st_name)

It <- diag(4)

N <- solve(It-matrix(x,4))

steadyStates(markov2)

### MULTINOMIAL TRANSITION ####

# Remove imposible transitions

Trans_df <- states_trans %>%
  filter(transition!="Boreal->Temperate" & transition!="Temperate->Boreal")

# remove forest where TYPEHUMUS not evaluated 
Trans_df <-  Trans_df[-which(is.na(Trans_df$DRAIN) | is.na(Trans_df$TYPEHUMUS)),]


#1. Scale variable ####
scale_tp <- scale(Trans_df$aTP)
scale_pp <- scale(Trans_df$aPP)
scale_slope_tp <- scale(Trans_df$slope_aTP)
scale_slope_pp <- scale(Trans_df$slope_aPP)
scale_harv <- scale(Trans_df$harvest_100)
scale_time <- scale(Trans_df$time_interv) 
scale_drain <- scale(Trans_df$DRAIN)
scale_age <- scale(Trans_df$age_mean)

var_to_scale <- c("time_interv",
                  "aTP", "aPP", "slope_aTP", "slope_aPP", 
                  "sTP", "sPP", "slope_sTP", "slope_sPP", 
                  "CMI", "slope_CMI",
                  "GS1", "GS2", "GSlength", "slope_GS1", "slope_GS2", "slope_GSlength",
                  "age_mean", "harvest_100", "nat_mort_100", 
                  "DRAIN", "EPMATORG", "PH_HUMUS")

Trans_df[, var_to_scale] <- scale(Trans_df[, var_to_scale])

Trans_df$nat_disturb <- as.factor(Trans_df$nat_disturb)

#3. Evaluate contribution of the coef: build all formula

coefs <-c("aTP", "aPP", "I(aTP^2)","I(aPP^2)", 
          "slope_aTP",
          "slope_aPP", 
          "CMI", "slope_CMI",
          "time_interv", "age_mean", 
          "TYPEHUMUS", "DRAIN", "ecoreg10",
          "harvest_100", "harvest_100:time_interv",
          "nat_mort_100", "nat_mort_100:time_interv") 


contrib_form <- as.character()
tested_coefs <- as.character()

for(i in 1:length(coefs)){
  temp_coefs <- coefs[-i]
  tested_coefs  <- append(tested_coefs, coefs[i])
  form <- paste0("To ~ ", paste(temp_coefs, collapse = "+"))
  contrib_form <- append(contrib_form, form)
}

# Add full model
full_formula <- paste0("To ~ ", paste(coefs, collapse = "+"))
null_formula <- "To ~ 1"

# Create DF of formula
formulas <- data.frame(tested = tested_coefs, formula = contrib_form, stringsAsFactors = FALSE)

## ---- multinom_trans1 ----

#5.1 Full Models and pseudo-R2
Trans_df$To <- relevel(Trans_df$To, ref = "Boreal")
B_multi_full <- multinom(full_formula, data = Trans_df, subset = From=="Boreal", trace = F)

Trans_df$To <- relevel(Trans_df$To, ref = "Mixed")
M_multi_full <- multinom(full_formula, data = Trans_df, subset = From=="Mixed", trace = F)

Trans_df$To <- relevel(Trans_df$To, ref = "Pioneer")
P_multi_full <- multinom(full_formula, data = Trans_df, subset = From=="Pioneer", trace = F)

Trans_df$To <- relevel(Trans_df$To, ref = "Temperate")
T_multi_full <- multinom(full_formula, data = Trans_df, subset = From=="Temperate", trace = F)

#5.2 Null Models and pseudo-R2
for(state in st_name) {
  st <- substring(state, 1, 1)
  # Null Models
  tmp.null <- nnet::multinom(null_formula, data=Trans_df, subset= From==state, trace = F)
  assign(paste0(st, '_multi_null'), tmp.null)
  # Compute McFadden pseudo-R2
  tmp.r2 <- 1 - (get(paste0(st, '_multi_full'))$deviance / tmp.null$deviance)
  print(assign(paste0(st, '_Rs'), tmp.r2))
}


#6. Compute contributions of each variable to the full model
ls_contrib <-list()

for(state in st_name){
  contrib <- as.numeric()
  st <- substring(state, 1, 1)
  for(i in 1:nrow(formulas)){
    multi <- nnet::multinom(formulas$formula[i], data = Trans_df, subset = From==state, trace = F)
    contrib <- append(contrib, (multi$AIC - get(paste0(st, '_multi_full'))$AIC))
  }
  ls_contrib[[state]] <- c(contrib, get(paste0(st, '_multi_full'))$AIC, get(paste0(st, '_Rs')))
}

# Summarize information in table
contrib.full <- do.call("rbind", ls_contrib)
colnames(contrib.full) <- c(gsub('[I(^)]', '', formulas$tested), "AIC", "R2")

contrib.full

#7. stepAIC
stepmod <- list()
for(state in st_name) {
  st <- substring(state, 1, 1)
  Trans_df$To <- relevel(Trans_df$To, ref = state)
  tmp.step <- MASS::stepAIC(get(paste0(st, '_multi_full')), trace = FALSE)
  stepmod[[state]] <- tmp.step
  # Pseudo R2
  tmp.r2 <- 1 - (tmp.step$deviance / get(paste0(st, '_multi_null'))$deviance)
  assign(paste0(st, '_step_Rs'), tmp.r2)
}

#8 variable contribution for reduced models

# ls_contrib <-list()
# 
# for(state in st_name){
#   contrib <- as.numeric()
#   st <- substring(state, 1, 1)
#   
#   contrib_form <- as.character()
#   tested_coefs <- as.character()
#   coefs_sel <- attr(terms(stepmod[[state]]), "term.labels")
#   for(i in 1:length(coefs_sel)){
#     temp_coefs <- coefs_sel[-i]
#     tested_coefs  <- append(tested_coefs, coefs_sel[i])
#     form <- paste0("To ~ ", paste(temp_coefs, collapse = "+"))
#     contrib_form <- append(contrib_form, form)
#   }
#   form.sel <- data.frame(tested = tested_coefs, formula = contrib_form, stringsAsFactors = FALSE)
#   
#   contrib.sel <- t(data.frame(row.names = coefs, deltaAIC = rep(NA, length(coefs))))
#   row.names(contrib.sel) <- state
#   for(i in 1:nrow(form.sel)){
#       l <- which(coefs %in% form.sel[i,1])
#       multi <- nnet::multinom(form.sel$formula[i], data = Trans_df, subset = From==state, trace = F)
#       contrib.sel[,l] <- multi$AIC - stepmod[[state]]$AIC
#   }
# 
#   ls_contrib[[state]] <- cbind(contrib.sel, 
#                                AIC = stepmod[[state]]$AIC, 
#                                R2 = get(paste0(st, '_step_Rs')))
# }
# 
# # Summarize information in table
# contrib.sel <- do.call("rbind", ls_contrib)
# colnames(contrib.sel) <- c(gsub('[I(^)]', '', formulas$tested), "AIC", "R2")
# contrib.sel

#### COEFPLOT ####

# Get coefficient and confidence interval
coef_stepmod <- lapply(stepmod, FUN = function(x) coefficients(x)[,-1])
CI_stepmod <- lapply(stepmod, FUN = function(x) confint(x)[-1,,])

# New coefficient names
newnames <- c("TP", "TP2", "∆TP", 
              "PP", "PP2", "∆PP", 
              "age", "harvest", "fire", "insect", 
              "time", "time2",
              "TP:PP", "TP:∆TP",
              "harvest:time", "TP:harvest")
newnames <- cbind.data.frame("coefs" = colnames(coefficients(M_multi_full))[-1], newnames)

# Plot of coefficient
quartz()
par(mfrow = c(2,2), mar = c(3,5,2,1))
for(st_from in st_name) {
  # Coefficient
  coef_st <- coef_stepmod[[st_from]]
  # Confidence interval
  CI_st <- CI_stepmod[[st_from]]
  rangeCI <- range(CI_st)
  # Number of variables
  nvar <- ncol(coef_st)
  # State_to
  nst <- row.names(coef_st)
  # Empty plot
  plot(1, type="n", xlab="", ylab="", 
       xlim = c(-3,3), #xlim=c(rangeCI[1],rangeCI[2]), 
       ylim=c(1, nvar), yaxt = "n", 
       main = paste0("From ", st_from))
  abline(v = 0, lty = 2, col = "grey25")
  abline(h = 1:nvar, lty = 3, col = "grey")
  # Labels
  coef_labs <- colnames(coef_st) #subset(newnames, coefs %in% colnames(coef_st))$newnames
  axis(side = 2, at = 1:nvar, labels = rev(coef_labs), las = 1, cex = 0.9)
  l = -0.12
  for(st_to in nst) {
    points(coef_st[st_to,], nvar:1+l, pch = 19, col = st_color[which(st_to == st_name)])
    for (i in 1:nvar) {
      yl <- nvar +1 - i
      lines(x=c(CI_st[i,1,st_to], CI_st[i,2,st_to]), y=c(yl+l, yl+l), 
            col = st_color[which(st_to == st_name)])
    }
    l = l + 0.12
    
  }
  
}

#### PREDICTION MAP ####

fit <- as.data.frame(fitted(M_multi_full2)) %>%
  mutate(st_pred = apply(., 1, function(x) names(which.max(x))))

xy_fit <- xy %>%
  right_join(cbind.data.frame(plot_id = subset(Trans_df, From=="Mixed")$plot_id, fit), 
             by = "plot_id")

plot(xy_fit["st_pred"], pal = alpha(st_color, .6), pch = 19, cex = .8)

mypal <- colorRampPalette(c("blue3", "grey90", "red3"))


quartz(height = 4, width = 7)
par(mfrow= c(2,2), mar=c(0,1.1,1,.1))
for(st in st_name) {
  levpal <- cut(xy_fit[[st]], 50)
  
  plot(st_geometry(ecoregion), border = "gray60",
       key.pos=NULL, main=NULL)
  plot(st_geometry(xy_fit), pch = 19, cex = .5, 
       col = alpha(mypal(50)[levpal],.5), add=T)
  mtext(paste0("Mixed -> ", st), 3, cex = .8, font = 2)
}
   

fit <- as.data.frame(fitted(T_multi_full)) %>%
  mutate(st_pred = apply(., 1, function(x) names(which.max(x))))

xy_fit <- xy %>%
  right_join(cbind.data.frame(plot_id = subset(Trans_df, From=="Temperate")$plot_id, fit), 
             by = "plot_id")
quartz(height = 4, width = 7)
par(mfrow= c(2,2), mar=c(0,1.1,1,.1))
for(st in st_name[-1]) {
  levpal <- cut(xy_fit[[st]], 50)
  
  plot(st_geometry(ecoregion), border = "gray60",
       key.pos=NULL, main=NULL)
  plot(st_geometry(xy_fit), pch = 19, cex = .5, 
       col = alpha(mypal(50)[levpal],.5), add=T)
  mtext(paste0("Temperate -> ", st), 3, cex = .8)
}



### PREDICTION PLOT ####

Trans_sf <- st_as_sf(Trans_df)
Trans_sf$nat_disturb <- as.ordered(Trans_sf$nat_disturb)
Trans_sf$harvest_intensity100 <- as.ordered(Trans_sf$harvest_intensity100)
Trans_sf <- Trans_sf %>% group_by(plot_id) %>%
  mutate_at(vars(TP,PP,slope_TP,slope_PP), mean) %>%
  mutate_at(vars(nat_disturb, harvest_intensity100), max)

Trans_sf <- Trans_sf[!duplicated(Trans_sf$plot_id),]
# Trans_sf <- filter(Trans_sf, From == 'Mixed')

##---------
## CHANGEMENT CLIMATIQUE ####
# 
# qTP <- quantile(Trans_df$slope_TP[which(Trans_df$From=="Mixed")], c(.05,.95))
# qPP <- quantile(Trans_df$slope_PP[which(Trans_df$From=="Mixed")], c(.1,.5,.9))
# 
# 
# eff_M_cc <- Effect(focal.predictors = c("slope_TP", "slope_PP"),
#                    mod = M_multi_full,
#                    xlevels = list(slope_TP = seq(qTP[1], qTP[2], length.out = 100),
#                                   slope_PP = c(qPP[1], qPP[2], qPP[3])))

eff_M_cc <- Effect(focal.predictors = c("slope_TP", "slope_PP"),
                   mod = M_multi_full,
                   xlevels = list(slope_TP = seq(-1.5, 1.1, length.out = 100),
                                  slope_PP = c(-.51, .27, 1.46)))

tp <- (eff_M_cc$model.matrix[,"slope_TP"] * attr(scale_slope_tp,"scaled:scale") + 
  attr(scale_slope_tp,"scaled:center")) *10
(c(-.51, .27, 1.46)* attr(scale_slope_pp,"scaled:scale") + 
  attr(scale_slope_pp,"scaled:center")) *10

effdf_M_cc <- cbind.data.frame(tp, 
                               pp = eff_M_cc$model.matrix[,"slope_PP"],
                               eff_M_cc$prob, 
                               eff_M_cc$lower.prob, 
                               eff_M_cc$upper.prob)

png("images/effect_cc0.png", res=300, 
    width = 7.2, height = 5, units = 'in', bg = "transparent", type='cairo')
# quartz(width = 7.2, height = 5)
layout(mat = matrix(c(1,1,2,1,1,3,4,4,0),3, byrow = T), heights = c(1,1, .4), 
       widths = c(1,1,.9))

par(mar=c(4,4.5,.1,.1), xaxs = "i")

plot.default(effdf_M_cc[,1:2], type = "n", ylim = c(0,1), 
             xaxt="n", yaxt="n", ylab = "", xlab = "")

axis(1, tcl = -0.3, col.ticks = "grey60", cex.axis=1.5)

axis(2, tcl = -0.3, col.ticks  = "grey60", las = 1, cex.axis=1.5)

# loop over "to" transition
for(st_to in st_name) {
  col.st <- alpha(st_color[which(st_name==st_to)], 1)
  
  tmp0 <- effdf_M_cc[which(effdf_M_cc$pp<0),] 
  tmp1 <- effdf_M_cc[which(effdf_M_cc$pp==.27),]
  tmp2 <- effdf_M_cc[which(effdf_M_cc$pp==1.46),]
  
  lines(tmp0[,c('tp', paste0("prob.", st_to))], col = col.st, lwd = 2, lty = 1)
  lines(tmp1[, c('tp', paste0("prob.", st_to))], col = col.st, lwd = 2, lty = 2)
  lines(tmp2[, c('tp', paste0("prob.", st_to))], col = col.st, lwd = 2, lty = 3)
}
legend("topright", legend = st_name[c(2,4,3,1)], bty = "n", col = "transparent", xpd=NA, 
       text.col = st_color[c(2,4,3,1)], text.font = 2, cex = 1.5, y.intersp = 1.4)

mtext("Temperature change (°C/decade)", side = 1, line = 3, cex = 1.1, font=2)
mtext("Probability of transition", side = 2, line = 3,  cex = 1.1, font =2, las = 0)

# map1

par(mar=c(5,0,0,0), xaxs = "r")
plot(st_geometry(ecoregion), border = "gray60",
     key.pos=NULL, main=NULL)
ordhar <- order(Trans_sf$slope_TP)
mypal <- colorRampPalette(c("blue3", "grey90", "red3"))
levpal <- cut(as.numeric(as.character(Trans_sf$slope_TP[ordhar])),50)
plot(st_geometry(Trans_sf)[ordhar], cex = .2, pch = 19, 
     col = mypal(50)[levpal], add = T)

bx <- par("usr")*.9
xs <- seq(bx[1],bx[2], len = 50)
for (i in 1:50) {
  polygon(c(xs[i], xs[i+1], xs[i+1], xs[i]), c(bx[3],bx[3],bx[3]+50000,bx[3]+50000),
          col = mypal(50)[i], border = NA)
}
mtext(side =1, "Temperature change", cex = .8, line = .5)



# map2
par(mar=c(5,0,0,0), xaxs = "r")
plot(st_geometry(ecoregion), border = "gray60",
     key.pos=NULL, main=NULL)
ordhar <- order(Trans_sf$slope_PP)
mypal <- colorRampPalette(c("blue", "grey90", "red3"))
levpal <- cut(as.numeric(as.character(Trans_sf$slope_PP[ordhar])),50)
plot(st_geometry(Trans_sf)[ordhar], cex = .2, pch = 19, 
     col = mypal(50)[levpal], add = T)

bx <- par("usr")*.9
xs <- seq(bx[1],bx[2], len = 50)
for (i in 1:50) {
  polygon(c(xs[i], xs[i+1], xs[i+1], xs[i]), c(bx[3],bx[3],bx[3]+50000,bx[3]+50000),
          col = mypal(50)[i], border = NA)
}
mtext(side =1, "Precipitation change", cex = .8, line = .5)



par(mar=c(0,5,0,.1))
plot0()
legend("center", 
       legend = c("Small decrease in precipitation (-1mm/decade)", 
                  "Moderate increase in precipitation (+11mm/decade)", 
                  "Large increase in precipitation (+30mm/decade)"), 
       bty = "n",  cex = 1.4, lty = c(1,2,3), lwd = 2)

dev.off()



##---------
## AGE ####

eff_M_age <- Effect(focal.predictors = "age_mean",
                mod = M_multi_full,
                xlevels = list(age_mean = seq(-1.4, 2, length.out = 100)))
age <- eff_M_age$model.matrix[,"age_mean"] * attr(scale_age,"scaled:scale") + 
  attr(scale_age,"scaled:center")

effdf_M_age <- cbind.data.frame(age, 
                             eff_M_age$prob, 
                             eff_M_age$lower.prob, 
                             eff_M_age$upper.prob)

png("images/effect_age.png", res=300, 
    width = 7.2, height = 5, units = 'in', bg = "transparent", type='cairo')
# quartz(width = 7.2, height = 5)
layout(mat = matrix(c(1,1,2,1,1,3,4,4,0),3, byrow = T), heights = c(1,1, .4), 
       widths = c(1,1,.9))

par(mar=c(4,4.5,.1,.1), xaxs="i")

plot.default(effdf_M_age[,1:2], type = "n", ylim = c(0,1), 
             xaxt="n", yaxt="n", ylab = "", xlab = "")

axis(1, tcl = -0.3, col.ticks = "grey60", cex.axis=1.5)

axis(2, tcl = -0.3, col.ticks  = "grey60", las = 1, cex.axis=1.5)

# loop over "to" transition
for(st_to in st_name) {
  col.st <- alpha(st_color[which(st_name==st_to)], 1)
  
  lines(effdf_M_age[,c('age', paste0("prob.", st_to))], col = col.st, lwd = 2)
  polygon(c(effdf_M_age[,'age'], 
            rev(effdf_M_age[,'age'])), 
          c(effdf_M_age[,paste0("L.prob.", st_to)], 
            rev(effdf_M_age[,paste0("U.prob.", st_to)])), 
          col = alpha(col.st, 0.1), border = NA)
}
mtext("Stand age", side = 1, line = 3,  cex = 1.1, font =2)
mtext("Probability of transition", side = 2, line = 3,  cex = 1.1, font =2, las = 0)

par(mar=c(0,0,0,0), xaxs="r")
plot0()

# legend 1
par(mar=c(5,0,0,0))
plot0()
legend("left", legend = st_name[c(2,4,3,1)], bty = "n", col = "transparent", xpd=NA, 
       text.col = st_color[c(2,4,3,1)], text.font = 2, cex = 1.5, y.intersp = 1.4)


par(mar=c(0,4,0,.1))
plot0()


dev.off()



##---------
## ECOREG ####

eff_M_reg <- Effect(focal.predictors = "ecoreg10",
                mod = M_multi_full)

region <- c("1", "2E", "2W", "3E", "3W", "4E", "4W", "5E", "5W", "6")
effdf_M_reg <- cbind.data.frame(region, 
                                eff_M_reg$prob, 
                                eff_M_reg$lower.prob, 
                                eff_M_reg$upper.prob)

# no mixed forests in the spruce domain
effdf_M_reg[10,-1] <- NA

reg_legend <- cbind(1:6, levels(Trans_df$ecoreg6))

png(paste0("images/effect_region.png"), res=300, 
    width = 7.2, height = 5, units = 'in', bg = "transparent", type='cairo')

# quartz(width = 7.2, height = 5)
layout(mat = matrix(c(1,1,2,1,1,3,4,4,0),3, byrow = T), heights = c(1,1, .4), 
       widths = c(1,1,.9))

par(mar=c(4,4.5,.1,.1))

plot.default(effdf_M_reg[,1:2], type = "n", ylim = c(0,1), xlim = c(1,10),
             xaxt="n", yaxt="n", ylab = "", xlab = "")

axis(1, at=1:10, labels = effdf_M_reg$region, tcl = -0.3, col.ticks = "grey60", cex.axis=1.5)

axis(2, tcl = -0.3, col.ticks  = "grey60", las = 1, cex.axis=1.5)



# loop over "to" transition
for(st_to in st_name) {
  col.st <- st_color[which(st_name==st_to)]
  
  points(effdf_M_reg[1:10,c(paste0("prob.", st_to))], type = "b",
         col = col.st, lwd = 1.5, pch =19, cex = 1.2)
  arrows(x0 = 1:10, 
         y0 = effdf_M_reg[,paste0("L.prob.", st_to)],
         y1 = effdf_M_reg[,paste0("U.prob.", st_to)], 
         angle = 90, code = 1, 
         length = 0, col = alpha(col.st, 0.2), lwd = 2.5)
  
}
mtext("Bioclimatic domains", side = 1, line = 3,  cex = 1.1, font =2)
mtext("Probability of transition", side = 2, line = 3,  cex = 1.1, font =2, las = 0)

# map
par(mar=c(5,0,0,0))
plot(st_geometry(ecoregion), border = "gray35", 
     col = rev(alpha(brewer.pal(6,"Greys"),.25))[ecoregion$SOUS_DOM6],
     key.pos=NULL, main=NULL)

# legend 1
par(mar=c(5,0,0,0))
plot0()
legend("left", legend = st_name[c(2,4,3,1)], bty = "n", col = "transparent", xpd=NA, 
       text.col = st_color[c(2,4,3,1)], text.font = 2, cex = 1.5, y.intersp = 1.4)

#legend 2
par(mar=c(0,6,0,0))
plot0()

legend("center", legend = reg_legend[,2], bty = "n", xpd=T,
       cex = 1.4, pt.cex =1.4, pch = reg_legend[,1],
       ncol = 2)
dev.off()



##---------
## HUMUS ####


eff_M_hum <- Effect(focal.predictors = "TYPEHUMUS",
                    mod = M_multi_full)

humus <- factor(c("Mull", "Moder", "Mor", "Peat", "Anmoor", "Organic", "No humus"), 
                levels = c("Mull", "Moder", "Mor", "Peat", "Anmoor", "Organic", "No humus"))

effdf_M_hum <- cbind.data.frame(humus, 
                                eff_M_hum$prob, 
                                eff_M_hum$lower.prob, 
                                eff_M_hum$upper.prob)

#remove no humus (only one occurence)
effdf_M_hum <- filter(effdf_M_hum, humus != "No humus")

png(paste0("images/effect_humus.png"), res=300, 
    width = 7.2, height = 5, units = 'in', bg = "transparent", type='cairo')
# quartz(width = 7.2, height = 5)
layout(mat = matrix(c(1,1,2,1,1,3,4,4,0),3, byrow = T), heights = c(1,1, .4), 
       widths = c(1,1,.9))

par(mar=c(4,4.5,.1,.1))

plot.default(effdf_M_hum[,1:2], type = "n", ylim = c(0,1), 
             xaxt="n", yaxt="n", ylab = "", xlab = "")

axis(1, at = 1:6, labels = F, tcl = -0.3, col.ticks = "grey60", cex.axis=1.5)

text(1:6, par("usr")[3] - 0.05, labels = effdf_M_hum$humus, srt = 35, adj=.9, 
     xpd = NA, cex =1.4)

axis(2, tcl = -0.3, col.ticks  = "grey60", las = 1, cex.axis=1.5)

# loop over "to" transition
for(st_to in st_name) {
  col.st <- st_color[which(st_name==st_to)]
  
  points(effdf_M_hum[,c('humus', paste0("prob.", st_to))], type = "b",
         col = col.st, pch =19, cex = 1.2)
  arrows(x0 = 1:10, 
         y0 = effdf_M_hum[,paste0("L.prob.", st_to)],
         y1 = effdf_M_hum[,paste0("U.prob.", st_to)], 
         angle = 90, code = 1, 
         length = 0, col = alpha(col.st, 0.2), lwd = 2.5)
  
}
mtext("Humus types", side = 1, line = 4.2,  cex = 1.1, font =2)
mtext("Probability of transition", side = 2, line = 3,  cex = 1.1, font =2, las = 0)

# map
par(mar=c(5,0,0,0))
plot(st_geometry(ecoregion), border = "gray60",
     key.pos=NULL, main=NULL)
tmp_sf <- Trans_sf[which(Trans_sf$TYPEHUMUS!="NA"),]
plot(st_geometry(tmp_sf), cex = .2, pch = 19, 
     col = brewer.pal(6, "Dark2")[tmp_sf$TYPEHUMUS], add = T)
legend("bottom", legend = effdf_M_hum$humus, col = brewer.pal(6, "Dark2"), pch = 19, 
       bty = "n", ncol = 2, inset = c(0,-.27), xpd = NA, cex = 1.1)


# legend 1
par(mar=c(5,0,0,0))
plot0()
legend("left", legend = st_name[c(2,4,3,1)], bty = "n", col = "transparent", xpd=NA, 
       text.col = st_color[c(2,4,3,1)], text.font = 2, cex = 1.5, y.intersp = 1.4)

#legend 2
par(mar=c(0,6,0,0))
plot0()


dev.off()

##---------
## TIME * NATURAL DISTURBANCES ####
plot(Trans_sf["nat_disturb"], pch = 19, cex = .5)

eff_M_nat <- Effect(focal.predictors = c("time_interv", "nat_disturb"),
                     mod = M_multi_full, 
                     xlevels = list(time_interv = seq(-1.34, 1.78, len = 100)))


nat_disturb1 <- eff_M_nat$model.matrix[,"nat_disturb1"]
nat_disturb2 <- eff_M_nat$model.matrix[,"nat_disturb2"]
time <- eff_M_nat$model.matrix[,"time_interv"] * attr(scale_time,"scaled:scale") +
  attr(scale_time,"scaled:center")

effdf_M_nat <- cbind.data.frame(time,
                                nat_disturb1,
                                nat_disturb2,
                                eff_M_nat$prob, 
                                eff_M_nat$lower.prob,
                                eff_M_nat$upper.prob)



png(paste0("images/effect_nat0.png"), res=300, 
    width = 7.2, height = 5, units = 'in', bg = "transparent", type='cairo')
# quartz(width = 7.2, height = 5)
layout(mat = matrix(c(1,1,2,1,1,3,4,4,0),3, byrow = T), heights = c(1,1, .4), 
       widths = c(1,1,.9))

par(mar=c(4,4.5,.1,.1), xaxs="i")

plot.default(effdf_M_nat[,1:2], type = "n", ylim = c(0,1), 
             xaxt="n", yaxt="n", ylab = "", xlab = "")

axis(1, tcl = -0.3, col.ticks = "grey60", cex.axis=1.5)

axis(2, tcl = -0.3, col.ticks  = "grey60", las = 1, cex.axis=1.5)

# loop over "to" transition
for(st_to in st_name) {
  col.st <- st_color[which(st_name==st_to)]
  tmp0 <- effdf_M_nat[which(effdf_M_nat$nat_disturb1==0 & effdf_M_nat$nat_disturb2==0),]
  tmp1 <- effdf_M_nat[which(effdf_M_nat$nat_disturb1==1),]
  tmp2 <- effdf_M_nat[which(effdf_M_nat$nat_disturb2==1),]
  
  lines(tmp0[,c('time', paste0("prob.", st_to))], col = col.st, lwd = 2, lty = 1)
  # polygon(c(tmp0[,'time'],
  #           rev(tmp0[,'time'])), 
  #         c(tmp0[,paste0("L.prob.", st_to)], 
  #           rev(tmp0[,paste0("U.prob.", st_to)])), 
  #         col = alpha(col.st, 0.08), border = NA)
  
  lines(tmp1[, c('time', paste0("prob.", st_to))], col = col.st, lwd = 2, lty = 2)
  lines(tmp2[, c('time', paste0("prob.", st_to))], col = col.st, lwd = 2, lty = 3)

}  
mtext("Time since first inventory (years)", side = 1, line = 3,  cex = 1.1, font =2)
mtext("Probability of transition", side = 2, line = 3,  cex = 1.1, font =2, las = 0)

# map
par(mar=c(5,0,0,0), xaxs="r")
plot(st_geometry(ecoregion), border = "gray60",
     key.pos=NULL, main=NULL)
ordnat <- order(Trans_sf$nat_disturb)
mypal <- c("#FFEDA0", "#EE870C", "#bc1800")
plot(st_geometry(Trans_sf)[ordnat], cex = .2, pch = 19, 
     col = mypal[Trans_sf$nat_disturb[ordnat]], add = T)
legend("bottom", legend = c("No disturbances", "Moderate disturbances", "Major disturbances"),
       col = mypal, pch = 19, cex = 1.1,
       bty = "n", inset = c(0,-.27), xpd = NA)


# legend 1
par(mar=c(5,0,0,0))
plot0()
legend("left", legend = st_name[c(2,4,3,1)], bty = "n", col = "transparent", xpd=NA, 
       text.col = st_color[c(2,4,3,1)], text.font = 2, cex = 1.5, y.intersp = 1.4)

#legend 2
par(mar=c(0,6,0,0))
plot0()

legend("center", 
       legend = c("No disturbances", "Moderate disturbances", "Major disturbances"), 
       bty = "n",  cex = 1.4, lty = c(1,2,3), lwd = 2)
dev.off()


##---------
## TIME * HARVEST ####

plot(Trans_sf["PP"], pch = 19, cex = .5)
eff_M_harv <- Effect(focal.predictors = c("time_interv", "harvest_100"),
                    mod = M_multi_full, 
                    xlevels = list(time_interv = seq(-1.34, 1.78, len = 100),
                                   harvest_100 = c(-0.3174316, 2.1, 4.46)))

# c(-0.3174316, 0.9, 3)* attr(scale_harv100,"scaled:scale") + attr(scale_harv100,"scaled:center")
harv <- eff_M_harv$model.matrix[,"harvest_100"]
time <- eff_M_harv$model.matrix[,"time_interv"] * attr(scale_time,"scaled:scale") +
  attr(scale_time,"scaled:center")

effdf_M_harv <- cbind.data.frame(time,
                                harv, 
                                eff_M_harv$prob, 
                                eff_M_harv$lower.prob,
                                eff_M_harv$upper.prob)

png(paste0("images/effect_harv0.png"), res=300, 
    width = 7.2, height = 5, units = 'in', bg = "transparent", type='cairo')

# quartz(width = 7.2, height = 5)
layout(mat = matrix(c(1,1,2,1,1,3,4,4,0),3, byrow = T), heights = c(1,1, .4), 
       widths = c(1,1,.9))

par(mar=c(4,4.5,.1,.1), xaxs = "i")

plot.default(effdf_M_harv[,1:2], type = "n", ylim = c(0,1), 
             xaxt="n", yaxt="n", ylab = "", xlab = "")

axis(1, tcl = -0.3, col.ticks = "grey60", cex.axis=1.5)

axis(2, tcl = -0.3, col.ticks  = "grey60", las = 1, cex.axis=1.5)

# loop over "to" transition
for(st_to in st_name) {
  col.st <- st_color[which(st_name==st_to)]
  tmp0 <- effdf_M_harv[which(effdf_M_harv$harv<0),]
  tmp1 <- effdf_M_harv[which(effdf_M_harv$harv==2.1),]
  tmp2 <- effdf_M_harv[which(effdf_M_harv$harv==4.46),]
  
  lines(tmp0[,c('time', paste0("prob.", st_to))], col = col.st, lwd = 2, lty = 1)
  # polygon(c(tmp0[,'time'],
  #           rev(tmp0[,'time'])),
  #         c(tmp0[,paste0("L.prob.", st_to)],
  #           rev(tmp0[,paste0("U.prob.", st_to)])),
  #         col = alpha(col.st, 0.08), border = NA)
  
  lines(tmp1[, c('time', paste0("prob.", st_to))], col = col.st, lwd = 2, lty = 2)
  lines(tmp2[, c('time', paste0("prob.", st_to))], col = col.st, lwd = 2, lty = 3)
  
}  
mtext("Time since first inventory (years)", side = 1, line = 3,  cex = 1.1, font =2)
mtext("Probability of transition", side = 2, line = 3,  cex = 1.1, font =2, las = 0)

# map
par(mar=c(5,0,0,0), xaxs = "r")
plot(st_geometry(ecoregion), border = "gray60",
     key.pos=NULL, main=NULL)
ordhar <- order(Trans_sf$harvest_intensity100)
mypal <- colorRampPalette(c("#FFEDA0","#FEB24C", "#b71901"))
levpal <- cut(as.numeric(as.character(Trans_sf$harvest_intensity100[ordhar])),50)
plot(st_geometry(Trans_sf)[ordhar], cex = .2, pch = 19, 
     col = mypal(50)[levpal], add = T)

bx <- par("usr")*.9
xs <- seq(bx[1],bx[2], len = 50)
for (i in 1:50) {
  polygon(c(xs[i], xs[i+1], xs[i+1], xs[i]), c(bx[3],bx[3],bx[3]+50000,bx[3]+50000),
          col = mypal(50)[i], border = NA)
}
mtext(side =1, "Harvest intensity", cex = .8, line = .5)


# legend 1
par(mar=c(5,0,0,0))
plot0()
legend("left", legend = st_name[c(2,4,3,1)], bty = "n", col = "transparent", xpd=NA, 
       text.col = st_color[c(2,4,3,1)], text.font = 2, cex = 1.5, y.intersp = 1.4)

#legend 2
par(mar=c(0,6,0,0))
plot0()
legend("center", 
       legend = c("Non harvested", "25% harvested", "75% harvested"), 
       bty = "n",  cex = 1.4, lty = c(1,2,3), lwd = 2)

dev.off()

## CONCLU ####
png(paste0("images/conclu_age_dist.png"), res=300, 
    width = 7.8, height = 3.7, units = 'in', bg = "transparent", type='cairo')

par(mfrow=c(1,2))
par(mar=c(3.5,3.5,.1,.1), xaxs="i")

plot.default(effdf_M_nat[,1:2], type = "n", ylim = c(0,1), 
             xaxt="n", yaxt="n", ylab = "", xlab = "")
axis(1, tcl = -0.3, col.ticks = "grey60", cex.axis=1)
axis(2, tcl = -0.3, col.ticks  = "grey60", las = 1, cex.axis=1)

# loop over "to" transition
for(st_to in st_name) {
  col.st <- st_color[which(st_name==st_to)]
  tmp0 <- effdf_M_nat[which(effdf_M_nat$nat_disturb1==0 & effdf_M_nat$nat_disturb2==0),]
  tmp1 <- effdf_M_nat[which(effdf_M_nat$nat_disturb1==1),]
  
  lines(tmp0[,c('time', paste0("prob.", st_to))], col = col.st, lwd = 2, lty = 1)
  lines(tmp1[, c('time', paste0("prob.", st_to))], col = col.st, lwd = 2, lty = 2)
}  
mtext("Time since first inventory (years)", side = 1, line = 2.5,  cex = 1.1, font =2)
mtext("Probability of transition", side = 2, line = 2.5,  cex = 1.1, font =2, las = 0)

legend("topright", legend = c("No disturbances", "Moderate disturbances"),
       lty = 1:2, cex = .8, bty = "n", lwd = 1.5)


plot.default(effdf_M_age[,1:2], type = "n", ylim = c(0,1), 
             xaxt="n", yaxt="n", ylab = "", xlab = "")
axis(1, tcl = -0.3, col.ticks = "grey60", cex.axis=1)
axis(2, tcl = -0.3, col.ticks  = "grey60", las = 1, cex.axis=1)
for(st_to in st_name) {
  col.st <- alpha(st_color[which(st_name==st_to)], 1)
  
  lines(effdf_M_age[,c('age', paste0("prob.", st_to))], col = col.st, lwd = 2)
  polygon(c(effdf_M_age[,'age'], 
            rev(effdf_M_age[,'age'])), 
          c(effdf_M_age[,paste0("L.prob.", st_to)], 
            rev(effdf_M_age[,paste0("U.prob.", st_to)])), 
          col = alpha(col.st, 0.1), border = NA)
}
mtext("Stand age", side = 1, line = 2.5,  cex = 1.1, font =2)
dev.off()


## methods ####
png(paste0("images/diag_trans_mixed2.png"), res=300, 
    width = 8, height = 4.5, units = 'in', bg = "transparent", type='cairo')
# quartz(width=8, height=4.5)
par(fig = c(0,1,0,1))
par(mar=c(3.5,3.5, 1,15), xaxs="i")
plot.default(effdf_M_cc[,1:2], type = "n", ylim = c(0,1), bty = "l",
             xaxt="n", yaxt="n", ylab = "", xlab = "")

axis(1, tcl = -0.3, col.ticks = "grey60", cex.axis=1)
axis(2, tcl = -0.3, col.ticks  = "grey60", las = 1, cex.axis=1)

# loop over "to" transition
for(st_to in st_name) {
  col.st <- alpha(st_color[which(st_name==st_to)], 1)
  tmp0 <- effdf_M_cc[which(effdf_M_cc$pp<0),]
  lines(tmp0[,c('tp', paste0("prob.", st_to))], col = col.st, lwd = 2, lty = 1)
}

mtext("Predictor", side = 1, line = 2.5,  cex = 1.1, font =2)
mtext("Probability of transition", side = 2, line = 2.5,  cex = 1.1, font =2, las = 0)

par(fig = c(0.5, 1, 0.35, 1), new = T)
par(mar=c(0,0,0,0), xaxs="r")
pos.box <- cbind (c(0.5, 0.2, 0.8, 0.5), 
                  c(0.8, 0.5, 0.5, 0.2))
lwd.mat <- matrix(c(2,2,2,2,
                    4,4,4,4, 
                    2,2,2,2,
                    2,2,2,2), 4, 4)
col.mat <- matrix(c("grey80","grey80","grey80","grey80",
                    "grey40","grey40","grey40","grey40",
                    "grey80","grey80","grey80","grey80",
                    "grey80","grey80","grey80","grey80"), 4, 4) 

plotmat(M.perc, pos = pos.box, curve = 0.07, name = st_name, lwd = lwd.mat, relsize=.98,
        box.cex = .9, cex.txt = 0, txt.col = "white", dtext = 0.15, txt.font = 2,
        box.lwd = 0.1, box.type = "rect", shadow.size = 0.005,
        box.prop = 0.4, box.size = 0.11, box.col = st_color,
        arr.length=.2, arr.width=.15,  arr.type ="curved",
        arr.col =col.mat, arr.lcol = col.mat,
        self.cex = 0.6, self.shifty = c(.07,0,0,-.07), self.shiftx = c(0,-.14,.14,0), 
        self.lwd = c(2,4,2,2))


dev.off()

#11. Predicting probability of transition with special interest for temperature and harvesting


## With CI

prob.df.all <- data.frame()
for(state in st_name){
  
  # predict probability of transition with newdata
  my.eff <- Effect(focal.predictors = c("time_interv","time2harvest", "TP"),
                   mod = stepmod[[state]],
                   xlevels = list(time2harvest = c(0, 0.6),
                                  TP = c(-0.9185971,  0.1046585,  0.9784375),
                                  time_interv = seq(-1.5, 2.5, length.out = 200)),
                   fixed.predictors = c(delta_TPy = mean(Trans_df$delta_TPy),
                                        PP = 0,
                                        delta_PPy = 0,
                                        age_mean = 0.5408876, #mean(Trans_df$age_mean),
                                        time2insect = 0, time2fire = 0))

  TP <- my.eff$model.matrix[,"TP"] * attr(scale_tp,"scaled:scale") + attr(scale_tp,"scaled:center")
  time_interv <- my.eff$model.matrix[,"time_interv"]*attr(scale_time,"scaled:scale") + attr(scale_time,"scaled:center")
  
  
  eff_df <- cbind.data.frame(TP, 
                             time = time_interv, 
                             harvest01 = my.eff$model.matrix[,"time2harvest"], 
                             my.eff$prob, 
                             my.eff$lower.prob, 
                             my.eff$upper.prob,
                             From = state)

  prob.df.all <- bind_rows(eff_df, prob.df.all)
  
}

prob.df.all <- prob.df.all %>% replace_na(replace=list(prob.Boreal=0, L.prob.Boreal=0, U.prob.Boreal=0,
                                                       prob.Temperate=0, L.prob.Temperate=0, U.prob.Temperate=0))



#12. Plot predictions

# setting layout
m <- matrix(c(1,1,1,0, 2:5, 6:9, 10:13, 14:17, 18,18,18,0), ncol = 4,byrow = T)
m <- cbind(c(0, rep(19,4), 0), m)

TPzone <- c("Warm (TP = 3.5)", "Mid (TP = 1.9)", "Cold (TP = 0.7)") 

col.tr <- factor(st_name, levels = st_name, labels = st_color)


quartz(width = 7, height = 7)
layout(m, heights= c(0.2, 1, 1, 1, 1, 0.4), widths = c(0.1, 1, 1, 1, 0.1))

# region title
par(mar=c(0,0,0,0), xaxs="i")
plot(0:3,0:3, type="n", axes = F)
for(i in 1:3) text(i-0.45, 1, TPzone[i], cex = 1.2, font = 2)

# loop over "from" transition
for(st_from in st_name) {
  tmp.from <- subset(prob.df.all, From == st_from)
  #loop over mean temperature
  for(tp in rev(unique(prob.df.all$TP))) {
    tmp.r <- tmp.from[which(tmp.from$TP == tp),]
    par(mar=c(1.8,1.8,0,0))
    plot(tmp.r[, c('time', paste0("prob.",st_from))], type = "n", ylim = c(0,1), 
         xaxt="n", yaxt="n", ylab = "", xlab = "", bty = "l")
    axis(1, labels = F, tcl = -0.3, col.ticks = "grey60")
    axis(2, labels = F, tcl = -0.3, col.ticks  = "grey60")
    if(st_from == "Temperate") axis(1, cex = 0.8, col = "grey60", tick=F, line = -.4)
    if(tp == rev(unique(prob.df.all$TP))[1]) axis(2, las = 1, cex = 0.8, col = "grey60", tick=F, line = -.4)
    # loop over "to" transition
    for(st_to in st_name) {
      col.st <- alpha(col.tr[which(st_name==st_to)], 1)
      
      tmp.h0 <- tmp.r[which(tmp.r$harvest01==0), ]
      lines(tmp.h0[,c('time', paste0("prob.",st_to))], col = col.st, lwd = 1.5)
      polygon(c(tmp.h0[,'time'], rev(tmp.h0[,'time'])), 
              c(tmp.h0[,paste0("L.prob.",st_to)], rev(tmp.h0[,paste0("U.prob.",st_to)])), 
              col = alpha(col.st,0.1), border = NA)
      
      tmp.h1 <- tmp.r[which(tmp.r$harvest01>0), ]
      lines(tmp.h1[,c('time', paste0("prob.",st_to))], col = col.st, lwd = 1.5, lty = 2)
      polygon(c(tmp.h1[,'time'], rev(tmp.h1[,'time'])), 
              c(tmp.h1[,paste0("L.prob.",st_to)], rev(tmp.h1[,paste0("U.prob.",st_to)])), 
              col = alpha(col.st,0.1), border = NA)
    }
    
  }
  par(mar=c(1.8,0,0,0))
  plot(1,1, type = "n", axes =F)
  text(1, 1, paste0("From ", st_from), cex = 1.1, font = 2, srt = 270)
}
par(mar=c(0,0,0,0))
plot.new()
mtext("Time (years since the first inventory)", side = 3, line = -1.5, adj = 0.5, cex = 0.9)
legend("bottom", legend = c(st_name, "non-harvested", "harvested"), 
       lty = c(1,1,1,1,1,2), lwd = 2, col = c(alpha(col.tr,1), "black", "black"), 
       bty = "n", ncol = 3)
plot.new()
mtext("Probability of transition", side = 2, line = -1.2,  cex = 0.9)





## Without CI ###

newdata <- data.frame(expand.grid(TP = c(-0.7518135, 0.4647507, 1.681315), 
                                  # quantile(Trans_df$TP, c(.1,.5,.9))
                                  delta_TPy = 0, #changing ∆TPy doesn't make difference in prediction
                                  PP = 0,
                                  delta_PPy = 0,
                                  harvest01 = c(0,1), 
                                  age_mean = 0.5408876, 
                                  outbreak2 = 0,
                                  burn2 = 0,
                                  windfall2 = 0,
                                  TYPEHUMUS = "MR",
                                  time_interv = seq(-2, 2, length.out = 200)))

prob.df.all <- data.frame()
for(state in st_name){
  
  # predict probability of transition with newdata
  prob <- as.data.frame(predict(stepmod[[state]], 
                                newdata = newdata, type = "probs", se.fit=T))
  
  TP <- newdata$TP * attr(scale_tp,"scaled:scale") + attr(scale_tp,"scaled:center")
  time_interv <- newdata$time_interv*attr(scale_time,"scaled:scale") + attr(scale_time,"scaled:center")
  
  prob.df <- cbind.data.frame(TP, 
                              time = time_interv, 
                              harvest01 = newdata$harvested, 
                              prob, From = state)
  prob.df.all <- bind_rows(prob.df, prob.df.all)
  
}

prob.df.all<- prob.df.all %>% replace_na(replace=list(Boreal=0, Temperate=0))


#12. Plot predictions

# setting layout
m <- matrix(c(1,1,1,0, 2:5, 6:9, 10:13, 14:17, 18,18,18,0), ncol = 4,byrow = T)
m <- cbind(c(0, rep(19,4), 0), m)

subzone <- c("Warm (TP = 3.5)", "Mid (TP = 1.9)", "Cold (TP = 0.7)") 

col.tr = factor(st_name, levels = st_name, 
                labels = st_color)


quartz(width = 7, height = 7)
layout(m, heights= c(0.2, 1, 1, 1, 1, 0.4), widths = c(0.1, 1, 1, 1, 0.1))

# region title
par(mar=c(0,0,0,0), xaxs="i")
plot(0:3,0:3, type="n", axes = F)
for(i in 1:3) text(i-0.45, 1, subzone[i], cex = 1.2, font = 2)

# loop over "from" transition
for(st_from in st_name) {
  tmp.from <- subset(prob.df.all, From == st_from)
  #loop over mean temperature
  for(tp in rev(unique(prob.df.all$TP))) {
    tmp.r <- tmp.from[which(tmp.from$TP == tp),]
    par(mar=c(1.8,1.8,0,0))
    plot(tmp.r[, c('time', st_from)], type = "n", ylim = c(0,1), 
         xaxt="n", yaxt="n", ylab = "", xlab = "", bty = "l")
    axis(1, labels = F, tcl = -0.3, col.ticks = "grey60")
    axis(2, labels = F, tcl = -0.3, col.ticks  = "grey60")
    if(st_from == "Temperate") axis(1, cex = 0.8, col = "grey60", tick=F, line = -.4)
    if(tp == rev(unique(prob.df.all$TP))[1]) axis(2, las = 1, cex = 0.8, col = "grey60", tick=F, line = -.4)
    # loop over "to" transition
    for(st_to in st_name) {
      col.st <- alpha(col.tr[which(st_name==st_to)], 1)
      lines(tmp.r[which(tmp.r$harvest01<0), c('time', st_to)], col = col.st, lwd = 1.5)
      lines(tmp.r[which(tmp.r$harvest01>0), c('time', st_to)], col = col.st, lty = 2, lwd = 1.5)
    }
    
  }
  par(mar=c(1.8,0,0,0))
  plot(1,1, type = "n", axes =F)
  text(1, 1, paste0("From ", st_from), cex = 1.1, font = 2, srt = 270)
}
par(mar=c(0,0,0,0))
plot.new()
mtext("Time (years since the first inventory)", side = 3, line = -1.5, adj = 0.5, cex = 0.9)
legend("bottom", legend = c(st_name, "non-harvested", "harvested"), 
       lty = c(1,1,1,1,1,2), lwd = 2, col = c(alpha(col.tr,1), "black", "black"), 
       bty = "n", ncol = 3)
plot.new()
mtext("Probability of transition", side = 2, line = -1.2,  cex = 0.9)



## ---- image_prob1 ----

pal.image <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))

figs_state <- matrix(c('Boreal','Pioneer',
                       'Boreal','Mixed',
                       "Mixed", "Pioneer",
                       'Mixed','Temperate',
                       'Pioneer','Temperate'),
                     ncol=2,nrow=5,byrow=TRUE)

# layout
nf <- matrix(c(1,2,0, 3:17), 6, 3, byrow = T)
nf <- cbind(c(0,rep(18,5)), nf)
nf <- rbind(nf, c(0,19,19,0))

quartz(height = 9, width = 5.2)
png("~/Desktop/pred_TSH_time.png",res=300, width = 5.2, height = 9, units = 'in', bg = "transparent",type='cairo')
layout(nf, widths = c(0.15,1,1,0.4), heights = c(0.1,1,1,1,1,1,0.15))
par(mar=c(0,0,0,0))

plot.new()
mtext("Non-harvested", side = 3, line = -1.1, adj = 0.5, cex = 0.8)
plot.new()
mtext("Harvested", side = 3, line = -1.1, adj = 0.5, cex = 0.8)

for(i in 1:nrow(figs_state)) {
  
  par(mar=c(2,1.8,0,0))
  for(h in c(0,0.6)) {
    fig_prob_time(figs_state[i,1],figs_state[i,2], harvest = h, axes = F)
    axis(1, labels = F, tcl = -0.3, col.ticks = "grey60")
    axis(2, labels = F, tcl = -0.3, col.ticks  = "grey60")
    if(h == 0) axis(2, las = 1, cex = 0.8, col = "grey60", line = -.4, tick=F)
    axis(1, cex = 0.8, col = "grey60", tick=F, line = -.4)
    box()
  }
  par(mar=c(0,0.3,0,0))
  plot.new()
  mtext(paste0("From ",figs_state[i,1], "\nTo ", figs_state[i,2]), side = 3, line = -2.4, adj = 0, cex = 0.7)
}
par(mar=c(0,0,0,0))
plot.new()
mtext("Annual Temperature", side = 2, line = -1.2, cex = 0.9)
plot.new()
mtext("Time (years)", side = 1, line = -1.2, cex = 0.9)

dev.off()


######################################################################
### MODEL WITH HARVEST AS A CONTINUOUS VARIABLE ####
######################################################################

#4. Evaluate contribution of the coef: build all formula

## ---- formulas2 ----

coefs <-c("TP" , "I(TP^2)", "delta_TP", "I(delta_TP^2)",
          "PP" , "I(PP^2)", "delta_PP", "I(delta_PP^2)", 
          "TP:PP", "TP:delta_TP", 
          "harvested", "stress_mecha", 
          "fire", "insect",
          "age_mean",
          "plantation", 
          "time")
contrib_form <- as.character()
tested_coefs <- as.character()

for(i in 1:length(coefs)){
  temp_coefs <- coefs[-i]
  tested_coefs  <- append(tested_coefs, coefs[i])
  form <- paste0("To ~ ", paste(temp_coefs, collapse = "+"))
  contrib_form <- append(contrib_form, form)
}
# Add full model
full_formula <- paste0("To ~ ", paste(coefs, collapse = "+"))
null_formula <- "To ~ 1"

# Create DF of formula
formulas <- data.frame(tested = tested_coefs, formula = contrib_form, stringsAsFactors = FALSE)

## ---- multinom_trans2 ----

#5.1 Full Models and pseudo-R2
Trans_df$To <- relevel(Trans_df$To, ref = "Boreal")
B_multi_full <- nnet::multinom(full_formula, data = Trans_df, subset = From=="Boreal", trace = F)

Trans_df$To <- relevel(Trans_df$To, ref = "Mixed")
M_multi_full <- nnet::multinom(full_formula, data = Trans_df, subset = From=="Mixed", trace = F)

Trans_df$To <- relevel(Trans_df$To, ref = "Pioneer")
P_multi_full <- nnet::multinom(full_formula, data = Trans_df, subset = From=="Pioneer", trace = F)

Trans_df$To <- relevel(Trans_df$To, ref = "Temperate")
T_multi_full <- nnet::multinom(full_formula, data = Trans_df, subset = From=="Temperate", trace = F)

#5.2 Null Models and pseudo-R2
for(state in st_name) {
  st <- substring(state, 1, 1)
  # Null Models
  tmp.null <- nnet::multinom(null_formula, data=Trans_df, subset= From==state, trace = F)
  assign(paste0(st, '_multi_null'), tmp.null)
  # Compute McFadden pseudo-R2
  tmp.r2 <- 1 - (get(paste0(st, '_multi_full'))$deviance / tmp.null$deviance)
  assign(paste0(st, '_Rs'), tmp.r2)
}


## ---- multinom_full_contrib2 ----

#6. Compute contributions of each variable to the full model
ls_contrib <-list()

for(state in st_name){
  contrib <- as.numeric()
  st <- substring(state, 1, 1)
  for(i in 1:nrow(formulas)){
    multi <- nnet::multinom(formulas$formula[i], data = Trans_df, subset = From==state, trace = F)
    contrib <- append(contrib, (multi$AIC - get(paste0(st, '_multi_full'))$AIC))
  }
  ls_contrib[[state]] <- c(contrib, get(paste0(st, '_multi_full'))$AIC, get(paste0(st, '_Rs')))
}

# Summarize information in table
contrib.full <- do.call("rbind", ls_contrib)
colnames(contrib.full) <- c(gsub('[I(^)]', '', formulas$tested), "AIC", "R2")


## ---- multinom_stepAIC2 ----

#7. stepAIC
stepmod <- list()
for(state in st_name) {
  st <- substring(state, 1, 1)
  Trans_df$To <- relevel(Trans_df$To, ref = state)
  tmp.step <- stepAIC(get(paste0(st, '_multi_full')), trace = FALSE)
  stepmod[[state]] <- tmp.step
  # Pseudo R2
  tmp.r2 <- 1 - (tmp.step$deviance / get(paste0(st, '_multi_null'))$deviance)
  assign(paste0(st, '_step_Rs'), tmp.r2)
}

## ---- multinom_reduce_contrib2 ----

#8 variable contribution for reduced models

ls_contrib <-list()

for(state in st_name){
  contrib <- as.numeric()
  st <- substring(state, 1, 1)
  
  contrib_form <- as.character()
  tested_coefs <- as.character()
  coefs_sel <- attr(terms(stepmod[[state]]), "term.labels")
  for(i in 1:length(coefs_sel)){
    temp_coefs <- coefs_sel[-i]
    tested_coefs  <- append(tested_coefs, coefs_sel[i])
    form <- paste0("To ~ ", paste(temp_coefs, collapse = "+"))
    contrib_form <- append(contrib_form, form)
  }
  form.sel <- data.frame(tested = tested_coefs, formula = contrib_form, stringsAsFactors = FALSE)
  
  contrib.sel <- t(data.frame(row.names = coefs, deltaAIC = rep(NA, length(coefs))))
  row.names(contrib.sel) <- state
  for(i in 1:nrow(form.sel)){
    l <- which(coefs %in% form.sel[i,1])
    multi <- nnet::multinom(form.sel$formula[i], data = Trans_df, subset = From==state, trace = F)
    contrib.sel[,l] <- multi$AIC - stepmod[[state]]$AIC
  }
  
  ls_contrib[[state]] <- cbind(contrib.sel, 
                               AIC = stepmod[[state]]$AIC, 
                               R2 = get(paste0(st, '_step_Rs')))
}

# Summarize information in table
contrib.sel <- do.call("rbind", ls_contrib)
colnames(contrib.sel) <- c(gsub('[I(^)]', '', formulas$tested), "AIC", "R2")

## ---- coeff_plot ----

#### COEFPLOT ####

# Get coefficient and confidence interval
coef_stepmod <- lapply(stepmod, FUN = function(x) coefficients(x)[,-1])
CI_stepmod <- lapply(stepmod, FUN = function(x) confint(x)[-1,,])

# New coefficient names
newnames <- c("TP", "TP2", "∆TP", "∆TP2", 
              "PP", "PP2", "∆PP", "∆PP2",
              "harvest", "stress", "fire", "insect", "age", "plantation", "time",
              "TP:PP", "TP:∆TP")

newnames <- cbind.data.frame("coefs" = colnames(coefficients(M_multi_full))[-1], newnames)

# Plot of coefficient
quartz()
par(mfrow = c(2,2), mar = c(3,5,2,1))
for(st_from in st_name) {
  # Coefficient
  coef_st <- coef_stepmod[[st_from]]
  # Confidence interval
  CI_st <- CI_stepmod[[st_from]]
  rangeCI <- range(CI_st)
  # Number of variables
  nvar <- ncol(coef_st)
  # State_to
  nst <- row.names(coef_st)
  # Empty plot
  plot(1, type="n", xlab="", ylab="", 
       xlim=c(rangeCI[1],rangeCI[2]), ylim=c(1, nvar), yaxt = "n", 
       main = paste0("From ", st_from))
  abline(v = 0, lty = 2, col = "grey25")
  abline(h = 1:nvar, lty = 3, col = "grey")
  # Labels
  coef_labs <- subset(newnames, coefs %in% colnames(coef_st))$newnames
  axis(side = 2, at = 1:nvar, labels = rev(coef_labs), las = 1, cex = 0.9)
  l = -0.12
  for(st_to in nst) {
    points(coef_st[st_to,], nvar:1+l, pch = 19, col = st_color[which(st_to == st_name)])
    for (i in 1:nvar) {
      yl <- nvar +1 - i
      lines(x=c(CI_st[i,1,st_to], CI_st[i,2,st_to]), y=c(yl+l, yl+l), 
            col = st_color[which(st_to == st_name)])
    }
    l = l + 0.12
    
  }
  
}
## p-value for coefficients
z <- summary(stepmod$Mixed)$coefficients/summary(stepmod$Mixed)$standard.errors
# 2-tailed Wald z tests to test significance of coefficients
p <- (1 - pnorm(abs(z), 0, 1)) * 2
p

afex::set_sum_contrasts() # use sum coding, necessary to make type III LR tests valid

car::Anova(stepmod$Mixed,type="III")

## ---- multinom_pred2 ----

#11. Predicting probability of transition with special interest for temperature and harvesting

newdata <- data.frame(expand.grid(TP = c(0.05, 0.35, 0.9), # from mid to warm temperature
                                  delta_TP = seq(-2.5, 2.5, length.out = 200),
                                  PP = 0,
                                  delta_PP = 0,
                                  harvested = seq(-0.22, 12, length.out = 10), 
                                  stress_mecha = 0,
                                  age_mean = median(Trans_df$age_mean),
                                  insect = as.factor(0), 
                                  plantation = as.factor(0), 
                                  fire = as.factor(0),
                                  time = 10))


prob.df.all <- data.frame()
for(state in st_name){
  
  # predict probability of transition with newdata
  prob <- as.data.frame(predict(stepmod[[state]], 
                                newdata = newdata, type = "probs", se.fit=T))
  
  TP <- newdata$TP * attr(scale_tp,"scaled:scale") + attr(scale_tp,"scaled:center")
  delta_TP <- newdata$delta_TP*attr(scale_delta_tp,"scaled:scale") + attr(scale_delta_tp,"scaled:center")
  harvest_unscale <- newdata$harvested*attr(scale_harv,"scaled:scale") + attr(scale_harv,"scaled:center")
  
  prob.df <- cbind.data.frame(TP, delta_TP, harvest=harvest_unscale, 
                              prob, From = state)
  prob.df.all <- bind_rows(prob.df, prob.df.all)
  
}

prob.df.all <- prob.df.all %>% replace_na(replace=list(Boreal=0, Temperate=0))

#12. Plot predictions

pal.h <- rev(brewer.pal(10, "RdYlBu"))

# setting layout
m <- matrix(c(1:4, 0, 5:24, 25,25,25,25,0), ncol = 5,byrow = T)
m <- cbind(c(0, rep(26,4), 0), m, c(0, rep(27,4), 0))

# Choose one temperature
tp <- rev(unique(prob.df.all$TP))[2]
prob.df.tp <- prob.df.all[which(prob.df.all$TP == tp),] 

harv_level <- round(unique(prob.df.all$harvest)*attr(scale_harv,"scaled:scale") + attr(scale_harv,"scaled:center"),0)

quartz(width = 8, height = 7)
layout(m, heights= c(0.2, 1, 1, 1, 1, 0.14), widths = c(0.15, 1, 1, 1, 1, 0.1, 0.6))

# region title
par(mar=c(0,0,0,0), xaxs="i")
for(i in 1:4) {
  plot(1,1, type = "n", axes =F)
  text(1, 1,paste0("To ", st_name[i]), cex = 1.1, font = 2)
}

# loop over "from" transition
for(st_from in st_name) {
  tmp.from <- subset(prob.df.tp, From == st_from)
  
  # loop over "to" transition
  par(mar=c(1.8,1.8,0,0))
  for(st_to in st_name) {
    
    plot(tmp.from[, c('delta_TP', st_to)], type = "n", ylim = c(0,1), 
         xaxt="n", yaxt="n", ylab = "", xlab = "", bty = "l")
    axis(1, labels = F, tcl = -0.3, col.ticks = "grey60")
    axis(2, labels = F, tcl = -0.3, col.ticks  = "grey60")
    if(st_from == "Temperate") {
      axis(1, cex = 0.8, col = "grey60", tick=F, line = -.4)
    }
    if(st_to == "Boreal") {
      axis(2, las = 1, cex = 0.8, col = "grey60", tick=F, line = -.4)
    }
    # loop over harvesting intensity
    
    for(h in unique(tmp.from$harvest)) {
      lines(tmp.from[which(tmp.from$harvest==h), c('delta_TP', st_to)], 
            col = pal.h[which(unique(tmp.from$harvest) == h)], lwd = 1.5)
    }
  }
  par(mar=c(1.8,0,0,0))
  plot(1,1, type = "n", axes =F)
  text(1, 1, paste0("From ", st_from), cex = 1.1, font = 2, srt = 270)
}

par(mar=c(0,0,0,0))
plot.new()
mtext("Temperature change (°C)", side = 3, line = -1.5, adj = 0.5, cex = 0.9)
plot.new()
mtext("Probability of transition", side = 2, line = -1.2,  cex = 0.9)
plot.new()
legend("center", title = "Harvesting intensity\nnb of trees cut",
       legend = harv_level, 
       lty = 1, lwd = 2, col = pal.h, bty = "n")


## ---- image_prob2 ----

pal.image <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))

figs_state <- matrix(c('Boreal','Pioneer',
                       'Boreal','Mixed',
                       "Mixed", "Pioneer",
                       'Mixed','Temperate',
                       'Pioneer','Temperate'),
                     ncol=2,nrow=5,byrow=TRUE)
#layout
nf <- matrix(c(1,2,3,0, 4:23), 6, 4, byrow = T)
nf <- cbind(c(0,rep(24,5)), nf)
nf <- rbind(nf, c(0,25,25,25,0))

quartz(height = 7.4, width = 6.6)

layout(nf, widths = c(0.15,1,1,1,0.45), heights = c(0.1,1,1,1,1,1,0.15))

par(mar=c(0,0,0,0))
plot.new()
mtext("No harvesting", side = 3, line = -1, adj = 0.5, cex = 0.8)
plot.new()
mtext("Moderate harvesting", side = 3, line = -1, adj = 0.5, cex = 0.8)
plot.new()
mtext("Intense harvesting", side = 3, line = -1, adj = 0.5, cex = 0.8)

for(i in 1:nrow(figs_state)) {
  par(mar=c(1.7,1.7,0,0))
  for(h in c(0, 6, 10)) {
    fig_prob_ima(figs_state[i,1], figs_state[i,2], harvest = h, axes = F)
    axis(1, labels = F, tcl = -0.3, col.ticks = "grey60")
    axis(2, labels = F, tcl = -0.3, col.ticks  = "grey60")
    if(h == 0) axis(2, las = 1, cex = 0.8, col = "grey60", line = -.4, tick=F)
    if(i == 5) axis(1, cex = 0.8, col = "grey60", tick=F, line = -.4)
    box()
  }
  par(mar=c(0,0.3,0,0))
  plot.new()
  mtext(paste0("From ",figs_state[i,1], "\nTo ", figs_state[i,2]), side = 3, line = -2.4, adj = 0, cex = 0.7)
}
par(mar=c(0,0,0,0))
plot.new()
mtext("Temperature change (°C)", side = 2, line = -1.2, cex = 0.9)
plot.new()
mtext("Mean annual temperature (°C)", side = 1, line = -1.2, cex = 0.9)

#### Modeling change in temperate species proportion ###

## ---- temperate_prop ----

state.prop <- sp.mat %>%
  mutate(Boreal = rowSums(dplyr::select(., boreal))) %>%
  mutate(Pioneer = rowSums(dplyr::select(., pioneer))) %>%
  mutate(Temperate = rowSums(dplyr::select(., temperate))) %>%
  mutate(Mixed = rowSums(dplyr::select(., mixed))) %>%
  dplyr::select(plot_id, year_measured, Boreal, Pioneer, Temperate, Mixed, TOTAL)

# Remove boreal plots
plot_boreal <- unique(subset(gr.prop, states == "Boreal", select = plot_id))

state.prop <- state.prop[which(!(state.prop$plot_id %in% plot_boreal$plot_id)),]

# Scale variable
bioclim[,c("TP","PP", "GSL")] <- scale(bioclim[,c("TP","PP", "GSL")])

state.prop <- state.prop %>%
  full_join(env_all_steps, by = c("plot_id", "year_measured")) %>%
  full_join(bioclim, by = c("plot_id", "year_measured")) %>%
  mutate(year_center = scale(year_measured)) # 1970 = year 0


mod1 <- lme4::glmer(cbind(Temperate, TOTAL - Temperate) ~ year_center + (1 + year_center|plot_id), 
                    data = state.prop, 
                    family = binomial,
                   control=lme4::glmerControl(optimizer="bobyqa"))
mod2 <- lme4::glmer(cbind(Temperate, TOTAL - Temperate) ~ TP + year_center + (1 + TP|plot_id), 
                    data = state.prop, 
                    family = binomial,
                    control=lme4::glmerControl(optimizer="bobyqa"))
mod3 <- lme4::glmer(cbind(Temperate, TOTAL - Temperate) ~ TP + PP + year_center + (1 + year_center|plot_id), 
                    data = state.prop, 
                    family = binomial,
                    control=lme4::glmerControl(optimizer="bobyqa"))
mod4 <- lme4::glmer(cbind(Temperate, TOTAL - Temperate) ~ TP*year_center + PP + harvested  + (year_center|plot_id), 
                    data = state.prop, 
                    family = binomial,
                    control=lme4::glmerControl(optimizer="bobyqa"))

mod5 <- lme4::glmer(cbind(Temperate, TOTAL - Temperate) ~ TP*year_center + PP*year_center + 
                      harvest01 + insect + age_mean +
                      (year_center|plot_id), 
                    data = state.prop, 
                    family = binomial, 
                    control=lme4::glmerControl(optimizer="bobyqa"))


anova(mod4, mod5)
library(effects)
my.eff <- Effect(c("year_center", "TP", "harvested"), mod5)
plot(my.eff)
plot(allEffects(mod5))

## ---- speed ----

#### SPEED ####

# PCA
pca <- rda(sp.hel)

spe.scores <- scores(pca, display="species", choices=c(1:3), scaling = 1) 
site.scores <- scores(pca, display="wa", choices=c(1:3), scaling = 1)

site.scores <- cbind.data.frame(sp.mat[,c("plot_id", "year_measured")], site.scores)


# Compute speed

delta_if <- as.data.frame(abs(aggregate(cbind(year_measured, PC1, PC2, PC3) ~ plot_id,
                                        site.scores, 
                                        FUN = function(x) last(x) - first(x))))
grTrans1 <- slice(group_by(grTrans, plot_id), 1)

distance_if <- sqrt(delta_if$PC1^2 + delta_if$PC2^2 + delta_if$PC3^2)

Speed_if <- cbind.data.frame("plot_id" = delta_if$plot_id, 
                             "rawSpeed" = distance_if/delta_if$year_measured,
                             "From" = Trans.map$From,
                             "To" = Trans.map$To,
                             "Transition" = paste0(substring(Trans.map$From, 1, 1), 
                                                   "-", substring(Trans.map$To, 1, 1)),
                             "Harvest" = as.factor(aggregate(harvest01 ~ plot_id, env_all_steps, 
                                                             FUN = function(x) ifelse(length(which(x==1))>0,1,0))[,2]))

Speed_if <- subset(Speed_if, Transition!="B-T" & Transition!="T-B")


ggplot(data = Speed_if, aes(x = rawSpeed, y = Transition, fill = Harvest)) +
  geom_density_ridges2(scale=1, alpha=0.3) +
  scale_fill_manual(values = c("#9C9797","#AB0505")) +
  labs(x = "Community speed", y = "State transition") +
  theme_ridges(font_size = 10) + 
  theme(axis.text=element_text(size=10), axis.title=element_text(size=11,face="bold"))

### MODEL ####

delta_clim4 <- slice(group_by(delta_clim, plot_id), n())

env_all_steps4 <- slice(group_by(env_all_steps, plot_id), n())

speed_df <- Speed_if %>%
  left_join(env_all_steps4, by = "plot_id") %>%
  left_join(delta_clim4, by = c("plot_id", "year1", "year2")) %>%
  ungroup(plot_id) %>%
  mutate(year_center = scale(year2)) # 1970 = year 0


mod1 <- lm(rawSpeed ~ year_center, 
           data = speed_df)
mod2 <- lm(rawSpeed ~ TP + year_center, 
           data = speed_df)
mod3 <- lm(rawSpeed ~ TP + PP + year_center,
           data = speed_df)
mod4 <- lm(rawSpeed ~ TP + PP + harvest01, 
           data = speed_df)

mod5 <- lm(rawSpeed ~ TP + delta_TP + Transition*harvest01 + 
             insect + fire + age_mean,
           data = speed_df)
mod6 <- lm(rawSpeed ~ TP + delta_TP + Transition*harvested +
             insect + fire + age_mean,
           data = speed_df)

anova(mod5, mod6)
AIC(mod4,mod5, mod6)
library(effects)
my.eff <- Effect(focal.predictors = c("Transition","harvested"), 
                 mod = mod6)
plot(my.eff)
plot(allEffects(mod6))

my.eff <- Effect(focal.predictors = c("Transition","harvest01"), 
                 mod = mod5,
                 # xlevels = list(TP=c(-0.9185971,  0.1046585,  0.9784375)),
                 fixed.predictors=list(given.values=c(TP = mean(speed_df$TP),
                                  delta_TP = mean(speed_df$delta_TP), 
                                  age_mean = median(speed_df$age_mean),
                                  insect1 = as.factor(0), 
                                  fire1 = as.factor(0))))

plot(my.eff)

str(my.eff)

eff_df <- data.frame(expand.grid(Transition = my.eff$variables$Transition$levels,
                                 harvested = c(-0.4,3,5,8,10)))


eff_df <- cbind.data.frame(eff_df, 
                           fit = my.eff$fit, 
                           lower = my.eff$lower, 
                           upper = my.eff$upper)

h.level <- c(-0.4,3,5,8,10)
col.h <- (brewer.pal(6, "YlOrRd"))[-1]
ntr <- length(my.eff$variables$Transition$levels)

quartz(width = 8, height = 4)
par(mar = c(4,4,1,1))
# Empty plot
plot(1, type="n", xlab="", ylab="", 
     xlim=c(1, ntr), ylim = c(min(eff_df$lower), max(eff_df$upper)), xaxt = "n", 
     cex.axis = 0.8, las =1)

l = -0.14

for(h in h.level) {
  tmp.h <- eff_df[which(eff_df$harvested==h),]
  
  points(1:ntr+l, tmp.h$fit, pch = 19, col = col.h[which(h.level==h)])
  
  for (i in 1:ntr) {
    lines(x = c(i+l, i+l), y = c(tmp.h[i, "lower"], tmp.h[i, "upper"]), 
          col = col.h[which(h.level==h)])
  }
  l = l + 0.07
}

axis(side = 1, at = 1:ntr, labels = my.eff$variables$Transition$levels, las = 1, cex = 0.8)
abline(v = 3.5, lty = 2)
abline(v = 7.5, lty = 2)
abline(v = 11.5, lty = 2)
mtext("Vitesse de changement", side = 2, line = 2.9, cex = 1.2)
mtext("Type de transitions", side = 1, line = 2.5, cex = 1.2)


## ---- latlon_shift ----

#### LATITUDINAL AND LONGITUDINAL SHIFT ####

coord_dist <- cbind.data.frame(Trans.map[,c("plot_id", "From", "To")], 
                               "Longitude" = st_coordinates(xy)[,1],
                               "Latitude" = st_coordinates(xy)[,2])

# ggplot(data = coord_dist, aes(x = Latitude, y=From)) +
#   geom_density_ridges2(scale=1.1,alpha=0.3) +
#   scale_fill_manual(values = c("#9C9797","#AB0505")) +
#   labs(x = "Latitude", y = "Group") +
#   theme_ridges()


latlim <- range(density(coord_dist$Latitude)$x)
lonlim <- range(density(coord_dist$Longitude)$x)
col1 <- alpha("#9C9797",0.6) #grey
border1 <- "#4f4e4e" #grey
col2 <- alpha("#AB0505", 0.2) #red
border2 <- "#300000" #red
  
#quartz()

nf <- matrix(1:8, ncol=2, byrow = T)
layout(nf)

par(mar=c(0,1,0,1),oma = c(4, 5, 0.5, 0.5), bty = "n")

for(reg in rev(levels(coord_dist$From))){
  tmp = subset(coord_dist, From == reg)
  tmp2 = subset(coord_dist, To == reg)
  
  # change in latitude density distribution
  d_lat <- density(tmp$Latitude, bw = 0.21)
  d_lat2 <- density(tmp2$Latitude, bw = 0.21)
  
  plot(d_lat, xlim = latlim, ylim = c(0,0.65), type = "n", 
       xaxt="n", yaxt="n", main = "")
  polygon(d_lat, col = col1, border = border1)
  polygon(d_lat2, col = col2, border = border2)
  
  axis(1, at = seq(44, 54, 2), labels = seq(44, 54, 2), outer = T, tick = F)
  abline(v = seq(44, 54, 2), xpd = T, col = alpha(c("#A6A6A6"), 0.5))
  abline(h = 0, xpd = F, col = alpha(c("#A6A6A6"), 0.5))
  
  mtext(reg, side = 2, las = 1, line = 0.1, at = 0.05)
  if(reg == "Temperate") mtext("Latitude", side = 1, line = 2, adj = 1)
  
  # change in longitute density distribution
  d_lon <- density(tmp$Longitude, bw = 0.83)
  d_lon2 <- density(tmp2$Longitude, bw = 0.83)
  
  plot(d_lon, xlim = lonlim, ylim = c(0,0.23), type = "n", 
       xaxt="n", yaxt="n", main = "")
  polygon(d_lon, col = col1, border = border1)
  polygon(d_lon2, col = col2, border = border2)
  
  axis(1, at = seq(-80, -60, 5), labels = seq(-80, -60, 5), outer = T, tick = F)
  abline(v = seq(-80, -60, 5), xpd = T, col = alpha(c("#A6A6A6"), 0.5))
  abline(h = 0, xpd = T, col = alpha(c("#A6A6A6"), 0.5))
  if(reg == "Temperate") mtext("Longitude", side = 1, line = 2, adj = 1)
  if(reg == "Boreal") legend("topright", legend = c("1970-1980", "2000-2010"), 
                             fill = c(col1, col2), border = c(border1, border2), bty = "n")

}

## ---- transition_coord_change ----
quant.xy.t1 <- aggregate(st_coordinates(xy), by = list(gr[l.t1]), FUN = function(x) quantile(x, .95))
quant.xy.t2 <- aggregate(st_coordinates(xy), by = list(gr[l.t2]), FUN = function(x) quantile(x, .95))

colnames(quant.xy.t1) <-  colnames(quant.xy.t2) <- c("Group.1", "longitude", "latitude")

quant.xy.t1 <- quant.xy.t1[order(mean.xy.t1$latitude),]
quant.xy.t2 <- quant.xy.t2[order(mean.xy.t2$latitude),]

quant.xy = cbind(quant.xy.t1, quant.xy.t2[,-1])
quant.xy$"∆Longitude" = quant.xy.t2$longitude - quant.xy.t1$longitude
quant.xy$"∆Latitude" = quant.xy.t2$latitude - quant.xy.t1$latitude

row.names(quant.xy) = quant.xy$Group.1 ; quant.xy = quant.xy[,-1]



# quartz()  
# par(mfcol=c(4,1), mar=c(3,3,4,1),oma = c(4, 4.2, 0.5, 0.5), xaxs="i", yaxs="i")
# 
# for(reg in rev(levels(Speed_if$From))){
#   tmp = subset(Speed_if, From == reg)
#   d <- density(tmp$rawSpeed)
#   plot(d,  main=reg)
#   polygon(d, col = alpha(tmp$col.tr,0.5), border="black")
#   tmp2 = subset(tmp, harvest=="coupés")
#   d2 <- density(tmp2$rawSpeed)
#   polygon(d2, col = alpha("#9C9797", 0.5))
# }
# mtext("Trajectory speed", side = 1, outer=T, line = 1)
# mtext("Density", side = 2, outer=T, line = 1)
# 

#### PCA SITE SCORES VS TIME ####

## ---- scores_vs_time ----
## Reg lin
# 
# plot(site.scores[,c("PC1", "PC2")],  type = 'n')
# points(site.scores[,c("PC1", "PC2")], col = alpha(col.gr, 0.5), pch = pch.gr, cex = 0.8)
# ordiellipse(site.scores[,c("PC1", "PC2")],  groups = gr, 
#             draw="polygon", alpha = 100,
#             col = levels(col.gr), border = "grey25", label =T)

site.scores.long <- gather(site.scores, PC, scores, PC1:PC2)

#quartz()
formula1 <- y ~ x
formula2 <- y ~ poly(x, 2, raw = TRUE)
ggplot(data = site.scores.long, aes(x = year_measured, y = scores, colour = harvest)) +
  geom_point(size= 0.5, alpha=0.3) +
  xlab("Year") +
  ylab("Score on PCA") + 
  geom_smooth(method=lm, formula = formula1) +
  scale_color_manual(values = c("#9C9797","#AB0505")) +
  theme_classic() + facet_grid(PC ~From) 


