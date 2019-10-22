### Plot best model ####

### PACKAGES ####
library(graphicsutils)
library(dplyr)
library(msm)
library(diagram)
library(latex2exp)

source('R/functions/markov.R')
### DATA ####

source('R/scripts/prep_data.R')

# Load msm results

load("res/msm_all75.rda")

msm_glb <- msm_all75[["msm_glb"]]


### Transition intensity diagram
qmat  <- t(coef(msm_glb)$baseline)
diag(qmat) <- 0 # remove self-transition
qmat <- round(qmat,4)

### letters

qmat_let <- c("q_{BM}", "q_{BP}",
          "q_{MB}", "q_{MP}", "q_{MT}",
          "q_{PB}", "q_{PM}", "q_{PT}",
          "q_{TM}", "q_{TP}")

pos.box <- cbind(c(0.5, 0.05, 0.95, 0.5), 
                 c(0.95, 0.5, 0.5, 0.05))

pdf("res/fig3_baseline.pdf", width = 4.7, height = 3.3)
#quartz(width = 4.7, height = 3.3)
par(mar=c(.5,0.5,.5,.5))
pm <- plotmat(qmat, pos = pos.box, curve = 0.07, name = states, 
              lwd = 1.2, relsize = .9,
              box.cex = 1, cex.txt = 0, txt.col = "white", 
              dtext = .55, txt.font = 2,
              box.lwd = 0.1, box.type = "rect", shadow.size = 0.005,
              box.prop = 0.35, box.size = 0.13, box.col = st_col,
              arr.length=qmat*10, arr.width=qmat*10,  arr.type ="triangle",
              arr.col ="grey40", arr.lcol = "grey40", 
              arr.lwd = qmat*250,
              self.cex = 0.6, self.lwd = 1.2, 
              self.shifty = c(.07,0,0,-.07), self.shiftx = c(0,-.14,.14,0),
              self.arrpos = c(-1,0,3,.5))
pm$arr$TextX[2] <- pm$arr$TextX[2]-.01
pm$arr$TextY[4] <- pm$arr$TextY[4]+.105
pm$arr$TextY[7] <- pm$arr$TextY[7]-.105
pm$arr$TextX[5] <- pm$arr$TextX[5]-.05
pm$arr$TextX[6] <- pm$arr$TextX[6]+.05

text(pm$arr$TextX, pm$arr$TextY, 
     TeX(paste(qmat_let, pm$arr$Value, sep = "=")), cex = .7)

dev.off()

qratio.msm(msm_glb, ind1 = c(2,1), ind2 = c(1,2))
qratio.msm(msm_glb, ind1 = c(2,3), ind2 = c(3,2))
# M-T vs T-M
qratio.msm(msm_glb, ind1 = c(2,4), ind2 = c(4,2))
# P-B vs B-P
qratio.msm(msm_glb, ind1 = c(3,1), ind2 = c(1,3))

###
qratio.msm(msm_glb, ind1 = c(2,1), ind2 = c(2,4)) # M-B << M-T
###
qratio.msm(msm_glb, ind1 = c(4,2), ind2 = c(2,1)) # T-M >> B-M

qratio.msm(msm_glb, ind1 = c(1,3), ind2 = c(3,1))
qratio.msm(msm_glb, ind1 = c(1,2), ind2 = c(2,1))
qratio.msm(msm_glb, ind1 = c(2,1), ind2 = c(1,2))



quartz()



pmat <- pmatrix.msm(msm_glb)
colSums(pmat)
diag(pmat) = 0
quartz()
transitionPlot(tr_intens, box_prop = cbind(rowSums(trans_nb)/sum(trans_nb), colSums(trans_nb)/sum(trans_nb)), fill_start_box = c("grey" , "black"))

(trans_perc <- trans_nb/rowSums(trans_nb))

trans_perc[1,4] <- 0
trans_perc[4,1] <- 0

quartz()
transitionPlot(trans_perc)

i = rowSums(trans_nb)/sum(trans_nb)



tr_intens  <- t(coef(msm0)$baseline)
diag(tr_intens) <- 0




