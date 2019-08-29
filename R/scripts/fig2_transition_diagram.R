### TRANSITION DIAGRAM ####

### PACKAGES ####

library(msm)
library(graphicsutils)
library(dplyr)
library(diagram)

### DATA ####

source('R/scripts/3_prep_trans.R')


(trans_nb <- statetable.msm(states_num, plot_id, data = states_ba))

rownames(trans_nb) <- colnames(trans_nb) <- states

kable(addmargins(trans_nb), format="markdown")



(trans_perc <- t(round(trans_nb/rowSums(trans_nb)*100, 1)))

trans_perc[1,4] <- 0
trans_perc[4,1] <- 0


qmat <- c("-sum(q[Bs], s!=B)","q[BM]", "q[BP]", 0,
          "q[MB]", "-sum(q[Ms], s!=M)", "q[MP]", "q[MT]",
          "q[PB]", "q[PM]", "-sum(q[Ps], s!=P)", "q[PT]",
          0, "q[TM]", "q[TP]", "-sum(q[Ts], s!=T)")
.expressions <- qmat
qmat_expres <- parse(text = .expressions)

pos.box <- cbind (c(0.5, 0.2, 0.8, 0.5), 
                  c(0.85, 0.5, 0.5, 0.15))

mat <- matrix(c(1,1,2,3),2, byrow = T)
 
pdf("res/fig2_trans_diagram.pdf", width=6.5, height=5)
#quartz(width=6.5, height=5)
layout(mat, heights = c(1,.6), widths = c(.75,1))
par(mar=c(0,3.7,0,3.7))

pm <- plotmat(trans_perc, pos = pos.box, curve = 0.07, name = states, 
              lwd = 1.2, relsize = .9,
              box.cex = 1, cex.txt = 0, txt.col = "white", dtext = .45, txt.font = 2,
              box.lwd = 0.1, box.type = "rect", shadow.size = 0.005,
              box.prop = 0.4, box.size = 0.1, box.col = st_col,
              arr.length=.12, arr.width=.12,  arr.type ="triangle",
              arr.col ="grey40", arr.lcol = "grey40", 
              self.cex = 0.6, self.lwd = 1.2, 
              self.shifty = c(.07,0,0,-.07), self.shiftx = c(0,-.14,.14,0),
              self.arrpos = c(-1,0,3,.5))
pm$arr$TextX[5] <- pm$arr$TextX[5]-.01
pm$arr$TextX[10] <- pm$arr$TextX[10]+.01
text(pm$arr$TextX, pm$arr$TextY, paste0(pm$arr$Value, "%"), cex = .9)

par(mar=c(0,0,0,0))
plot0(xlim = c(0,9), ylim = c(1,9), xpd = NA, yaxs = "i")
text(.5, 5, expression(paste(bold("Q")," =")))
lines(c(1,1), c(1.5,8.5))
lines(c(9,9), c(1.5,8.5))
text(rep(c(2,4,6,8),4),rep(c(8,6,4,2),ea=4), qmat_expres, cex = .9)

par(mar=c(0,0,0,.1))
plot0(0:5)
text(2.5,4, expression(atop(q[rs]==q[rs.0]%*%exp(beta[rs.1]%*%climate+beta[rs.2]%*%soil+beta[rs.3]%*%disturbances), "for"~r!=s~and~s!=Pioneer)), cex = .9, adj = 0.5)

text(2.5,2, expression(atop(q[rs]==q[rs.0]%*%exp(beta[rs.3]%*%disturbances), "for"~s==Pioneer)), cex = .9, adj = .5)


dev.off()




