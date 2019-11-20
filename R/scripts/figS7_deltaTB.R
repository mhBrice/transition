### PLOT CHANGE IN BA OF TEMPERATE AND BOREAL STATES VS TRANSITION TYPES ####

### PACKAGES ####
library(graphicsutils)
library(dplyr)
library(msm)
library(sf)
library(data.table)


### DATA ####

source('R/functions/prep_data.R')



states_trans <- states_ba %>% 
  select(plot_id, year1 = year_measured, from = states_ba, 
         Temperate1 = Temperate, Boreal1 = Boreal, Pioneer1 = Pioneer) %>% 
  data.table()

cols = c("year1","from", "Temperate1", "Boreal1", "Pioneer1")
anscols = c("year2", "to", "Temperate2", "Boreal2", "Pioneer2")
states_trans[, (anscols) := shift(.SD, 1, NA, "lead"), .SDcols=cols, by=plot_id]

states_trans <- states_trans %>% 
  mutate(t = year2 - year1, trans = paste0(from,"-",to),
         deltaT = Temperate2 - Temperate1,
         deltaB = Boreal2 - Boreal1,
         deltaP = Pioneer2 - Pioneer1) %>% 
  filter(!is.na(to)) %>% select(plot_id, t, from, to, trans, deltaT:deltaP)

states_trans=states_trans[-which(states_trans$trans=="Boreal-Temperate"),]
states_trans=states_trans[-which(states_trans$trans=="Temperate-Boreal"),]
bp_trans <- rbind(cbind.data.frame(type = "T", trans = states_trans$trans, 
                  delta = states_trans$deltaT),
                  cbind.data.frame(type = "B", trans = states_trans$trans, 
                                   delta = states_trans$deltaB))


bp <- barplot(matrix(1:28, 2), beside=T, space = c(0,1), 
              xaxs = "i", plot = T)

pdf("res/figSupp_deltaTB.pdf", width = 7, height = 7)
par(mar=c(4.3,10,2,.1))
plot0(xlim =  c(-50,50), ylim = rev(range(bp)+c(.8,-.8)), frame.plot=TRUE)
abline(v = 0, lty = 2, col = alpha("black", .8), lwd = 1.2)
abline(h = mean(as.vector(bp)[6:7]), col = "grey", lwd = 1.2)
abline(h = mean(as.vector(bp)[14:15]), col = "grey", lwd = 1.2)
abline(h = mean(as.vector(bp)[22:23]), col = "grey", lwd = 1.2)
boxplot(delta ~ type+trans, data = bp_trans, outline = FALSE, add = TRUE,
        las=1, names = rep("", 28), ann = FALSE, axes = FALSE, horizontal = TRUE, 
        col = rep(st_col[c(4,1)],times=14), at = as.vector(bp))

axis(1, cex.axis = .9, tcl = -.4)
axis(3, cex.axis = .9, tcl = -.4)
mtext("Change in state basal area", 1, line = 2.1, font = 2)
legend(-25, 46.5, legend = c("Boreal", "Temperate"), fill = st_col[c(1,4)], 
       xpd = NA, bty = "n", horiz = TRUE, cex = .85)
par(las = 1, cex = .8)
mtext("From Boreal", 2, at = bp[1], line = 4, 
      font = 2, padj = 0, cex = .9)
mtext("From Mixed", 2, at = bp[7], line = 4, 
      font = 2, padj = 0, cex = .9)
mtext("From Pioneer", 2, at = bp[15], line = 4, 
      font = 2, padj = 0, cex = .9)
mtext("From Temperate", 2, at = bp[23], line = 4, 
      font = 2, padj = 0, cex = .9)
axis(2, at = apply(bp, 2, mean), paste('To', sub(".*-", "",levels(bp_trans$trans))))

xx <- as.factor(paste0(bp_trans$trans,bp_trans$type))
mylevels <- levels(xx)
levelProportions <- table(xx)/nrow(bp_trans)
for(i in 1:length(mylevels)){
  
  thislevel <- mylevels[i]
  thisvalues <- bp_trans[xx==thislevel, "delta"]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  vec = bp[2:1,]
  myjitter <- jitter(rep(vec[i], length(thisvalues)), amount=levelProportions[i]/2)
  points(thisvalues, myjitter, pch=20, col=rgb(0,0,0,.1), cex = .5) 
  
}
dev.off()

