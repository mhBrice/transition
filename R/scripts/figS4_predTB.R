
### PACKAGES & FUNCTIONS ####

source("R/functions/packages.R")

### DATA ####

source('R/functions/prep_data.R')

pal <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))

states_ba1 <- states_ba %>% group_by(ID_PE) %>% slice(n())

## GLM with quasibinomial family for proportion ####

states_ba1$prop = states_ba1$Temperate/states_ba1$tot_TB
glmBT <- glm(prop ~ poly(sTP, 3) + poly(sCMI, 3) + sTP:sCMI, data = states_ba1, 
             family = quasibinomial)


summary(glmBT)

(pseudoR2 <- 1 - (glmBT$deviance/glmBT$null))

## Model prediction ###

tp <- seq(quantile(states_ba1$sTP, .1), quantile(states_ba1$sTP, .9), len = 100)
cmi <- seq(quantile(states_ba1$sCMI, .1), quantile(states_ba1$sCMI, .9), len = 100)

predBT <- predict(glmBT, newdata = data.frame(expand.grid(sTP = tp, sCMI = cmi)), 
                 type = "response")


tp <- tp * sc_sTP[2] + sc_sTP[1]
cmi <- cmi * sc_sCMI[2] + sc_sCMI[1]

pal_br <- seq(0, 1, .01)


### Figure S ####
pdf("res/figS_predBT.pdf", width = 4, height = 4)

#quartz(width = 4, height = 4)
par(mar = c(2.8,3.5,2.8,.5))
image(x = cmi, y = tp, 
      z = t(matrix(predBT, ncol = length(cmi), nrow = length(tp))),
      xlab = "", ylab = "", col = pal(100), breaks = pal_br, 
      las = 1, axes = F)
box2()
axis(1, labels = FALSE, tcl = -.4)
axis(1, tick = FALSE, line = -.2, cex.axis = 0.8)
axis(2, labels = FALSE, tcl = -.4)
axis(2, tick = FALSE, line = -.2, cex.axis = 0.8, las = 1)
# points(x= states_ba1$sCMI* sc_sCMI[2] + sc_sCMI[1], y = states_ba1$sTP*sc_sTP[2] + sc_sTP[1],
#        cex = .1, pch = 20, col = alpha("grey15",.3), xpd = FALSE)
contour(x = cmi, y = tp, 
        z = t(matrix(predBT, ncol = length(cmi), nrow = length(tp))),
        cex = 1, add = TRUE)

mtext("Predicted proportion of temperate trees", side = 3, font = 2, line = 1.5)
mtext("Temperate/(Boreal + Temperate)", side = 3, line = .5)

mtext("Temperature", side = 2, cex = 1, font = 2, line = 2.5)
mtext("CMI", side = 1, cex = 1, font = 2, line = 1.8)
dev.off()






xy <- st_read("data/plot_xy32198_nov2019.gpkg") %>%
  filter(ID_PE %in% states_ba$ID_PE)

map = eigenmap(x=st_coordinates(xy), boundaries = c(0,1000))
plot(map)
mca1 = MCA(Y = as.matrix(states_ba1[,c("Temperate", "Boreal")]), 
           X = as.matrix(states_ba1[,c("sTP", "sCMI")]), 
           emobj = map)

mca1_partest <- test.cdp(mca1)
mca1_partest
summary(mca1_partest)
par(mar = c(6,4,2,4))
plot(mca1_partest, las = 3)
mca1_pertest <- permute.cdp(mca1, permute=999)




