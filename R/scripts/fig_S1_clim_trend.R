### FIGURE S1. Climate trend ####

# Data
source("R/scripts/3_prep_trans.R")

bioclim_all <- readRDS("data/bioclim_corrected.RDS")
bioclim_all$year <- as.numeric(bioclim_all$year)
bioclim_all <- bioclim_all %>% subset(ID_PE %in% states_ba$ID_PE & year >= 1960) %>%
  left_join(states_ba)


mean_TP <- aggregate(bioclim_all$sg_15, by = list(bioclim_all$year), mean)
mean_CMI <- aggregate(bioclim_all$cmi, by = list(bioclim_all$year), function(x) mean(x, na.rm=T))


lm_cmi <- lm(x ~ Group.1, data = mean_CMI)
summary(lm_cmi)
lm_TP <- lm(x ~ Group.1, data = mean_TP)
summary(lm_TP)


pdf("res/figS1_clim_trend.pdf",
    width = 3.5, height = 5.5)
# quartz(width = 3.5, height = 5.5)
par(mfrow = c(2,1), mar = c(1,4,.5,.5), oma = c(2.5,0,0,0))

plot(x ~ Group.1, data = mean_TP, type = "l", las = 1,
     xlab = "", ylab = "", cex.axis=.8,
     col = "grey45", axes = F, frame.plot = TRUE)
axis(2, cex.axis=.8, las = 1)
axis(1, labels = FALSE)
abline(lm_TP, lwd = 1.2)
mtext(paste("Slope =", round(lm_TP$coef[2],3), "°C/year"),
      3, line = -1, at = 1963, adj = 0,
      cex = .75)
mtext("p-value < 0.001",
      3, line = -2, at = 1963, adj = 0,
      cex = .75)

mtext("Growing season temperature", 2, line = 2.7, cex = 0.8, font = 2)


plot(x ~ Group.1, data = mean_CMI, type = "l", las = 1,
     xlab = "", ylab = "", cex.axis=.8,
     col = "grey45", ylim = c(45,77))
abline(lm_cmi, lwd = 1.2, lty = 2)
mtext(paste("Slope =", round(lm_cmi$coef[2],3), "cm/year"),
      3, line = -1, at = 1977, adj = 0,
      cex = .75)
mtext("p-value = 0.54",
      3, line = -2, at = 1977, adj = 0,
      cex = .75)
mtext("Year", 1, line = 2.5, cex = 0.8, font = 2)
mtext("Annual Climate Moisture Index", 2, line = 2.7, cex = 0.8, font = 2)

dev.off()


plot(x~Group.1, data = mean_CMI, type = "l")
lines(x, y, col = "red")
