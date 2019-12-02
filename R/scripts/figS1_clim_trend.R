### FIGURE S1. Climate trend ####


### PACKAGES & FUNCTIONS ####

source("R/functions/packages.R")

### DATA ####

source("R/functions/prep_data.R")

bioclim_all <- readRDS("data/bioclim_ally_nov2019.RDS")
bioclim_all$year <- as.numeric(bioclim_all$year)
bioclim_all <- bioclim_all %>% 
  subset(ID_PE %in% states_ba$ID_PE)

bioclim_all <- bioclim_all %>% ungroup() %>%
  mutate(sCMI = rowSums(select(., "cmi_05":"cmi_09"))) 

mean_TP <- aggregate(bioclim_all$sg_15, by = list(bioclim_all$year), mean)
mean_CMI <- aggregate(bioclim_all$sCMI, by = list(bioclim_all$year), function(x) mean(x, na.rm=T))


lm_cmi <- lm(x ~ Group.1, data = mean_CMI)
summary(lm_cmi)
lm_TP <- lm(x ~ Group.1, data = mean_TP)
summary(lm_TP)


pdf("res/figS1_clim_trend.pdf",
    width = 3.5, height = 5.5)
# quartz(width = 3.5, height = 5.5)
par(mfrow = c(2,1), mar = c(1,4,.5,.5), oma = c(2.5,0,0,0))

plot(x ~ Group.1, data = mean_TP, type = "l", las = 1,
     xlab = "", ylab = "", cex.axis=.8, ylim = c(10.9,14.1), xaxs = "i", yaxs = "i",
     col = "grey45", axes = F, frame.plot = TRUE)
axis(2, cex.axis=.75, las = 1)
axis(1, labels = FALSE)
abline(lm_TP, lwd = 1.2)
mtext(paste("Slope =", round(lm_TP$coef[2],3), "°C/year"),
      3, line = -1, at = 1960, adj = 0,
      cex = .75)
mtext("p-value < 0.001",
      3, line = -2, at = 1960, adj = 0,
      cex = .75)

mtext("Mean temperature\nduring the growing season", 2, line = 2.3, cex = 0.8, font = 2)


plot(x ~ Group.1, data = mean_CMI, type = "l", las = 1,
     xlab = "", ylab = "", cex.axis=.75, xaxs = "i", yaxs = "i",
     col = "grey45", ylim = c(0,30))
abline(lm_cmi, lwd = 1.2, lty = 2)
mtext(paste("Slope =", round(lm_cmi$coef[2],3), "cm/year"),
      3, line = -1, at = 1960, adj = 0,
      cex = .75)
mtext("p-value = 0.83",
      3, line = -2, at = 1960, adj = 0,
      cex = .75)
mtext("Year", 1, line = 2.5, cex = 0.8, font = 2)
mtext("Climate Moisture Index\nfrom May to September", 2, line = 2.3, cex = 0.8, font = 2)

dev.off()


