### PLOT STEADY-STEADY PROPORTION ####

# Figure 6

### PACKAGES & FUNCTIONS ####

source("R/functions/packages.R")

source('R/functions/markov.R')

### DATA ####

source('R/functions/prep_data.R')

# Load msm results

load("res/msm_all75.rda")

msm_glb <- msm_all75[["msm_glb"]]

### DATA PREP ####

# Mean environmental conditions at the ecotone

ll_ecotone <- which(states_ba$ecoreg3=="Mixed")

mu_ecotone <- as.list(apply(states_ba[ll_ecotone, c("sTP", "sCMI", "DRAIN", "PH_HUMUS")],
                            2, mean))

# Temperature boundaries of the ecotone

tp_ecotone <- quantile(states_ba$sTP[ll_ecotone], c(.2, .8))

# Temperature gradient

tp_grad <- seq(-1.9, 1.6, len = 50)

# Expand data frame with temperature gradient and disturbances levels

dfl <- expand.grid(sTP = tp_grad, logging = c(0, 1, 2), 
                  sCMI = mu_ecotone[[2]],
                  DRAIN = mu_ecotone[[3]],
                  PH_HUMUS = mu_ecotone[[4]])

dfn <- expand.grid(sTP = tp_grad, natural = c(0, 1, 2), 
                   sCMI = mu_ecotone[[2]],
                   DRAIN = mu_ecotone[[3]],
                   PH_HUMUS = mu_ecotone[[4]])

## BARPLOT

# Initial observed state proportion

init <- states_ba %>% filter(ecoreg3 == "Mixed") %>%
  group_by(ID_PE) %>% 
  arrange(year_measured) %>% 
  slice(1) %>% filter(year_measured<15)
init <- table(init$states_ba)
init <- init/sum(init)

dld <- list(c(mu_ecotone), 
           c(mu_ecotone, natural1 = 1), 
           c(mu_ecotone, logging1 = 1),
           c(mu_ecotone, natural2 = 1), 
           c(mu_ecotone, logging2 = 1))

# Q matrix

qmats <- lapply(dld, function(x)
  qmatrix.msm(msm_glb, covariates = x, ci = "none"))

# Steady states

SS <- lapply(qmats, function(x) steady_state(qmat = x))
SS <- rbind(init, do.call(rbind, SS))



### FIGURE 6. STEADY STATE ALONG TP GRADIENT ####


mat <- matrix(c(4,4,4,5,1:3,5), 4)

pdf("res/fig5_steady.pdf", width = 3.5, height = 7)
#quartz(width = 3.5, height = 7.2)

layout(mat, widths = c(0.07, 1), heights = c(1,1,1,0.2))
par(mar = c(2,2,1.4,0.2), oma = c(0,0,0.1,0))

barplot_index(index = SS, ylim = c(0,1))
mtext(paste0("(", letters[1],")"), 3, adj = 0, line = .2, font = 2, cex = .8)

# Steady state for natural disturbances

x=plot_ss(mod = msm_glb, df = dfn, tp_ecotone = tp_ecotone, 
               dist = "natural",
               unscale = sc_sTP, axes = 2, 
               xlab = NULL, ylab = NULL, main = "Natural")
mtext(paste0("(", letters[2],")"), 3, adj = 0, line = 0.2, font = 2, cex = .8)

# Steady state for logging

x=plot_ss(mod = msm_glb, df = dfl, tp_ecotone = tp_ecotone, 
               dist = "logging",
               unscale = sc_sTP,
               ylab = NULL, xlab = "Mean temperature of the growing season (°C)",
               main = "Logging")
mtext(paste0("(", letters[3],")"), 3, adj = 0, line = 0.2, font = 2, cex = .8)


# Axis

par(mar = c(1,0,1.5,0))
plot0()
mtext("State proportion at equilibrium", 2, line = -1.5, cex =.85)

# Legend

par(mar = c(0,2.2,0,0.2))
plot0()
legend("bottom", legend = c("No or minor", "Moderate", "Major"), cex = 1.1,
       col = c("grey75","grey45","grey15"), lty = 1:3, lwd = 1.5, 
       seg.len = 2, x.intersp = 0.9,
       horiz = TRUE, bty = "n")

dev.off()




