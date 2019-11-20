### FIGURE supp. 2 DISTURBANCE FREQUENCY #####

### PACKAGES ####
library(dplyr)
library(scales)
library(sf)
library(graphicsutils)
library(RColorBrewer)

### DATA ####
source('R/functions/prep_data.R')

env_all <- readRDS("data/env_all.RDS") %>%
  filter(ID_PE %in% states_ba$ID_PE) %>%
  mutate_at(vars(logging:natural), factor)

env_all <- env_all %>% group_by(plot_id) %>%
  arrange(year_measured) %>% 
  mutate_at(vars(ORIGINE, PERTURB), lead)

### Frequency of disturbances ####

disturb_summ <- env_all %>% 
  group_by(ID_PE) %>% 
  summarise("Minor" = all(natural==0 & logging==0),
            "Moderate natural" = any(natural==1),
            "Major natural" = any(natural==2),
            "Moderate logging" = any(logging==1),
            "Major logging" = any(logging==2)) 

disturb_summ <- apply(disturb_summ[,-1], 2, sum)


### BARPLOT ####

pdf("res/figS2_disturb_bars.pdf", width = 4, height = 4)
#quartz(width = 4, height = 4)
par(mar = c(3,4,.5,.5))
bp <- barplot(disturb_summ, las = 1, border = NA, axisnames = FALSE, width = .8,
              cex.axis = .75, ylab = "Number of forest plots", cex.lab = .9)
text(bp, -200, c("Minor", rep(c("Moderate", "Major"), 2)), xpd = NA, cex = .85,
     adj = c(.5, 1))
text(mean(bp[2:3]), -700, "Natural", xpd = NA, cex = .9, font = 2)
text(mean(bp[4:5]), -700, "Logging", xpd = NA, cex = .9, font = 2)
lines(x=bp[2:3]+c(-.2,.2), y = c(-450,-450), lwd = 2, xpd = NA)
lines(x=bp[4:5]+c(-.2,.2), y = c(-450,-450), lwd = 2, xpd = NA)
dev.off()


### WAFFLE PLOTS ####

log2 <- table(env_all$ORIGINE, env_all$logging)[,3][c("CBA", "CBT", "CPR", "CT")]
nat2 <- table(env_all$ORIGINE, env_all$natural)[,3][c("BR", "ES", "CHT","DT")]

log1 <- table(env_all$PERTURB, env_all$logging)[,2][c("CAM","CB","CD","CDL","CE","CJ","CP","EPC")]
nat1 <- table(env_all$PERTURB, env_all$natural)[,2][c("BRP", "EL", "CHP", "VEP", "DP")]

disturb_types <- list("No or minor" = disturb_summ[1],
                      "Moderate natural" = nat1,
                      "Major natural" = nat2,
                      "Moderate logging" = log1,
                      "Major logging" = log2)

disturb_types <- lapply(disturb_types, function(x) sort(x, decreasing = TRUE))
lapply(disturb_types,sum)

disturb_code <- list(c("Outbreak", "Windfall", "Decline", "Burn", "Ice storm"),
                  c("Burn", "Outbreak", "Windfall", "Decline"),
                  c("Partial cut", "Commercial thinning", "Strip cut", 
                    "Selection cut", "Diameter-limit cut", 
                    "Partial cut + outbreak", "Improvement cut", "Checkerboard cut"),
                  c("Clearcut", "Cut with protection\nof regeneration",
                  "Final strip cut", "Strip cut"))



my_waffle <- function(x, nrows, ncols = 65, 
                      pal = "Dark2", main = NULL, lgd = TRUE, ...) {
  nd <- length(x)
  if(nd != 1) {
    cols <- brewer.pal(length(x), name = pal)
  } else {
    cols <- brewer.pal(3, name = pal)[1]
  }
  xx <- rep(cols, times = x)
  lx <- length(xx)

  m <- matrix(nrow = nrows, ncol = ncols)#ncol = (lx %/% rows) + (lx %% rows != 0))
  m[1:length(xx)] <- xx
  
  o <- cbind(c(row(m)), c(col(m))) + 1
  plot0(xlim = c(0, max(o[, 2]) + 1), ylim = c(0, max(o[, 1]) + 1),
        asp = 1, xaxs = 'i', yaxs = 'i')
  
  mtext(main, 3, line = .85, cex = .8, font = 2)
  mtext(paste("n =", sum(x)), 3, line = 0, cex = .65)
  usr <- par("usr")
  rect(o[, 2], o[, 1], o[, 2] + .85, o[, 1] + .85, col = c(m), border = NA)
  
  if(lgd) {
  lgd_pos <- which(is.na(m), arr.ind = T)[1,2]
  legend(lgd_pos+1, mean(usr[3:4]), yjust = 0.5, cex = .8,
         legend = paste0(names(x), " ", round(x/sum(x)*100, 2), "%"), 
         fill = cols, border = NA,
         bty = "n", x.intersp = .7, y.intersp = 1,
         inset = c(0, -1), xpd = NA)
  }
  invisible(list(m = m, o = o))
}



lgd = c(FALSE, TRUE, TRUE, TRUE, TRUE)

png("res/figS2_waffle.png", width = 7.4, height = 7.3, units = "in", res = 300)
#quartz(width = 7.4, height = 7.3)
layout(matrix(c(1,1,2:5), 3, byrow = TRUE))
par(oma=c(0,0,0,1))
for(i in 1:length(disturb_types)) {
  d <- disturb_types[[i]]
  if(i!=1) names(d) <- disturb_code[[i-1]]
  ncols <- ifelse(i==1, 100, 65)

  par(mar = c(0,0,1.8,0))
  if(i==1) par(mar = c(0,0,1.8,22)) 
  if(i %in% c(2,4)) par(mar = c(0,0,1.8,5.8))
  
  my_waffle(d, nrows = 50, ncols = ncols,
            main = names(disturb_types[i]), lgd = lgd[i])

}

dev.off()
