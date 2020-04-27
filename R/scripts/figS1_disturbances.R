### FIGURE S1. DISTURBANCE FREQUENCY #####

### PACKAGES & FUNCTIONS ####

source("R/functions/packages.R")

source("R/functions/my_waffle.R")

### DATA ####

source('R/functions/prep_data.R')

env_all <- readRDS("data/env_all.RDS") %>%
  filter(ID_PE %in% states_ba$ID_PE) %>%
  mutate_at(vars(logging:natural), factor)

env_all <- env_all %>% group_by(ID_PE) %>%
  arrange(year_measured) %>% 
  mutate_at(vars(ORIGINE, PERTURB), lead)

env_all1 <- env_all %>% 
  group_by(ID_PE) %>%
  arrange(year_measured) %>%
  slice(1)

reg <- aggregate(env_all1$ID_PE_MES, by = list(env_all1$ecoreg3), length)

### Frequency of disturbances ####

disturb_summ <- env_all %>% 
  group_by(ID_PE, ecoreg3) %>% 
  summarise("Minor" = all(natural==0 & logging==0),
            "Moderate natural" = any(natural==1),
            "Major natural" = any(natural==2),
            "Moderate logging" = any(logging==1),
            "Major logging" = any(logging==2)) 

disturb_summ <- aggregate(disturb_summ[,-c(1:2)], by=list(disturb_summ$ecoreg3), sum)



### Frequency of the 21 original disturbance types ####

log2 <- table(env_all$ecoreg3, env_all$ORIGINE, env_all$logging)[,,3][,c("CBA", "CBT", "CPR", "CT")]
nat2 <- table(env_all$ecoreg3, env_all$ORIGINE, env_all$natural)[,,3][,c("BR", "ES", "CHT","DT")]

log1 <- table(env_all$ecoreg3, env_all$PERTURB, env_all$logging)[,,2][,c("CAM","CB","CD","CDL","CE","CJ","CP","EPC")]
nat1 <- table(env_all$ecoreg3, env_all$PERTURB, env_all$natural)[,,2][,c("BRP", "EL", "CHP", "VEP", "DP")]

disturb_types <- list("No or minor" = disturb_summ[,2],
                      "Moderate natural" = nat1,
                      "Major natural" = nat2,
                      "Moderate logging" = log1,
                      "Major logging" = log2)


disturb_types <- lapply(disturb_types, sort_col)
lapply(disturb_types,sum)

disturb_code <- list(c("Outbreak", "Windfall", "Decline", "Burn", "Ice storm"),
                     c("Burn", "Outbreak", "Windfall", "Decline"),
                     c("Partial cut", "Commercial thinning", "Strip cut", 
                       "Selection cut", "Diameter-limit cut", 
                       "Partial cut + outbreak", "Improvement cut", "Checkerboard cut"),
                     c("Clearcut", "Cut with protection\nof regeneration",
                       "Final strip cut", "Strip cut"))


disturb_col <- list("No or minor" = "#c7cbce",
                    "Moderate natural" = c("#3c784f", "#89a8ad", "#68566e", "#a83f02", "#3c7dba"),
                    "Major natural" = c("#3c784f","#a83f02", "#89a8ad", "#68566e"),
                    "Moderate logging" = c("#4d8585", gpuPalette('insileco', 13)[3:9]),
                    "Major logging" = c("#4d8585", gpuPalette('insileco', 13)[c(8,4,4)]))


### FIGURE S1. WAFFLE PLOTS ####

png("res/figS1_waffle.png", width = 5.5, height = 6.5, units = "in", res = 600)

#quartz(width = 5.5, height = 6.5)
par(oma=c(0,0,0,0))
layout(matrix(c(1:25,0,26:28,0), 6, byrow = TRUE), 
       heights = c(1, rep(.72, 4), .155), 
       widths = c(.8, rep(.9,3), 1))

for(i in 1:length(disturb_types)) {
  d <- disturb_types[[i]]
  if(i!=1) colnames(d) <- disturb_code[[i-1]]
  nrows <- ifelse(i==1, 60, 42)
 
  par(mar = c(0,0,.2,.1))
  plot0()
  mtext(names(disturb_types[i]), 3, line = -2, 
        cex = .65, font = 2, xpd = NA)
  mtext(paste0("n = ", sum(d)), 3, line = -3, cex = .55)
  
  for(r in 1:3) {
    if(is.null(dim(d))) { 
      dr <- d[r] 
    }else { 
      dr <- d[r,]
      }
    my_waffle(dr, nrows = nrows, ncols = 40, cols = disturb_col[[i]],
              main = "", lgd = FALSE)

  }
 
  plot0()
  if(i>1)
    legend(-1, 0, xjust = 0, yjust = .5, cex = .8,
           legend = rev(names(dr)), 
           fill = rev(disturb_col[[i]]), border = NA,
           bty = "n", x.intersp = .7, y.intersp = 1,
           inset = c(0, -1), xpd = NA)
  
}
par(mar = c(0,0,0,0))
plot0()
mtext("Boreal zone", 3, line = -.9, cex = 0.75, font = 2)
mtext(paste(reg[1,2], "forest plots"), 3, line = - 1.9, cex = .6)
plot0()
mtext("Ecotonal zone", 3, line = -.9, cex = 0.75, font = 2)
mtext(paste(reg[2,2], "forest plots"), 3, line = - 1.9, cex = .6)
plot0()
mtext("Temperate zone", 3, line = -.9, cex = 0.75, font = 2)
mtext(paste(reg[3,2], "forest plots"), 3, line = - 1.9, cex = .6)

dev.off()
