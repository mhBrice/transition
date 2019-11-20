### TABLE SUPP. HR + CI ####

### PACKAGES ####
library(graphicsutils)
library(dplyr)
library(msm)

library(knitr)
library(kableExtra)

### DATA ####

source('R/functions/prep_data.R')

# Load msm results

load("res/msm_all75.rda")

msm_glb <- msm_all75[["msm_glb"]]

### FORMAT DATA FRAME ####

print_msm <- print(msm_glb)
coef_names <- rownames(print_msm)
est <- lapply(print_msm, function(x) round(x, 3))
est[-c(1:3)] <- lapply(est[-c(1:3)], function(x) round(x, 2))


signif <- function(coef){
  x <- ifelse(coef[,2] <= 1 & coef[,3] <= 1 | coef[,2] >= 1 & coef[,3] >= 1, TRUE, FALSE)
  x[is.na(x)] <- FALSE
  x
} 

est_signif <- lapply(est[-c(1:4)], signif)
est_signif

msm_res <- cbind.data.frame(Transitions = coef_names,
                            Baseline = paste0(est$base.Estimate, "\n(", est$base.L, ', ', est$base.U, ")"),
                            Temperature = paste0(est$sTP[,1], "\n(", est$sTP[,2], ', ', est$sTP[,3], ")"),
                            CMI = paste0(est$sCMI[,1], "\n(", est$sCMI[,2], ', ', est$sCMI[,3], ")"),
                            Drainage = paste0(est$DRAIN[,1], "\n(", est$DRAIN[,2], ', ', est$DRAIN[,3], ")"),
                            pH = paste0(est$PH_HUMUS[,1], "\n(", est$PH_HUMUS[,2], ', ', est$PH_HUMUS[,3], ")"),
                            Natural1 = paste0(est$natural1[,1], "\n(", est$natural1[,2], ', ', est$natural1[,3], ")"),
                            Natural2 = paste0(est$natural2[,1], "\n(", est$natural2[,2], ', ', est$natural2[,3], ")"),
                            Logging1 = paste0(est$logging1[,1], "\n(", est$logging1[,2], ', ', est$logging1[,3], ")"),
                            Logging2 = paste0(est$logging2[,1], "\n(", est$logging2[,2], ', ', est$logging2[,3], ")"))

msm_res <- msm_res %>% mutate_if(is.factor, as.character)
msm_res[msm_res == "NA\n(NA, NA)"] <- " "

msm_res[msm_res == "1\n(NA, NA)"] <- "1.000"


### KABLE ####

options(knitr.kable.NA = '')

msm_res %>% 
  mutate_at(-1, function(x) linebreak(x, align = "c")) %>%
  mutate(Temperature = cell_spec(Temperature, "latex", background = ifelse(est_signif$sTP, "#cccccc", "white"), escape = FALSE)) %>%
  mutate(CMI = cell_spec(CMI, "latex", background = ifelse(est_signif$sCMI, "#cccccc", "white"), escape = FALSE)) %>%
  mutate(CMI = cell_spec(CMI, "latex", background = ifelse(est_signif$sCMI, "#cccccc", "white"), escape = FALSE)) %>%
  mutate(Drainage = cell_spec(Drainage, "latex", background = ifelse(est_signif$DRAIN, "#cccccc", "white"), escape = FALSE)) %>%
  mutate(pH = cell_spec(pH, "latex", background = ifelse(est_signif$PH_HUMUS, "#cccccc", "white"), escape = FALSE)) %>%
  mutate(Natural1 = cell_spec(Natural1, "latex", background = ifelse(est_signif$natural1, "#cccccc", "white"), escape = FALSE)) %>%
  mutate(Natural2 = cell_spec(Natural2, "latex", background = ifelse(est_signif$natural2, "#cccccc", "white"), escape = FALSE)) %>%
  mutate(Logging1 = cell_spec(Logging1, "latex", background = ifelse(est_signif$logging1, "#cccccc", "white"), escape = FALSE)) %>%
  mutate(Logging2 = cell_spec(Logging2, "latex", background = ifelse(est_signif$logging2, "#cccccc", "white"), escape = FALSE)) %>%
  kable(format = "latex", booktabs = T, linesep = "", escape = FALSE, align = "lccccccccc") %>%
  kable_styling(font_size = 6.8) %>%
  column_spec(1, bold = T) %>%
  column_spec(2, width = "1.9cm") %>%
  column_spec(3:10, width = "1.8cm") %>%
  row_spec(0, bold = TRUE)

