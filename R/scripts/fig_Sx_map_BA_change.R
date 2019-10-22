### FIGURE SUPP - MAP OF TEMPERATE CHANGE #####

### PACKAGES ####

library(graphicsutils)
library(dplyr)
library(sf)

source('R/functions/plot_map.R')

### DATA ####

source('R/functions/prep_data.R')

# Compute change in basal area for Boreal and Temperate species ####

delta_T <- states_ba %>%
  group_by(plot_id) %>%
  filter(sum(Temperate)>0) %>% 
  mutate(delta_T = (lead(Temperate) - Temperate)/(lead(year_measured) - year_measured)) %>%
  summarise(delta_T = mean(delta_T, na.rm = TRUE)) %>% 
  select(plot_id, delta_T) %>% 
  left_join(xy, by = "plot_id")

delta_B <- states_ba %>%
  group_by(plot_id) %>%
  filter(sum(Boreal)>0) %>% 
  mutate(delta_B = (lead(Boreal) - Boreal)/(lead(year_measured) - year_measured)) %>%
  summarise(delta_B = mean(delta_B, na.rm = TRUE)) %>% 
  select(plot_id, delta_B) %>% 
  left_join(xy, by = "plot_id")



# Graphical parameters

n_T <- nrow(delta_T)
inc_T <- round(length(which(delta_T$delta_T>0))/n_T*100, 2)
dec_T <- round(length(which(delta_T$delta_T<=0))/n_T*100, 2)

n_B <- nrow(delta_B)
inc_B <- round(length(which(delta_B$delta_B > 0))/n_B*100, 2)
dec_B <- round(length(which(delta_B$delta_B <= 0))/n_B*100, 2)

# Point size
cex_Ti <- delta_T$delta_T/2
cex_Ti[cex_Ti < 0.1] <- 0.1
cex_Td <- (-delta_T$delta_T)/2
cex_Td[cex_Td < 0.1] <- 0.1

cex_Bi <- delta_B$delta_B/2
cex_Bi[cex_Bi < 0.1] <- 0.1
cex_Bd <- (-delta_B$delta_B)/2
cex_Bd[cex_Bd < 0.1] <- 0.1

# Color of points
bg_Ti <- ifelse(delta_T$delta_T > 0, "#2C7BB6", "transparent")
bg_Td <- ifelse(delta_T$delta_T <= 0, "#D7191C", "transparent")

bg_Bi <- ifelse(delta_B$delta_B > 0, "#2C7BB6", "transparent")
bg_Bd <- ifelse(delta_B$delta_B <= 0, "#D7191C", "transparent")

col_Ti <- ifelse(delta_T$delta_T > 0, "#0B4775", "transparent")
col_Td <- ifelse(delta_T$delta_T <= 0, "#8C0808", "transparent")

col_Bi <- ifelse(delta_B$delta_B > 0, "#0B4775", "transparent")
col_Bd <- ifelse(delta_B$delta_B <= 0, "#8C0808", "transparent")


### Plot map of change in boreal and temperate


pdf("res/figSx_map_BA_change.pdf", width = 7.5, height = 4.7)
#quartz(width = 7.5, height = 4.7)
par(mfcol = c(2,2), mar = c(.5,1,.3,.3), oma = c(1,1,1,0))

### Temperate

### Increase
plot_map(ecoregion, xy_pts = st_as_sf(delta_T), 
         pch = 21, col = col_Ti, bg = bg_Ti, cex = cex_Ti, 
         axes = 2)
text(-62, 45.2, labels = paste0(inc_T,"% of plots increase"), cex=0.85)
mtext("Temperate", 3, font = 2)

### Decrease
plot_map(ecoregion, xy_pts = st_as_sf(delta_T), 
         pch = 21, col = col_Td, bg = bg_Td, cex = cex_Td)
text(-62, 45.2, labels = paste0(dec_T,"% of plots decrease"), cex=0.85)


### Temperate
### Increase
plot_map(ecoregion, xy_pts = st_as_sf(delta_B), 
         pch = 21, col = col_Bi, bg = bg_Bi, cex = cex_Bi, 
         axes = NULL)
text(-62, 45.2, labels = paste0(inc_B,"% of plots increase"), cex=0.85)
mtext("Boreal", 3, font = 2)

### Decrease
plot_map(ecoregion, xy_pts = st_as_sf(delta_B), 
         pch = 21, col = col_Bd, bg = bg_Bd, cex = cex_Bd, 
         axes = 1)
text(-62, 45.2, labels = paste0(dec_B,"% of plots decrease"), cex=0.85)

dev.off()

