### DEMOGRAPHICAL PROCESSES & COMMUNITY TRANSITION ####

### PACKAGES ####

source("R/functions/packages.R")

### DATA ####

source('R/functions/prep_data.R')

tree_data <- readRDS("data/tree_data_nov2019.RDS") %>%
  filter(ID_PE_MES %in% states_ba$ID_PE_MES) %>% 
  arrange(ID_PE_MES)

tree_code <- read.csv2("data/ref_spCode.csv")


#### SPECIES OF INTEREST #####

# most abundant species in BA
MySpecies <- c("ACERUB", "ACESAC", "BETALL", "FAGGRA", "THUOCC", "PINSTR", 
               "ABIBAL", "PICGLA", "PICMAR", "PINBAN",
               "BETPAP", "POPTRE")

MySpecies <- c("ACERUB", "ACESAC", "BETALL", "ABIBAL", "PICMAR", "BETPAP", "POPTRE")


### FUNCTIONS ####


### COMPUTE BASAL AREA PER HECTARE ####

ba_fun <- function(dhp) sum(pi*(dhp/(2 * 1000))^2, na.rm = T)*(10000/399.7312)

bai_fun <- function(dhp) (pi*(dhp/(2 * 1000))^2)*(10000/399.7312)

### SHAPE DEMOGRAPHIC DATA ####

### RECRUITS ####

# Count
tree_data <- tree_data %>% mutate(recru = ifelse(ETAT %in% c(40:48), 1, 0))

recru_mat <- reshape2::dcast(tree_data, 
                             ID_PE + ID_PE_MES + year_measured ~ sp_code, 
                             fun.aggregate = sum,
                             fill = 0,
                             value.var = "recru") %>%
  select(ID_PE, ID_PE_MES, year_measured, MySpecies)

# Basal area
tree_data <- tree_data %>% mutate(recru_ba = recru*DHP)
recru_ba <- reshape2::dcast(tree_data, 
                            ID_PE + ID_PE_MES + year_measured ~ sp_code, 
                            fun.aggregate = ba_fun,
                            fill = 0,
                            value.var = "recru_ba") %>%
  select(ID_PE, ID_PE_MES, year_measured, MySpecies)

### DEAD ####

# Count    
tree_data <- tree_data %>% mutate(dead = ifelse(state=="dead", 1, 0))

dead_mat <- reshape2::dcast(tree_data, 
                            ID_PE + ID_PE_MES + year_measured ~ sp_code, 
                            fun.aggregate = sum, 
                            fill = 0, 
                            value.var = "dead") %>%
  select(ID_PE, ID_PE_MES, year_measured, MySpecies)

# Basal area 
# Some dead trees have DHP mesures others don't

tree_data <- tree_data %>%  mutate(dead_ba = dead*DHP)

tree_data <- tree_data %>% 
  group_by(ID_ARBRE) %>% 
  arrange(year_measured) %>% 
  mutate(dead_ba = ifelse(is.na(dead_ba), dead*lag(DHP), dead_ba)) %>% 
  ungroup() %>% 
  arrange(ID_PE_MES)

dead_ba <- reshape2::dcast(tree_data, 
                           ID_PE + ID_PE_MES + year_measured ~ sp_code, 
                           fun.aggregate = ba_fun,
                           fill = 0,
                           value.var = "dead_ba") %>%
  select(ID_PE, ID_PE_MES, year_measured, MySpecies)


### HARVESTED ####

# Count
tree_data <- tree_data %>% mutate(harv = ifelse(state=="harvested", 1, 0))
harv_mat <- reshape2::dcast(tree_data, 
                            ID_PE + ID_PE_MES + year_measured ~ sp_code, 
                            fun.aggregate = sum,
                            fill = 0, 
                            value.var = "harv") %>%
  select(ID_PE, ID_PE_MES, year_measured, MySpecies)


# Basal area 

tree_data <- tree_data %>% mutate(harv_ba = harv*DHP)

tree_data <- tree_data %>% 
  group_by(ID_ARBRE) %>% 
  arrange(year_measured) %>% 
  mutate(harv_ba = ifelse(is.na(harv_ba), harv*lag(DHP), harv_ba)) %>% 
  ungroup() %>% 
  arrange(ID_PE_MES)

harv_ba <- reshape2::dcast(tree_data, 
                           ID_PE + ID_PE_MES + year_measured ~ sp_code, 
                           fun.aggregate = ba_fun,
                           fill = 0,
                           value.var = "harv_ba") %>%
  select(ID_PE, ID_PE_MES, year_measured, MySpecies)

### GROWTH ####

tree_data <- tree_data %>% mutate(BA = bai_fun(DHP))

tree_data <- tree_data %>% 
  group_by(ID_ARBRE) %>% 
  arrange(year_measured) %>% 
  mutate(G = ifelse(state=="alive", BA - lag(BA), 0)) %>% 
  ungroup() %>% 
  arrange(ID_PE_MES)

growth_ba <- reshape2::dcast(tree_data, 
                             ID_PE + ID_PE_MES + year_measured ~ sp_code, 
                             fun.aggregate = function(x) sum(x, na.rm = TRUE),
                             fill = 0,
                             value.var = "G")  %>%
  select(ID_PE, ID_PE_MES, year_measured, MySpecies)


### Lists ####

colnames(recru_ba) <- c("ID_PE", "ID_PE_MES", "year_measured", paste0(MySpecies, "_R"))
colnames(growth_ba) <- c("ID_PE", "ID_PE_MES", "year_measured", paste0(MySpecies, "_G"))
colnames(dead_ba) <- c("ID_PE", "ID_PE_MES", "year_measured", paste0(MySpecies, "_D"))
colnames(harv_ba) <- c("ID_PE", "ID_PE_MES", "year_measured", paste0(MySpecies, "_H"))

states_demo <- states_ba %>% 
  select(ID_PE, ID_PE_MES, year_measured, states_ba, logging, natural) %>%
  mutate(year_measured = year_measured+1970) %>%
  left_join(recru_ba, by = c("ID_PE", "ID_PE_MES", "year_measured")) %>%
  left_join(growth_ba, by = c("ID_PE", "ID_PE_MES", "year_measured")) %>%
  left_join(dead_ba, by = c("ID_PE", "ID_PE_MES", "year_measured")) %>%
  left_join(harv_ba, by = c("ID_PE", "ID_PE_MES", "year_measured"))
  
states_demo <- states_demo %>% 
  rename(year2 = year_measured, to = states_ba) %>%
  select(ID_PE, ID_PE_MES, year2, to, everything()) %>% data.table()

cols = c("year2","to", "logging", "natural")
anscols = c("year1", "from", "logging", "natural")
states_demo[, (anscols) := shift(.SD, 1, NA, "lag"), .SDcols=cols, by=ID_PE]

states_demo <- states_demo %>% 
  mutate(t = year2 - year1, trans = paste0(from,"-",to)) %>% 
  filter(!is.na(from)) %>% 
  select(ID_PE, ID_PE_MES, year1, year2, t, from, to, trans, everything())


#########
### PLOT #####

sp_name <- tree_code$complete.name[match(MySpecies,tree_code$spCode)]


indic_demo <- indval(states_demo[,-c(1:10)], states_demo$trans)

# Table of the significant indicator species
gr <- indic_demo$maxcls[indic_demo$pval <= 0.05]
iv <- indic_demo$indcls[indic_demo$pval <= 0.05]
pv <- indic_demo$pval[indic_demo$pval <= 0.05]
fr <- apply(states_demo[,11:38] > 0, 2, sum)[indic_demo$pval <= 0.05]
fidg <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
fidg <- fidg[order(fidg$group, -fidg$indval),]
fidg$group <- factor(fidg$group, 1:16, levels(as.factor(states_demo$trans)))


tab=t(indic_demo$indval)

tab[tab>.2] = .2

tab <- tab[c(1,2,4,3, 5,6,8,7, 13,14,16,15, 9,10,12,11), ]

pal0 <- c("#fcfaf0", "#afc24e", "#9eb625", "#045579", "#0b2d40")
pal <- colorRampPalette(pal0)(30)
seqx = seq(0, 1, len = nrow(tab))
seqy = seq(0, 1, len = ncol(tab))



png("res/fig7_demo_trans.png", width = 7, height = 6, units = "in", res = 300)
#quartz(width = 7, height = 6)
par(mar=c(2.5, 8, .2, .2))
image(tab[,28:1], col = (pal), axes = F)

text(seqx, rep(-.04, 16), gsub("[^A-Z]", "", row.names(tab)), 
     xpd = NA, cex = .8, font = 2)
mtext("Transitions", 1, line = 1.3, font = 2)

text(rep(-.04, 28), rev(seqy), gsub("_.", "", sp_name), xpd = NA, 
     cex = .8, adj = 1, font = 3)
text(-.32, mean(seqy[1:7]), "Logging", srt = 90, xpd = NA, font = 2)
text(-.32, mean(seqy[8:14]), "Mortality", srt = 90, xpd = NA, font = 2)
text(-.32, mean(seqy[15:21]), "Growth", srt = 90, xpd = NA, font = 2)
text(-.32, mean(seqy[22:28]), "Recruitment", srt = 90, xpd = NA, font = 2)

arrows(x0 = mean(seqx[4:5]), y0 = -0.05, y1 = 1.02, length = 0, xpd = NA)
arrows(x0 = mean(seqx[8:9]), y0 = -0.05, y1 = 1.02, length = 0, xpd = NA)
arrows(x0 = mean(seqx[12:13]), y0 = -0.05, y1 = 1.02, length = 0, xpd = NA)
# polygon(x = c(c(mean(seqx[7:8]),mean(seqx[8:9])),
#                 rev(c(mean(seqx[7:8]),mean(seqx[8:9])))),
#         y = c(1.02,1.02, -.02,-.02), xpd = NA, border = "red3", lwd = 2)

arrows(y0 = mean(seqy[7:8]), x0 = -0.33, x1 = 1.04, length = 0, xpd = NA)
arrows(y0 = mean(seqy[14:15]), x0 = -0.33, x1 = 1.04, length = 0, xpd = NA)
arrows(y0 = mean(seqy[21:22]), x0 = -0.33, x1 = 1.04, length = 0, xpd = NA)

dev.off()