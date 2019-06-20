#### Model of community transition ####

## ---- packages ----

### PACKAGES ####
require(dplyr)
require(reshape2)
require(data.table)
#require(tibble)
require(tidyr)

#require(labdsv)
require(cluster)
#require(factoextra)
#require(adespatial)
require(vegan)
#require(fmsb)
#require(ROCR)
#require(coefplot2)

require(nnet)
require(MASS)
require(effects)
require(sf)
#require(nngeo)

require(ggplot2)
#require(ggridges)
require(RColorBrewer)
require(scales)
require(diagram)


fct_path <- "~/Documents/GitHub/Doctorat/community_trajectory/functions/"
source(paste0(fct_path, "modelTransition.R"))
#source(paste0(fct_path, "RegTreeTransition.R"))

## ---- data ----

### DATA ####

# Species
path_data <- "~/Documents/GitHub/Doctorat/data/"
sp.mat.abun <- readRDS(paste0(path_data, "species/sp_mat_abun.RDS"))


# Spatial
ecoregion <- st_read(paste0(path_data, "map/ecoregion_simple.gpkg"), quiet = T)
xy_epsg4269 <- st_read(paste0(path_data, "map/xy_epsg4269.gpkg"), quiet = T)

# plot info
ecoreg_df <- readRDS(paste0(path_data,"map/ecoreg_df.RDS"))
plotInfo <- readRDS(paste0(path_data,"map/plotInfo.RDS"))
mortalityTable <- readRDS(paste0(path_data,"map/mortalityTable.RDS"))

# Climate
bioclim <- readRDS(paste0(path_data,"env/bioclim5_mat.RDS"))

# disturbance
plot_disturb <- readRDS(paste0(path_data,"env/plot_disturb.RDS"))


## ---- format_sp ----

### FORMATTING SPECIES ####

### Keep only plots that have been sampled first before 1980 and last after 2000

sp.mat.abun <- sp.mat.abun %>% 
  dplyr::select(-id_pep_mes) %>%
  filter(!(PICABI>0 | PINSYL>0)) %>% # should test it, since plantation is an explanatory variable 
  mutate(TOTAL = rowSums(dplyr::select(., -plot_id, -year_measured))) %>%
  group_by(plot_id) %>% 
  filter(first(year_measured)<=1980 & last(year_measured)>=2000) %>% # before 1980 & after 2000
  filter(first(TOTAL)>0 & last(TOTAL)>0) # remove first and last rows that = 0


# Removing rare species (abundance <110)
sp_keep <- names(which(colSums(sp.mat.abun) > 110))

sp.mat.abun <- sp.mat.abun %>% dplyr::select(sp_keep) 

MySpecies <- sp_keep[!sp_keep %in% c("plot_id", "year_measured", "TOTAL")]

sp.mat.abun <- merge(sp.mat.abun, ecoreg_df, by="plot_id", sort=F)
sp.mat.abun <- merge(sp.mat.abun, plotInfo, by=c("plot_id", "year_measured"), sort=F) %>%
  dplyr::select(-ecoreg11, -ecoreg10)

sp.mat <- sp.mat.abun

### FORMATTING SPATIAL ####
xy <- subset(xy_epsg4269, plot_id %in% sp.mat$plot_id)

st_crs(xy) <- 4269


plotInfo <- subset(plotInfo, plot_id %in% sp.mat$plot_id)

### FORMATTING CLIMATE ####

bioclim <- bioclim %>% 
  right_join(sp.mat[, c("plot_id", "year_measured")]) %>%
  dplyr::select(plot_id, year_measured, bio_01, bio_12, sg_03) %>%
  rename(TP = bio_01, PP = bio_12, GSL = sg_03)


delta_clim <- bioclim %>%
  group_by(plot_id) %>%
  mutate(year1 = year_measured) %>%
  mutate(year2 = lead(year_measured, 1L)) %>%
  mutate(year_measured = NULL) %>%
  mutate(TP = lead(TP, 1L) - TP, PP = lead(PP, 1L) - PP, GSL = lead(GSL, 1L) - GSL) %>%
  filter(!is.na(year2))

### FORMATTING DISTURBANCES ####
plotInfo2 <- sp.mat[, c("plot_id", "year_measured")] %>%
  group_by(plot_id) %>%
  mutate(year1 = year_measured) %>%
  mutate(year2 = lead(year_measured, 1L)) %>%
  mutate(year_measured = NULL, harvest = NULL) %>%
  filter(!is.na(year2))

plotInfo2 <- data.table(plotInfo2)
setkey(plotInfo2, plot_id, year1, year2)

# disturbance
plot_disturb <- readRDS(paste0(path_data,"env/plot_disturb.RDS"))

plot_disturb <- plot_disturb %>% 
  subset(plot_id %in% sp.mat$plot_id)

# fire
fire <- plot_disturb %>% 
  dplyr::select(plot_id, fire_year) %>% replace_na(replace=list(fire_year=0))

fire <- data.table(fire)
fire[, fire_year2 := fire_year]

fire <- foverlaps(fire, plotInfo2, type = "within", by.x = names(fire), nomatch = 0L)

fire <- fire %>% dplyr::select(plot_id, year2) %>% 
  mutate(fire = 1) %>% 
  rename(year_measured = year2) %>%
  distinct()


# insect

insect <- plot_disturb %>% 
  dplyr::select(plot_id, insect_year) %>% replace_na(replace=list(insect_year=0))

insect <- data.table(insect)
insect[, insect_year2 := insect_year]

insect <- foverlaps(insect, plotInfo2, type = "within", by.x = names(insect), nomatch = 0L)

insect <- insect %>% dplyr::select(plot_id, year2) %>% 
  mutate(insect = 1) %>% 
  rename(year_measured = year2) %>%
  distinct()


# harvest
harvest <- mortalityTable %>% 
  subset(plot_id %in% sp.mat$plot_id) %>% 
  mutate(harvested = ifelse(harvested > 10, 1, 0)) %>%
  select_("plot_id", "year_measured", "harvested")

# plantation

plantation <- plot_disturb %>% 
  filter(startsWith(origine, "P")) %>% 
  dplyr::select(plot_id, an_origine)

plantation <- data.table(plantation)
plantation[, an_origine2 := an_origine]

plantation <- foverlaps(plantation, plotInfo2, type = "within", by.x = names(plantation), nomatch = 0L)

plantation <- plantation %>% dplyr::select(plot_id, year2) %>% 
  mutate(plantation = 1) %>% 
  rename(year_measured = year2) %>%
  distinct()


env_all_steps <- harvest %>%
  full_join(plantation,  by = c("plot_id", "year_measured")) %>% 
  full_join(fire, by = c("plot_id", "year_measured")) %>%
  full_join(insect,  by = c("plot_id", "year_measured")) %>%
  replace_na(replace=list(plantation=0, fire=0, insect=0)) %>%
  mutate_at(vars(harvested, plantation, fire, insect), funs(factor))


### DISTANCE MATRIX & TRANSFORMATION ####

# compute Hellinger transformation
sp.hel <- decostand(sp.mat[,MySpecies], "hel")


# 4 CLUSTERS
grClara_res <- clara(sp.hel, k=4, samples = 100, sampsize = 1000, pamLike = T)

### GROUP MEMBERSHIP ####
grClara <- grClara_res$clustering

## identify medoid
sp.hel.medoid <- grClara_res$medoids

gr.name.medoid <- paste0("cl_", colnames(sp.hel.medoid[, apply(sp.hel.medoid, 1, which.max)]))

gr.name <-c("Mixed", "Temperate","Pioneer", "Boreal")

grClara <- factor(grClara, levels = 1:4, 
                  labels = gr.name)

grClara_res$clustering <- grClara

## ---- transition_df ----

### TRANSITION DF ####

gr.foret <- as.factor(grClara) 

l.t1 <- sp.mat %>% mutate(id = row_number()) %>% group_by(plot_id) %>% slice(1)
l.t2 <- sp.mat %>% mutate(id = row_number()) %>% group_by(plot_id) %>% slice(n()-1)


tr <- addmargins(table(gr.foret[l.t1$id], gr.foret[l.t2$id], dnn=c("From", "To")))

# disturbance based on change in nb of individuals
disturb <- as.data.frame(aggregate(TOTAL ~ plot_id, sp.mat, FUN = function(x) diff(x)))

disturb$fac <- ifelse(disturb$TOTAL < (-15), "disturb", "undisturb")


Trans.df <- cbind.data.frame(sp.mat[, c("plot_id", "year_measured")], gr.foret)

Trans.df <- Trans.df %>%
  group_by(plot_id) %>%
  mutate(year1 = first(year_measured)) %>%
  mutate(year2 = lead(year_measured, 1L)) %>%
  mutate(From = first(gr.foret)) %>%
  mutate(To = lead(gr.foret, 1L)) %>%
  mutate(year_measured = NULL, gr.foret = NULL) %>%
  filter(!is.na(year2))


Trans.df <- mutate(Trans.df, transition = paste0(From, "->", To))

delta_clim[,c("TP", "PP")] <- scale(delta_clim[,c("TP", "PP")])
Trans.df = merge(Trans.df, delta_clim, by = c("plot_id", "year1", "year2"))
Trans.df = merge(Trans.df, env_all_steps, by.x = c("plot_id", "year2"), by.y = c("plot_id", "year_measured"))
Trans.df2 = merge(Trans.df, cbind.data.frame(plot_id=xy$plot_id, st_coordinates(xy)), by = "plot_id")
#2. Add time interval
Trans.df$time <- Trans.df$year2 - Trans.df$year1

coefs <-c("TP" , "I(TP^2)" ,"I(TP^3)"  , "PP" , "I(PP^2)" ,"I(PP^3)", "TP:PP","harvested","insect","fire","plantation", "X", "Y", "time")
contrib_form <- as.character()
tested_coefs <- as.character()

for(i in 1:length(coefs)){
  temp_coefs <- coefs[-i]
  tested_coefs  <- append(tested_coefs ,coefs[i])
  form <- paste0("To ~ ",paste(temp_coefs,collapse="+"))
  contrib_form <- append(contrib_form,form)
}
# Add full model
full_formula <- paste0("To ~ ",paste(coefs,collapse="+"))
null_formula <- "To ~ 1"

# Create DF of formula
formulas <- data.frame(tested=tested_coefs,formula=contrib_form,stringsAsFactors=FALSE)

Trans.df2$To <- relevel(Trans.df2$To, ref = "Temperate")

#5. Full Models
T_multi_full <- nnet::multinom(full_formula, data=Trans.df, subset= From=='Temperate')
M_multi_full <- nnet::multinom(full_formula, data=Trans.df, subset= From=='Mixed')
R_multi_full <- nnet::multinom(full_formula, data=Trans.df, subset= From=='Pioneer')
B_multi_full <- nnet::multinom(full_formula, data=Trans.df, subset= From=='Boreal')

#6. Null Models
T_multi_null <- nnet::multinom(null_formula,data=Trans.df, subset= From=='Temperate')
M_multi_null <- nnet::multinom(null_formula,data=Trans.df, subset= From=='Mixed')
R_multi_null <- nnet::multinom(null_formula,data=Trans.df, subset= From=='Pioneer')
B_multi_null <- nnet::multinom(null_formula,data=Trans.df, subset= From=='Boreal')

#7. Compute McFadden pseudo-R2
T_Rs <- 1 - (T_multi_full$deviance / T_multi_null$deviance)
M_Rs <- 1 - (M_multi_full$deviance / M_multi_null$deviance)
R_Rs <- 1 - (R_multi_full$deviance / R_multi_null$deviance)
B_Rs <- 1 - (B_multi_full$deviance / B_multi_null$deviance)

#8. Compute contributions to the full model
states <- c("M", "B", "T", "R")
ls_contrib <-list()

for (state in 1:length(gr.name)){
  contrib <- as.numeric()
  for(i in 1:nrow(formulas)){
    multi <- nnet::multinom(formulas$formula[i], data = Trans.df, subset = From==gr.name[state])
    contrib <- append(contrib, (multi$AIC - get(paste0(states[state], '_multi_full'))$AIC))
  }
  ls_contrib[[state]] <- contrib
}

#9. Summarize information in table and export in latex
contrib <- do.call("rbind", ls_contrib)
colnames(contrib) <- formulas$tested
rownames(contrib) <- states
contrib <- data.frame(contrib, AIC=c(T_multi_full$AIC, M_multi_full$AIC, R_multi_full$AIC, B_multi_full$AIC), R2=c(T_Rs,M_Rs,R_Rs,B_Rs))


# stepAIC
# graph fixer tout sauf TP et harvested
x = stepAIC(M_multi_full)

pp = 0
tp = seq(-3,3, length.out = 100)
harvested = c(0,1)
Y = c(46.7, 48.96, 51.16)
X = mean(Trans.df$X)

time=10
newdata <- data.frame(expand.grid(PP= pp, TP =tp, harvested = as.factor(harvested), insect= as.factor(insect),Y=Y, X=X, time = time))

prob = as.data.frame(predict(x, newdata = newdata, type = "probs", se.fit = T))

prob2 = cbind(newdata[,c("TP", "harvested", "Y")], prob)


lpp <- melt(prob2, id.vars = c("TP", "harvested", "Y"), value.name = "probability")
ggplot(lpp, aes(x = TP, y = probability, colour = harvested)) + 
  geom_line() + 
  facet_grid(variable ~.)
