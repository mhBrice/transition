### DATA FORMATTING FOR COMMUNITY TRANSITION ####

### PACKAGES ####
library(dplyr)
library(reshape2)
library(sf)
library(data.table)


### DATA ####

# Environmental

env_data <- readRDS("data/env_data_may2018.RDS")

site_plantation <- env_data$ID_PE[which(env_data$ORIGINE2=="plantation")]

env_data <- env_data %>%
  group_by(ID_PE) %>% arrange(year_measured) %>%
  filter(n()>1) %>% # remove plots that were sample only once
  filter(!(ID_PE %in% site_plantation)) 

# Ecoregion

ecoreg_df <- readRDS("data/ecoreg_df.RDS") %>% filter(ID_PE %in% env_data$ID_PE)

# Spatial

xy <- st_read("data/plot_xy32198_may2018.gpkg") %>% filter(ID_PE %in% env_data$ID_PE)
st_crs(xy) <- 32198

# Climate

bioclim10 <- readRDS("data/bioclim10_mat.RDS") %>% filter(ID_PE %in% env_data$ID_PE)
bioclim_all <- readRDS("data/bioclim_corrected.RDS") %>% filter(ID_PE %in% env_data$ID_PE)



###########################
### FORMATTING CLIMATE ####
###########################

bioclim_all <- bioclim_all %>% 
  mutate(year = as.numeric(year)) %>%
  subset(year >= 1960)



# sg_15 <- mean temperature for period 3 (growing season)
# sg_06 <- total precipitation for period 3 (growing season)
# cmi <- climate moisture index
# bio_06 <- min temperature of coldest period

# sg_1 <- julian day number of start of the growing season
# sg_2 <- julian day number at end of the growing season
# sg_3 <- number of days of the growing season

# Join to other climate variables

bioclim <- bioclim10 %>% 
  select(ID_PE, ID_PE_MES, plot_id, year_measured,
         bio_01, bio_12,
         sg_15, sg_06, bio_06, cmi,
         sg_01, sg_02, sg_03) %>%
  rename(aTP = bio_01, aPP = bio_12,
         sTP = sg_15, sPP = sg_06, TP_min = bio_06, CMI = cmi,
         GS1 = sg_01, GS2 = sg_02, GSlength = sg_03)

### CLIMATE ANOMALIES ####

bioclim_all <- bioclim_all %>% 
  group_by(ID_PE) %>% 
  mutate(TP05mean6090 = mean(subset(bio_05, year <= 1990))) %>% 
  mutate(TP05sd6090 = sd(subset(bio_05, year <= 1990))) %>% 
  mutate(TP06mean6090 = mean(subset(bio_06, year <= 1990))) %>% 
  mutate(TP06sd6090 = sd(subset(bio_06, year <= 1990))) %>% 
  mutate(CMImean6090 = mean(subset(cmi, year <= 1990))) %>%
  mutate(CMIsd6090 = sd(subset(cmi, year <= 1990))) %>%
  mutate(year2 = year)

x_trans <- select(env_data, ID_PE:year_measured) %>%
  group_by(plot_id) %>%
  mutate(year = year_measured) %>%
  mutate(year2 = dplyr::lead(year_measured, 1L)) %>%
  mutate(year_measured=NULL) %>%
  filter(!is.na(year2)) %>% data.table



setkey(x_trans, ID_PE, year, year2)


bioclim_trans <- foverlaps(data.table(bioclim_all), x_trans, type="within", nomatch=NULL)

bioclim_trans <- bioclim_trans %>%
  group_by(ID_PE_MES) %>%
  mutate(CMI_xmin = ifelse(any(cmi < CMImean6090-2*CMIsd6090), "1", "0"),
         TP_xmax = ifelse(any(bio_05 > TP05mean6090+2*TP05sd6090), "1", "0"),
         TP_xmin = ifelse(any(bio_06 < TP06mean6090-2*TP06sd6090), "1", "0")) %>%
  select(ID_PE, ID_PE_MES, plot_id, year, year2,
         TP_xmax, TP_xmin, CMI_xmin) %>%
  distinct()

bioclim2 <- bioclim %>% 
  left_join(bioclim_trans, 
            by = c('ID_PE', "ID_PE_MES", "plot_id", "year_measured" = "year"))

bioclim3 <- bioclim2 %>% 
  group_by(ID_PE) %>% 
  arrange(year_measured) %>%
  mutate_at(vars(TP_xmax:CMI_xmin), zoo::na.locf) %>% select(-year2)

################################
### FORMATTING DISTURBANCES ####
################################

# for msm disturbance between t1 and t2 must be reported at t1

# Disturbance

env_data <- env_data %>%
  mutate(major_disturb = ifelse(is.na(ORIGINE2), 0, 2),
         minor_disturb = ifelse(is.na(PERTURB2), 0, 1)) %>%
  mutate(disturb = major_disturb + minor_disturb) %>%
  mutate(major_logging = ifelse(is.na(ORIGINE2) |ORIGINE2!="logging", 0, 2),
         minor_logging = ifelse(is.na(PERTURB2) |PERTURB2!="partial_logging", 0, 1)) %>%
  mutate(logging = major_logging + minor_logging) %>%
  mutate(major_natural = ifelse(is.na(ORIGINE2) | ORIGINE2=="logging", 0, 2),
         minor_natural = ifelse(is.na(PERTURB2) | PERTURB2=="partial_logging", 0, 1)) %>%
  mutate(natural = major_natural + minor_natural) %>%
  mutate_at(vars(disturb,logging,natural), funs(ifelse(.>2, 2, .))) %>%
  select(-c(starts_with('major_'), starts_with('minor'), ORIGINE2, PERTURB2)) 

env_data <- env_data %>% group_by(plot_id) %>%
  arrange(year_measured) %>%
  mutate(old_disturb = first(disturb), old_logging = first(logging), old_natural = first(natural)) %>%
  mutate(disturb_lag = disturb, logging_lag = logging, natural_lag = natural) %>% 
  mutate_at(vars(disturb, logging, natural), lead) %>% 
  tidyr::replace_na(list(disturb = 0, logging = 0, natural = 0, 
                         disturb_lag = 0, logging_lag = 0, natural_lag = 0))

################################
### SOIL ####
################################

# Drainage

env_data$DRAIN <- factor(env_data$CL_DRAI2,
                          levels = c("excessif", "rapide", "bon", "modere",
                                     "imparfait", "mauvais", "tres_mauvais", "complexe"),
                          ordered = T)
levels(env_data$DRAIN) <- c("1","1", "2", "3","4","5","6","4")
env_data$DRAIN <- as.numeric(env_data$DRAIN)


# Humus
env_data$TYPEHUMUS <- factor(env_data$TYPEHUMUS,
                             levels = c("MU", "MD", "MR", "TO", "AN", "SO", "NA"))

env_data <- env_data %>% group_by(plot_id) %>%
  mutate_at(vars(TYPEHUMUS:CL_DRAI2, DRAIN), last) %>%
  select(ID_PE:POURCPIERR, DRAIN, age_mean, disturb:natural_lag, ORIGINE, PERTURB) %>%
  ungroup()


################################
### PREVALENCE ####
################################

# Get the 100 nearest neighbors
# neigh <- nngeo::st_nn(xy, xy, k = 100, maxdist = 50000)
# 
# names(neigh) <- xy$plot_id
# neigh <- lapply(neigh, FUN = function(x) xy$plot_id[x])
# 
# which(lapply(neigh, length)<2)
# 
# neigh_B <- lapply(neigh, 
#                   FUN = function(x) mean(filter(states_ba, plot_id %in% x)$Boreal))
# 
# neigh_T <- lapply(neigh, 
#                   FUN = function(x) mean(filter(states_ba, plot_id %in% x)$Temperate))
# 
# neigh <- cbind.data.frame(plot_id = xy$plot_id,
#                           neigh_T = unlist(neigh_T), 
#                           neigh_B = unlist(neigh_B))


# Join all environmental data
env_all <- env_data %>%
  #left_join(neigh, by = "plot_id") %>%
  left_join(bioclim3, by = c("ID_PE", "plot_id", "ID_PE_MES", "year_measured")) %>%
  left_join(ecoreg_df, by = c("ID_PE", "plot_id")) %>%
  ungroup() %>%
  arrange(plot_id, year_measured)


saveRDS(env_all, "data/env_all.RDS")



