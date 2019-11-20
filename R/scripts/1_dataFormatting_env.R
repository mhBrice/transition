### DATA FORMATTING FOR COMMUNITY TRANSITION ####

### PACKAGES ####
library(dplyr)
library(reshape2)
library(sf)
library(data.table)


### DATA ####

# Environmental

env_data <- readRDS("data/env_data_nov2019.RDS")

site_plantation <- env_data$ID_PE[which(env_data$ORIGINE2 %in% c("plantation", "wasteland"))]

env_data <- env_data %>%
  group_by(ID_PE) %>% arrange(year_measured) %>%
  filter(n()>1) %>% # remove plots that were sample only once
  filter(!(ID_PE %in% site_plantation))  %>% ungroup()

# Ecoregion

ecoreg_df <- readRDS("data/ecoreg_df_nov19.RDS") %>% 
  filter(ID_PE %in% env_data$ID_PE)

# Spatial

xy <- st_read("data/plot_xy32198_nov2019.gpkg") %>% 
  filter(ID_PE %in% env_data$ID_PE)
st_crs(xy) <- 32198

# Climate

# 10-year average
bioclim10 <- readRDS("data/bioclim_roll10_nov2019.RDS") %>% 
  filter(ID_PE %in% env_data$ID_PE) %>% 
  left_join(env_data[,1:4], by = c("ID_PE", "year" = "year_measured")) %>%
  ungroup()

###########################
### FORMATTING CLIMATE ####
###########################
# sg_15 <- mean temperature for period 3 (growing season)
# sg_12 <- Annual mean temperature
# sg_03 <- number of days of the growing season
# sg_06 <- total precipitation for period 3
# cmi_sum <- annual climate moisture index
# sum(cmi_05:cmi_09) <- climate moisture index for summer months

# Join to other climate variables

bioclim <- bioclim10 %>% 
  mutate(sCMI = rowSums(select(., "cmi_05":"cmi_09"))) %>% 
  select(ID_PE, ID_PE_MES, plot_id, year_measured = year,
         sTP = sg_15, sCMI)


################################
### FORMATTING DISTURBANCES ####
################################

# for msm, disturbance between t1 and t2 must be reported at t1

# Disturbance

env_data <- env_data %>%
  mutate(major_logging = 
           ifelse(is.na(ORIGINE2) | ORIGINE2!="logging", 0, 2),
         minor_logging = 
           ifelse(is.na(PERTURB2) | PERTURB2!="partial_logging", 0, 1)) %>%
  mutate(logging = major_logging + minor_logging) %>%
  mutate(major_natural = 
           ifelse(is.na(ORIGINE2) | ORIGINE2=="logging", 0, 2),
         minor_natural = 
           ifelse(is.na(PERTURB2) | PERTURB2=="partial_logging", 0, 1)) %>%
  mutate(natural = major_natural + minor_natural) %>%
  mutate_at(vars(logging, natural), ~ifelse(.>2, 2, .)) %>%
  select(-c(starts_with('major_'), starts_with('minor'), ORIGINE2, PERTURB2)) 

env_data <- env_data %>% group_by(plot_id) %>%
  arrange(year_measured) %>% 
  mutate_at(vars(logging, natural), lead) %>% 
  tidyr::replace_na(list(logging = 0, natural = 0))

################################
### SOIL ####
################################

# Drainage

env_data$DRAIN <- factor(env_data$CL_DRAI2,
                          levels = c("excessif", "rapide", 
                                     "bon", "modere",
                                     "imparfait", "mauvais", 
                                     "tres_mauvais", "complexe"),
                          ordered = T)
levels(env_data$DRAIN) <- c("1","1", "2", "3","4","5","6","4")
env_data$DRAIN <- as.numeric(env_data$DRAIN)


# Humus
env_data$TYPEHUMUS <- factor(env_data$TYPEHUMUS,
                             levels = c("MU", "MD", "MR", 
                                        "TO", "AN", "SO", "NA"))

env_data <- env_data %>% group_by(plot_id) %>%
  mutate_at(vars(TYPEHUMUS:CL_DRAI2, DRAIN), last) %>%
  select(ID_PE:POURCPIERR, DRAIN, logging:natural, ORIGINE, PERTURB) %>%
  ungroup()


# Join all environmental data
env_all <- env_data %>%
  left_join(bioclim, by = c("ID_PE", "plot_id", "ID_PE_MES", "year_measured")) %>%
  left_join(ecoreg_df, by = c("ID_PE", "plot_id")) %>%
  ungroup() %>%
  arrange(plot_id, year_measured)


saveRDS(env_all, "data/env_all.RDS")



