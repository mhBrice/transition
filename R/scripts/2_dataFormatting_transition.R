### DATA FORMATTING FOR COMMUNITY TRANSITION ####

### PACKAGES ####
require(dplyr)
require(reshape2)

### DATA ####

# Species
sp_mat <- readRDS("data/sp_mat_abun_jul2018.RDS")
sp_ba <- readRDS("data/sp_mat_ba_jul2018.RDS")

# Environmental data
env_all <- readRDS("data/env_all.RDS")

### FORMATTING SPECIES ####

### Keep only plots that have been sampled first before 1980 and last after 2000
id_plantation <- sp_mat$ID_PE[which(sp_mat$PICABI>0 | sp_mat$PINSYL > 0)]

sp_mat <- sp_mat %>%
  filter(ID_PE %in% env_all$ID_PE) %>% # keep same site as env_all
  filter(!(ID_PE %in% id_plantation)) %>% # remove remaining site with planted tree species
  group_by(plot_id) %>% arrange(year_measured) %>%
  filter(n()>1)  # remove plots that were sample only once



# Species selection
# Combine Sorbus americana and Sorbus decora as forest technicians are not able to distinguish these species
# Remove Malus sp (introduced, difficult to interpret)
# Keep all other species
# conservative treshold because 1. all species contribute to beta diversity and 2. rare species doesn't affect our analysis...
sp_mat <- sp_mat %>%
  mutate(SORSP = SORAME + SORDEC) %>%
  select(-MALSP, -SORAME, -SORDEC)

MySpecies <- names(which(colSums(sp_mat[,-c(1:4)])>1))

sp_mat <- sp_mat %>%
  select(ID_PE, ID_PE_MES, plot_id, year_measured, MySpecies)

sp_mat$TOTAL <- rowSums(sp_mat[,MySpecies])

## Basal area
sp_ba <- sp_ba %>%
  mutate(SORSP = SORAME + SORDEC) %>%
  select(-MALSP, -SORAME, -SORDEC)

sp_ba <- sp_ba %>% filter(ID_PE_MES %in% sp_mat$ID_PE_MES) %>%
  select(ID_PE, ID_PE_MES, plot_id, year_measured, MySpecies)

sp_ba$TOTAL <- rowSums(sp_ba[,MySpecies])


### DEFINING STATES BASED ON DOMINANT SPECIES ####

pioneer <- c("BETPAP", "BETPOP", "CRASP",
             "POPBAL", "POPDEL", "POPGRA", "POPTRE", "PRUPEN", "SALSP", "SORSP")

temperate <- c("ACENEG", "ACENIG", "ACEPEN", "ACERIN", "ACERUB", "ACESAC", "ACESPI", 
               "AMESP", "BETALL",
               "CARCAR", "CARCOR", "FAGGRA", "FRAAME", "FRANIG", "FRAPEN", "JUGCIN",
               "OSTVIR",
               "PICRUB", "PINRES", "PINRIG", "PINSTR",
               "PRUSER", "PRUVIR",
               "QUEALB", "QUEBIC", "QUEMAC", "QUERUB", "THUOCC",
               "TILAME", "TSUCAN", "ULMAME", "ULMRUB", "ULMTHO")

boreal <- c("ABIBAL","ALNRUG", "LARLAR", "PICGLA", "PICMAR", "PINBAN")



### SPECIES GROUP BASED ON ABUNDANCE ####

states_ab <- sp_mat %>%
  ungroup() %>%
  mutate(Boreal = rowSums(select(., boreal))) %>%
  mutate(Pioneer = rowSums(select(., pioneer))) %>%
  mutate(Temperate = rowSums(select(., temperate))) %>%
  tidyr::replace_na(replace = list(Temperate=0, Pioneer=0, Boreal=0)) %>%
  select(ID_PE:year_measured, Temperate, Pioneer, Boreal) %>%
  mutate(tot_TB = Temperate + Boreal) %>%
  mutate(states_ab = as.factor(case_when((Pioneer > tot_TB) ~ "Pioneer",
                               (Temperate/tot_TB >= .25 & Boreal/tot_TB >= .25) ~ "Mixed",
                               (Boreal/tot_TB > .75) ~ "Boreal",
                               (Temperate/tot_TB > .75) ~ "Temperate")))

### SPECIES GROUP BASED ON BASAL AREA ####

states_ba <- sp_ba %>%
  mutate(Boreal = rowSums(select(., boreal))) %>%
  mutate(Pioneer = rowSums(select(., pioneer))) %>%
  mutate(Temperate = rowSums(select(., temperate))) %>%
  tidyr::replace_na(replace = list(Temperate=0, Pioneer=0, Boreal=0)) %>%
  select(ID_PE:year_measured, Temperate, Pioneer, Boreal) %>%
  mutate(tot_TB = Temperate + Boreal) %>%
  mutate(states_ba75 = as.factor(case_when((Pioneer > tot_TB) ~ "Pioneer",
                                (Temperate/tot_TB >= .25 & Boreal/tot_TB >= .25) ~ "Mixed",
                                (Boreal/tot_TB > .75) ~ "Boreal",
                                (Temperate/tot_TB > .75) ~ "Temperate"))) %>%
  mutate(states_ba85 = as.factor(case_when((Pioneer > tot_TB) ~ "Pioneer",
                                         (Temperate/tot_TB >= .15 & Boreal/tot_TB >= .15) ~ "Mixed",
                                         (Boreal/tot_TB > .85) ~ "Boreal",
                                         (Temperate/tot_TB > .85) ~ "Temperate"))) %>%
  mutate(states_ba5 = as.factor(case_when((Pioneer > tot_TB) ~ "Pioneer",
                                         (Temperate/tot_TB >= .5 & Temperate/tot_TB != 1) ~ "Mixed_t",
                                         (Boreal/tot_TB >= .5 & Boreal/tot_TB != 1) ~ "Mixed_b",
                                         (Boreal/tot_TB == 1) ~ "Boreal",
                                         (Temperate/tot_TB == 1) ~ "Temperate"))) %>%
  mutate(states_ba6 = as.factor(case_when((Pioneer > tot_TB) ~ "Pioneer",
                                            (Temperate/tot_TB >= .25 & Boreal/tot_TB >= .25) ~ "Mixed",
                                          (Temperate/tot_TB >= .75 & Temperate/tot_TB != 1) ~ "Temperate_M",
                                          (Boreal/tot_TB >= .75 & Boreal/tot_TB != 1) ~ "Boreal_M",
                                          (Boreal/tot_TB == 1) ~ "Boreal",
                                          (Temperate/tot_TB == 1) ~ "Temperate")))


# Sites with BA less than  5m2/ha or less than 5 trees were considered pioneer whatever the species present

states_ba$states_ba75[which(sp_ba$TOTAL <= 5)] <- "Pioneer"
states_ba$states_ba85[which(sp_ba$TOTAL <= 5)] <- "Pioneer"
states_ba$states_ba5[which(sp_ba$TOTAL <= 5)] <- "Pioneer"
states_ba$states_ba6[which(sp_ba$TOTAL <= 5)] <- "Pioneer"


states_ab$states_ab[which(sp_mat$TOTAL <= 5)] <- "Pioneer"

# Order

states_ba <- states_ba %>%
  arrange(plot_id, year_measured) %>%
  distinct(plot_id, year_measured, .keep_all = T)

states_ab <- states_ab %>%
  arrange(plot_id, year_measured) %>%
  distinct(plot_id, year_measured, .keep_all = T)

# In msm, observed states should be numeric variables or factors called "1","2",...
states_ba$states_num75 <- as.factor(states_ba$states_ba75)
levels(states_ba$states_num75) <- c("1","2","3","4")

states_ba$states_num85 <- as.factor(states_ba$states_ba85)
levels(states_ba$states_num85) <- c("1","2","3","4")

states_ba$states_num5 <- as.factor(states_ba$states_ba5)
levels(states_ba$states_num5) <- c("1","2","3","4", "5")

states_ba$states_ba6 <- factor(states_ba$states_ba6, 
                               levels = c("Boreal", "Boreal_M", "Mixed", "Pioneer", "Temperate_M", "Temperate"))
states_ba$states_num6 <- states_ba$states_ba6
levels(states_ba$states_num6) <- c("1","2","3","4", "5","6")

states_ab$states_num <- as.factor(states_ab$states_ab)
levels(states_ab$states_num) <- c("1","2","3","4")

saveRDS(states_ba, "data/states_ba.RDS")
saveRDS(states_ab, "data/states_ab.RDS")

#Conifer- dominated cover types (>75% conifer, based on canopy coverage), which were labelled as “softwood”, “swamp softwood”, or “black spruce” on the 1930 map, were merged here into the “conifer” type. The 1930 map also included “deciduous” (>75% deciduous), “mixed” (>25% of both deciduous and conifers), and “no cover” (recently disturbed and naturally nonwooded areas) cover types. Boucher et al 2006



# Join state to environmental data

states_envba <- states_ba %>%
  left_join(env_all, by = c("ID_PE", "ID_PE_MES", "plot_id", "year_measured")) %>%
  arrange(plot_id, year_measured) %>%
  mutate_at(vars(disturb:natural_lag), factor)

states_envab <- states_ab %>%
  left_join(env_all, by = c("ID_PE", "ID_PE_MES", "plot_id", "year_measured")) %>%
  arrange(plot_id, year_measured) %>%
  mutate_at(vars(disturb:natural_lag), factor)



saveRDS(states_envba, "data/states_envba.RDS")
saveRDS(states_envab, "data/states_envab.RDS")

