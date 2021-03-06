### DATA FORMATTING FOR COMMUNITY TRANSITION ####

### PACKAGES ####

source("R/functions/packages.R")

### DATA ####

# Species
sp_ba <- readRDS("data/sp_mat_ba_nov2019.RDS")

# Environmental data
env_all <- readRDS("data/env_all.RDS")

### PLOT SELECTION ####


id_plantation <- sp_ba$ID_PE[which(sp_ba$PICABI > 0 | sp_ba$PINSYL > 0)]

sp_ba <- sp_ba %>%
  filter(ID_PE %in% env_all$ID_PE) %>% # keep same site as env_all
  filter(!(ID_PE %in% id_plantation)) %>% # remove remaining site with planted tree species
  group_by(ID_PE) %>% arrange(year_measured) %>%
  filter(n() > 1) %>% ungroup() # remove plots that were sample only once



# Species selection
# Combine Sorbus americana and Sorbus decora as forest technicians are not able to distinguish these species
# Remove Malus sp (introduced, difficult to interpret)
# Remove shrub species
# Keep all other species
# conservative treshold because rare species doesn't affect our analysis...
sp_ba <- sp_ba %>%
  mutate(SORSP = SORAME + SORDEC) %>%
  select(-MALSP, -SORAME, -SORDEC, -ALNCRI, -CRASP, -ALNRUG, -PRUVIR)

MySpecies <- sort(names(which(colSums(sp_ba[,-c(1:4)])>0)))

sp_ba <- sp_ba %>%
  select(ID_PE, ID_PE_MES, plot_id, year_measured, all_of(MySpecies))

sp_ba$TOTAL <- rowSums(sp_ba[, MySpecies])


### ASSIGN SPECIES TO GROUP ####

pioneer <- c("BETPAP", "BETPOP",
             "POPBAL", "POPDEL", "POPGRA", "POPTRE",
             "PRUPEN", "SALSP", "SORSP")

temperate <- c("ACENEG", "ACENIG", "ACEPEN", "ACERIN", "ACERUB", "ACESAC", "ACESPI",
               "AMESP", "BETALL",
               "CARCAR", "CARCOR", "FAGGRA",
               "FRAAME", "FRANIG", "FRAPEN",
               "JUGCIN", "OSTVIR",
               "PICRUB", "PINRES", "PINRIG", "PINSTR",
               "PRUSER",
               "QUEALB", "QUEBIC", "QUEMAC", "QUERUB",
               "THUOCC", "TILAME", "TSUCAN",
               "ULMAME", "ULMRUB", "ULMTHO")

boreal <- c("ABIBAL", "LARLAR", "PICGLA", "PICMAR", "PINBAN")

c(pioneer,temperate,boreal) %in% MySpecies
MySpecies %in% c(pioneer,temperate,boreal)

### DEFINING STATES BASED ON DOMINANT SPECIES GROUP ####
#"Conifer- dominated cover types (>75% conifer, based on canopy coverage), which were labelled as “softwood”, “swamp softwood”, or “black spruce” on the 1930 map, were merged here into the “conifer” type. The 1930 map also included “deciduous” (>75% deciduous), “mixed” (>25% of both deciduous and conifers), and “no cover” (recently disturbed and naturally nonwooded areas) cover types." (Boucher et al 2006)


states_ba <- sp_ba %>%
  mutate(Boreal = rowSums(select(., boreal))) %>%
  mutate(Pioneer = rowSums(select(., pioneer))) %>%
  mutate(Temperate = rowSums(select(., temperate))) %>%
  tidyr::replace_na(replace = list(Temperate=0, Pioneer=0, Boreal=0)) %>%
  select(ID_PE:year_measured, Temperate, Pioneer, Boreal) %>%
  mutate(tot_TB = Temperate + Boreal) %>%
  mutate(states_ba = as.factor(case_when((Pioneer > tot_TB) ~ "Pioneer",
                                (Temperate/tot_TB >= .25 & Boreal/tot_TB >= .25) ~ "Mixed",
                                (Boreal/tot_TB > .75) ~ "Boreal",
                                (Temperate/tot_TB > .75) ~ "Temperate"))) %>%
  mutate(states_ba85 = as.factor(case_when((Pioneer > tot_TB) ~ "Pioneer",
                                         (Temperate/tot_TB >= .15 & Boreal/tot_TB >= .15) ~ "Mixed",
                                         (Boreal/tot_TB > .85) ~ "Boreal",
                                         (Temperate/tot_TB > .85) ~ "Temperate")))


# Sites with BA less than  5m2/ha were considered pioneer whatever the species present

states_ba$states_ba[which(sp_ba$TOTAL <= 5)] <- "Pioneer"
states_ba$states_ba85[which(sp_ba$TOTAL <= 5)] <- "Pioneer"


# Order

states_ba <- states_ba %>%
  arrange(ID_PE, year_measured) %>%
  distinct(ID_PE, year_measured, .keep_all = TRUE)

# In msm, observed states should be numeric variables or factors called "1","2",...
states_ba$states_num <- as.factor(states_ba$states_ba)
levels(states_ba$states_num) <- c("1", "2", "3", "4")

states_ba$states_num85 <- as.factor(states_ba$states_ba85)
levels(states_ba$states_num85) <- c("1", "2", "3", "4")




### ENVIRONMENT ####

# Join state to environmental data

states_envba <- states_ba %>%
  left_join(env_all, by = c("ID_PE", "ID_PE_MES", "plot_id", "year_measured")) %>%
  arrange(ID_PE, year_measured) %>%
  mutate_at(vars(logging:natural), factor)

### SELECT & SUBSET ####

states_envba <- states_envba %>%
  select(ID_PE:states_num85,
         "sTP", "sCMI",
         "natural", "logging",
         "DRAIN","PH_HUMUS", "ecoreg3", "ecoreg6") %>%
  na.omit() %>%
  filter(complete.cases("DRAIN")) %>%
  group_by(ID_PE) %>% filter(n() > 1) %>% ungroup()


saveRDS(states_envba, "data/states_envba.RDS")
