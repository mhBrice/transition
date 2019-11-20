### Data prep ####

# Open requiered data and scale environmental variables

# STATE + ENV

states_ba <- readRDS("data/states_envba.RDS") %>% arrange(plot_id, year_measured)

# Scale variables 


states_ba <- states_ba %>% 
  mutate(year_measured = year_measured - 1970) %>% 
  select(ID_PE:states_num90, 
         "sTP", "sCMI", 
         "natural", "logging", 
         "DRAIN","PH_HUMUS", "ecoreg3", "ecoreg6") %>% 
  na.omit() %>%
  group_by(ID_PE) %>% filter(n() > 1) %>% ungroup() %>% 
  mutate_at(vars("sTP", "sCMI", "DRAIN", "PH_HUMUS"), scale) %>%
  mutate_at(vars(logging:natural), factor)

sc_sTP <- c(attr(states_ba$sTP, "scaled:center"), attr(states_ba$sTP, "scaled:scale"))

states_ba <- states_ba %>% mutate_if(is.matrix, as.vector)


# STATE NAMES + COLORS

states <- c("Boreal", "Mixed", "Pioneer", "Temperate")

st_col <- c("#158282", "#A1BD93", "#FEAC19", "#D43650")


# ECOREGION MAP
ecoregion <- st_read("data/ecoregion_simple.gpkg")
ecoregion$SOUS_DOM6 <- factor(ecoregion$SOUS_DOM6, c("Sugar maple-bitternut hickory",
                                                     "Sugar maple-basswood",
                                                     "Sugar maple-yellow birch",
                                                     "Balsam fir-yellow birch",
                                                     "Balsam fir-white birch",
                                                     "Spruce-moss"))

# COORDINATES
xy <- st_read("data/plot_xy32198_nov2019.gpkg") %>%
  st_transform(st_crs(ecoregion)) %>% 
  filter(ID_PE %in% states_ba$ID_PE)

