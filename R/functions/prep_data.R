### Data prep ####

# Open requiered data and scale environmental variables

# STATE + ENV

states_ba <- readRDS("data/states_envba.RDS") %>% arrange(plot_id, year_measured)

# Scale variables 


states_ba <- states_ba %>% 
  mutate(year_measured = year_measured - 1970) %>% 
  mutate_at(vars("sTP", "CMI", "DRAIN", "PH_HUMUS"), scale)

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
xy <- st_read("data/plot_xy32198_may2018.gpkg") %>%
  st_transform(st_crs(ecoregion)) %>% 
  filter(plot_id %in% states_ba$plot_id)

