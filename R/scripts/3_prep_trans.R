### Data prep ####

### DATA ####

states_ba <- readRDS("data/states_envba.RDS") %>% arrange(plot_id, year_measured)

states <- c("Boreal", "Mixed", "Pioneer", "Temperate")

st_col <- c("#158282", "#A1BD93","#FEAC19", "#D43650")


# Select and subset

states_ba <- states_ba %>% 
  select(ID_PE:states_num6, 
         "sTP", "CMI",
         "natural", "logging", 
         "DRAIN","PH_HUMUS", "ecoreg3", "ecoreg6") %>% 
  na.omit() %>% mutate_if(is.matrix, as.vector) %>% 
  filter(complete.cases("DRAIN")) %>%
  group_by(ID_PE) %>% filter(n() > 1) %>% ungroup()

states_ba$states_num = states_ba$states_num75

states_ba$states = states_ba$states_ba75

# Scale variables 

#### \!/ scale after removing lines!! \!/ ####

states_ba <- states_ba %>% 
  mutate(year_measured = year_measured - 1970) %>% 
  mutate_at(vars("sTP", "CMI", "DRAIN", "PH_HUMUS"), scale)

sc_sTP <- c(attr(states_ba$sTP, "scaled:center"), attr(states_ba$sTP, "scaled:scale"))

states_ba <- states_ba %>% mutate_if(is.matrix, as.vector)

ecoregion <- st_read("data/ecoregion_simple.gpkg")
ecoregion$SOUS_DOM6 <- factor(ecoregion$SOUS_DOM6, c("Sugar maple-bitternut hickory",
                                                     "Sugar maple-basswood",
                                                     "Sugar maple-yellow birch",
                                                     "Balsam fir-yellow birch",
                                                     "Balsam fir-white birch",
                                                     "Spruce-moss"))
xy <- st_read("data/plot_xy32198_may2018.gpkg") %>%
  st_transform(st_crs(ecoregion)) %>% 
  filter(plot_id %in% states_ba$plot_id)

