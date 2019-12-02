### Data prep ####

# 1. STATE + ENV

states_ba <- readRDS("data/states_envba.RDS") %>% arrange(plot_id, year_measured)

# 2. SCALE VARIABLES 

states_ba <- states_ba %>% 
  mutate(year_measured = year_measured - 1970) %>% 
  select(ID_PE:states_num85, 
         "sTP", "sCMI", 
         "natural", "logging", 
         "DRAIN","PH_HUMUS", "ecoreg3", "ecoreg6") %>% 
  na.omit() %>%
  group_by(ID_PE) %>% filter(n() > 1) %>% ungroup() %>% 
  mutate_at(vars("sTP", "sCMI", "DRAIN", "PH_HUMUS"), scale) %>%
  mutate_at(vars(logging:natural), factor)

sc_sTP <- c(attr(states_ba$sTP, "scaled:center"), attr(states_ba$sTP, "scaled:scale"))

states_ba <- states_ba %>% mutate_if(is.matrix, as.vector)


# 3. STATE NAMES + COLORS

states <- c("Boreal", "Mixed", "Pioneer", "Temperate")

st_col <- c("#158282", "#A1BD93", "#FEAC19", "#D43650")


