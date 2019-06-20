### Data prep ####

### DATA ####

states_ba <- readRDS("data/states_envba.RDS") %>% arrange(plot_id, year_measured)

states <- c("Boreal", "Mixed", "Pioneer", "Temperate")

st_col <- c("#158282", "#A1BD93","#FEAC19", "#D43650")

# Scale variables

states_ba <- states_ba %>% 
  mutate(year_measured = year_measured - 1970) %>% 
  mutate_at(vars(EPMATORG:age_mean, aTP:GSlength), scale) 

# Select and subset

states_ba <- states_ba %>% 
  select(ID_PE:states_num6, 
         "sTP", "GSlength", "sPP", "CMI",
         "TP_xmax","TP_xmin", "CMI_xmin",
         "natural", "logging", "disturb", 
         "natural_lag", "logging_lag", "disturb_lag", 
         "old_natural", "old_logging", "old_disturb", 
         "TYPEHUMUS", "PH_HUMUS", "DRAIN", "ecoreg3") %>% 
  na.omit() %>% mutate_if(is.matrix, as.vector) %>%
  filter(TYPEHUMUS %in% c("MU", "MD", "MR", "SO", "TO")) %>% # Subset plots for humus type
  mutate(TYPEHUMUS = droplevels(TYPEHUMUS)) #%>%
  #group_by(ID_PE) %>% filter(n() > 1)


# State transition matrix to compute transition matrix property

states_trans <- states_ba %>%
  group_by(plot_id) %>%
  mutate(year1 = year_measured + 1970) %>%
  mutate(year2 = lead(year_measured + 1970, 1L)) %>%
  mutate(time_interv = year2 - year1) %>%
  mutate(From = states_ba) %>%
  mutate(To = lead(states_ba, 1L)) %>%
  mutate(transition = paste0(From, "->", To)) %>%
  filter(!is.na(year2)) %>%
  select(ID_PE:plot_id, year1:transition, disturb, natural, logging)
