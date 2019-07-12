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
         "natural", "logging", "disturb", 
         "natural_lag", "logging_lag", "disturb_lag", 
         "TYPEHUMUS",  "DRAIN", "ecoreg3") %>% 
  na.omit() %>% mutate_if(is.matrix, as.vector) %>%
  filter(TYPEHUMUS %in% c("MU", "MD", "MR", "SO", "TO")) %>% # Subset plots for humus type
  mutate(TYPEHUMUS = droplevels(TYPEHUMUS)) #%>%
  #group_by(ID_PE) %>% filter(n() > 1)

# Remove test data (fifth inventory)
train_dat <- filter(states_ba, year_measured+1970 <= 2015) %>%
  group_by(ID_PE) %>% filter(n() > 1)
test_dat <- filter(states_ba, year_measured+1970 >= 2015) %>%
  filter(ID_PE %in% train_dat$ID_PE)
 

# State transition matrix to compute transition matrix property

# states_trans <- states_ba %>%
#   group_by(plot_id) %>%
#   mutate(year1 = year_measured + 1970) %>%
#   mutate(year2 = lead(year_measured, 1L) + 1970) %>%
#   mutate(time_interv = year2 - year1) %>%
#   mutate(From = states_ba) %>%
#   mutate(To = lead(states_ba, 1L)) %>%
#   mutate(transition = paste0(From, "->", To)) %>%
#   filter(!is.na(year2)) %>%
#   select(ID_PE:plot_id, year1:transition, disturb, natural, logging)
