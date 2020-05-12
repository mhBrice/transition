### Analysis pipeline ####

### Load packages
source("R/functions/packages.R")

### Data formatting for environmental data
source("R/scripts/1_dataFormatting_env.R")

### Data formatting for community transition
source("R/scripts/2_dataFormatting_transition.R")

### Continuous time markov model
source("R/scripts/3_model_msm.R")

### Model validation
source("R/scripts/4_msm_valid.R")

### Map
# Figure 1
source("R/scripts/fig1_map.R")

### Conceptual figure
# Figure 2
source("R/scripts/fig2_transition_diagram.R")

### Plot model baseline intensities
# Figure 3
source("R/scripts/fig3_baseline.R")

### Plot best model
# Table 2
# Figure 4
# Figure S6
# Figure S7
source("R/scripts/5_plot_model.R")



### Plot steady-state proportion
# Figure 5
source("R/scripts/6_plot_steady.R")

### Plot transient measures
# Figure 6
# Figure S8
# Figure S9
source("R/scripts/7_plot_transient.R")

### Demographical processes & community transition
# Figure 7
source("R/scripts/8_demography.R")

### Supplementary figures and tables
# Figures S1, S2, S3, S4
# Tables S1, S4
source("R/scripts/fig_S1_disturbances.R")
source("R/scripts/fig_S2_time_interval.R")
source("R/scripts/fig_S3_clim_trend.R")
source("R/scripts/fig_S4_predTB.R")
source("R/scripts/tableS1_speciesGr.R")
source("R/scripts/tableS4_model.R")
