# Paper: Moderate disturbances accelerate forest transition dynamics under climate change at the temperate-boreal ecotone

This repository includes the data and R scripts to reproduce the analyses and figures found in the article *Moderate disturbances accelerate forest transition dynamics under climate change at the temperate-boreal ecotone* by Brice, Vissault, Vieira, Gravel, Legendre and Fortin submitted to Global Ecology and Biogeography.

## Installation

The analyses were carried out with [R (a free software environment for statistical computing and graphics)](https://www.r-project.org/) and require the installation of a recent version of it.

Analyses were reproduced in MacOSX Catalina.

<details>
R version 3.6.1 (2019-07-05)
Platform: x86_64-apple-darwin19.0.0 (64-bit)
Running under: macOS Catalina 10.15.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /usr/local/Cellar/openblas/0.3.7/lib/libopenblasp-r0.3.7.dylib

locale:
[1] en_CA.UTF-8/en_CA.UTF-8/en_CA.UTF-8/C/en_CA.UTF-8/en_CA.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] scoring_0.6            latex2exp_0.4.0        kableExtra_1.1.0      
 [4] knitr_1.25             sf_0.8-0               reshape2_1.4.3        
 [7] dplyr_0.8.3            data.table_1.12.6      RColorBrewer_1.1-2    
[10] scales_1.0.0           graphicsutils_1.4.9000 diagram_1.6.4         
[13] shape_1.4.4            msm_1.6.7             

loaded via a namespace (and not attached):
 [1] tidyselect_0.2.5   xfun_0.10          purrr_0.3.3        splines_3.6.1     
 [5] lattice_0.20-38    colorspace_1.4-1   vctrs_0.2.0        expm_0.999-4      
 [9] viridisLite_0.3.0  htmltools_0.4.0    survival_2.44-1.1  rlang_0.4.1       
[13] e1071_1.7-2        pillar_1.4.2       glue_1.3.1         DBI_1.0.0         
[17] plyr_1.8.4         stringr_1.4.0      munsell_0.5.0      rvest_0.3.5       
[21] mvtnorm_1.0-11     evaluate_0.14      class_7.3-15       Rcpp_1.0.2        
[25] KernSmooth_2.23-15 readr_1.3.1        backports_1.1.5    classInt_0.4-2    
[29] webshot_0.5.1      hms_0.5.2          digest_0.6.22      stringi_1.4.3     
[33] grid_3.6.1         tools_3.6.1        magrittr_1.5       tibble_2.1.3      
[37] crayon_1.3.4       pkgconfig_2.0.3    zeallot_0.1.0      Matrix_1.2-17     
[41] xml2_1.2.2         assertthat_0.2.1   rmarkdown_1.16     httr_1.4.1        
[45] rstudioapi_0.10    R6_2.4.0           units_0.6-5        compiler_3.6.1
</details>

The following packages must be installed:

- msm
- scoring
- data.table
- dplyr
- reshape2
- sf
- knitr
- kableExtra
- scales
- RColorBrewer
- diagram
- [graphicsutils](https://github.com/inSileco/graphicsutils) (not on CRAN)

Below are the R commands to install them all:

```R
install.packages(
  "msm", "scoring",
  "data.table", "dplyr",  "reshape2",  "sf",
  "scales", "RColorBrewer", "diagram",
  "knitr", "kableExtra", "remotes"
)
remotes::install_github("inSileco/graphicsutils")
```


## Guidelines

Several steps required for data preparation were performed using scripts
available in the [Quebec_data
repository](https://github.com/mhBrice/Quebec_data). The rest of the data
manipulations are performed directly in this repository:

1. You can reproduce all the steps to tidy data:

```R
source("R/scripts/1_dataFormatting_env.R")
source("R/scripts/2_dataFormatting_transition.R")
```

2. Then, run the multi-state model:

```R
source("R/scripts/3_model_msm.R")
```

3. The cross-validation of the candidate models takes a fair amount of time, and this step is not necessary for the other scripts to work. but you may run:

```R
source("R/scripts/4_msm_valid.R")
```

4. Then, figures of the main text (as well as fig. SX) are obtained as follow:

```R
source("R/scripts/fig1_map.R")
source("R/scripts/fig2_transition_diagram.R")
source("R/scripts/fig3_baseline.R")
source("R/scripts/5_plot_model.R")
source("R/scripts/6_plot_steady.R")
source("R/scripts/7_plot_transient.R")
```

5. Finally, run the following scripts to obtain supplementary figures and tables S1 and S4:

```R
source("R/scripts/fig_S1_clim_trend.R")
source("R/scripts/fig_S2_disturbances.R")
source("R/scripts/figS7_deltaTB.R")

source("R/scripts/tableS1_speciesGr.R")
source("R/scripts/tableS4_model.R")
```

Figures and tables are saved in [res](https://github.com/mhBrice/transition/tree/master/res).