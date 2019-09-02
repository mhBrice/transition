---
title: Supplementary Information
geometry: margin=1in
header-includes:
    - \usepackage{setspace}
    - \setstretch{1,5}
    - \usepackage{lscape}
---

# Supplementary tables

*Add a table with species and their group*

Table S1. List of R packages used.

| Packages | Main functions | Uses                                 | References             |
|----------|----------------|--------------------------------------|------------------------|
| msm      | msm            | Multi-state Markov models in continuous time | @jackson_multi-state_2011 |
|          | lrtest.msm     | Likelihood ratio test                | |
|          | pmatrix.msm    | Transition probability matrix        | |
|          | hazard.msm     | Calculate tables of hazard ratios for covariates on transition intensities  | |
|          | pnext.msm      | Probability of each state being next | |
|          | sojourn.msm    | Mean sojourn times from a multi-state model | |
| sf       |                | Manipulation and mapping of spatial data   | @pebesma_simple_2018 |
| pROC     | multiclass.roc | Compute multi-class AUC              | @robin_proc_2011 |
| scoring  | logscore       | Compute logarithmic score            | @merkle_choosing_2013 |

\pagebreak

Table S2. Frequency of all observed transitions between the four forest states during the study period. Transitions are from rows to columns.

|From/To   | Boreal| Mixed| Pioneer| Temperate| Total|
|:---------|------:|-----:|-------:|---------:|-----:|
|Boreal    |   9632|   210|    1171|        17| 11030|
|Mixed     |    131|  2121|     188|       656|  3096|
|Pioneer   |   1484|   345|    4839|       281|  6949|
|Temperate |     11|   383|     215|      6019|  6628|
|Total     |  11258|  3059|    6413|      6973| 27703|

Table S3. Frequency of observed transitions between the four forest states from the first to the final survey of each plot. Transitions are from rows to columns.

|          | Boreal| Mixed| Pioneer| Temperate| Total at $t_0$|
|:---------|------:|-----:|-------:|---------:|-----:|
|Boreal    |   3637|   181|     716|        33|  4567|
|Mixed     |     72|   528|      97|       446|  1143|
|Pioneer   |   1064|   211|    1285|       152|  2712|
|Temperate |     22|   186|      78|      1680|  1966|
|Total at $t_{final}$ |   4795|  1106|    2176|      2311| 10388|


Table S4. Table of risk ratios + CI

<!--**Table S3**. Frequency of transitions between the four forest states from the first to the last survey of each plot.
...-->

# Supplementary figures

![Spatial distribution of observed change in tree basal area for temperate (left) and boreal (right) species. Forest plots in blue (top) have seen an increase in the basal area of one of the two groups, whereas forest plots in red (bottom) have seen a decrease in basal area.](res/figSupp_map_BA_change.pdf)

<!-- add prevalence plot-->

![Probability of transition between forest states through time for natural disturbances (left) and logging (right) and for different levels (line type) as predicted by the best multi-state model. All other covariates are fixed at the average conditions found in the ecotone, i.e. the balsam fir-yellow birch domain.](res/figSupp_proba.pdf)

![Changes in forest state proportion at equilibrium (a-d) as well as change in the position of the Boreal-Temperate transition (e) along the temperature (latitudinal) gradient for different disturbance scenarios . Proportion of Boreal (blue) and Temperate + Mixed forests (red) for increasing frequency of natural disturbances (a,c) and logging (b,d). Increasing frequency of disturbances is illustrated as increasing color intensity (from pale, low intensity, to dark, high intensity). The circles at the top of each plot (a-d) indicate the position of the boundary between dominance of Boreal forests and dominance of Temperate and Mixed forests (i.e. the advancing front). The change in the position of these circles along increasing frequency of each disturbance type and intensity are illustrated in (e). All other covariates are fixed at the average conditions found in the ecotone, i.e. the balsam fir-yellow birch domain, to focus solely on the effect of disturbances along the temperature gradient. The colors at the top of the plots approximate the position of the bioclimatic domains along the temperature gradient.](res/figSupp_SS_dfreq.pdf)

![State contribution to forest turnover (see Fig. 7a,b in main text) along the temperature (latitudinal) gradient for different disturbance scenarios: minor (solid), moderate (dashed) and major (dotted) disturbances for both natural (a,c,e) and logging (b,d,f). All other covariates are fixed at the average conditions found in the ecotone, i.e. the balsam fir-yellow birch domain, to focus solely on the effect of disturbances along the temperature gradient. The turnover time of a state (or sojourn time) measures the time spent in this state before transitioning to the next. Long turnover time can translate to large resistance. Here, at any point along the gradient, state turnover time is scaled by the steady state distribution and the sum of all scaled state turnover gives the the turnover time of the transition matrix. The colors at the top of the plots approximate the position of the bioclimatic domains.](res/figSupp_contrib2turnover.pdf)


![State contribution to forest entropy (see Fig. 7c,d in main text) along the temperature (latitudinal) gradient for different disturbance scenarios: minor (solid), moderate (dashed) and major (dotted) disturbances for both natural (a,c,e) and logging (b,d,f). All other covariates are fixed at the average conditions found in the ecotone, i.e. the balsam fir-yellow birch domain, to focus solely on the effect of disturbances along the temperature gradient. The entropy of a state measures the incertitude of its next transition. Here, at any point along the gradient, state entropy is scaled by the steady state distribution and the sum of all scaled state entropy gives the the entropy of the transition matrix. The colors at the top of the plots approximate the position of the bioclimatic domains.](res/figSupp_contrib2entropy.pdf)


# References
