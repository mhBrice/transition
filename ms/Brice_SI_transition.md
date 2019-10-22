---
title: Supplementary Information
geometry: margin=1in
header-includes:
    - \usepackage{setspace}
    - \setstretch{1,5}
    - \usepackage{lscape}
    - \newcommand{\blandscape}{\begin{landscape}}
    - \newcommand{\elandscape}{\end{landscape}}
    - \usepackage{caption}
    - \usepackage{makecell}
    - \usepackage{colortbl}
    - \usepackage{xcolor}
    - \usepackage{float}
    - \floatplacement{figure}{H}
    - \renewcommand{\figurename}{Fig.}
    - \renewcommand\thefigure{S\arabic{figure}}
    - \renewcommand\thetable{S\arabic{table}}
---

# Supplementary tables

**Table S1**. List of species included in the analyses and their corresponding group. The species groups were defined using their trait values and knowledge of species ecology [see @brice_disturbances_2019 for details].

\begin{longtable}{>{\em}ll}
\toprule
\textbf{Species name} & \textbf{Vernacular name}\\
\midrule
\endfirsthead
\multicolumn{2}{@{}l}{\textit{(continued)}}\\
\toprule
\textbf{Species name} & \textbf{Vernacular name}\\
\midrule
\endhead
\
\endfoot
\bottomrule
\endlastfoot
\addlinespace[0.3em]
\multicolumn{2}{l}{\textbf{Boreal}}\\
\hspace{1em}Abies balsamea & Balsam fir\\
\hspace{1em}Picea glauca & White spruce\\
\hspace{1em}Picea mariana & Black spruce\\
\hspace{1em}Pinus banksiana & Jack pine\\
\hspace{1em}Larix laricina & Tamarack\\
\hspace{1em}Alnus incana & Speckled alder\\
\addlinespace[0.3em]
\multicolumn{2}{l}{\textbf{Pioneer}}\\
\hspace{1em}Betula papyrifera & White birch\\
\hspace{1em}Betula populifolia & Grey birch\\
\hspace{1em}Populus tremuloides & Trembling aspen\\
\hspace{1em}Populus deltoides & Cottonwood\\
\hspace{1em}Populus balsamifera & Balsam poplar\\
\hspace{1em}Populus grandidentata & Large tooth aspen\\
\hspace{1em}Salix sp. & Willow\\
\hspace{1em}Prunus pensylvanica & Pin cherry\\
\hspace{1em}Crataegus sp. & Hawthorn\\
\hspace{1em}Sorbus sp. & Mountain-ash\\
\addlinespace[0.3em]
\multicolumn{2}{l}{\textbf{Temperate}}\\
\hspace{1em}Picea rubens & Red spruce\\
\hspace{1em}Pinus resinosa & Red pine\\
\hspace{1em}Pinus strobus & Eastern white pine\\
\hspace{1em}Tsuga canadensis & Eastern hemlock\\
\hspace{1em}Ulmus americana & American elm\\
\hspace{1em}Ulmus rubra & Red elm\\
\hspace{1em}Ulmus thomasii & Rock elm\\
\hspace{1em}Carya cordiformis & Bitternut hickory\\
\hspace{1em}Juglans cinerea & Butternut\\
\hspace{1em}Quercus macrocarpa & Bur oak\\
\hspace{1em}Quercus alba & White oak\\
\hspace{1em}Quercus bicolor & Swamp white oak\\
\hspace{1em}Quercus rubra & Red oak\\
\hspace{1em}Fagus grandifolia & American beech\\
\hspace{1em}Betula alleghaniensis & Yellow birch\\
\hspace{1em}Carpinus caroliniana & Blue beech\\
\hspace{1em}Ostrya virginiana & Ironwood\\
\hspace{1em}Tilia americana & Basswood\\
\hspace{1em}Prunus serotina & Black cherry\\
\hspace{1em}Acer rubrum & Red maple\\
\hspace{1em}Acer saccharum & Sugar maple\\
\hspace{1em}Acer pensylvanicum & Striped maple\\
\hspace{1em}Acer saccharinum & Silver maple\\
\hspace{1em}Acer spicatum & Mountain maple\\
\hspace{1em}Fraxinus pennsylvanica & Red ash\\
\hspace{1em}Fraxinus americana & White ash\\
\hspace{1em}Fraxinus nigra & Black ash\\
\hspace{1em}Thuja occidentalis & White cedar\\
\hspace{1em}Amelanchier sp. & Serviceberry\\
\hspace{1em}Acer negundo & Manitoba maple\\
\hspace{1em}Acer nigrum & Black maple\\
\hspace{1em}Pinus rigida & Pitch pine\\
\hspace{1em}Prunus virginiana & Chokecherry\\*
\end{longtable}

\pagebreak

**Table S2**. List of R packages used.

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

**Table S3**. Frequency of all observed transitions between the four forest states during the study period. Transitions are from rows to columns.

|From/To   | Boreal| Mixed| Pioneer| Temperate| Total|
|:---------|------:|-----:|-------:|---------:|-----:|
|Boreal    |   9632|   210|    1171|        17| 11030|
|Mixed     |    131|  2121|     188|       656|  3096|
|Pioneer   |   1484|   345|    4839|       281|  6949|
|Temperate |     11|   383|     215|      6019|  6628|
|Total     |  11258|  3059|    6413|      6973| 27703|

<!--
Table S3. Frequency of observed transitions between the four forest states from the first to the final survey of each plot. Transitions are from rows to columns.

|          | Boreal| Mixed| Pioneer| Temperate| Total at $t_0$|
|:---------|------:|-----:|-------:|---------:|-----:|
|Boreal    |   3637|   181|     716|        33|  4567|
|Mixed     |     72|   528|      97|       446|  1143|
|Pioneer   |   1064|   211|    1285|       152|  2712|
|Temperate |     22|   186|      78|      1680|  1966|
|Total at $t_{final}$ |   4795|  1106|    2176|      2311| 10388|
-->
\pagebreak
\blandscape

**Table S4**. Table of risk ratios + CI

\begin{table}[H]
\centering\begingroup\fontsize{6.8}{8.8}\selectfont

\begin{tabular}{>{\bfseries}l>{\centering\arraybackslash}p{1.95cm}>{\centering\arraybackslash}p{1.95cm}>{\centering\arraybackslash}p{1.95cm}>{\centering\arraybackslash}p{1.95cm}>{\centering\arraybackslash}p{1.95cm}>{\centering\arraybackslash}p{1.95cm}>{\centering\arraybackslash}p{1.95cm}>{\centering\arraybackslash}p{1.95cm}>{\centering\arraybackslash}p{1.95cm}}
\toprule
\textbf{Transitions} & \textbf{Baseline} & \textbf{Temperature} & \textbf{CMI} & \textbf{Drainage} & \textbf{pH} & \textbf{Natural1} & \textbf{Natural2} & \textbf{Logging1} & \textbf{Logging2}\\
\midrule
Boreal - Boreal & -0.008\newline (-0.009 ,-0.007) & \cellcolor{white}{ } & \cellcolor{white}{\cellcolor{white}{ }} & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ }\\
Boreal - Mixed & 0.002\newline (0.001 ,0.003) & \cellcolor[HTML]{cccccc}{8.546\newline (6.406 ,11.401)} & \cellcolor[HTML]{cccccc}{\cellcolor[HTML]{cccccc}{1.47\newline (1.199 ,1.802)}} & \cellcolor[HTML]{cccccc}{0.714\newline (0.613 ,0.832)} & \cellcolor{white}{1.023\newline (0.867 ,1.207)} & \cellcolor[HTML]{cccccc}{2.803\newline (2.014 ,3.902)} & \cellcolor{white}{2.739\newline (0.999 ,7.506)} & \cellcolor[HTML]{cccccc}{2.784\newline (1.628 ,4.762)} & \cellcolor{white}{1.333\newline (0.028 ,63.96)}\\
Boreal - Pioneer & 0.006\newline (0.005 ,0.007) & \cellcolor{white}{1.000} & \cellcolor{white}{\cellcolor{white}{1.000}} & \cellcolor{white}{1.000} & \cellcolor{white}{1.000} & \cellcolor[HTML]{cccccc}{5.202\newline (4.208 ,6.431)} & \cellcolor[HTML]{cccccc}{29.474\newline (23.455 ,37.037)} & \cellcolor[HTML]{cccccc}{11.067\newline (7.631 ,16.05)} & \cellcolor[HTML]{cccccc}{164.842\newline (98.705 ,275.293)}\\
Mixed - Boreal & 0.005\newline (0.003 ,0.008) & \cellcolor{white}{0.784\newline (0.485 ,1.266)} & \cellcolor{white}{\cellcolor{white}{1.291\newline (0.989 ,1.686)}} & \cellcolor{white}{0.99\newline (0.783 ,1.251)} & \cellcolor{white}{0.792\newline (0.606 ,1.034)} & \cellcolor{white}{0.845\newline (0.41 ,1.741)} & \cellcolor{white}{1.056\newline (0.041 ,27.276)} & \cellcolor{white}{1.553\newline (0.842 ,2.863)} & \cellcolor{white}{0.669\newline (0.009 ,50.155)}\\
Mixed - Mixed & -0.034\newline (-0.039 ,-0.03) & \cellcolor{white}{ } & \cellcolor{white}{\cellcolor{white}{ }} & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ }\\
Mixed - Pioneer & 0.005\newline (0.004 ,0.007) & \cellcolor{white}{1.000} & \cellcolor{white}{\cellcolor{white}{1.000}} & \cellcolor{white}{1.000} & \cellcolor{white}{1.000} & \cellcolor[HTML]{cccccc}{1.939\newline (1.042 ,3.609)} & \cellcolor[HTML]{cccccc}{9.714\newline (2.933 ,32.17)} & \cellcolor[HTML]{cccccc}{2.27\newline (1.117 ,4.613)} & \cellcolor[HTML]{cccccc}{27.499\newline (16.906 ,44.73)}\\
Mixed - Temperate & 0.024\newline (0.021 ,0.027) & \cellcolor{white}{0.921\newline (0.756 ,1.122)} & \cellcolor[HTML]{cccccc}{\cellcolor[HTML]{cccccc}{0.785\newline (0.694 ,0.887)}} & \cellcolor{white}{0.968\newline (0.885 ,1.058)} & \cellcolor{white}{0.953\newline (0.879 ,1.034)} & \cellcolor[HTML]{cccccc}{2.401\newline (1.931 ,2.984)} & \cellcolor[HTML]{cccccc}{4.507\newline (2.249 ,9.032)} & \cellcolor[HTML]{cccccc}{3.225\newline (2.635 ,3.946)} & \cellcolor[HTML]{cccccc}{5.32\newline (3.309 ,8.553)}\\
Pioneer - Boreal & 0.028\newline (0.027 ,0.03) & \cellcolor{white}{0.985\newline (0.915 ,1.059)} & \cellcolor[HTML]{cccccc}{\cellcolor[HTML]{cccccc}{1.518\newline (1.428 ,1.615)}} & \cellcolor{white}{0.998\newline (0.949 ,1.05)} & \cellcolor[HTML]{cccccc}{0.934\newline (0.877 ,0.995)} & \cellcolor{white}{0.987\newline (0.803 ,1.212)} & \cellcolor[HTML]{cccccc}{0.728\newline (0.547 ,0.97)} & \cellcolor[HTML]{cccccc}{1.625\newline (1.277 ,2.068)} & \cellcolor[HTML]{cccccc}{2.05\newline (1.227 ,3.425)}\\
Pioneer - Mixed & 0.004\newline (0.003 ,0.005) & \cellcolor[HTML]{cccccc}{4.243\newline (3.399 ,5.297)} & \cellcolor[HTML]{cccccc}{\cellcolor[HTML]{cccccc}{1.62\newline (1.369 ,1.916)}} & \cellcolor{white}{0.959\newline (0.842 ,1.093)} & \cellcolor{white}{0.934\newline (0.827 ,1.055)} & \cellcolor{white}{1.234\newline (0.778 ,1.957)} & \cellcolor{white}{0.864\newline (0.365 ,2.046)} & \cellcolor[HTML]{cccccc}{3.267\newline (2.256 ,4.73)} & \cellcolor{white}{1.117\newline (0.599 ,2.083)}\\
Pioneer - Pioneer & -0.034\newline (-0.036 ,-0.032) & \cellcolor{white}{ } & \cellcolor{white}{\cellcolor{white}{ }} & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ }\\
Pioneer - Temperate & 0.001\newline (0.001 ,0.002) & \cellcolor[HTML]{cccccc}{13.666\newline (9.914 ,18.839)} & \cellcolor[HTML]{cccccc}{\cellcolor[HTML]{cccccc}{2.523\newline (2.032 ,3.132)}} & \cellcolor[HTML]{cccccc}{0.832\newline (0.709 ,0.976)} & \cellcolor{white}{0.989\newline (0.878 ,1.114)} & \cellcolor{white}{0.134\newline (0.007 ,2.534)} & \cellcolor{white}{0.776\newline (0.212 ,2.837)} & \cellcolor{white}{0.914\newline (0.409 ,2.046)} & \cellcolor{white}{0.382\newline (0.12 ,1.216)}\\
Temperate - Mixed & 0.014\newline (0.011 ,0.017) & \cellcolor[HTML]{cccccc}{0.466\newline (0.361 ,0.601)} & \cellcolor{white}{\cellcolor{white}{0.889\newline (0.769 ,1.029)}} & \cellcolor[HTML]{cccccc}{1.319\newline (1.159 ,1.502)} & \cellcolor[HTML]{cccccc}{0.8\newline (0.709 ,0.901)} & \cellcolor{white}{1.248\newline (0.763 ,2.042)} & \cellcolor{white}{2.054\newline (0.583 ,7.232)} & \cellcolor{white}{1.284\newline (0.968 ,1.705)} & \cellcolor[HTML]{cccccc}{2.796\newline (1.263 ,6.193)}\\
Temperate - Pioneer & 0.001\newline (0.001 ,0.002) & \cellcolor{white}{1.000} & \cellcolor{white}{\cellcolor{white}{1.000}} & \cellcolor{white}{1.000} & \cellcolor{white}{1.000} & \cellcolor{white}{1.362\newline (0.238 ,7.787)} & \cellcolor[HTML]{cccccc}{30.08\newline (8.99 ,100.652)} & \cellcolor[HTML]{cccccc}{3.942\newline (2.181 ,7.126)} & \cellcolor[HTML]{cccccc}{122.836\newline (78.146 ,193.084)}\\
Temperate - Temperate & -0.015\newline (-0.018 ,-0.013) & \cellcolor{white}{ } & \cellcolor{white}{\cellcolor{white}{ }} & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ }\\
\bottomrule
\end{tabular}
\endgroup{}
\end{table}

\elandscape

\pagebreak



# Supplementary figures

![Temporal trends in growing season temperatures (top) and annual climate moisture index (bottom). Grey lines represent averaged climate values across the 10,388 studied forest plots. Straight black lines show the fitted least-squared linear regression lines.](res/figS1_clim_trend.pdf)


\pagebreak

![Spatial distribution of observed change in tree basal area for temperate (left) and boreal (right) species. Forest plots in blue (top) have seen an increase in the basal area of one of the two groups, whereas forest plots in red (bottom) have seen a decrease in basal area.](res/figSupp_map_BA_change.pdf)


\pagebreak

![Probability of transition between forest states through time for natural disturbances (left) and logging (right) and for different levels (line type) as predicted by the best multi-state model. All other covariates are fixed at the average conditions found in the ecotone, i.e. the balsam fir-yellow birch domain.](res/figSupp_proba.pdf)

\pagebreak

![Changes in forest state proportion at equilibrium (a-d) as well as change in the position of the Boreal-Temperate transition (e) along the temperature (latitudinal) gradient for different disturbance scenarios . Proportion of Boreal (blue) and Temperate + Mixed forests (red) for increasing frequency of natural disturbances (a,c) and logging (b,d). Increasing frequency of disturbances is illustrated as increasing color intensity (from pale, low intensity, to dark, high intensity). The circles at the top of each plot (a-d) indicate the position of the boundary between dominance of Boreal forests and dominance of Temperate and Mixed forests (i.e. the advancing front). The change in the position of these circles along increasing frequency of each disturbance type and intensity are illustrated in (e). All other covariates are fixed at the average conditions found in the ecotone, i.e. the balsam fir-yellow birch domain, to focus solely on the effect of disturbances along the temperature gradient. The colors at the top of the plots approximate the position of the bioclimatic domains along the temperature gradient.](res/figSupp_SS_dfreq.pdf)

\pagebreak

![State contribution to forest turnover (see Fig. 7a,b in main text) along the temperature (latitudinal) gradient for different disturbance scenarios: minor (solid), moderate (dashed) and major (dotted) disturbances for both natural (a,c,e) and logging (b,d,f). All other covariates are fixed at the average conditions found in the ecotone, i.e. the balsam fir-yellow birch domain, to focus solely on the effect of disturbances along the temperature gradient. The turnover time of a state (or sojourn time) measures the time spent in this state before transitioning to the next. Long turnover time can translate to large resistance. Here, at any point along the gradient, state turnover time is scaled by the steady state distribution and the sum of all scaled state turnover gives the the turnover time of the transition matrix. The colors at the top of the plots approximate the position of the bioclimatic domains.](res/figSupp_contrib2turnover.pdf)

\pagebreak

![State contribution to forest entropy (see Fig. 7c,d in main text) along the temperature (latitudinal) gradient for different disturbance scenarios: minor (solid), moderate (dashed) and major (dotted) disturbances for both natural (a,c,e) and logging (b,d,f). All other covariates are fixed at the average conditions found in the ecotone, i.e. the balsam fir-yellow birch domain, to focus solely on the effect of disturbances along the temperature gradient. The entropy of a state measures the incertitude of its next transition. Here, at any point along the gradient, state entropy is scaled by the steady state distribution and the sum of all scaled state entropy gives the the entropy of the transition matrix. The colors at the top of the plots approximate the position of the bioclimatic domains.](res/figSupp_contrib2entropy.pdf)

\pagebreak


# Supplementary methods

### Performance of candidate models

We fitted all models on the full data sets but also used cross-validation to
estimate the predictive performance on held-out data. We used two statistics,
the area under the receiver operating characteristic (ROC) curve (AUC) and the
logarithmic scoring rule (LS), to assess the agreement between the observed
state and the models’ predictions. The AUC is a popular performance metric for
binary classifiers that measures the probability that a randomly drawn
member of state $s$ has a lower estimated probability of belonging to state $r$
than a randomly drawn member of state $r$. The AUC ranges from 0 to 1, where a
score of 1 indicates perfect discrimination, while a score of 0.5 is as good as
random. @hand_simple_2001 has extended the AUC method to multi-class problems.
For any pair of states $r$ and $s$, we can compute $\hat{A}(r|s)$, the
probability that a randomly drawn member of state $s$ has a lower estimated
probability of belonging to state $r$ than a randomly drawn member of state $r$.
We can measure the discrimination rate between all pairs of states by computing
the pairwise AUC:

$$\hat{A}(r,s) = [\hat{A}(r|s)+ \hat{A}(s|r)]/2$$

Averaging the pairwise AUC gives the overall multi-class AUC (hereafter mAUC) of
the model:

$$mAUC = \frac{2}{c(c-1)} \sum_{r<s} \hat{A}(r,s)$$

The LS was proposed by Good (1952) and is often used in weather forecasts
[@gneiting_strictly_2007]. While AUC is a function of different classification
thresholds, LS measures the degree to which predicted probabilities are close to
the observed outcomes. We computed a global score for each model:

$$LS = \frac{1}{N} \sum_{i=1}^N -log(P(S_i = s_i))$$

where $S_i$ is the random variable describing the state of the forest in the
$i^{th}$ plot and $s_i$ is the observed state. So, LS only depends upon the
predicted probability of the realised state and not on the probabilities
assigned to the other possible states. The score is very sensitive to incorrect
predictions: if a model predicted the observed state with a probability of 100%,
the score for that plot would be 0, while if a probability of zero was assigned
to the observed state, the score would go to infinity. Hence, this sensitivity
emphasises the differences between model predictions and strongly penalises a
model that only gives high probabilities to self-transitions.

To assess the quality of prediction for the four states individually, we
computed LS for each state $r$ where we summed the predicted probabilities $P(S_i
= r)$ if the observed state is indeed $r$ and $1 - P(S_i = r)$ otherwise.


We evaluated and compared the predictive performance of our five models using
the overall mAUC and the pairwise AUCs, as well as the overall LS and the
state-specific LS. These metrics were estimated using stratified *K*-fold
cross-validation [@burnham_model_2002]. We first stratified the data set by
bioclimatic domains to ensure that each fold was representative of the plot
geographical distribution and randomly split the data set in *k*=10 folds. The
cross-validation process was repeated *k* times, during which *k* − 1 folds were
used to train the models and the remaining fold was used to validate the model
predictions against the observed state transitions. The cross-validated
performance metrics were then averaged for each model.


### Results

Model evaluation using 10-fold
cross-validation revealed that including climate and disturbances improved
overall model predictive performances, while soil variables had a
negligible effect (*Fig. Supp*). All models were good at distinguishing Boreal from
Temperate (high pairwise AUC). Soil variables slightly help to predict Mixed and
Temperate states. Including climate variables help to distinguish Mixed from the
other states, while including disturbances help to distinguish Pioneer from the
other states, especially Boreal.

Table Sx. Comparisons of the five candidate multi-state models. The number of
parameters used in each model corresponds to the number of modelled transitions
(10) $\times$ the number of covariates. The $\Delta$AIC is the difference
between the Akaike information criterion of each model (AIC~m~) and the minimum
of AIC among all the models (AIC~min~): $\Delta$AIC = AIC~m~ – AIC~min~.
Multi-class area under the curve (mAUC) were obtained
through 10-fold cross-validation. Higher mAUC indicate better model
predictive performance. The best model is the one in bold with $\Delta$AIC = 0.

|             |Covariates        | Nb of parameters| -2 Log-likelihood|$\Delta$AIC|LR test|  mAUC|     LS|
|:------------|:-----------------|----------------:|-----------------:|---------:|-------:|-----:|-------:|
|Baseline     |Intercept         |               10|           32032.3|    6132.6|< 0.001 | 0.899|   0.578|
|Soil         |Drainage, pH      |               24|           31886.9|    6015.2|< 0.001 | 0.906|   0.576|
|Climate      |Temperature, CMI  |               24|           30438.7|    4566.9|< 0.001 | 0.921|   0.550|
|Disturbances |Natural, Logging  |               50|           27341.3|    1521.6|< 0.001 | 0.925|   0.495|
|**Full**     |**All**           |           **78**|       **25763.7**|   **0.0**|< 0.001 |**0.940**|**0.468**|


![Mean pairwise AUC (areas under the receiver operating characteristic curves) obtained through 10-fold cross-validation. Higher values indicate a better capacity to discriminate between the four forest states: (B)oreal, (M)ixed, (P)ioneer and (T)emperate. The overall mAUC of each model is given next to the legend.](res/figSupp_cv_auc.pdf)

**AND/OR**

![Mean state-specific logarithmic skill score where each model including covariates is compared to the baseline model. Values were obtained through 10-fold cross-validation. Higher values indicate a larger improvement (predicted probabilities are closer to the observed outcomes) compared to the baseline model. The overall logarithmic score of each model is given next to the legend.](res/figSupp_cv_LSS.pdf)

**OR**

![Mean state-specific logarithmic score obtained through 10-fold cross-validation. Lower values indicate that the predicted probabilities are closer to the observed outcomes. The overall logarithmic score of each model is given next to the legend.](res/figSupp_cv_LS.pdf)


# References
