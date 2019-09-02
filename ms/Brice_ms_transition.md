---
title: Moderate disturbances accelerate forest transition dynamics under climate change
documentclass: article
font: 12pt
papersize: a4paper
geometry: margin=1in
header-includes:
    - \usepackage{setspace}
    - \setstretch{1,5}
    - \usepackage{lineno}
    - \linenumbers


---

## List of authors  {-}

Marie-Hélène Brice ^1,2^
Steve Vissault ^2,3^
Willian Vieira ^2,3^
Dominique Gravel ^2,3^
Pierre Legendre ^1,2^
Marie-Josée Fortin ^4^


## Institutional affiliations {-}

1. Département de Sciences Biologiques, Université de Montréal, Montréal, Québec,  Canada.
2. Québec Centre for Biodiversity Sciences, McGill University, Montréal, Québec, Canada.
3. Département de biologie, Université de Sherbrooke, Sherbrooke, Québec, Canada.
4. Department of Ecology and Evolutionary Biology, University of Toronto, Toronto, Ontario, Canada.

## Contact Information {-}

Marie-Hélène Brice

- email: marie-helene.brice@umontreal.ca

## Running title
?

# Abstract [draft]

Several temperate tree species in North America are expected to expand their
distribution northward, where the boreal forest is already established and
dominated by conifers. However, these transitions from boreal to mixed and from
mixed to temperate forests could be hampered by non-climatic factors, such as
unsuitable soil conditions and competition by resident species.

In this paper, we describe the state transition dynamics of Quebec's forest
communities in recent decades and identify the factors and processes that
promote or prevent these state transitions. Specifically, we ask (1) Is forest
transition dynamic affected by recent climate change? (2) How do disturbances
and soil characteristics influence the transition probability among forest
states under climate change? Can natural or anthropogenic disturbances
accelerate climate related transitions? And, conversely, can soil
characteristics constrain these transitions? (3) How do different disturbance
type and intensity influence the transient dynamics and equilibrium distribution
of forest states?

To answer these questions, we analyzed the transition dynamics in over 10,000
forest inventory plots of southern Quebec (1970-2016). Using a continuous time
Markov multi-state model, we investigated the relationship between the
transition probabilities between the different states of the forest community
(temperate, boreal, mixed and pioneer) and variables related to climate, local
soil conditions and disturbances. We also describe the long-term equilibrium
state distribution under different disturbance scenarios, as well as the
transient dynamics using complementary measures: the time to reach equilibrium,
the turnover time and the entropy of the modeled system.

The transition probabilities are mainly related by natural and anthropogenic
disturbances and secondarily to climatic variables, whereas soil characteristics
(drainage and pH) were not a strong constraint. Moreover, moderate disturbances
increased the probability of transition from Mixed to Temperate forests and thus
accelerate long-term climate-induced transitions. [equilibrium... transient dynamic...]

Moderate disturbances were found to catalyze rapid transitions in forest communities and accelerate shift in response to recent climate change.
In ecotones, areas of ecological tension, disturbances may provide opportunities
for some migrating species to establish in otherwise competitive environments...

## Keywords

Climate change,
Disturbances,
Forest,
Continuous-time Markov model,
Multi-state model,
Québec,
Temperate-boreal ecotone,
Transition probabilities,
Resilience

<!-- fathom
upheaval = major change, catastrophe
Perhaps-->

\pagebreak

# Introduction

## *P1. Forest dynamics under CC*

Global climate warming is forcing species to move [@parmesan_globally_2003]. In
the temperate-boreal forest ecotone, several temperate deciduous tree species
are slowly migrating northward, colonizing conifer dominated forests
[@fisichelli_temperate_2014; @evans_borealtemperate_2017;
@boisvert-marsh_shifting_2014; @sittaro_tree_2017]. As climate warms and tips
the balance in favor of temperate over boreal species, forests may transition
from coniferous to mixedwood and from mixedwood to temperate deciduous
[@boulanger_climate_2019; @price_anticipating_2013]. While boreal forest dynamic
is characterized by broad-scale disturbances, mainly fires and insect outbreaks,
and slow decomposition of an acidic and nutrient poor litter, temperate forest
dynamic is characterized by small scale canopy gaps and rapid decomposition of a
rich litter [@goldblum_deciduous_2010]. Hence, as ecological processes strongly
differ among those biomes, climate-induced range shifts not only impact species
distributions, but also alter the structure of communities, microclimates,
biogeochemical cycles and ecosystem functioning and thus trigger a “regime
shift” [@scheffer_catastrophic_2001].

While tree climate niches are expected to shift northward by several hundred
kilometers by the end of the century [@mckenney_potential_2007], multiple
studies indicate that tree migration will not keep pace with global warming
[e.g. @zhu_failure_2012; @woodall_assessing_2013; @sittaro_tree_2017;
@renwick_temporal_2015; @talluto_extinction_2017; @vissault_biogeographie_2016].
The slow response of forests to environmental changes is not surprising given
that trees cannot walk, live for several decades to centuries and disperse over
very short distance. Hence, if forests are undisturbed, the speed of transitions
between forest biomes will be mainly limited by the natural turnover rate of the
resident species as well as the dispersal and establishment rates of the
migrating species [@neilson_transient_1993]. However, as the disequilibrium
between climate conditions and forest composition grows larger, alternative
stable states become possible and forests may fail to return to their previous
state following a disturbance [@johnstone_changing_2016].


## *P2. Disturbances*

Both gradual and abrupt transitions from one biome to another will likely take
place concurrently, but, as forests are increasingly subjected to both pervasive
climatic stresses and direct human-related disturbances, catastrophic
transitions are likely to play a dominant role in driving the climate shift in
biomes. Established communities are generally though to be resistant to
competitive displacement and resilient to commonly experienced disturbances
[@grondin_have_2018; @seidl_disturbance_2014]. Resistance can be defined as the
ability of an ecosystem to persist through time following a disturbance, whereas
resilience is the capacity to recover its pre-disturbance composition
[@gunderson_ecological_2000; @holling_resilience_1973]. Because the ability to
persist and recover following a pulse disturbance (e.g., fire or logging) can be
altered by a press disturbance (e.g., climate change), investigating resistance
and resilience can provide precious insights into their interacting effects.

However, as climate change weaken forest resilience, disturbances can provide
niche opportunities for migrating species and community composition can shift
abruptly to species that are better suited to current conditions
[@johnstone_changing_2016; @renwick_temporal_2015; @turner_disturbance_2010].
For example, canopy gaps have been shown to locally facilitate establishment of
temperate species in mixed forests of Ontario [@leithead_northward_2010]. In
Alaska, white spruce (*Picea glauca*) is invading black spruce (*Picea mariana*)
stands following fire and permafrost degradation [@wirth_white_2008]. Similarly,
moderate disturbances favored the increase of warm-adapted species and led to a
broad-scale community thermophilization of forests in Québec
[@brice_disturbances_2019]. However, the coupling of warming and disturbances
may depend on the intensity and type (natural or anthropogenic) of disturbances
[@johnstone_changes_2010]. Hence, in some cases, disturbances can promote
invasions by early successional species which then displace long-lived
shade-tolerant. For instance, clearcutting has been found to favor the expansion
of trembling aspen (*Populus tremuloides*) in mixed and boreal stands of Québec
[@laquerre_augmentation_2009; @grondin_have_2018] and Alberta
[@landhausser_disturbance_2010]. In contrast, other simulation studies have
found that, although disturbances can influence some types of transitions, they
are unlikely to drive extensive biome shifts in the coming decades
[@vanderwel_how_2014; @liang_how_2018]. As a result, disturbances could
accelerate shifts to an increasingly deciduous‐dominated landscape, but
empirical evidence of this process are necessary to understand how various
intensities and types of disturbances may catalyze or hinder transitions between
certain forest states.


*other ref: @xu_importance_2012; @danneyrolles_stronger_2019; @johnstone_non-equilibrium_2003*


## *P3. Soil*

The northward migration of temperate species might be further constrained by
their capacity to colonize different types of soil [@lafleur_response_2010;
@bennett_plant-soil_2017; @brown_non-climatic_2014]. Indeed, soils of cold
boreal forests generally have lower pH, lower microbial activity and slower
decomposition rate of organic matter than warmer southern temperate forests
[@goldblum_deciduous_2010]. These local and regional variations in soil
properties (e.g., quality of drainage, availability of nutrients, pH,
mycorrhizae) are expected to slow down or inhibit the establishment of temperate
trees into the boreal forest. For instance, transplant experimental studies have
shown that seedlings of sugar maple (*Acer saccharum*) in conifer-dominated
stands were negatively affected by seed predators and fungal pathogens
[@brown_non-climatic_2014] as well as by soil acidity through reduced foliar
nutrition [@collin_conifer_2017]. However, @kellman_sugar_2004 found that, after
initial high mortality due to seed predation, survival of *Acer saccharum*
seedling in boreal stands was high, even superior to that in the temperate
stands, potentially because of increased light availability. Hence, it has been
suggested that soil properties in boreal forests may not be a major impediment
to the migration of temperate species showing broad ecological tolerance
[@lafleur_response_2010; @barras_supply_1998; @kellman_sugar_2004]. Nonetheless,
suboptimal soil conditions could delay forest transition under climate change
[@brown_non-climatic_2014] and,  moreover, waterlogged conditions and deep
organic matter accumulation may prevent the regeneration of several species
[@lafleur_response_2010]. While experimental studies provide valuable insights
on the potential role of soils at local scales, we need to test the generality
of such constraints, or the lack thereof, on long term forest dynamic, across
species and across scales to better anticipate future biome transition.

<!--
@johnstone_changes_2010 "An important insight that emerges from this broad‐scale study is that the potential for fire to drive shifts in successional trajectories is contingent on landscape factors such as site moisture. Consequently, we should expect that the resilience of black spruce forests to changing climate and fire regime will not be uniform across the landscape and that drier spruce forests may have the greatest potential to switch to deciduous‐dominated forests under future environmental change."-->


## *P4. Markov multi-state models*

One approach to investigating the process of biome shifts in response to climate
change is to model transition probabilities of forest plots among states based
upon the knowledge of their current state, as well as their current
environmental characteristics. Given the unequivocal distinction between
temperate and boreal forests, the dynamic of tree communities at ecotone can be
adequately characterized using discrete functional and successional states,
namely Boreal, Mixed, Temperate and Pioneer [@vissault_biogeographie_2016] and
thus can be formalized as a multi-state Markov model
[@jackson_multi-state_2018]. Markov models provide a useful framework for
modeling changes of state over time using longitudinal data. In epidemiology,
for example, the models are often used to describe the progression of diseases
[@van_den_hout_multi-state_2016]. In ecology, the models have been used to study
processes such as ecological succession [@runkle_gap_1981; @hill_markov_2004],
metapopulation dynamics [@hanski_metapopulation_2003; @moilanen_patch_1999],
landcover changes [@yang_land_2014; @muller_markov_1994], or stage class
transitions [@caswell_matrix_2008].

Despite the simplicity of a four-state transition model, this modeling framework
allows to answer questions related to the coarse-scale dynamics of biome shift.
Moreover, without explicitly modeling them, multi-state models represent
mechanisms [@wootton_prediction_2001]. For example, transitions to pioneer
reflect disturbance, transitions from pioneer reflect colonization, dispersal
and recruitment limitation and transitions between the other states reflect
competition. In addition to their simplicity, there exists numerous
well‐established properties of Markov transition matrices [@hill_markov_2004].
Transition matrices can be estimated from the model output and their properties
can then be compared under different scenarios to further explore the mechanisms
of forest dynamics. For instance, the steady state distribution can be derived
from a transition matrix and allow to infer the long-term forest composition
[REF]. Multi-state models can not only be used to explicitly test hypotheses
about equilibrium, but they can also be used to describe state changes during
dynamic periods, i.e. transient dynamics [@boulangeat_transient_2018]. As most
ecosystems never reach equilibrium [@holling_resilience_1973], transient
dynamics, especially in forests, may persist over very long time periods and are
now thus recognized as crucial part of ecosystem dynamics [REF]. The time of
convergence to reach the steady state distribution describes how fast the system
converges; the turnover time indicates how fast the transitions occur and
provides insights about the stability and the resistance of a forest state,
while the entropy reveals the predictability of the transitions. Contrasting
empirically derived transition matrices and their properties among disturbance
scenarios can shed new light on forest dynamics under climate change and may
even provide insights on management measures.




## *P5. Objectives*

In this paper, we investigate the response of forests to recent climate warming
by estimating the transition probability among four community states, boreal,
mixed, temperate and pioneer. Specifically, we address these questions: 1. Is
forest transition dynamic affected by recent climate change? 2. How do
disturbances and soil characteristics influence the transition probability among
forest states under climate change? Can natural or anthropogenic
disturbances catalyze climate related transitions? And, conversely, can soil
characteristics constrain these transitions? 3. How do different disturbance
type and intensity influence the transient dynamics and equilibrium distribution
of forest states?


We expect that the probability of self transitions will be the highest, i.e.
most forests will not change states, because of tree slow demography. However,
climate warming should promote more transitions from boreal to mixed forests and
from mixed to temperate than the reverse. We also anticipate that natural and
anthropogenic disturbances will further favor these climate-related transitions
while soil characteristics will slow them down. We apply a time-continuous
multi-state model to the dynamics of forest communities to estimate state
transition probabilities and evaluate the influence of environmental covariates
on these transitions. Using the result from our multi-state model, we
investigate the impact of disturbances on forest equilibrium and transient
dynamics under recent climate change using several measures: equilibrium state
distribution, time to converge to equilibrium, turnover time and entropy.

<!--"Our ability to predict future ecosystem responses to climate change at time scales relevant to society depends on improving our understanding of the mechanisms, such as disturbance‐dependent migration that may cause species migration rates to lag behind their potential climate limits. "-->


# Methods

## Study area and forest inventory data

To investigate large-scale transition dynamics in forest communities, we used
forest inventory plots in Quebec, Canada, which have been sampled approximately
every ten years since 1970 and ongoing by the Ministère des forêts, de la Faune
et des Parcs [@mffp_placettes-echantillons_2016]. The study area extends from
approximately 45° to 52°N of latitude (ca. 795 000 km^2^) and covers six
bioclimatic domains and three different forest subzones (Fig. 1). The mixed
forest (from 47°N to 48°N) marks the transition between the hardwood forest to
the south, which is dominated by *Acer saccharum*, and the boreal forest to the
north, which is dominated by *Abies balsamea*  and *Picea mariana*.

We selected all inventory plots that were sampled at least twice as well as the
ones where soil covariates were available. We disregarded plots that were
subjected to active reforestation during the study period because we were
interested in transition dynamics resulting from natural recolonization. This
yielded a total of 10,388 plots analyzed (Fig. 1). The time intervals between
plot surveys varied from 4 to 43 years, with a mean time of 11 years ($\sigma$ = 3.85).

![Locations of the 10,388 forest inventory plots in meridional Québec, Canada. Colors delimit the six bioclimatic domains. The two southernmost domains (red) are here combined. The number of forest plots in each domain is written in parentheses. The balsam fir-yellow birch domain is the ecotonal zone between hardwood and boreal forests.](res/fig1_region.pdf)

## Community states

We classified the forest inventory plots into four community states (Boreal,
Mixed, Temperate and Pioneer) using species basal area and composition at each
sampling date. We first assigned each studied species as boreal, temperate or
pioneer according to their functional traits [see Table SX in @brice_disturbances_2019]. For
each plot, we computed the total basal area of each species group and then
classified the plot following the @mffp_placettes-echantillons_2016 definitions
to one of the four states; Boreal (boreal species represent >75% of the plot
basal area), Temperate (temperate species represent >75% of the plot basal
area), Mixed (temperate and boreal species both occupy between >25% and <75% of
the plot basal area) and Pioneer (pioneer species represent >75% of the plot
basal area or plot total basal area <5m^2^/ha). We analyzed state transitions
between each consecutive plot survey. Based on this classification, from the
38,091 observations (plots x number of years measured), we observed 27,703 state
transitions (Fig. 2).


## Environmental variables

The annual past climatic conditions, covering a period from 1960 to 2013, were
extracted from a 2km^2^ (60 arc sec) resolution grid for the entire study area
using the ANUSPLIN climate modelling software
[http://cfs.nrcan.gc.ca/projects/3/8; @mckenney_customized_2011]. Plot locations
were intercepted with two bioclimatic variables hypothesized to influence
tree establishment, survival and growth: the mean temperature during the
growing season and the annual climate moisture index (CMI),
which is the difference between annual precipitation and potential
evapotranspiration (Table 1). To reduce the effect of inter-annual climate
variability, each climate variable was averaged over a 10-year period prior to
the plot measurement. Over the past four decades, growing season temperature
have increased by 0.14 °C/decade, while CMI has decreased by 1.2 cm/decade<!-- (Fig.
SX)-->.

We also collected information pertaining to natural and anthropogenic
disturbances that have affected the forest plots during the study period (Table
1<!--, Fig. Sx-->). At each plot, the type of disturbances (21 types) and their level
of intensity (moderate or major) were recorded [see Table S2 in @brice_disturbances_2019;
@mffp_placettes-echantillons_2016]. The MFFP defined major disturbances as
events that have eliminated more than 75% of the tree basal area, whereas
moderate disturbances have eliminated between 25% and 75%. For our multi-state
model, we differentiated two main types of disturbances: natural disturbances
and harvest, with three levels of intensity each (minor, moderate or major).

Finally, at each plot, several edaphic characteristics were recorded
[@mffp_placettes-echantillons_2016]. Of the available variables, we selected
drainage and pH because they largely affects nutrient availability, soil structural
properties and vegetation development [REF]. These two variables also capture
most of the variance in soil characteristics in plots across Quebec and were
orthogonal in a PCA (not shown).

Table 1. Description of the covariates used in the multi-state models.

|Covariate name   |Covariate description                                       |
|:----------------|:-----------------------------------------------------------|
|**Climate**      |                                                            |
|Temp             |Mean temperature during growing season, 10-year average prior to first measurement (°C). |
|CMI              |Mean annual Climate Moisture Index, 10-year average prior to first measurement (cm). |
|**Soil**         |                                                            |
|pH               |pH of of the surface horizon                                |
|Drainage         |7 classes of soil drainage, which range from excessive to very poor, that were treated as numeric.|
| **Disturbances**|                                                            |
|Logging          |Tree harvesting, including clearcutting, selection cutting, shelterwood cutting, seed-tree cutting, etc. None or minor (0), moderate (1) or major (2). |
|Natural          |Natural disturbances, including forest fires, insect outbreaks, windfall, etc. None (0), moderate (1) or major (2). |


## Analysis

*NB: very detailed... where can we cut?*

*NB2: Lot of mathematical equations with different letters... verify if everything matches!*

### Continuous-time multi-state Markov model

We derived our modeling framework from methods widely used in survival
analysis and disease progression model [@jackson_multi-state_2018;
@van_den_hout_multi-state_2016]. Similarly to @vissault_biogeographie_2016, we
formalized forest dynamics as a four-state model, but here we used a
continuous-time multi-state model [@jackson_multi-state_2018] in which
transitions among states depend on the previous state, time interval, climate,
disturbances and soil characteristics (Fig. 2).

In ecology, Markov models are often built using discrete time steps. However,
because (1) time interval between surveys are irregular, (2) for each time
interval, multiple transitions are possible and (3) the exact times of state
changes are unobserved (i.e. observations are interval-censored), a
continuous-time Markov model, in which time is treated as continuous, is
preferable. Given observations at fixed time intervals, a homogeneous
continuous-time Markov chain is a special case of a discrete-time Markov chain.

For states $r,s \in {B, M, P, T}$ and time $h,t ≥ 0$, transition probabilities
($p_{rs}$) are defined as the probability that a plot in state $r$ at time $h$
is in state $s$ at time $h + t$ and can be denoted by:

$P_{r,s}(h, t) = P(S_{h+t} = s | S_{h} = r)$.

The Markov process is assumed to be time homogeneous, meaning that the
transition probabilities are constant over time (i.e. independent of $t$, but
dependent of the time interval), hence $P(S_{h+t} = s | S_{h} = r) = P(S_{t} = s
| S_{0} = r)$. However, this assumption can be relaxed (see below). In a
four-state transition model, the transition probability matrix $P(t)$ is a 4
$\times$ 4 matrix, where the rows are the current state and the columns the future
state, containing the transition probabilities $p_{rs}(t)$ for a specified time
interval. For a time-homogeneous model, $P(t)$ can be solved by taking the
matrix exponential of the intensity or generator matrix $Q$ scaled by the time
interval:

$P(t) = e^{tQ}$

The intensity matrix $Q$ contains transition intensities $q_{r,s}$ which
represent the instantaneous risk of moving from state $r$ to state $s$:

$q_{r,s} = \lim_{\Delta \to 0} \frac{P(Y_{t+\Delta} = s | Y_t = r)}{\Delta}$, on
off-diagonal elements.

$q_{r,r} = − \sum_{s \neq r}{q_{rs}}$, on diagonal elements.

Transition-specific hazard regression models can be defined for those $r,s \in
S$ between which a direct transition is possible according to the specified
multi-state process (Fig. 2). The intensities $q_{r,s}$ can be modeled as a
combination of a baseline hazard $q_{rs.0}$ with a vector of explanatory
variables $x(t)$ and a vector of coefficients $\beta_{rs}$ :

$q_{rs}(t) = q_{rs}(t|x(t)) = q_{rs.0}(t)exp(\beta_{rs}'x(t))$,

The definition of the log-linear regression hazard model allows to fit a
site-specific and time-dependent covariate vector $x(t)$ to transition
intensities. Time-dependent covariates, such as climate and disturbances, are
assumed to be piecewise-constant, i.e. the hazard is constant within a specified
time interval $[h, h+t]$ and depends on the covariate value at $h$, but is
allowed to change between the intervals. The time homogeneity assumption is thus
relaxed by the inclusion of time-dependent covariates in the model.

Estimation of model parameters can be obtained by maximizing the log-likelihood
function using the transition probability matrix. The contribution of the plot
$i$ at time $j$ to the likelihood is given by:

$LL_i(\theta | s,x) = \prod\limits_{j=1}^{J} P(S_j=s_j|S_{j-1}=s_{j-1},\theta,x)$,

where $\theta$ is the vector with all the model parameters, $x$ denotes the
vector with the covariate values, and $s$ denotes the observed state trajectory
$s_1,...,s_J$ at times $t_1,...,t_J$. The full likelihood function is the
product of contributions for all $N$ plots:

$LL(\theta) = \prod\limits_{i=1}^{N} LL_i(\theta | s,x)$,

### Definition of candidate models

It is important to consider which transitions can realistically occur in
continuous time. Because the states are defined based on the proportion of each
species group, it is assumed that in order for a site to travel from one state
to a non-adjacent state, the plot also has to travel through the intermediate
states. Thus, in this model, we assumed that an instantaneous transition from
Boreal to Temperate and from Temperate to Boreal is impossible (necessary
transition through Mixed), however all states can transition directly to Pioneer
when disturbed (Fig. 2).

We built five different models: one baseline model with intercept only, one for
each subgroup of covariates independently (climate, soil and disturbances) and
one full model that is a combination of all the covariates (Table 1). Because we
estimate multiple state transitions in a single model (all $q_{rs}$ in Fig. 2),
the number of parameters increase rapidly with the number of covariates (number of
modeled transitions (here 10) $\times$ number of covariates). Thus, to reduce the
number of parameters, we hypothesized that transitions from any state to Pioneer
were only determined by disturbances while climate and soil variables should not
directly influence these transitions. All quantitative variables were
standardized ($\mu$ = 0, $\sigma$ = 1) prior to running the models.

![Multi-state transition diagram (a), intensity matrix (b) and equations of our
full model (c). Directional arrows in (a) depict the allowed transitions between
states. The numbers represent the percentage of observed transitions between
states (nb~rs~/nb~r.~ x 100). Instantaneous transition from Boreal to Temperate
and vice versa are considered impossible in the model (hence the absence of
arrows in the diagram and the zeros in the Q matrix), however rare transitions
from Boreal to Temperate and from Temperate to Boreal were observed in the data
(less than 0.2%). All transitions from any states to Pioneer were modeled only
as dependent of disturbances.](res/fig2_trans_diagram.pdf)


### Evaluation of candidate models

We first evaluated the goodness-of-fit of each model containing covariates
(climate, soil, disturbances and full) against the baseline model using
likelihood ratio tests [@jackson_multi-state_2011], which test if the addition
of one or more new parameters significantly increases the likelihood of the
model. The statistical significance of individual covariates in
the presence of the other was determined by comparing the full model with the
correspondingly reduced model using likelihood ratio tests.

### Prediction performance of candidate models

We fitted all models on the full data sets but also used cross-validation to
estimate the predictive performance on held-out data. We used two statistics,
the area under the receiver operating characteristic (ROC) curve (AUC) and the
logarithmic scoring rule (LS), to assess the agreement between the observed
state and the models' predictions. The AUC is a popular performance metric for
binary classifiers that measures the probability of correct ranking of a random
positive-negative pair. The AUC ranges from 0 to 1, where a score of 1 indicates
perfect discrimination, while a score of 0.5 is as good as random.
@hand_simple_2001 has extended the AUC method to multi-class problems, here
denoted mAUC. For any pair of states $r$ and $s$, we can compute $\hat{A}(r|s)$,
the probability that a randomly drawn member of state $s$ has a lower estimated
probability of belonging to state $r$ than a randomly drawn member of state $r$.
We can measure the discrimination rate between all pairs of states by computing
the pairwise AUC:

$$\hat{A}(r,s) = [\hat{A}(r|s)+ \hat{A}(s|r)]/2$$

The mAUC can be obtained by averaging the pairwise AUC:

$$mAUC = \frac{2}{c(c-1)} \sum_{r<s} \hat{A}(r,s)$$


The LS was proposed by Good (1952) and is often used in weather forecasts
[@gneiting_strictly_2007]. While AUC is a function of different classification
thresholds, LS measures the degree to which predicted probabilities are close to
the observed outcomes. We computed a global score for each model:

$$LS = \frac{1}{N} \sum_{i=1}^N -log(P(S_i = s_i))$$

where $S_i$ is the random variable describing the state of the forest in the
$i^{th}$ plot and $s_i$ is the observed state. So, LS only depends upon the
predicted probability of the realized state and not on the probabilities
assigned to the other possible states. The score is very sensitive to incorrect
predictions: if a model predicted the observed state with a probability of 100%,
the score for that plot is 0, while if a probability of zero was assigned to the
observed state, the score goes to infinity. Hence, this sensitivity allows to
emphasize the differences between model predictions and strongly penalize a
model that only gives high probabilities to self transitions.

To assess the quality of prediction for the four states individually, we also
computed LS for each state $r$ where we summed the predicted probability $P(S_i
= r)$ if the observed state is indeed $r$ and $1 - P(S_i = r)$ otherwise.

We evaluated and compared the predictive performance of our five models using the
overall mAUC and the pairwise AUCs, as well as the overall LS and the state-specific
LS. These metrics were estimated using stratified K-fold cross-validation [REF].
We first stratified the data set by bioclimatic domain to ensure that each fold
was representative of the plot geographical distribution and randomly split the
data set in k=10 folds. The cross-validation process was repeated k times,
during which k − 1 folds were used to train the models and the remaining fold
was used to validate the model predictions against the observed state
transitions. The cross-validated performance metrics were then averaged for each
model.

### Effects of covariates on the transition dynamic

Using the result from the best model, we compared the estimated hazard ratios to
investigate the influence of environmental covariates on transition dynamic. We
also computed the predicted 10-year transition probabilities of forest plots
under different disturbance scenarios, while keeping all other covariates at
their average found in the ecotonal zone (i.e. the balsam fir-yellow birch domain),
to facilitate visual interpretation of the impacts of disturbances on the matrix
structure.

To further understand how disturbances can modify the forest dynamic under
climate change, we characterized different properties of the transition dynamic
and compared them among levels and types of disturbances along the latitudinal
temperature gradient. An extensive literature describes the multiple properties
of discrete time Markov transition matrix [@caswell_matrix_2001;
@hill_markov_2004] and these properties can also be measured for continuous time
Markov models. We chose 4 informative and complementary properties that fully
characterized the dynamic of our modeled system: (1) the steady state
distribution, which corresponds to the long-term equilibrium state distribution
of the system; (2) the time of convergence to the steady state; (3) the turnover
time, which measures the speed of transient successional changes; and (4) the
entropy, which captures the predictability of transitions.

First, we estimated the steady state distribution, $\pi$. For a regular
Markov process, any initial state distribution $s(0)$ converges to the same
equilibrium as $t$ tends toward infinity:

$\displaystyle \lim_{t \to \infty} s(0)P(t)=\pi$

The vector of equilibrium $\pi$ can be obtained by taking the left
eigenvector of $Q$ with eigenvalue 0, normalized to sum to 1 or by taking
the dominant eigenvector of $P$ with eigenvalue 1, normalized to sum to 1
[REF].

Then, the convergence rate to the equilibrium distribution can be measured
using the damping ratio [@hill_markov_2004]:

$\rho = \lambda_{P1} / \lambda_{P2}$ or $\rho = exp(\lambda_{Q1} - \lambda_{Q2})$,

where $\lambda_{P1}$ and $\lambda_{P2}$ are the largest and second-largest
eigenvalues of $P$ ($\lambda_{P1}$ = 1 for stochastic $P$) and $\lambda_{Q1}$
and $\lambda_{Q2}$ are the largest and second-largest eigenvalues of $Q$
($\lambda_{Q1}$ = 0 for stochastic $Q$). The half-life to equilibrium is given
by:

$t_{1/2} = log(2)/log(\rho)$

We also measured the turnover time in each forest state, also called the mean
sojourn time in multi-state models, which corresponds to the time spent in one
state before transitioning to the next. The turnover time can be estimated by
$Turnover_{r} =−1/q_{rr}$, where $q_{rr}$ is the r^{th} entry on the
diagonal of the estimated generator matrix. The turnover of the whole system is
given by the average of each state turnover time over the steady state
distribution:

$Turnover = \displaystyle -\sum_r{\pi_r \times Turnover_{r}}$

Finally, @hill_markov_2004 suggested using the entropy of a discrete-time
transition matrix as an index of the predictability of successional changes. It
measures how uncertain we are about the next new state of a site knowing its
current state. For a continuous-time process, the entropy can be measured using
the jump matrix [@spencer_continuous-time_2005]. The jump matrix contains the
probabilities that the next state after state $r$ is state $s$:

$j_{rs} = −q_{rs}/q_{rr}$,

where $qrs$ is the transition intensity from $r$ to $s$. The entropy of state
$s$ is then:

$H(j_{.s}) = \displaystyle -\sum_r{j_{rs} \times log(j_{rs})}$

The normalized entropy of the whole system is the average of the entropies
over the steady state, divided by $H_{max} = log(n_{state} = 4)$:

$\text{Entropy} = \displaystyle \frac{-\sum_r{\pi_r \times H(j_{.s})}}{H_{max}}$

Values of entropy closer to zero indicated a more deterministic transition
dynamic and values closer to one indicated a more random dynamic.

All analyses were performed using the R programming language version 3.5.1
[@r_core_team_r_2018]. The list of R packages that have been used throughout the
analysis is provided in the Supporting Information (Table S1). All the data used
in the study, in addition to R scripts to reproduce the analyses and the
figures, will be made available online on Github unpon manuscript acceptance.


# Results

*Note: standardize verb tense + add letters to panel figures*

## Transition dynamics

In this study, there are a total of 38,091 surveys that were recorded for the
10,388 forest plots. At the beginning of the surveys, there were 4567 forest
plots assigned as Boreal, 1143 as Mixed, 2712 as Pioneer, and 1966 as Temperate.
At the end, there were 4795 Boreal plots, 1106 Mixed plots, 2176 Pioneer plots
and 2311 Temperate plots (Fig. 2; Table S2). A large fraction of Mixed forests
transitioned to Temperate forests (21.2%) but few did the opposite (5.8%). There
were many transitions between Boreal and Pioneer, but more Pioneer recovered to
Boreal than the reverse. Temperate and Boreal forests were generally more stable
(large fraction of forests did not transition) than Mixed and Pioneer forests
(Fig. 2).

## Model evaluation

Overall, the full model including climate, soil and disturbance variables had
the best fit and predictive performance (Fig. 3; Table 2). All variable subsets
improved significantly the likelihood of the model (all likelihood ratio tests
were highly significant, p << 0.001; Table 2). Model evaluation using 10-fold
cross-validation revealed that including climate and disturbances improved
overall model predictive performance (mAUC and LS), while soil variables had a
negligible effect (Fig. 3). All models were good at distinguishing Boreal from
Temperate (high pairwise AUC). Soil variables slightly help to predict Mixed and
Temperate states. Including climate variables help to distinguish Mixed from the
other states, while including disturbances help to distinguish Pioneer from the
other states, especially Boreal. Hereafter, all inferences about transition
probability parameters were derived from the full model.

![Mean pairwise AUC (areas under the receiver operating characteristic curves) obtained through 10-fold cross-validation. Higher values indicate a better capacity to discriminate between the four forest states: (B)oreal, (M)ixed, (P)ioneer and (T)emperate. The overall mAUC of each model is given next to the legend.](res/fig3_cv_auc.pdf)

**AND/OR**

![Mean state-specific logarithmic skill score where each model including covariates is compared to the baseline model. Values were obtained through 10-fold cross-validation. Higher values indicate a larger improvement (predicted probabilities are closer to the observed outcomes) compared to the baseline model. The overall logarithmic score of each model is given next to the legend.](res/fig3_cv_LSS.pdf)

**OR**

![Mean state-specific logarithmic score obtained through 10-fold cross-validation. Lower values indicate that the predicted probabilities are closer to the observed outcomes. The overall logarithmic score of each model is given next to the legend.](res/fig3_cv_LS.pdf)

Table 2. Comparisons of the five candidate multi-state models. The number of
parameters used in each model corresponds to the number of modeled transitions
(10) $\times$ the number of covariates. The $\Delta$AIC is the difference
between the Akaike information criterion of each model (AIC~m~) and the minimum
of AIC among all the models (AIC~min~): $\Delta$AIC = AIC~m~ – AIC~min~.
Multi-class area under the curve (mAUC) and logarithmic score were obtained
through 10-fold cross-validation. Higher mAUC and lower LS indicate better model
predictive performance. The best model is the one in bold with $\Delta$AIC = 0.

|             |Covariates        | Nb of parameters| -2 Log-likelihood|$\Delta$AIC|LR test|  mAUC|     LS|
|:------------|:-----------------|----------------:|-----------------:|---------:|-------:|-----:|-------:|
|Baseline     |Intercept         |               10|           32032.3|    6132.6|< 0.001 | 0.899|   0.578|
|Climate      |Temperature, CMI  |               24|           30438.7|    4566.9|< 0.001 | 0.921|   0.550|
|Soil         |Drainage, pH      |               24|           31886.9|    6015.2|< 0.001 | 0.906|   0.576|
|Disturbances |Natural, Logging  |               50|           27341.3|    1521.6|< 0.001 | 0.925|   0.495|
|**Full**     |**All**           |           **78**|       **25763.7**|   **0.0**|< 0.001 |**0.940**|**0.468**|


\pagebreak

## Effect of covariates on transition intensity

The full multi-state model allows to reveal interesting relationships between
probabilities of forest state transition and environmental covariates. All
transitions to Pioneer were highly influenced by disturbances (Fig. 4). As
expected, major disturbances exert stronger effects than moderate disturbances
(for both natural and logging), but logging had stronger effects for both levels
of intensity. For example, the risk of transition from Boreal to Pioneer has
surged up to 165 times higher for plots that suffered major logging (logging 2)
compared to that which were not disturbed (minor). Disturbances of all types and
intensities favored transitions from Mixed to Temperate forests. Major
disturbances increase the risk of transition from Mixed to Temperate by ca. 5
times (Hazard Ratio (HR) = 4.51 and 5.32, for natural and logging,
respectively). Moderate disturbances (natural and logging) also favored
transitions from Boreal to Mixed (HR = 2.80 and 2.78, respectively), while major
disturbances had no significant effect.

Climate variables also had a significant influence on most transitions (Fig. 4).
Warmer summer temperature (higher sTP) and higher humidity (higher CMI) favored
the transitions from Boreal to Mixed as well as from Pioneer to Mixed and
Pioneer to Temperate. Interestingly, warmer temperature did not significantly
influence the risk of transition from Mixed to Temperate and higher CMI had a
negative effect.

Although less important (Fig. 3), soil variables also impacted state transitions
(Fig. 4). Holding the other covariates constant, poorer drainage (more humid)
decreased the instantaneous risk of transition from Boreal to Mixed by 29% and
from Pioneer to Temperate by 17%, but increased the risk of transition from
Temperate to Mixed by 32% (HR = 0.71, 0.83 and 1.32, respectively). Higher pH
(acidic soil) had a considerable negative effect on the transitions from
Temperate to Mixed and a smaller one on transitions from Pioneer to Boreal (HR =
0.80 and 0.93, respectively). These changes in risk ratio associated to soil
variables appear almost irrelevant compare to the effect of disturbances, but a
slight increase in drainage can dampen the positive effect of disturbances.
For instance, under moderate natural disturbances, the estimated risk of
transition from Boreal to Mixed is 2.12 (1.04-4.33) at moderate drainage but
decreases to 1.50 (0.71-3.20) when increasing drainage by 1 point.


![Hazard ratios (HR) and 95% confidence intervals as estimated from the best multi-state transition model. The y-axis is in log scale. The HR of predictors are interpretable as multiplicative effects on the hazard, where values above 1 (in blue) indicate that the predictor is associated with a greater risk of state transition, while values below 1 (in red) indicate a lower risk of transition. Predictors different from 1 are colored in dark blue or red. Numbers following disturbance predictors indicate their levels of intensity: 1 = moderate and 2 = major. ](res/fig4_HR.pdf)


\pagebreak

Disturbances completely altered the structure of the 10-year transition
probability matrix (Fig. 5). The largest values across most matrices were
generally associated with self transitions (matrix diagonal), meaning that the
vast majority of forest plots (characterized by the average environmental
conditions of ecotone) remain in the same state after 10 years. For undisturbed
forest plots (minor), the self transitions are very strong but transitions from
Pioneer to Boreal, from Mixed to Temperate, and from Temperate to Mixed were
not trivial. At moderate disturbances, probabilities of self transitions
decrease, while transitions from Boreal to Pioneer, from Mixed to Temperate
increase the most. Transitions from Mixed and Temperate to Pioneer do not
increase much at moderate disturbances, likely because such disturbances were
less frequent and less severe than in Boreal forests. The difference between
natural disturbances and logging emerges only at major disturbances. For both
types of major disturbances, the probabilities in the third column, transitions
to Pioneer, showed a great increase compare to moderate disturbances, but these
values exploded in severely logged transition matrix, exceeding self
transitions. Interestingly, the estimated probability of Mixed to Temperate
remain quite high at major disturbances.


![Predicted change in 10-year transition probabilities for different disturbance types and levels. All other covariates are fixed at the average conditions found in the ecotone, i.e. the balsam fir-yellow birch domain. Letters correspond to the four forest states: (B)oreal, (M)ixed, (P)ioneer and (T)emperate. Numbers are the modeled transition probabilities from rows to columns and darker color highlights stronger transitions.](res/fig5_pmatrix.pdf)

<!--
- Moderate disturbances strongly increase the probability of transition from Boreal to Pioneer, and even more so with major disturbances. This effect seems long-lasting.

- Strong effect of major harvesting on the probability of transition from Temperate to Pioneer. Weaker for major natural disturbances.

- Without disturbances (Minor), the transition probability from Mixed to Temperate increase constantly through time.

- Moderate disturbances further increase the transition probability from Mixed to Temperate through time.

- However, major disturbances increase the transition probability from Mixed to Temperate comparatively to moderate on the short term, but on the long term it decreases the transition probability. This effect is weak for major natural disturbances, but very pronounced for major harvesting.
-->


\pagebreak

## Effect of disturbances on long term equilibrium

The steady state forest dominance changes as expected along the temperature
gradient (Fig. 6); the Boreal state dominates at low temperature (high latitude)
and Temperate and Mixed states dominate at high temperature (low latitude), with
a transition boundary located at a growing season temperature of about 12.75°C,
which falls in the ecotonal balsam fir-yellow birch domain. When moderate
disturbances were included, the steady state proportion of Temperate and Mixed
forests increases and the Boreal-Temperate boundary occurred at lower
temperatures, hence further north in the balsam fir-white birch domain (Fig. 6).
Moderate natural disturbances and logging had a similar positive effect on the
proportion of Temperate and Mixed. For example, at 12.58°C, northern end of
balsam fir-yellow birch domain, the steady state proportion of Temperate and
Mixed almost doubles with moderate disturbances (minor: 36%; moderate natural:
58%; moderate logging: 60%). However, because there was a larger reduction in
the proportion of Boreal to the benefit of Pioneer, the displacement of the
boundary was slightly larger for logging than natural disturbances (displaced at
12.05°C for logging and at 12.18°C for natural disturbances; Fig. 6). Finally,
for major disturbances, the proportion of Boreal forests at steady state
collapsed while that of Temperate and Mixed forests also decreased to a lesser
extent. The landscape is now dominated by Pioneer forests. The boundary modestly
moved north with major natural disturbances (12.60°C), while it retreated to the
south with major logging (13.05°C). Simulating change in frequency of
disturbances (instead of all or nothing scenarios) shows that the
Boreal-Temperate boundary slowly moved north with increasing moderate
disturbance frequency (natural or logging), but it stagnates with increasing
frequency of major natural disturbances and recedes south with increasing
frequency of major logging (Fig. S3).

![Changes in forest state proportion at equilibrium along the temperature (latitudinal) gradient for different disturbance scenarios. Proportion of Boreal (blue) and Temperate + Mixed forests (red) for minor (solid), moderate (dashed) and major (dotted) disturbances. All other covariates are fixed at the average conditions found in the ecotone, i.e. the balsam fir-yellow birch domain, to focus solely on the effect of disturbances along the temperature gradient. The white (minor), grey (moderate) and black (major) circles indicate the position of the boundary between dominance of Boreal forests and dominance of Temperate and Mixed forests (i.e. the advancing front). The colors at the top of the plots approximate the position of the bioclimatic domains along the temperature gradient.](res/fig6_SS_gradient.pdf)

## Effect of disturbances on transient dynamics

Natural disturbances and logging affected forest transient dynamics, with greater
impacts for higher disturbance intensity (Fig. 7). In the minor disturbance
scenario, turnover time was generally longer at low temperature, indicating a
slow transition dynamic in forests of northern latitudes (Fig. 7a,b). The
turnover time then rapidly declines to reach a minimum at 13.20°C, between the
sugar maple-yellow birch and the balsam fir-yellow birch domains, and goes back
up after this point. This trough, where transition dynamics are the fastest, is
located just a little south of the point of transition between Boreal and
Temperate dominance found in Figure 6. Major disturbances accelerate transition
dynamics all along the temperature gradient, while moderate disturbances also
decrease turnover time but more strongly in the boreal domains (balsam fir-white
birch and spruce-moss domains; Fig. 7a,b). These spatial patterns reflect the
turnover time of the dominant state at each point along the temperature gradient
(Fig. S4).

At minor disturbances, the entropy of the system generally increases from north
to south and peaked at 12.56°C, at the southern end of the balsam white-birch
domain (Fig. 7c,d). This peak illustrates where the transition dynamic is most
uncertain (transition to all states are possible at this point), while it is
very predictable in northern boreal forests (Boreal stays Boreal until it
transitions to Pioneer later on). The peak can be mainly attributed to the
entropy of the Boreal state and the generally high values at high temperature
can be principally attributed to the Temperate state (Fig. S5). This
latitudinal pattern of entropy is modified by disturbances. Moderate natural
disturbances decrease the entropy throughout the gradient, but especially where
was the peak (Fig. 7c). With moderate logging, the peak disappears and entropy
increases monotonically from north to south (Fig. 7d). When major disturbances
are included, wether natural or logging, the peak of entropy is displaced to the
south (Fig. 7c,d) where it is dominated by the entropy of the Pioneer state
(Fig. S5).

The half-life to equilibrium is the longest at 12.24°C, in the balsam fir-white
birch domain, while it is the fastest in the southernmost latitudes (Fig. 7e,f).
Interestingly, the peak for half-life closely matches the peak for entropy, but
matches no feature of the turnover time curve.  Moderate disturbances flatten
and shift this peak to the north and the effect of moderate logging (Fig. 7f) is
stronger than natural disturbances (Fig. 7e). Hence, in the balsam fir-white
birch, the half-life to reach equilibrium distribution is reduced almost by half
by moderate logging. With major disturbances, forests all along the temperature
gradient can reach very quickly their steady state distribution (maximum of
about 8 years for major logging and 25 years for major natural disturbances).

![Changes in the characteristics of the forest transient dynamics along the temperature (latitudinal) gradient for different disturbance scenarios: minor (solid), moderate (dashed) and major (dotted) disturbances for both natural (a,c,e) and logging (b,d,f). All other covariates are fixed at the average conditions found in the ecotone, i.e. the balsam fir-yellow birch domain, to focus solely on the effect of disturbances along the temperature gradient. The turnover of the whole system (i.e. whole transition matrix) (a,b) corresponds the time spent in a state before transitioning to the next and is given by the average of each state turnover time over the steady state distribution. The entropy of the whole system (c,d) corresponds to the incertitude of the next transition and is given by the average of each state entropy over the steady state distribution. The half-life to equilibrium (e,f) is the time taken to reach 50% of the steady state distribution, i.e. when the first eigenvalue to become twice as great as the contribution of the second eigenvalue. The colors at the top of the plots approximate the position of the bioclimatic domains.](res/fig7_index_gradient.pdf)


# Discussion

## P1. Main results

Our study reveals that disturbances are likely to accelerate forest response to
climate change by promoting transitions from boreal to mixed forests and,
particularly, from mixed to temperate forests. Our analysis of the equilibrium
further highlights that the long term forest dynamics under moderate
disturbances favors an increase proportion of temperate forests and a northward
shift of the boreal-temperate ecotone. Disturbances also modified the forest
transient dynamics, accelerating both the turnover and convergence time and
making the dynamics more predictable. In accordance with the hypothesis
formulated by previous studies [@johnstone_changes_2010;
@johnstone_changing_2016; @vissault_biogeographie_2016;
@brice_disturbances_2019], our findings demonstrate that moderate disturbances
catalyze transition to alternate, temperate-dominated forest state and, as a
result, promote regime shifts.



## P2. Is forest transition dynamic affected by recent climate change?

Our results seem to support the hypothesis that recent climate warming
influences transition dynamics among forest states at the boreal-temperate
ecotone. The higher number of transitions from mixed to temperate is consistent
with the expectation of a northward shift in range of temperate trees species
into the mixed and boreal forests. Indeed, the warming trends of the last
decades (fig supp?) have been shown to improve growth and reproductive rates of
temperate species, but to reduce growth of boreal species
[@reich_geographic_2015; @fisichelli_temperate_2014;
@boisvertmarsh_divergent_2019; @goldblum_tree_2005], thus providing a
competitive advantage of temperate over boreal species.  

Alternatively, the increased transition to temperate species may be a response
to other pressures. Comparisons of pre-settlement and present-day forested
landscape of North America have demonstrated an important deciduous encroachment
in response to historical human activities [@danneyrolles_stronger_2019;
@terrail_reorganization_2019; @boucher_logging-induced_2006]. Similarly,
@johnstone_non-equilibrium_2003 suggests that the northern range expansion of
lodgepole pine following fire is not related to current climate change, but was
rather a continued migration initiated in the early Holocene. However,
historical legacies and climate change are likely mutually non-exclusive
explanations. Simulations by @boulanger_climate_2019 showed that the future
climate-induced expansion in temperate species to the detriment of boreal
species would amplify the already ongoing trend since preindustrial times.


<!--Previous studies have shown that successional pathways were strongly affected by disturbances and climate change [].
 Forest transition dynamics are slow. Even if our model overestimates transition, the changes are still very slow. Once again, migration lags.
Nevertheless, forest dynamics appear to be responsive to CC.-->

## P3. Environmental catalysts and inhibitors of forest state transition

Our study demonstrated that moderate disturbances favor climate-related
transitions, while major disturbances merely promote pioneer states.
Disturbances directly remove trees which lead to immediate and substantial
changes in forest composition and successional trajectories. Without climate
change, forests are expected to be resilient to normally experienced
disturbances and should thus return to their previous states. However, climate
change alters the conditions that initially supported the persistence of a
forest state, making them more vulnerable to disturbances
[@johnstone_changing_2016]. Hence, moderate disturbances remove resident species
and reduce competition for light and nutrient, which likely facilitate
colonization and establishment by opportunistic temperate species under warmer
conditions [@brice_disturbances_2019; @leithead_northward_2010;
@landhausser_disturbance_2010]. In contrast, severe disturbances in the study
area, primarily clearcutting but also large fires, create openings of very large
extent which likely favor early-successional species that can disperse seed over
long distances, such as *Populus sp* and *Betula sp*
[@landhausser_disturbance_2010].


Compared to the catalyzing effect of disturbances, soil characteristics do not
appear as a large impediment to state transition but have the potential to slow
down transitions. Poor drainage constrained climate-related transitions, from
Boreal to Mixed states, but not from Mixed to Temperate. This indicates that
temperate species can readily colonize soils found in mixedwoods, but may have
more difficulty to colonize hydric boreal soils. Very poor drainage, often
associated with peatland and thick organic layer, is indeed considered as one of
the only major edaphic limit for the regeneration of temperate species
[@lafleur_response_2010]. Numerous studies found that *Acer saccharum*
regenerate well across the ecotone because of its large tolerance to various
soil conditions [@barras_supply_1998; @goldblum_age_2002; @kellman_sugar_2004;
@fisichelli_temperate_2014; @collin_can_2018]. At their northern range limit,
*A. saccharum* and *A. rubrum*, the species contributing most to compositional
changes [@brice_disturbances_2019], are hypothesized to be mostly limited by
soil temperature [@barras_supply_1998; @goldblum_age_2002].

Moreover, disturbances may counteract any effect of soil properties. Indeed,
disturbances, such as logging and fire, often remove the surface organic layers
and expose mineral soil and can, consequently, provide an appropriate seedbed
for temperate species recruitment [@archambault_fifty_2006;
@landhausser_disturbance_2010]. In combination with climate warming,
disturbances may also facilitate temperate migration by increasing understory
air and soil temperatures [@stevens_forest_2015; @de_frenne_microclimate_2013].

<!--
Our model also predicts faster transition rates from Mixed to Temperate than from Boreal to Mixed. The inertia of boreal forests was previously observed [@vissault_biogeographie_2016].

- Why M-T > B-M? Dispersal, soil

- Major disturbances only promoted invasions by pioneer species.
  - Contrast with @johnstone_changes_2010 where increasing fire intensity favor northern range expansion of lodgepole pine.

@brown_non-climatic_2014 found that seed predation might be more important than soil quality for the regeneration of *A. saccharum* along an altitudinal gradient.

cold air drainage is preventing sugar maples from expanding their range downslope


-->



## P4. Equilibrium and potential range limits

Our model highlights the potential role of disturbances in controlling the
position of the boreal-temperate boundary as well as the proportion of temperate
and boreal biomes at equilibrium. As a result of the increased replacement of
Mixed by Temperate states and a decline of Boreal to Pioneer states, the
equilibrium boreal-temperate boundary shifts northward with moderate
disturbances. While our results should not be interpreted as predictions of the
future, they are useful to highlight the direction of the forest dynamic. Our
results support the simulations of @boulanger_climate_2019 where harvesting
under future climate warming was projected to promote further invasions of
pioneer species, such as *Populus*, and temperate species, such as *Acer* and
*Fagus*, in mixedwoods of Québec. In contrast, based on their simulations,
@liang_how_2018 and @vanderwel_how_2014 concluded that logging would accelerate
the expansion of pioneer forests, but have little or no effect on extensive
biome shifts over the next century in eastern United States. These apparently
conflicting results could be due to the contrasting tree species responses to
disturbance. Disturbances may facilitate the range expansion of some species but
hinder that of others depending on their functional traits
[@matthews_disturbance_2013; @aubin_traits_2016]. For instance, because of its
positive response to past [@danneyrolles_stronger_2019], recent
[@brice_disturbances_2019] and future [@boulanger_climate_2019] disturbances in
Québec, *Acer rubrum* is likely to play a disproportional role in the temperate
biome shift.

- *Can we really compare?
- in all disturbance scenarios, boreal forests are losing ground primarily to pioneer forests

- Similar to @liang_how_2018 and @vanderwel_how_2014 We also predict that major disturbances will only promote pioneer states to the detriment of boreal states...

Other results:
Simulations by Vieira showed that plantation and enrichment planting shifted northward the boreal-temperate range limits, but not harvesting. Thinning increased the transition from mixed to temperate stands, but did not have any effect on range limits shift.

@boulangeat_transient_2018 even predicted
a southern shift of the boreal-temperate boundary when including the effects of
herbivores.*



## P5. Transient dynamics

*Help! Not sure how to interpret these results!*

Beyond their impacts on steady state distribution, our results suggest that
disturbances may have a substantial influence on forest transient dynamics. In
the continuous boreal zone (spruce-moss domain), forests dominated by *Picea
mariana* are usually characterized by a dynamic of stand self replacement with
minimal compositional changes across disturbance cycles
[@goldblum_deciduous_2010]. Consistent with this dynamic, the turnover time of
undisturbed northern boreal forests is very long and the entropy very low. The
turnover becomes very rapid with disturbances, but the entropy remains low
indicating that the dynamic is still as predictable (back and forth transitions
between boreal and pioneer states) and that there is no shift dynamic associated
with disturbances. Hence, while resistant at low disturbances, boreal forests
loose their resistance when moderately disturbed but remain resilient as they
return to their previous boreal state. Under major disturbances, boreal forests
collapse to pioneer state and reach this new equilibrium swiftly (short
half-life). This interpretation agrees with previous studies suggesting that
boreal forests can easily shift into an alternative treeless state in response
to severe or repeated disturbances [@payette_shift_2003;
@sanchez-pinillos_resistance_2019].

In contrast, the ecotone is characterized by a rapid turnover and a high entropy
indicating abrupt compositional shift which can go in all directions, hence this
zone is neither resistant nor resilient. Compared to northern boreal forests,
the short turnover time implies a low persistence, hence a low resistance, of
the forest states in this region even under minor disturbances. This result
corroborate the predictions made by @vissault_biogeographie_2016, where mixed
forests undergo a swift conversion to temperate forests in the next decades
whereas boreal forests present a large amount of inertia. The dynamic of the
ecotone appear unstable because it is caught between the two stable states, i.e.
boreal to the north and temperate to the south. Under moderate disturbances, the
probability of transitioning to Temperate increases to the detriment of the
other possible states (Fig. 6), hence the entropy is decreased and the dynamic
becomes more predictable. Such a clear directional shift strongly indicates a
non‐equilibrium dynamic in this region. Although turnover is fast, half-life to
equilibrium is long because a forest may not move in the right direction and may
undergo multiple transitions.


- *implication for resilience and resistance
  - long turnover = high persistence = high resistance?
  - high entropy = low stability/resilience?

- Compare to @boulangeat_transient_2018 + Vieira
  -In contrast to our results, @boulangeat_transient_2018 showed that disturbances
  by browsers "reduced the asymptotic resilience of the system, decreasing the
  rate at which the new equilibrium was approached." Effect of browsers are likely less dramatic than tree logging (here moderate disturbances remove at least 25% of the basal area).

- active migration zone more sensitive

- at this scale there appear to be no difference between natural and anthropogenic disturbances.*






## P6. Limitations of the study

<!--Although we tested reasonable parameters, we did not attempt to include every ecological process affecting the forests of boreal-temperate ecotone. Therefore, our simulations are not predictions, rather plausible futures where we can explore the effects of alternative assumptions.-->

Our goal was not to make predictions about the future state of Québec forests,
but rather to explore how disturbances and soils may interact with recent
climate change and modify transition dynamics. For this reason, it is
inappropriate to regard the state distribution at equilibrium or the half-life
to convergence as specific predictions about the extent and timing of biome
shifts. Indeed, we did not simulate future climate change, we have only
considered recent climate change that was observed during the surveyed period.
Moreover, we did not include any dispersal limitations, although it is known to
affect tree migration [@pearson_climate_2006]. Because detailed temporal data on
tree composition around each plot is not available, neighborhood composition
should have been approximated using climate data, but we considered that it was
redundant as they are already used in our multi-state model. Dispersal
limitations could potentially reduce the effect of disturbances on long term
forest distribution, but should not influence our empirically derived results.

Finally, the definition of states can affect to some extent the results. A
higher threshold to define the boreal and temperate states (e.g., >90% instead
of >75% of dominance of boreal and temperate, respectively) would change the
transition probability but the direction of the dynamic would remain the same.


## P7. Ecological and management implications

<!--natural and anthropogenic disturbances, combine with the
underlying effect of climate change, will likely modify forest dynamics by providing
establishment opportunities for some migrating species in otherwise competitive
environments, thus triggering state shifts in forest ecosystems.

Inertia, persistence –Ability of a living system to survive moderate disturbances

Resilience –Ability of a living system to be restored through secondary succession after a moderate disturbance

Tipping point –Any additional stress can cause the system to change in an abrupt
-->

The relationships that we demonstrate between forest state transitions and
disturbances provide a strong empirical basis for predicting the types of
changes in forest dynamics that are likely to unfold in the coming century. A
shift in dominant forest cover from conifer to deciduous broadleaf species not
only entail changes in tree species diversity and composition, but a complete
transformation of forest dynamics and functions. In the long term, this regime
shift could locally increase tree diversity and productivity [REF], increase
carbon sequestration [@thurner_carbon_2014; as long as mortality is limited
@nrdc_pandoras_2018], modify disturbance regime [reduced flammability of
broadleaf species, @terrier_potential_2013 and reduced sensitivity to current
outbreak-prone pest, @mffp_insectes_2018], alter soil microbial activities [e.g.
@laganiere_how_2010] and affect wildlife distribution [e.g.
@mizel_rapidly_2016]. However, such regime shift also have large repercussions
on forest management strategies in area where silvicultural practices are
tailored to regional disturbance regimes and rely on natural regeneration.

Our study reveals the potential of moderate disturbances to facilitate
climate-related transitions and thus suggests that alternative silvicultural
practices could be used to reduced tree migration lags. But even well planned
logging may not benefit all species equally [@brice_disturbances_2019] and may
interact will other natural disturbances to exacerbate tree mortality and
compromise forest resilience in the long run [@buma_disturbance_2011]. In
addition to alternative silvicultural strategies, some studies suggest that
temperate species could be planted farther north to speed up succession
[@iverson_tree-species_2013; @duveneck_measuring_2016; Vieira et al.]. However,
whether promoting temperate tree migration is "desirable" or not depends on the
decisions we make regarding forest management. Do we try to maintain historical
conditions, let nature takes its course, or actively move species to
climatically suitable locations outside their current ranges
[@frelich_will_2010]? In Québec, ecosystem-based forest management seek to
maintain the composition and structure of a reference state, define by the
preindustrial forest conditions [@pinna_amenagement_2009]. Yet,
@boulanger_climate_2019 showed that such management would fail to restore or
approach historical forest conditions under future climate change. In the
context of ongoing climate change, our study reinforces that forest management
have to consider the present system state in relationship to its transient
dynamic as well as its likely trajectory. But even then, simulations by
@duveneck_measuring_2016 suggest that alternative management strategies,
including modified silviculture and climate suitable planting, will have limited
ability to increase resilience and resistance of forest under climate change.
Therefore, in order to insure long term forest resilience at the
boreal-temperate ecotone, adaptative management does not appear sufficient and
drastic reduction of greenhouse gas emission is necessary to limit global
warming [@ipcc_climate_2014].





<!--
However, important questions still need to be answered, such as: How will forest respond to multiple interacting disturbances? And how will these rapid ecosystem changes feedback to ecosystem processes and global cycles? [@turner_disturbance_2010].

Moreover, the fact that, once disturbed, many forest are transitioning to new states suggest the existence of alternative stable state. Management practices might not be adapted to transient dynamics.

- But disturbance-related changes may exacerbate tree mortality and compromise forest resilience in the long run + not all species will benefit from disturbances...
the inherent differences be-tween fire and logging can alter the composition of plant communities
The climate of the future will likely lead to higher mortality among mature trees, because of the greater frequency of droughts, fires, forest-leveling windstorms, and outbreaks of native and exotic insect pests and diseases.

- Stop climate change! The global scientific community [@ipcc_climate_2014], as well as Canadian and other national governments, recognize that climate change is driven by human activity.
Global climate change will alter natural disturbance regimes because many disturbances have a significant climate forcing. What will happen when disturbance regimes change?
-->

<!--

It remains unclear if such transition is gradual and abrupt
[@clements_indicators_2018; @johnstone_changing_2016].
forests may gradually
transition or abruptly jump from one stable state to another
[@clements_indicators_2018; @johnstone_changing_2016]


In some cases, climate change may push forests past critical thresholds such that, upon perturbation, they undergo drastic changes in community composition and ecosystem properties (‘catastrophic shift’) and fail to return to their previous state.

Moreover, as species move, they alter themselves the forest dynamics
and processes and, through a positive feedback loop [@bergeron_fire_2004],
forests may not be able to return to their previous state.

disturbance that fall outside the range of normal background levels.

As a result, conditions that support the persistence of mature forests may not be amenable to forest regeneration.
-->

## References
