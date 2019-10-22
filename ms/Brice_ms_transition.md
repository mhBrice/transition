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

# Abstract [max 300]

Several temperate tree species are expected to migrate northward and colonise
boreal forests in response to climate change. Tree migrations could lead to
transitions in forest types, but these could be hampered by unsuitable soil
conditions and competition by resident species.

Here, we model the state transition dynamics of Quebec's forests in
recent decades to identify the environmental conditions that promote or prevent
these transitions. We also investigate how different disturbance types and
intensities impact the potential long-term equilibrium as well as the
short-term transient dynamics.

We analysed over 10,000 forest inventory plots, sampled from 1970 to 2016 in
meridional Québec, Canada. We used a continuous time Markov multi-state model to
quantify the transition probabilities between forest states (temperate, boreal,
mixed, pioneer) in relation to climate, soil conditions and disturbances.
We described the equilibrium and transient dynamics under different
disturbance scenarios, using complementary properties of Markov transition
matrices.

Although the majority of forests persist in the same state, most transitions
were conversions from mixed to temperate stands, as well as back and forth
between boreal and pioneer forests. Transition probabilities were mainly driven
by natural and anthropogenic disturbances and secondarily by climate, whereas
soils exerted only minor constraints. Moderate disturbances increased the
probability of transition from boreal to mixed and from mixed to temperate
forests. Boreal forests were characterised by a great inertia and predictable
dynamics of stand self-replacement with no climate-induced changes when
disturbed. In contrast, mixed forests at the ecotone presented rapid and
unpredictable transitions at low disturbances, but a clear directional shift
toward temperate forests when disturbed. This led to an increased landscape
proportion of temperate forests and a northward shift of the boreal-temperate
boundary at equilibrium. Hence, moderate disturbances could catalyse rapid
forest transitions and accelerate biome shifts by reducing competition and
providing establishment opportunities for temperate species.


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

Global climate warming is forcing species to move [@parmesan_globally_2003].
Several temperate deciduous tree species at the temperate-boreal forest ecotone
are slowly migrating northward, colonising conifer dominated forests
[@fisichelli_temperate_2014; @evans_borealtemperate_2017;
@boisvert-marsh_shifting_2014; @sittaro_tree_2017]. As climate warms up and tips
the balance in favour of temperate over boreal species, forests are expected to
transition from coniferous to mixedwood and from mixedwood to temperate
deciduous [@boulanger_climate_2019; @price_anticipating_2013]. While boreal
forest dynamics are characterised by broad-scale disturbances, mainly fires and
insect outbreaks, and slow decomposition of an acidic and nutrient-poor litter,
temperate forest dynamics are characterised by small-scale canopy gaps and rapid
decomposition of a rich litter [@goldblum_deciduous_2010]. Hence, as ecological
processes strongly differ among these biomes, climate-induced range shifts not
only impact species distributions, but also alter the structure of communities,
microclimates, biogeochemical cycles and ecosystem functioning and might trigger
a "regime shift" [@scheffer_catastrophic_2001].

Favourable climatic conditions are expected to shift northward by several
hundred kilometres by the end of the century [@mckenney_potential_2007], but
many studies indicate that tree migration will not keep pace with global warming
[e.g. @zhu_failure_2012; @woodall_assessing_2013; @sittaro_tree_2017;
@renwick_temporal_2015; @talluto_extinction_2017; @vissault_biogeographie_2016].
Indeed, warmer climate may improve reproduction, recruitment, survival, as well
as growth of temperate species at their northern range limit
[@goldblum_tree_2005; @boisvertmarsh_divergent_2019;
@graignic_geographical_2014], while reducing growth and increasing mortality of
boreal species at their southern range limit [@goldblum_tree_2005;
@peng_drought-induced_2011]. These changes in immigration and extinction rates
are driving recent species range shifts [@talluto_extinction_2017], which would,
in theory, continue until forests attain a steady-state distribution, i.e. the
asymptotic proportion of forest states, at a new equilibrium with climate.
However, because trees are long-lived species that disperse over very short
distances, immigration and extinction events in response to environmental
changes are often delayed and forest ecosystems rarely reach their equilibrium
[@talluto_extinction_2017]. If forests are undisturbed, the speed of transitions
between forest biomes will be mainly limited by the persistence and turnover
rate of resident species as well as the dispersal and establishment rates of
migrating species [@neilson_transient_1993], resulting in transient dynamics
that may last a very long time [@hastings_transient_2018;
@jackson_balancing_2010; @talluto_extinction_2017]. While the study of
equilibrium highlights the potential long-term direction of forest dynamics
under some environmental conditions, transients dynamics reveal the short-term
response, which is more relevant to the changes likely to occur in the 21st
century and within the time frame of management plans
[@hastings_transient_2018]. Combining asymptotic and transient dynamics can thus
expand the understanding of forest response at different time scales.

 <!--

@jackson_balancing_2010:
- biodiversity balance is the difference between pre- and post-forcing equilibrium biodiversity

- Immigration and extinction processes can interact to accelerate or delay each other. If persistence of incumbents slows down establishment of potential competitors, then a negative feedback ensues whereby both extinction and immigration are slowed down. The dynamic reverses, however, if successful immigrants exert competitive pressure on incumbents, leading to acceleration of extinction and immigration.
- Recent climate change has also triggered immigration and extinction processes across the globe [6,49], and current climate forecasts indicate that more ecological change will occur over the coming decades [50–52]


Transient behaviour is essential to study and to consider as it directly influences the long run process. Transient dynamics are
now thus recognised as crucial parts of ecosystem dynamics
-->

<!--

Les paragraphes qui suivent sont intéressants et pertinents mais mélangent deux éléments qui sont importants de distinguer : l’effet sur la dynamique (la vitesse de remplacement) et sur l’équilibre. Il faut clairement, en amont, établir la différence entre ces deux effets et décrire un peu mieux l’écologie qui sous-tend ces différences.

Par exemple, *les perturbations accélèrent les changements parce que avec plus de ressources, la croissance est plus rapide. Les perturbations détruisent aussi la biomasse des résidents et créent des opportunités.*

Les perturbations permettent aussi le maintient d’espèces moins tolérantes à l’ombre, ainsi elles ont un effet sur la compétitivité des différentes espèces.

-->

## *P2. Disturbances*

Both gradual invasions and abrupt transitions from one biome to another will
likely take place concurrently. But, given that forests are increasingly
subjected to pervasive climatic stresses and direct human disturbances,
catastrophic transitions are likely to play a dominant role in driving the
climate shift in biomes. With long transients, the disequilibrium between
climate conditions and forest composition will grow larger, which could give
rise to alternative stable states. In this context, forests may fail to return
to their previous states following a disturbance [@johnstone_changing_2016],
resulting in sudden abrupt transitions. Indeed, as climate change slowly
modifies the competitive balance among species, disturbances destroy the
resident community in whole or in part, thus providing establishment
opportunities for migrating species and making resources available for a fast
growth. Consequently, community composition can shift abruptly to species that
are better suited to current conditions [@johnstone_changing_2016;
@renwick_temporal_2015; @turner_disturbance_2010]. For example, canopy gaps have
been shown to locally facilitate establishment of temperate species in mixed
forests of Ontario [@leithead_northward_2010]. In Alaska, white spruce (*Picea
glauca*) is invading black spruce (*Picea mariana*) stands following fire and
permafrost degradation [@wirth_white_2008]. Similarly, moderate disturbances
favoured the increase of warm-adapted species and led to a broad-scale community
thermophilization of forests in Québec [@brice_disturbances_2019].

The coupling of warming and disturbances may however depend on the intensity and
type (natural or anthropogenic) of disturbances [@johnstone_changes_2010]. In
some cases, disturbances may not accelerate biome shifts but instead promote
invasions by early successional species which then displace long-lived
shade-tolerant species. For instance, clearcutting has been found to favour the
expansion of trembling aspen (*Populus tremuloides*) in mixed and boreal stands
of Québec [@laquerre_augmentation_2009; @grondin_have_2018] and Alberta
[@landhausser_disturbance_2010]. While some of the examples above suggest that
disturbances have the potential to catalyse shifts to an increasingly
deciduous-dominated landscape, other simulation studies have concluded that they
are unlikely to drive extensive biome shifts in the coming decades
[@vanderwel_how_2014; @liang_how_2018]. Therefore, more empirical evidence of
this process is essential to determine its importance in broad-scale biome shifts
and to disentangle the role of various intensities and types of disturbances.

## *P3. Soil*

The northward migration of temperate species might be contingent on
their capacity to colonise different types of soil [@lafleur_response_2010;
@bennett_plant-soil_2017; @brown_non-climatic_2014]. For instance, soils of cold
boreal forests generally have lower pH, lower microbial activity and slower
decomposition rates of organic matter than warmer southern temperate forests
[@goldblum_deciduous_2010]. These local and regional variations in soil
properties are expected to slow down or inhibit the establishment of temperate
trees into the boreal forest. For instance, transplant experimental studies have
shown that seedlings of sugar maple (*Acer saccharum*) in conifer-dominated
stands were negatively affected by seed predators and fungal pathogens
[@brown_non-climatic_2014] as well as by soil acidity through reduced foliar
nutrition [@collin_conifer_2017]. However, @kellman_sugar_2004 found that, after
initial high mortality due to seed predation, survival of *Acer saccharum*
seedlings in boreal stands was high, even superior to that in the temperate
stands, potentially because of increased light availability. Hence, it has been
suggested that soil properties in boreal forests may not be a major impediment
to the migration of temperate species showing broad ecological tolerance
[@lafleur_response_2010; @barras_supply_1998; @kellman_sugar_2004]. Nonetheless,
suboptimal soil conditions could delay forest transition under climate change
[@brown_non-climatic_2014]. While experimental studies provide valuable insights
on the potential role of soils at local scales, we need to test the generality
of such constraints, or the lack thereof, on long-term forest dynamics, across
species and scales, to better anticipate future biome transitions.


## *P4. Markov multi-state models*

One approach to investigating biome shifts in response to climate change is to
model transitions of forest plots among states as a stochastic process
influenced by their current state, as well as their current environmental
characteristics. Given the unequivocal distinction between temperate and boreal
forests, the dynamic of tree communities at the boreal-temperate ecotone can be
adequately characterised using discrete functional and successional states,
namely boreal, mixed, temperate and pioneer [@vissault_biogeographie_2016] and
thus can be formalised as a multi-state Markov model
[@jackson_multi-state_2018]. *Markov models provide a useful framework for
modelling changes of state over time using longitudinal data. In epidemiology,
for example, these models are often used to describe the progression of diseases
[@van_den_hout_multi-state_2016]. In ecology, they have been used to study
processes such as forest succession [@runkle_gap_1981; @waggoner_transition_1970;
@lienard_data-intensive_2015], metapopulation dynamics
[@hanski_metapopulation_2003; @moilanen_patch_1999], landcover changes
[@yang_land_2014; @muller_markov_1994], or stage class transitions
[@caswell_matrix_2008].* [**Move to methods?**]

Representation of forest dynamics with Markov chains allows to explore
ecological mechanisms [@wootton_prediction_2001] underlying biome shifts. For
example, transitions to pioneer reflect disturbance, transitions from pioneer
reflect colonisation, dispersal and recruitment limitation and transitions
between the other states reflect competitive exclusion. In addition, multi-state
models can be used to investigate biome shifts from the perspective of both
transient dynamics and long-term equilibrium. Markov transition matrices can be
estimated from the model output and their well-established properties can then
be compared under different scenarios [@hill_markov_2004;
@boulangeat_transient_2018]. For instance, the steady state distribution can be
derived from a transition matrix and allows to infer the potential long-term
forest composition under some environmental conditions, providing insights about
the direction of the current forest dynamic [@waggoner_transition_1970;
@hill_markov_2004]. Moreover, transient periods can be described using the time
of convergence to reach the steady state distribution which measures the length
of the transient period; the turnover time which indicates how fast the
transitions occur and informs about the persistence of forest states; and the
entropy reveals the uncertainty about the next transition. Contrasting
empirically derived transition matrices and their properties among disturbance
scenarios can shed new light on forest dynamics under climate change and may
even provide insights on management measures.



## *P5. Objectives*

In this paper, we investigate the response of forests to recent climate warming
by estimating the transition probabilities among four forest states: boreal,
mixed, temperate and pioneer. Specifically, we address the following questions:
(1) Are forest transition dynamics affected by recent climate change? (2) How do
disturbances and soil characteristics influence the transition probabilities
among forest states under climate change? (3) Do different disturbance types and
intensities impact the potential long-term equilibrium distribution of forest
states? And (4) how do different disturbance types and intensities influence
the short-term transient dynamics?


We expect that the probability of self-transitions will be the highest, i.e.
most forests will not change states, because of the tree slow demography.
However, climate warming should promote colonisation by temperate species and
competitive exclusion of boreal species, resulting in higher transition
probabilities from boreal to mixed and from mixed to temperate than the reverse.
The most conspicuous effect of natural and anthropogenic disturbances is
expected to be the destruction of trees in place, which should provoke
transitions from every states to pioneer. Nevertheless, we also anticipate that
disturbances will catalyse the climate-related transitions (boreal-mixed and
mixed-temperate). In contrast, soil characteristics of coniferous forests (low
pH and poor drainage) should slow down the colonisation by temperate trees.
Together, we expect that these effects on the transient forest dynamics will
result in an increased proportion of temperate forests at equilibrium relative
to the current state distribution. We apply a continuous-time Markov multi-state model
to the dynamics of forest communities to estimate state transition probabilities
and evaluate the influence of environmental covariates on these transitions.
Using results from our multi-state model, we investigate the impact of
disturbances on forest equilibrium and transient dynamics under recent climate
change using several measures: equilibrium state distribution, time to converge
to equilibrium, turnover time and entropy.

# Methods

## Study area and forest inventory data


We used forest inventory plots in Quebec, Canada, to investigate large-scale
transition dynamics in forest communities. Plots have been sampled approximately
every ten years since 1970 and ongoing by the *Ministère des forêts, de la Faune
et des Parcs* [@mffp_placettes-echantillons_2016] in order to monitor changes in
forest productivity and growth. The study area extends from approximately 45° to
52°N of latitude (ca. 795 000 km^2^). It covers six bioclimatic domains (Fig. 1)
and three different forest subzones: the mixed forest (from 47°N to 48°N) marks
the transition between the hardwood forest to the south, which is dominated by
*Acer saccharum*, and the boreal forest to the north, which is dominated by
*Abies balsamea*  and *Picea mariana*.

We selected all inventory plots that had been sampled at least twice as well as
the ones where soil covariates were available. We disregarded plots that were
subjected to active reforestation during the study period because we were
interested in transition dynamics resulting from natural recolonisation. This
yielded a total of 10,388 plots analysed (Fig. 1). The time intervals between
plot surveys varied from 4 to 43 years, with a mean interval of 11 years
($\sigma$ = 3.85).

![Locations of the 10,388 forest inventory plots in meridional Québec, Canada. Colours delimit the six bioclimatic domains. The two southernmost domains (red) are here combined. The number of forest plots in each domain is written in parentheses. The balsam fir-yellow birch domain is the ecotonal zone between hardwood and boreal forests.](res/fig1_region.pdf)

## Forest states

We classified the forest inventory plots into four forest states (Boreal (B), Mixed (M),
Temperate (T) and Pioneer (P)) using species basal area and composition at each sampling
date. We first assigned each studied species to a state: boreal, temperate or
pioneer according to their functional traits [see Table SX in
@brice_disturbances_2019]. For each plot, we computed the total basal area of
each species group and then classified the plot to one of the four states
similar to the @mffp_placettes-echantillons_2016 definitions; Boreal (boreal
species represent >75% of the plot basal area), Temperate (temperate species
represent >75% of the plot basal area), Mixed (when temperate and boreal species
both occupy between >25% and <75% of the plot basal area), and Pioneer (pioneer
species represent >75% of the plot basal area or plot total basal area
<5m^2^/ha). We analysed state transitions between each consecutive plot survey.
Based on this classification, for the 38,091 observations (plots $\times$ number
of years measured), we observed 27,703 state transitions (Fig. 2).

It should be noted that the definition of forest states can affect the results
to some extent. A higher threshold to define the boreal and temperate states
(e.g., >90% instead of >75% of dominance of boreal and temperate, respectively)
can influence the transition probabilities, but the direction of the dynamics
would remain the same.



## Environmental variables

The past annual climatic conditions, covering a period from 1960 to 2013, were
extracted from a 2-km^2^ (60 arc sec) resolution grid for the entire study area
using the ANUSPLIN climate modelling software
[http://cfs.nrcan.gc.ca/projects/3/8; @mckenney_customized_2011]. Plot locations
were intercepted with two bioclimatic variables hypothesised to influence tree
establishment, survival and growth: the mean temperature during the growing
season and the annual climate moisture index (CMI), which is the difference
between annual precipitation and potential evapotranspiration (Table 1). To
reduce the effect of inter-annual climate variability, each climate variable was
averaged over a 10-year period prior to the plot measurement. During the past
five decades, growing season temperatures have increased by 0.16 °C/decade in the plots,
while CMI have shown no clear trends (Fig. S1).

We also collected information pertaining to natural and anthropogenic
disturbances that have affected the forest plots during the study period (Table
1; **Fig. Sx**). At each plot, the type of disturbances (21 types) and their
level of intensity (moderate or major) were recorded during field surveys [see Table S2
in @brice_disturbances_2019; @mffp_placettes-echantillons_2016]. The MFFP
defined major disturbances as events that have eliminated more than 75% of the
total tree basal area, whereas moderate disturbances have eliminated between 25%
and 75%. For our multi-state model, we differentiated two main types of
disturbances: natural disturbances and harvest, with three levels of intensity
each (minor, moderate or major).

Finally, at each plot, several edaphic characteristics were recorded
[@mffp_placettes-echantillons_2016]. Of the available variables, we selected
drainage and pH because they largely affect nutrient availability, soil structural
properties and vegetation development [@tan_environmental_2009]. These two variables also capture
most of the variance in soil characteristics in plots across Quebec and were
orthogonal in a PCA (not shown).

Climate and disturbances were included as time-varying covariates, while soil
was considered as static. Climate variables at time $t$ were used to model
transitions during the interval $t$ and $t + \Delta$t$. Disturbances that
occurred during the interval $t$ and $t + \Delta$t$ were used to model
transitions during the same time period.

Although we tested a well-founded set of covariates, we did not attempt to
include every ecological process affecting the forests of boreal-temperate
ecotone. Notably, because our goal was not to make projections about the future
state of Québec forests, but rather to explore how disturbances and soils have
modified transition dynamics under recent climate change, we did not include any
neighbourhood effect. A neighbourhood effect, which approximate for propagule
availability, is known to affect tree migration [@pearson_climate_2006], but
would be important mostly when making projections, i.e. when climate and
distribution will be decoupled. In contrast, when studying recent forest
dynamics, the neighbourhood composition is very strongly correlated with
climate, which is already included in our multi-state model.



Table 1. Description of the covariates used in the multi-state models.

|Covariate name   |Covariate description                                       |
|:----------------|:-----------------------------------------------------------|
|**Climate**      |                                                            |
|Temp             |Mean temperature during growing season, 10-year average prior to plot measurement (°C). |
|CMI              |Mean annual Climate Moisture Index, 10-year average prior to plot measurement (cm). |
|**Soil**         |                                                            |
|pH               |pH of of the surface horizon                                |
|Drainage         |7 classes of soil drainage, which range from excessive to very poor, that were treated as numeric.|
| **Disturbances**|                                                            |
|Logging          |Tree harvesting, including clearcutting, selection cutting, shelterwood cutting, seed-tree cutting, etc. None or minor (0), moderate (1) or major (2). |
|Natural          |Natural disturbances, including forest fires, insect outbreaks, windfall, etc. None (0), moderate (1) or major (2). |


## Analysis


### Continuous-time multi-state Markov model

We formalised forest dynamics with a continuous-time multi-state model
[@jackson_multi-state_2018] in which transitions among states depend upon the
current state, time interval, climate, disturbances and soil characteristics
(Fig. 2). Fitting methods were based on survival analysis and disease
progression models [@jackson_multi-state_2018; @van_den_hout_multi-state_2016].

Markov models are often built using discrete time steps. However, because (1)
time intervals between surveys are irregular, (2) for each time interval,
multiple transitions are possible, and (3) the exact times of state changes are
unobserved (i.e. observations are interval-censored), a continuous-time Markov
model, in which time is treated as continuous, is preferable
[@van_den_hout_multi-state_2016].

For states $r,s \in {B, M, P, T}$ and time $t, \Delta{t} ≥ 0$, transition probabilities
($p_{rs}$) are defined as the probability that a plot in state $r$ at time $t$
is in state $s$ at time $t + \Delta{t}$ and can be written as:

$$P_{r,s}(t + \Delta{t}) = P(S_{t + \Delta{t}} = s | S_{t} = r)$$

In a four-state transition model, the transition probability matrix $P(t +
\Delta{t})$, hereafter simplified to $P(t)$, is a 4 $\times$ 4 matrix, where the
rows are the current states and the columns the future states, containing the
transition probabilities $p_{rs}(t)$ for a specified time interval. For a
time-homogeneous model, $P(t)$ is solved by taking the matrix exponential of
the intensity matrix $Q$ scaled by the time interval:

$$P(t) = e^{tQ}$$

The intensity matrix $Q$ contains transition intensities $q_{r,s}$ which
represent the instantaneous risk of moving from state $r$ to state $s$:

$q_{r,s} = \lim_{\Delta \to 0} \frac{P(Y_{t+\Delta} = s | Y_t = r)}{\Delta}$, on
off-diagonal elements,

$q_{r,r} = − \sum_{s \neq r}{q_{rs}}$, on diagonal elements.

We can define transition-specific hazard regression models for those states $r,s
\in {B, M, P, T}$ between which a direct transition is possible according to the
specified multi-state process (Fig. 2). The intensities $q_{r,s}$ can be
modelled as a product of a baseline hazard $q_{rs.0}$ and a log-linear effect of
the explanatory variables $x(t)$ and their coefficients $\beta_{rs}$:

$q_{rs}(t) = q_{rs}(t|x(t)) = q_{rs.0}(t)exp(\beta_{rs}'x(t))$.

In this model, $q_{rs.0}(t)$ is a baseline hazard function that describes the
risk for a reference plot $i$ with environment $x_i(t) = 0$, and
$exp(\beta_{rs}'x(t))$ is the relative increase or decrease in risk associated
with the set of characteristics $x_i(t)$. This model allows one to include the
effect of time-dependent covariates on transition intensities and therefore to
relax the time homogeneity assumption of Markov models. Time-dependent
covariates, such as climate and disturbances, are assumed to be
piecewise-constant, i.e. the hazard is constant within a specified time interval
$[t, t + \Delta{t}]$ and depends on the covariate value at $t$, but is allowed
to change between the intervals.

Estimation of model parameters can be obtained by maximising the log-likelihood
function using the transition probability matrix. The contribution of the plot
$i$ at time $j$ to the likelihood is given by:

$LL_i(\theta | s,x) = \prod\limits_{j=1}^{J} P(S_j=s_j|S_{j-1}=s_{j-1},\theta,x)$,

where $\theta$ is the vector with all the model parameters, $x$ denotes the
vector with the covariate values, and $s$ denotes the observed state trajectory
$s_1,...,s_J$ at times $t_1,...,t_J$. The full likelihood function is the
product of contributions for all $N$ plots:

$LL(\theta) = \prod\limits_{i=1}^{N} LL_i(\theta | s,x)$.

### Definition of candidate models

Because the states are defined based on stand composition, we assumed that an
instantaneous transition from Boreal to Temperate and from Temperate to Boreal
was impossible (there is a necessary transition through the Mixed state).
However, all states can transition directly to Pioneer when disturbed (Fig. 2).

We built five different models: one baseline model with intercept only, one for
each subgroup of covariates independently (climate, soil and disturbances), and
one full model, which is a combination of all the covariates (Table 1). Because we are
estimating multiple state transitions in a single model (all $q_{rs}$ in Fig. 2),
the number of parameters increases rapidly with the number of covariates (number of
modelled transitions (here 10) $\times$ number of covariates). Thus, to reduce the
number of parameters, we hypothesised that transitions from any state to Pioneer
were only determined by disturbances, while climate and soil variables should not
directly influence these transitions. All quantitative variables were
standardised ($\mu$ = 0, $\sigma$ = 1) prior to running the models.

![Multi-state transition diagram (a), intensity matrix (b) and equations of our
full model (c). Directional arrows in (a) depict the allowed transitions between
states. The numbers represent the percentage of observed transitions between
states (nb~rs~/nb~r.~ $\times$ 100). Instantaneous transition from Boreal to
Temperate and vice versa are considered impossible in the model (hence the
absence of arrows in the diagram and the zeros in the Q matrix), however rare
transitions from Boreal to Temperate and from Temperate to Boreal were observed
in the data (less than 0.2%). The Q matrix (b) contains the instantaneous risk
to move from (row) to (column) one of the four states: (B)oreal, (M)ixed,
(P)ioneer and (T)emperate (in that order). All transitions from any states to
Pioneer were modelled only as dependent on disturbances
(c).](res/fig2_trans_diagram.pdf)


### Evaluation of candidate models

We first evaluated the goodness-of-fit of each model containing covariates
(climate, soil, disturbances and full) against the baseline model using
likelihood ratio tests [@jackson_multi-state_2011], which test if the addition
of one or more new parameters significantly increases the likelihood of the
model. The statistical significance of individual covariates in the presence of
the other was determined by comparing the full model with the correspondingly
reduced model using likelihood ratio tests. We also compare and rank the models
using the Akaike's information criterion [AIC; @burnham_model_2002]. The model
that minimises AIC is the best compromise between parsimony and goodness of fit.


### Baseline transition intensities

To evaluate whether forest transition dynamics were affected by recent climate
change (objective 1), we used the baseline hazards estimated by our best model
($q_{rs.0}$) as indicators of the underlying forest response. For each pair of
states, the baseline hazard describes the risk to make the transition for a mean
forest plot (all covariates set to 0).

### Effects of covariates on transition probabilities

To investigate the influence of environmental covariates on transition dynamics
(objective 2), we compared the estimated hazard ratios derived from our best
model ($exp(\beta_{rs})$). We also computed the predicted 10-year transition
probabilities of forest plots under different disturbance scenarios, while
keeping all other covariates at the average found in the ecotonal zone (i.e. the
balsam fir-yellow birch domain), to facilitate visual interpretation of the
impacts of disturbances on the transition matrix structure.

### Effects of disturbances on transient dynamics and equilibrium

To further understand how disturbances can modify the long-term equilibrium
(objective 3) and the forest transient dynamics (objective 4) under climate
change, we computed different properties on the Markov transition matrix and
compared them among levels and types of disturbances along the latitudinal
temperature gradient. An extensive literature describes the multiple properties
of discrete time Markov transition matrix [@caswell_matrix_2001;
@hill_markov_2004] and the latter can also be measured for continuous time
Markov models. We chose four informative and complementary properties that fully
characterise both the short and long time scale dynamics of our modelled system:
(1) the steady-state distribution, which corresponds to the potential long-term
proportion of forest states at equilibrium; (2) the time of convergence towards
the steady-state distribution, which measures the length of the transient
period; (3) the turnover time, which measures the rate of transient successional
changes; and (4) the entropy, which captures the incertitude of the transitions.

First, to measure the potential direction of forest dynamics under a given
scenario, we estimated the steady state distribution, $\pi$. For a regular
Markov process, any initial state distribution $s(0)$ converges to the same
equilibrium as $t$ tends toward infinity:

$\displaystyle \lim_{t \to \infty} s(0)P(t)=\pi$

The vector of equilibrium $\pi$ can be obtained by taking the left eigenvector
of $Q$, which has an eigenvalue of 0, normalised to sum to 1, or the
left eigenvector of $P$, which has an eigenvalue of 1, normalised to sum to
1 [@norris_1997].

Then, the convergence rate to the equilibrium distribution can be measured
using the damping ratio [@hill_markov_2004]:

$\rho = \lambda_{P1} / \lambda_{P2}$ or $\rho = exp(\lambda_{Q1} - \lambda_{Q2})$,

where $\lambda_{P1}$ and $\lambda_{P2}$ are the largest and second-largest
eigenvalues of $P$ ($\lambda_{P1}$ = 1 for stochastic $P$), whereas $\lambda_{Q1}$
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
transition matrix as an index of the incertitude of successional changes. It
measures how uncertain we are about the next new state of a site knowing its
current state. For a continuous-time process, the entropy can be measured using
the jump matrix [@spencer_continuous-time_2005]. The jump matrix contains the
probabilities that the next state after state $r$ is state $s$:

$j_{rs} = −q_{rs}/q_{rr}$,

where $qrs$ is the transition intensity from $r$ to $s$. The entropy of state
$s$ is then:

$H(j_{.s}) = \displaystyle -\sum_r{j_{rs} \times log(j_{rs})}$

The normalised entropy of the whole system is the average of the entropies
over the steady state, divided by $H_{max} = log(n_{state} = 4)$:

$\text{Entropy} = \displaystyle \frac{-\sum_r{\pi_r \times H(j_{.s})}}{H_{max}}$

Values of entropy closer to zero indicated more deterministic transition
dynamics whereas values closer to one indicate more random dynamics.

All analyses were performed using the R programming language version 3.5.1
[@r_core_team_r_2018]. The list of R packages that were used throughout the
analysis is provided in the Supporting Information (Table S1). All data used
in the study, in addition to R scripts that reproduced the analyses and
figures, will be made available online on Github upon manuscript acceptance.


# Results

<!--retirer des résultats que tu ne discutes pas dans les sections: Effect of covariates on transition intensity & Effect of disturbances on transient dynamics

Même que dans la section précédente, il faut essayer de marquer les transitions en faisant appel aux questions/prédictions de l’intro. L’étude des chaines de markov est très technique et difficile à suivre pour une majorité de lecteurs, il faut marquer la logique pour que celui qui ne comprend pas le modèle puisse au moins suivre le développement de l’histoire.-->

In this study, there were a total of 38,091 surveys that were recorded for the
10,388 forest plots. A large fraction of Mixed forests transitioned to Temperate
forests (21.2%) but few did the opposite (5.8%). There were many transitions
from Boreal to Pioneer (10.6%), but even more from Pioneer to Boreal (21.4%).
Temperate and Boreal forests were generally more persistent (90.8 and 87.3%,
respectively, did not transition during the study period) than Mixed and Pioneer
forests (68.5 and 69.6%, respectively; Fig. 2). [*keep or remove this paragraph*]

## Model evaluation

Overall, the full model, which include climate, soil and disturbance variables,
had the best fit and predictive performance for the observed data (Table 2; *Fig. SX*). The
second-best model was the disturbance model, but was far behind with a
difference in AIC of more than 1500 units from the full model (Table 2). All
variable subsets improved significantly the likelihood of the model (all
likelihood ratio tests were highly significant, p << 0.001; Table 2). Moreover,
a 10-fold cross-validation was also performed to evaluate model performance
(see Supplementary Materials for methods*) and revealed that including climate and disturbances
improved overall model predictive performances, while soil variables had a
negligible effect (*Fig. SX*). Hereafter, all inferences about transition
probability parameters were derived from the full model.


Table 2. Comparisons of the five candidate multi-state models. The number of
parameters used in each model corresponds to the number of modelled transitions
(10) $\times$ the number of covariates. The $\Delta$AIC is the difference
between the Akaike's information criterion of each model (AIC~m~) and the
minimum of AIC among all the models (AIC~min~): $\Delta$AIC = AIC~m~ – AIC~min~.
Models are ordered in terms of their $\Delta$AIC. The best model is the one in
bold with $\Delta$AIC = 0.

|             |Covariates        | Nb of parameters| -2 Log-likelihood|$\Delta$AIC|LR test|  
|:------------|:-----------------|----------------:|-----------------:|---------:|-------:|
|Baseline     |Intercept         |               10|           32032.3|    6132.6|< 0.001 |
|Soil         |Drainage, pH      |               24|           31886.9|    6015.2|< 0.001 |
|Climate      |Temperature, CMI  |               24|           30438.7|    4566.9|< 0.001 |
|Disturbances |Natural, Logging  |               50|           27341.3|    1521.6|< 0.001 |
|**Full**     |**All**           |           **78**|       **25763.7**|   **0.0**|**< 0.001**|

\pagebreak

## Baseline transition intensities

The baseline transition intensities of the full model provides insight about the
background rate of forest changes (Fig. 3). The forest dynamics over all the
study area were largely dominated by transitions between Mixed and Temperate
(q~MT~ = 0.0238 and q~TM~ = 0.014) and transitions Pioneer to Boreal states
(q~PB~ = 0.0284; Fig. 3). Mixed forests were 1.7 times more likely to transition
to Temperate than the reverse (q~MT~ / q~TM~ = 1.7), indicating that temperate
species have been abundantly colonising mixedwoods, outcompeting boreal species,
during the study period. For Boreal forests, regeneration from Pioneer to Boreal
was 4.8 times more likely than the resetting of succession by disturbances
(q~PB~ / q~BP~ = 4.8).

![Baseline transition intensities as estimated from the best multi-state transition model... instantaneous risk of moving from one state to another when all covariates are set to 0 (hence, the mean of standardized covariates and the no disturbance level).](res/fig3_baseline.pdf)

## Effect of covariates on transition probabilities

The full multi-state model allows one to reveal interesting relationships
between probabilities of forest state transitions and environmental covariates.
All transitions to Pioneer were highly influenced by disturbances (Fig. 4,
Table S4). As could be expected, major disturbances exert stronger effects
than moderate disturbances (for both natural and logging), but logging had
stronger effects than natural for both levels of intensity. For example, the
risk of transition from Boreal to Pioneer has surged up to 165 times higher for
plots that suffered major logging (logging 2) compared to that which were not
disturbed (minor). Disturbances of all types and intensities favoured
transitions from Mixed to Temperate forests. Major disturbances increased the
risk of transition from Mixed to Temperate by ca. 5 times (Hazard Ratio (HR) =
4.51 and 5.32, for natural and logging, respectively). Moderate disturbances
(natural and logging) also favoured transitions from Boreal to Mixed (HR = 2.80
and 2.78, respectively), while major disturbances had no significant effect.

Climate variables also had a significant influence on most transitions (Fig. 4).
Warmer summer temperature (higher sTP) and higher humidity (higher CMI) favoured
transitions from Boreal to Mixed as well as from Pioneer to Mixed and
Pioneer to Temperate. Interestingly, warmer temperature did not significantly
influence the risk of transition from Mixed to Temperate and higher CMI had a
negative effect.

Although less important (*Fig. SX*), soil variables also impacted state
transitions (Fig. 4, Table S4). Holding the other covariates constant, poorer
drainage (more humid) decreased the instantaneous risk of transition from Boreal
to Mixed by 29% and from Pioneer to Temperate by 17%, but increased the risk of
transition from Temperate to Mixed by 32% (HR = 0.71, 0.83 and 1.32,
respectively). Higher pH (acidic soil) had a considerable negative effect on the
transitions from Temperate to Mixed and a smaller one on transitions from
Pioneer to Boreal (HR = 0.80 and 0.93, respectively). *These changes in risk
ratio associated to soil variables appear almost irrelevant compare to the
effect of disturbances, but a slight increase in drainage can dampen the
positive effect of disturbances. For instance, under moderate natural
disturbances, the estimated risk of transition from Boreal to Mixed is 2.12
(1.04-4.33) at moderate drainage but decreases to 1.50 (0.71-3.20) when
increasing drainage by 1 point.* [*remove?*]


![Hazard ratios (HR) and 95% confidence intervals as estimated from the best multi-state transition model. Each plot represent the estimated HR for transitions from row to column state, e.g. the plot on the first row, second column shows the HR for the Boreal to Mixed transition. The ordinate is in log scale. The HR of predictors are interpretable as multiplicative effects on the hazard, where values above 1 (in blue) indicate that the predictor is associated with a greater risk of state transition, whereas values below 1 (in red) indicate a lower risk of transition. Predictors statistically different from 1 are coloured in dark blue or red. Numbers following disturbance predictors indicate their levels of intensity: 1 = moderate and 2 = major.](res/fig4_HR.pdf)


\pagebreak

Disturbances completely altered the structure of the 10-year transition
probability matrix (Fig. 5). The largest values across most matrices were
generally associated with self-transitions (matrix diagonal), meaning that the
vast majority of forest plots (characterised by the average environmental
conditions found in the ecotone) remained in the same state after 10 years. For
undisturbed forest plots (minor), the self-transitions were very strong but
transitions from Pioneer to Boreal, from Mixed to Temperate, and from Temperate
to Mixed were not trivial. At moderate disturbances, probabilities of
self-transitions decreased, while transitions from Boreal to Pioneer, and from
Mixed to Temperate increased the most. Transitions from Mixed and Temperate to
Pioneer did not increase much at moderate disturbances, likely because such
disturbances were less frequent and less severe than in Boreal forests. The
difference between natural disturbances and logging emerges only at major
disturbances. For both types of major disturbances, the probabilities in the
third column, transitions to Pioneer, showed a great increase compared to
moderate disturbances, but these values exploded in the severely logged
transition matrix, exceeding self-transitions. Interestingly, the estimated
probability of Mixed to Temperate transition remained quite high at major
disturbances.


![Predicted change in 10-year transition probabilities for different disturbance types and levels. All other covariates are fixed at the average conditions found in the ecotone, i.e. the balsam fir-yellow birch domain. Letters correspond to the four forest states: (B)oreal, (M)ixed, (P)ioneer and (T)emperate. Numbers are the modelled transition probabilities from rows to columns and darker colour highlights stronger transitions.](res/fig5_pmatrix.pdf)


\pagebreak

## Effect of disturbances on long-term equilibrium

The steady state forest dominance changes as expected along the temperature
gradient (Fig. 6); the Boreal state dominates at low temperature (high latitude)
and Temperate and Mixed states dominate at high temperature (low latitude), with
a transition boundary located at a growing season temperature of about 12.75°C,
which falls in the ecotonal balsam fir-yellow birch domain (Fig. 6). When moderate
disturbances were included, the steady state proportion of Temperate and Mixed
forests increases and the Boreal-Temperate boundary occurred at lower
temperatures, hence further north in the balsam fir-white birch domain (Fig. 6).
Moderate natural disturbances and logging had a similar positive effect on the
proportion of Temperate and Mixed states. For example, at 12.58°C, which approximates the northern limit of
the balsam fir-yellow birch domain, the steady state proportion of Temperate and
Mixed almost doubles with moderate disturbances (minor: 36%; moderate natural:
58%; moderate logging: 60%). However, because there was a larger reduction in
the proportion of Boreal to the benefit of Pioneer, the displacement of the
boundary was slightly larger for logging than for natural disturbances (displaced at
12.05°C for logging and at 12.18°C for natural disturbances; Fig. 6). Finally,
for major disturbances, the proportion of Boreal forests at steady state
collapsed while that of Temperate and Mixed forests also decreased but to a lesser
extent. The landscape is now dominated by Pioneer forests. The boundary modestly
moved north with major natural disturbances (12.60°C), while it retreated to the
south with major logging (13.05°C). Simulating change in frequency of
disturbances (instead of all or nothing scenarios) shows that the
Boreal-Temperate boundary slowly moved north with increasing moderate
disturbance frequency (natural or logging), but it stagnates with increasing
frequency of major natural disturbances and recedes south with increasing
frequency of major logging (Fig. S4).

![Changes in forest state proportion at equilibrium along the temperature (latitudinal) gradient for different disturbance scenarios. Proportion of Boreal (blue) and Temperate + Mixed forests (red) for minor (solid), moderate (dashed) and major (dotted) disturbances. All other covariates are fixed at the average conditions found in the ecotone, i.e. the balsam fir-yellow birch domain, to focus solely on the effect of disturbances along the temperature gradient. The white (minor disturbances), grey (moderate) and black (major) circles indicate the positions of the boundary between dominance of Boreal forests and dominance of Temperate and Mixed forests (i.e. the advancing front) while the grey (moderate) and black (major) arrows show how disturbances move the boundary. The colours at the top of the plots approximate the positions of the bioclimatic domains along the temperature gradient.](res/fig6_SS_gradient.pdf)


## Effect of disturbances on transient dynamics

Natural disturbances and logging affected forest transient dynamics, with greater
impacts for higher disturbance intensity (Fig. 7). In the minor disturbance
scenario, turnover time was generally longer at low temperature, indicating
slow transition dynamics in forests of northern latitudes (Fig. 6a,b). The
turnover time then rapidly declined to reach a minimum at 13.20°C, between the
sugar maple-yellow birch and the balsam fir-yellow birch domains, and went back
up after this point. This trough, where transition dynamics are the fastest, is
located just a little south of the boundary between the Boreal and
Temperate dominances found in Figure 6. Major disturbances accelerate transition
dynamics all along the temperature gradient, while moderate disturbances also
decrease turnover time but more strongly in the boreal domains (balsam fir-white
birch and spruce-moss domains; Fig. 7a,b). These spatial patterns reflect the
turnover time of the dominant state at each point along the temperature gradient
(Fig. S5).

At minor disturbances, the entropy of the system generally increased from north
to south and peaked at 12.56°C, at the southern end of the balsam white-birch
domain (Fig. 7c,d). This peak illustrates where the transition dynamics are most
uncertain (transition to all states are possible at this point), while it is
very predictable in northern boreal forests (Boreal stays Boreal until it
transitions to Pioneer later on). The peak can be mainly attributed to the
entropy of the Boreal state at the ecotone, and the generally high values at low latitudes
can be principally attributed to the Temperate state (Fig. S6). This
latitudinal pattern of entropy is modified by disturbances. Moderate natural
disturbances decreased the entropy throughout the gradient, but especially where
was the peak (Fig. 7c). With moderate logging, the peak disappeared, and entropy
increased monotonically from north to south (Fig. 7d). When major disturbances
are included, whether natural or logging, the peak of entropy was displaced to the
south (Fig. 7c,d) where it is dominated by the entropy of the Pioneer state
(Fig. S6).

The half-life to equilibrium is the longest at 12.24°C, in the balsam fir-white
birch domain, while it is the fastest in the southernmost latitudes (Fig. 7e,f).
Interestingly, the peak for half-life closely matches the peak for entropy but
matches no feature of the turnover time curve. Moderate disturbances flatten
and shift this peak to the north and the effect of moderate logging (Fig. 7f) is
stronger than natural disturbances (Fig. 7e). Hence, in the balsam fir-white
birch, the half-life to reach equilibrium distribution is reduced almost by half
by moderate logging. With major disturbances, forests all along the temperature
gradient can reach very quickly their steady state distribution (maximum of
about 8 years for major logging and 25 years for major natural disturbances).

![Changes in the characteristics of the forest transient dynamics along the temperature (latitudinal) gradient for different disturbance scenarios: minor (solid), moderate (dashed) and major (dotted) disturbances for both natural (a,c,e) and logging (b,d,f). All other covariates are fixed at the average conditions found in the ecotone, i.e. the balsam fir-yellow birch domain, to focus solely on the effect of disturbances along the temperature gradient. The turnover of the whole system (i.e. whole transition matrix) (a,b) corresponds the time spent in a state before transitioning to the next and is given by the average of each state turnover time over the steady state distribution. The entropy of the whole system (c,d) corresponds to the uncertainty of the next transition and is given by the average of each state entropy over the steady state distribution. The half-life to equilibrium (e,f) is the time taken to reach 50% of the steady state distribution, i.e. when the first eigenvalue becomes twice as large as the contribution of the second eigenvalue. The colours at the top of the plots approximate the positions of the bioclimatic domains.](res/fig7_index_gradient.pdf)


# Discussion

Our study reveals that disturbances are likely to accelerate forest response to
climate change by promoting transitions from boreal to mixed forests and,
particularly, from mixed to temperate forests. Our analysis of the equilibrium
further highlights that the long-term forest dynamics under moderate
disturbances favours an increase proportion of temperate forests and a northward
shift of the boreal-temperate ecotone. Disturbances also modified the forest
transient dynamics, accelerating both the turnover and convergence time and
making the dynamics more predictable. In accordance with the hypothesis
formulated by previous studies [@johnstone_changes_2010;
@johnstone_changing_2016; @vissault_biogeographie_2016;
@brice_disturbances_2019], our findings show that moderate disturbances catalyse
transitions to alternate, temperate-dominated forest state and could therefore
promote regime shifts.



## Are forest transition dynamics affected by recent climate change?

Our results support the hypothesis that recent climate warming influences
transition dynamics among forest states at the boreal-temperate ecotone. The
high baseline transition intensity as well as the high number of observed
transitions from mixed to temperate is consistent with the expectation of a
northward shift in range of temperate trees species into the mixed and boreal
forests. Indeed, the warming trends of the last decades (Fig. S1) have been
shown to increase growth and reproductive rates of temperate species and reduce
growth of boreal species [@reich_geographic_2015; @fisichelli_temperate_2014;
@boisvertmarsh_divergent_2019; @goldblum_tree_2005], thus providing a
competitive advantage of temperate over boreal species.  

Alternatively, the increased transition to temperate species may be a response
to other pressures. Comparisons of pre-settlement and present-day forested
landscapes of North America have highlighted an important deciduous encroachment
in response to historical human activities [@danneyrolles_stronger_2019;
@terrail_reorganization_2019; @boucher_logging-induced_2006]. Similarly,
@johnstone_non-equilibrium_2003 suggested that the northern range expansion of
lodgepole pine following fire is not related to current climate change but was
rather a continued migration initiated in the early Holocene. However,
historical legacies and climate change are likely mutually non-exclusive
explanations. Simulations by @boulanger_climate_2019 showed that the future
climate-induced expansion in temperate species to the detriment of boreal
species would amplify the already ongoing trend since preindustrial times.


<!--Previous studies have shown that successional pathways were strongly affected by disturbances and climate change [].
 Forest transition dynamics are slow. Even if our model overestimates transition, the changes are still very slow. Once again, migration lags.
Nevertheless, forest dynamics appear to be responsive to CC.-->

## Disturbances catalyse forest state transition

Our study highlighted that moderate disturbances favour climate-related
transitions, whereas major disturbances merely promote pioneer states.
Disturbances directly remove trees, which leads to immediate and substantial
changes in forest composition and successional trajectories. Without climate
change, forests are expected to be resilient to normally experienced
disturbances and should thus return to their preceding states. However, climate
change alters the conditions that initially supported the persistence of a
given state, making forests more vulnerable to disturbances
[@johnstone_changing_2016]. Hence, when moderate disturbances remove resident species
and reduce competition for light and nutrient, it likely facilitates
colonisation and establishment by opportunistic temperate species under warmer
conditions [@brice_disturbances_2019; @leithead_northward_2010;
@landhausser_disturbance_2010]. In contrast, severe disturbances in the study
area, primarily clearcutting but also large fires, create openings of very large
extent which likely favour early-successional species that can disperse seed over
long distances, such as *Populus sp* and *Betula sp*
[@landhausser_disturbance_2010].


Compared to the catalysing effect of disturbances, soil characteristics do not
appear to represent a large impediment to state transitions but have the potential to slow
down transitions. Poor drainage constrained climate-related transitions, from
Boreal to Mixed states, but not from Mixed to Temperate. This indicates that
temperate species can readily colonise soils found in mixedwoods but may have
more difficulty to colonise hydric boreal soils. Very poor drainage, often
associated with peatland and thick organic layer, is indeed considered to be one of
the only major edaphic limits for the regeneration of temperate species
[@lafleur_response_2010]. Several studies found that *Acer saccharum*
regenerates well across the ecotone because of its large tolerance to various
soil conditions [@barras_supply_1998; @goldblum_age_2002; @kellman_sugar_2004;
@fisichelli_temperate_2014; @collin_can_2018]. At their northern range limit,
*A. saccharum* and *A. rubrum*, the species contributing most to compositional
changes [@brice_disturbances_2019], are hypothesised to be mostly limited by
cold soil temperature [@barras_supply_1998; @goldblum_age_2002].


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
<!--
SOIL:
@johnstone_changes_2010 "An important insight that emerges from this broad‐scale study is that the potential for fire to drive shifts in successional trajectories is contingent on landscape factors such as site moisture. Consequently, we should expect that the resilience of black spruce forests to changing climate and fire regime will not be uniform across the landscape and that drier spruce forests may have the greatest potential to switch to deciduous‐dominated forests under future environmental change."-->



## Changes in potential long-term equilibrium and range limits

Our model highlights the potential role of disturbances in influencing the
position of the boreal-temperate boundary as well as the proportion of temperate
and boreal biomes at equilibrium. As a result of the increased replacement of
Mixed by Temperate states and a decline of Boreal to Pioneer states, the
equilibrium boreal-temperate boundary shifts northward with moderate
disturbances. While our results should not be interpreted as predictions for the
future, they are useful to highlight the direction of forest dynamics under
different disturbance scenarios. Our results support the simulations of
@boulanger_climate_2019 where harvesting under future climate warming was
projected to promote further invasions of pioneer species, such as *Populus*,
and temperate species, such as *Acer* and *Fagus*, in mixedwoods of Québec. In
contrast, based on their simulations, @liang_how_2018 and @vanderwel_how_2014
concluded that logging would primarily accelerate the expansion of pioneer
forests, but would have little or no effect on extensive biome shifts over the
next century in eastern United States. These apparently conflicting results
could be due to the contrasting tree species responses to disturbances.
Disturbances may facilitate the range expansion of some species but hinder that
of others depending on their functional traits [@matthews_disturbance_2013;
@aubin_traits_2016]. For instance, because of its positive response to past
[@danneyrolles_stronger_2019], recent [@brice_disturbances_2019] and future
[@boulanger_climate_2019] disturbances in Québec, *Acer rubrum* is likely to
play a disproportional role in the temperate biome shift.

<!--
- *Can we really compare?
- in all disturbance scenarios, boreal forests are losing ground primarily to pioneer forests

- Similar to @liang_how_2018 and @vanderwel_how_2014 We also predict that major disturbances will only promote pioneer states to the detriment of boreal states...

Other results:
Simulations by Vieira showed that plantation and enrichment planting shifted northward the boreal-temperate range limits, but not harvesting. Thinning increased the transition from mixed to temperate stands, but did not have any effect on range limits shift.

@boulangeat_transient_2018 even predicted
a southern shift of the boreal-temperate boundary when including the effects of
herbivores.*
-->


## Transient dynamics

[*add short definition of resilience and resistance*]

Beyond their impacts on steady-state distributions, our results suggest that
disturbances may have a substantial influence on forest transient dynamics. In
the continuous boreal zone (spruce-moss domain), forests dominated by *Picea
mariana* are usually characterised by dynamics of stand self-replacement with
minimal compositional changes across disturbance cycles
[@goldblum_deciduous_2010]. Consistent with these dynamics, the turnover time of
undisturbed northern boreal forests was very long and the entropy very low. The
turnover became very rapid with disturbances, but the entropy remained low,
indicating that the dynamics was still as predictable (back and forth transitions
between boreal and pioneer states) and that there was no directional shift associated
with disturbances. *Hence, while resistant at low disturbances, boreal forests
loose their resistance when moderately disturbed, and remain resilient as they
return to their previous boreal state.* Under major disturbances, boreal forests
collapsed to pioneer state and reached this new equilibrium swiftly (short
half-life). This interpretation agrees with previous studies suggesting that
boreal forests can easily shift into an alternative treeless state in response
to severe or repeated disturbances [@payette_shift_2003;
@sanchez-pinillos_resistance_2019].

In contrast, the ecotone is characterised by a rapid turnover and a high entropy
indicating abrupt compositional shift which can go in all directions; *hence this
zone is neither resistant nor resilient*. Compared to northern boreal forests,
the short turnover time implies a low persistence, hence a low resistance, of
the forest states in this region even under minor disturbances. This result
corroborate the predictions made by @vissault_biogeographie_2016, where mixed
forests would undergo a swift conversion to temperate forests in the next decades
whereas boreal forests would present a large inertia. The dynamics of the
ecotone appear unstable because it is caught between two stable states, i.e.
boreal to the north and temperate to the south. Under moderate disturbances, the
probability of transitioning to Temperate increases to the detriment of the
other possible states (Fig. 6), hence the entropy is decreased, and the dynamics
become more predictable. Such a clear directional shift strongly indicates
non-equilibrium dynamics in this region. Although turnover is fast, half-time to
equilibrium is long because a forest may not move in the right direction and may
undergo multiple transitions.

<!--
- *implication for resilience and resistance
  - long turnover = high persistence = high resistance?
  - high entropy = low stability/resilience?

- Compare to @boulangeat_transient_2018 + Vieira
  -In contrast to our results, @boulangeat_transient_2018 showed that disturbances
  by browsers "reduced the asymptotic resilience of the system, decreasing the
  rate at which the new equilibrium was approached." Effect of browsers are likely less dramatic than tree logging (here moderate disturbances remove at least 25% of the basal area).

- active migration zone more sensitive

- at this scale there appear to be no difference between natural and anthropogenic disturbances.*


clarifier définition de résilience

MOVE TO DISCUSSION:
relaxation time = time needed to reach a novel equilibrium following perturbations (Diamond, 1972; Hylander & Ehrlén, 2013).
Established communities are generally thought to be resistant to
competitive displacement and resilient to commonly experienced disturbances
[@grondin_have_2018; @seidl_disturbance_2014]. Resistance can be defined as the
ability of an ecosystem to persist through time following a disturbance, whereas
*resilience is the capacity to recover its pre-disturbance composition*
[@gunderson_ecological_2000; @holling_resilience_1973]. Because the ability to
persist and recover following a pulse disturbance (e.g., fire or logging) can be
affected by a press disturbance (e.g., climate change), investigating resistance
and resilience can provide precious insights into their interacting effects.
-->



## Ecological and management implications

<!--natural and anthropogenic disturbances, combine with the
underlying effect of climate change, will likely modify forest dynamics by providing
establishment opportunities for some migrating species in otherwise competitive
environments, thus triggering state shifts in forest ecosystems.

Inertia, persistence –Ability of a living system to survive moderate disturbances

Resilience –Ability of a living system to be restored through secondary succession after a moderate disturbance

Tipping point –Any additional stress can cause the system to change in an abrupt
-->

The relationships that we revealed between forest state transitions and
disturbances provide a strong empirical basis for predicting the types of
changes in forest dynamics that are likely to unfold in the 21st century. A
shift in dominant forest cover from conifer to deciduous broadleaf species not
only entails changes in tree species diversity and composition, but a complete
transformation of forest dynamics and functions. In the long term, this regime
shift could locally increase tree diversity [@berteaux_cc-bio_2010] and productivity [REF], increase
carbon sequestration [@thurner_carbon_2014; as long as mortality is limited
@nrdc_pandoras_2018], modify disturbance regimes [reduced flammability of
broadleaf species @terrier_potential_2013, and reduced sensitivity to current
outbreak-prone pest @mffp_insectes_2018], alter soil microbial activity
[@laganiere_how_2010] and affect wildlife distribution [@mizel_rapidly_2016].

Such regime shifts will also have large repercussions on forest management
strategies in area where silvicultural practices are tailored to regional
disturbance regimes and rely on natural regeneration. In Québec, ecosystem-based
forest management seeks to maintain the composition and structure of a reference
state, defined by the preindustrial forest conditions [@pinna_amenagement_2009].
Yet, @boulanger_climate_2019 showed that such management would fail to restore
or approach historical forest conditions under future climate change. Trying to
maintain a historical state might not be possible and our results reinforce the
idea that forest management should consider the present system state in relation
to their transient dynamics as well as its most likely trajectory. Our study
also reveals the potential of moderate disturbances to catalyse climate-related
transitions. This suggests that planned logging practices could be used to
reduce extinction debt and colonisation credit, and thus tree migration lags. In
addition to modified silviculture, temperate species could be planted farther
north outside their current range to speed up their migration
[@iverson_tree-species_2013; @duveneck_measuring_2016; Vieira et al.]. However,
there is still plenty of uncertainty and many important unanswered questions
that need to be addressed. For instance, will multiple interacting disturbances
exacerbate tree mortality? Which species will be able to benefit from the
opportunities created by canopy openings? And how will these rapid transitions
feedback to ecosystem processes and global cycles? The rate of recent climate
change already outpace tree migration capacity [@sittaro_tree_2017], and even more so the
scientific capacity to understand and mitigate its consequences
[@jackson_balancing_2010]. Therefore, in order to insure long-term forest health
at the boreal-temperate ecotone, we need to simultaneously limit global warming
through drastic reduction of greenhouse gas emission and intensify research
effort to develop effective adaptation strategies.


<!--
Moreover, the fact that, once disturbed, many forest are transitioning to new states suggest the existence of alternative stable state. Management practices might not be adapted to transient dynamics.

- But disturbance-related changes may exacerbate tree mortality and compromise forest resilience in the long run + not all species will benefit from disturbances...
the inherent differences be-tween fire and logging can alter the composition of plant communities
The climate of the future will likely lead to higher mortality among mature trees, because of the greater frequency of droughts, fires, forest-leveling windstorms, and outbreaks of native and exotic insect pests and diseases.

- Stop climate change! The global scientific community [@ipcc_climate_2014], as well as Canadian and other national governments, recognize that climate change is driven by human activity.
Global climate change will alter natural disturbance regimes because many disturbances have a significant climate forcing. What will happen when disturbance regimes change?
-->

<!--

Moreover, as species move, they alter themselves the forest dynamics
and processes and, through a positive feedback loop [@bergeron_fire_2004],
forests may not be able to return to their previous state.

disturbance that fall outside the range of normal background levels.

As a result, conditions that support the persistence of mature forests may not be amenable to forest regeneration.
-->

## References
