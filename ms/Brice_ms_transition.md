---
title: Moderate disturbances accelerate forest transition dynamics under climate change in the temperate-boreal ecotone of eastern North America
documentclass: article
font: 12pt
indent: 2m
papersize: a4paper
geometry: margin=1in
header-includes:
    - \usepackage{setspace}
    - \doublespacing
    - \usepackage{lineno}
    - \linenumbers
    - \usepackage{float}
    - \floatplacement{figure}{H}
    - \usepackage{caption}
    - \renewcommand{\figurename}{\bfseries Figure}

---

**Running title:**
Forest transitions under climate change

# Abstract

Several temperate tree species are expected to migrate northward and colonise
boreal forests in response to climate change. Tree migrations could lead to
transitions in forest types, and these could be influenced by other non-climatic
factors, such as disturbances and soil conditions.

Using decades of forest inventory data, we model the state transition dynamics
of Québec's forests to identify the environmental conditions that promote or
prevent forest transitions. We further investigate how different disturbance
types and intensities impact forests' short-term transient dynamics and
long-term equilibrium.

We analysed over 10,000 forest inventory plots, sampled from 1970 to 2018 in
meridional Québec, Canada. We used a continuous-time multi-state Markov model to
quantify the transition probabilities between forest states (temperate, boreal,
mixed, pioneer) **in relation to climate (mean temperature and climate moisture
index during the growing season), soil conditions (pH and drainage) and
disturbances (levels of severity of natural disturbances and logging)**. We
described transient dynamics and equilibrium under different disturbance
scenarios, using properties of Markov transition matrices.

Whereas the larger portion of the forest plots remained in the same state during
the study period, the most common transitions were from mixed to temperate
states, as well as from pioneer to boreal forests. In our study, variations in
transition probabilities were mainly driven by natural and anthropogenic
disturbances and secondarily by climate, whereas soil characteristics exerted
relatively minor constraints. While major disturbances only promoted transitions
to the pioneer state, moderate disturbances increased the probability of
transition from mixed to temperate states. **Long-term projections of our model
under the actual environmental conditions indicates that moderate disturbances
would promote a northward shift of the temperate forests.** Moreover, using the
properties of Markov transition matrices, we found that disturbances reduce
turnover and convergence time for all transitions, leading to an acceleration of
forest dynamics. Hence, moderate disturbances could catalyse rapid forest
transitions and accelerate broad-scale biome shifts.




## Keywords

Climate change,
Natural disturbances and logging,
Continuous-time multi-state Markov model,
Québec,
Temperate-boreal ecotone,
Transition probabilities,
Equilibrium,
Transient dynamics


\pagebreak

# Introduction

Global climate warming has led to upward (altitudinal) or poleward (latitudinal)
migration of species across the globe [@parmesan_globally_2003;
@walther_ecological_2002; @chen_rapid_2011]. In ecotones, where transition
between vegetation biomes occurs, these shifts in species distributions entail
far reaching consequences for forest ecosystems [@goldblum_deciduous_2010;
@evans_borealtemperate_2017; @bright_climate_2014]. In some cases,
climate-induced shifts in tree species distributions might trigger a "regime
shift" [@scheffer_catastrophic_2001] and transform treeless tundra into boreal
forests [@danby_variability_2007; @kapralov_changes_2006], tropical forests into
savanna [@hirota_global_2011] or coniferous forests into deciduous forests
[@boulanger_climate_2019; @price_anticipating_2013]. As ecological processes may
strongly differ among these biomes, this reorganisation of biodiversity not only
impacts local species composition [@williams_novel_2007], but also alters the
functional and structural characteristics of communities
[@esquivel-muelbert_compositional_2018; @danneyrolles_stronger_2019;
@scheffer_thresholds_2012], hence feedbacks to microclimates, biogeochemical
cycles and ecosystem functioning [@anderson_biophysical_2011]. However, these
large transitions in forest types are still poorly understood notably because
ecotones are not solely controlled by **regional** climate but by many other
**landscape- and local-scale** factors that could accelerate or slow-down these
changes.

Range shift dynamics ultimately arise from change in **fine-scale demographic
processes** (e.g., recruitment, growth, mortality) that determine where a tree
species can establish and persist [@godsoe_integrating_2017; @schurr_how_2012].
Whereas range expansion depends upon dispersal and establishment of new
individuals, range contraction is the result of declining vitality and mortality
[@jump_altitude-for-latitude_2009]. In the temperate-boreal forest ecotone,
recent climate warming has indeed been shown to improve recruitment, survival
and growth of some temperate tree species at their northern range limit
[@fisichelli_temperate_2014; @boisvertmarsh_divergent_2019; @sittaro_tree_2017;
@grundmann_impact_2011; @bolte_understory_2014], leading to range expansion at
their leading edge. At the same time, boreal species were competitively
disadvantaged by slower growth and larger increase in mortality associated to
heat and drought stress [@peng_drought-induced_2011; @goldblum_tree_2005;
@grundmann_impact_2011; @bolte_understory_2014]. Hence, as climate warms and tips
the balance in favour of temperate over boreal species, forests at the ecotone
are expected to transition from coniferous to mixedwood and from mixedwood to
temperate deciduous [@boulanger_climate_2019; @price_anticipating_2013;
@chen_modeling_2002; @lindner_climate_2010].

The reported shifts in species distributions are, however, much slower than the
rate of climate change [@sittaro_tree_2017; @talluto_extinction_2017]. Such lags
in species responses are hypothesised as primarily due to demographic
constraints [@svenning_disequilibrium_2013; @renwick_temporal_2015]. Because
trees are long-lived species that disperse over very short distances,
colonisation and extinction events in response to environmental changes are
often delayed, such that forests are rarely in equilibrium with their
environment [@talluto_extinction_2017]. Hence, if forests are undisturbed,
transition rates between forest types following natural succession pathways will
be mainly limited by the persistence and turnover of resident species
[@loehle_forest_2000; @bouchard_tree_2019] as well as the dispersal and
establishment rates of migrating species [@neilson_transient_1993], resulting in
large disequilibrium and transient dynamics that may last a very long time
[@hastings_transient_2018; @jackson_balancing_2010; @talluto_extinction_2017].

Disturbance **events, such as fire and harvesting,** also affect demographic
processes and increase turnover, and are thus likely to influence forest
responses to climate change [@anderson-teixeira_altered_2013;
@serra-diaz_disturbance_2015; @boulanger_climate_2019; @bolte_understory_2014].
Indeed, as **global trends** in climate warming slowly modifies the competitive balance
among species, pulse disturbances remove the resident community in whole or in
part, thus providing establishment opportunities for migrating species and
making resources available for a fast growth. Consequently, following a
disturbance, forest composition may shift to species that are better suited to
current conditions and fail to return to its previous state
[@johnstone_changing_2016; @renwick_temporal_2015; @turner_disturbance_2010].
For example, canopy gaps have been shown to locally facilitate establishment of
temperate species in mixed forests of Ontario [@leithead_northward_2010]. In a
nature reserve of Scandinavia, @bolte_understory_2014 showed that Norway spruce
(*Picea abies*) was particularly sensitive to the combination disturbances and
warming which benefited the growth of European beech (*Fagus sylvatica*). In a
previous study, we showed that moderate disturbances (i.e., disturbances that
removed between 25-75% of the tree basal area) have favoured the increase of
warm-adapted species and led to a broad-scale community thermophilization of
forests at the temperate-boreal ecotone in Québec [@brice_disturbances_2019].
These results call for further investigation on how disturbances may affect
forest dynamics under recent climate change and whether **their effects can
scale-up to** trigger punctuated and episodic shifts in forest types.


**Cross-scale interactions between landscape disturbances and global warming could
drive abrupt transitions between forest types [@peters_crossscale_2007; @allen_interactions_2007]. Given
that forests are increasingly subject to human disturbances, such non-linear
processes** could play a key role in driving the climate shift in biomes. Some
simulation studies have however concluded that they are unlikely to drive
extensive biome shifts in the coming decades [@vanderwel_how_2014;
@liang_how_2018]. In some cases, disturbances may only promote early
successional species, which can then displace long-lived shade-tolerant species.
For instance, clearcutting has been found to favour the expansion of a pioneer
species, the trembling aspen (*Populus tremuloides*), in mixed and boreal stands
of North America [@laquerre_augmentation_2009; @grondin_have_2018;
@landhausser_disturbance_2010]. This suggests that the effect of disturbances on
forest dynamics may depend on their intensity and type (natural or
anthropogenic). Indeed, logging strongly differs from natural disturbances in
severity, frequency, selectivity and spatial extent
[@boucher_logging-induced_2006; @schulte_homogenization_2007], which could alter
successional pathways. Therefore, more empirical evidence is essential to
disentangle the role of various intensities and types of disturbances in
broad-scale biome shifts.


The northward migration of temperate species may nevertheless be contingent on
their capacity to colonise different types of soil [@lafleur_response_2010;
@bennett_plant-soil_2017; @brown_non-climatic_2014; @carteron_soil_2020]. Soils
of cold boreal forests generally have lower pH, lower microbial activity and
slower decomposition rates of organic matter than warmer southern temperate
forest soils [@goldblum_deciduous_2010]. These local and regional variations in
soil properties are expected to slow down or inhibit the establishment of
temperate trees into the boreal forest. For instance, transplant experimental
studies have shown that seedlings of sugar maple (*Acer saccharum*) in boreal
soils were negatively affected by soil biotic and abiotic conditions
[@carteron_soil_2020; @brown_non-climatic_2014]. In contrast,
@kellman_sugar_2004 found a higher survival of *Acer saccharum* seedlings in
boreal stands than in hardwood stands, potentially because of better light
availability. Hence, it has been suggested that soil properties in boreal
forests may not be a major impediment to the migration of temperate species
showing broad ecological tolerance [@lafleur_response_2010; @barras_supply_1998;
@kellman_sugar_2004]. Nonetheless, suboptimal soil conditions under a boreal
canopy could delay forest transitions under climate change
[@solarik_priority_2019]. While experimental studies provide valuable knowledge
on the role of soils at local scales, the importance of such constraints on
long-term forest dynamics should be evaluated at regional scale and across
species to better anticipate future biome transitions.


One approach to investigating biome shifts in response to climate change is to
model transition probabilities between forest states using a Markov chain
approach. Given the unequivocal distinction between temperate and boreal
forests, the dynamics of tree communities at the temperate-boreal ecotone of
North America can be adequately characterised using discrete **ecological** and
successional states, namely boreal (stands dominated by boreal coniferous
species), mixed (mixed stands of coniferous and deciduous species), temperate
(stands dominated by temperate deciduous species) and pioneer (stands dominated
by early successional species, **which can be found any disturbed habitats
across the latitudinal gradient**) [@vissault_biogeographie_2016]. Using such
classification, the forest dynamics thus can be formalised as a multi-state
Markov model, where transitions among states are represented by a stochastic
process influenced by their current state and environmental characteristics of
interest [@jackson_multi-state_2018]. The Markov framework has been previously
used to study forest succession [@runkle_gap_1981; @waggoner_transition_1970;
@lienard_data-intensive_2015] notably because it is based on a straightforward
definition of transitions between various forest states and provides a simple
mechanistic interpretation of the estimated transition probabilities. This
method thereby offers the possibility of exploiting the full complexity and
temporal depth of forest inventory data, while buffering the idiosyncrasies of
species responses [@strigul_modelling_2012; @lienard_modelling_2016].


Representation of forest dynamics with Markov chains allows the exploration of
**fine-scale** ecological mechanisms [@wootton_prediction_2001] underlying
**broad-scale** biome shifts. For example, transitions to pioneer reflect
disturbance, transitions from pioneer reflect colonisation, dispersal and
recruitment limitation, and transitions between the other states reflect
competitive exclusion. In addition, multi-state models can be used to
investigate biome shifts from the perspective of both transient dynamics and
long-term equilibrium. Markov transition matrices can be estimated from the
model output and their well-established properties can then be compared under
different scenarios [@hill_markov_2004; @boulangeat_transient_2018]. For
instance, the equilibrium or steady-state distribution can be derived from a
transition matrix and used to infer the potential long-term forest composition
under some environmental conditions [@scheffer_catastrophic_2001], providing
insights about the direction of current forest dynamics
[@waggoner_transition_1970; @hill_markov_2004]. Moreover, transient periods can
also be described: the time of convergence to equilibrium measures the length of
the transient period; the turnover time indicates how fast the transitions occur
and informs about the persistence of forest states; and the entropy reveals the
uncertainty about the next transition. Contrasting empirically derived
transition matrices and their properties among disturbance scenarios can shed
new light on forest dynamics under climate change and may even provide hints
about management measures.


Here, we investigate how **regional-scale** forest dynamics is influenced by disturbances and soil
conditions under recent climate warming. In particular, we ask the following
questions: (1) How recent forest transitions dynamics vary with climate, soil
and disturbances? (2) Do different disturbance types and intensities impact
the potential long-term equilibrium distribution of forest states? (3) How do
different disturbance types and intensities influence the short-term transient
dynamics under climate change? And (4) what is the relative importance of tree
demographic processes underlying the transition dynamics? We answer those
question by estimating the influence of environmental covariates on transition
probabilities among four forest states (boreal, mixed, temperate and pioneer)
using a continuous-time Markov multi-state model. Using results from our model,
we then examine the impact of disturbances on forest equilibrium and transient
dynamics by comparing different complementary matrix properties.

We expect that climate warming should promote colonisation by temperate species
into mixed and boreal forests and competitive exclusion of boreal species,
resulting in higher transition probabilities from boreal to mixed and from mixed
to temperate, rather than the reverse. The most conspicuous effect of
disturbances is expected to be the destruction of trees in place, which should
provoke transitions from other states to pioneer. Nevertheless, we also
anticipate that disturbances will favour climate-related transitions
(boreal-mixed and mixed-temperate), whereas soil characteristics of coniferous
forests (low pH and poor drainage) should slow down colonisation by temperate
trees. Disturbances should also accelerate the transient dynamics by shortening
turnover and convergence times. Together, these effects on transitions should
influence the steady-state distribution by promoting an increase in the
proportion of temperate forests in the long run.

# Methods

## Study area and forest inventory data

We used forest inventory plots in Québec, Canada, to investigate broad-scale
transition dynamics in forest communities. Permanent plots have been sampled
approximately every ten years from 1970 to 2018 (and ongoing) by the *Ministère
des forêts, de la Faune et des Parcs* [@mffp_placettes-echantillons_2016] in
order to monitor changes in forest productivity and growth. The study area
extends from approximately 45° to 52° North latitude (ca. 795 000 km^2^). It
covers six bioclimatic domains (Fig. 1) and three different vegetation zones; the
mixed forest, which corresponds to the balsam fir-yellow birch domain (from 47°N
to 48°N; hereafter, the ecotone), marks the transition between the hardwood
forest to the south, dominated by *Acer saccharum*, and the boreal forest to the
north, dominated by *Abies balsamea*  and *Picea mariana*. Plots were randomly
positioned across these three zones with a decreasing sampling intensity
northward [@mffp_reseaux_2014].

The natural disturbance regimes vary considerably along the latitudinal gradient
of the study area, with fires in the northern boreal forests, spruce budworm
outbreaks in the mixedwood forests, and small windthrows and treefall gaps in
the southernmost deciduous forests [Fig. S1; @goldblum_deciduous_2010].
Anthropogenic disturbances are not homogeneously distributed either; clearcuts
are more frequent in northern regions, while in southern regions partial cuts
are more common [Fig. S1; @boucher_logging_2009].


We first selected all inventory plots that had been sampled at least twice. We
then disregarded plots that were subjected to active reforestation (i.e.,
plantation) during the study period because we were interested in transition
dynamics resulting from natural recolonisation processes. Finally, we kept plots
for which soil covariates were available. This yielded a total of 11,058 plots
analysed (Fig. 1). The time intervals between plot surveys varied from 3 to 39
years, with a mean interval of 11 years (sd = 3.45; Fig. S2).


## Forest states

We classified the forest inventory plots into four forest states using species
basal area and composition at each sampling date. We first assigned each
studied species to a group based on their traits and their distribution [Table
S1; see @brice_disturbances_2019 for details]: boreal species are mostly
coniferous trees with a northern distribution; temperate species are mostly
deciduous trees with a southern distribution; and pioneer species have low shade
tolerance and are generally found in any disturbed habitats. For each plot, we
computed the total basal area of each species group and then classified the plot
to one of the four states similar to the @mffp_placettes-echantillons_2016
definitions; Boreal (boreal species representing >75% of the plot basal area),
Temperate (temperate species representing >75% of the plot basal area), Mixed
(when temperate and boreal species both occupy between >25% and <75% of the plot
basal area), and Pioneer (when the basal area of pioneer species is superior to
that of boreal and temperate species or when plot total basal area <5m^2^/ha).
We analysed state transitions between consecutive plot surveys. Based on this
classification, for the 42,633 observations (plots $\times$ number of years
measured), we recorded 31,690 state transitions, including self-transitions
(Fig. 2, Table S2).

The definitions of forest states can affect the results to some extent. A higher
threshold to define the boreal and temperate states (e.g., >85% instead of >75%
of dominance of boreal and temperate, respectively) influences the transition
probabilities, but the direction of the dynamics remains the same (see
comparison between Tables S3 and S4).



## Environmental variables

Annual climatic conditions, covering a period from 1960 to 2018, were
extracted from a 2-km^2^ (60 arc sec) resolution grid for the entire study area
using the ANUSPLIN climate modelling software
[http://cfs.nrcan.gc.ca/projects/3/8; @mckenney_customized_2011]. Plot locations
were intercepted with two bioclimatic variables hypothesised to influence tree
establishment, survival and growth: the mean temperature during the growing
season and the climate moisture index (CMI; difference between precipitation and
potential evapotranspiration) from May to September (Table 1). To reduce the
effect of inter-annual climate variability, each climate variable was averaged
over a 10-year period prior to the plot measurement. During the 1950 until the
present day, growing season temperatures have increased by 0.17 °C/decade in the
plots, while CMI have shown no trends (Fig. S3).

We also collected information pertaining to natural and anthropogenic
disturbances that have affected the forest plots during the study period (Table
1; Fig. S1). At each plot, the type of disturbances (21 types) and their level
of severity were recorded during field surveys [see Fig. S1 for details;
@mffp_placettes-echantillons_2016]. For our multi-state model, we differentiated
two main types of disturbances: natural disturbances and logging, with three
levels of severity each (0, no or minor; 1, moderate; 2, major). The MFFP defined major
disturbances as events that have resulted in a loss of more than 75% of the
total tree basal area, whereas moderate disturbances have caused between 25 and
75% of loss. When the loss in basal area is less than 25%, it is considered to be
minor.

Finally, at each plot, several edaphic characteristics were recorded
[@mffp_placettes-echantillons_2016]. We selected drainage and pH because they
largely affect nutrient availability, soil structural properties and vegetation
development [@tan_environmental_2009], and also because they captured most of
the variance in soil characteristics in our plots.

Climate and disturbances were included as time-varying explanatory variables
(often called covariates in survival models), while soil variables were
considered as static. Climate variables at time $t$ were used to model
transitions during the interval $t$ and $t + \Delta t$. Disturbances that
occurred during the interval $t$ and $t + \Delta t$ were used to model
transitions during the same time period.

Note that we solely focused on a parsimonious set of variables that allowed us
to determine how climate, disturbances and soils influence transition dynamics.
We decided not to include an index of propagule availability, even though it is
known to affect tree range shifts [@pearson_climate_2006], as forest
composition is already very strongly correlated with our climate covariates
[Fig. S4; @goldblum_deciduous_2010; @vissault_biogeographie_2016;
@paquette_effect_2011]. Our model is therefore well-suited for our
research goals; however it is not designed to make future range shift
projections.



## Analysis


### Continuous-time multi-state Markov model

We formalised forest dynamics with a continuous-time multi-state model
[@jackson_multi-state_2018; @van_den_hout_multi-state_2016] in which transitions
among states depend upon the current state, time interval, climate, disturbances
and soil characteristics (Fig. 2). This type of model takes into account the
fact that (1) time intervals between surveys were irregular, (2) multiple
transitions were possible during an interval, and (3) the exact moments of
transitions were not observed (i.e. observations are interval-censored)
[@van_den_hout_multi-state_2016; @logofet_mathematics_2000].

In a four-state transition model in continuous time, the Markov process is
governed by a 4 $\times$ 4 transition intensity matrix, $Q$, where rows are the
current states and columns are the future states (Fig. 2b). For each state $r,s
\in {B, M, P, T}$, the transition intensity ($q_{rs}$) represents the
instantaneous risk that a plot transitions from state $r$ to state $s$.
Because the states were defined based on stand basal area, instantaneous
transitions from Boreal to Temperate ($q_{BT}$) or from Temperate to Boreal
($q_{TB}$) were impossible without disturbance; there is a necessary transition
through Mixed or Pioneer. For this reason and the fact that these transitions
were very rare in the data, we fixed $q_{BT}$ and $q_{TB}$ at 0 (Fig. 2b).
However, all states can transition directly to Pioneer when disturbed (Fig. 2).

The intensities $q_{r,s}$ can be modelled as follows:

$$q_{rs}(t|x(t)) = q_{rs.0}(t)exp(\beta_{rs}'x(t))\text{,}$$

where $x(t)$ is the matrix of explanatory variables (surveys as rows, covariates
as columns), $\beta_{rs}$ are coefficients to be estimated, and $q_{rs.0}(t)$ is
a baseline hazard that describes the risk when environment $x(t) = 0$. Hence,
$exp(\beta_{rs}'x(t))$ is the relative increase or decrease in risk associated
with a set of characteristics $x(t)$. In this model, time-dependent variables,
such as climate and disturbances, are assumed to be piecewise-constant, i.e.,
the hazard is constant within a time interval $[t, t + \Delta{t}]$ and depends
on the variable value at $t$, but can change between the intervals. The
inclusion of time-dependent variables in the model allows one to fit a
non-homogeneous Markov process. Estimation of model parameters were obtained by
maximising the log-likelihood (see Supplementary Methods for details).

We built five models: one baseline model that solely includes the $q_{rs.0}$,
one model for each category of covariates independently (climate, soil and
disturbances), and one full model, which combines all covariates (Table 1).
Because multiple state transitions are estimated in a single model (all $q_{rs}$
in Fig. 2b), the number of parameters increases rapidly with the number of
covariates (number of modelled transitions (here 10) $\times$ (number of
covariates + 1)). Thus, to reduce the number of parameters, we assumed that
transitions from any state to Pioneer were only determined by disturbances,
while climate and soil variables should not directly influence these
transitions. All quantitative variables were standardised ($\mu$ = 0, $\sigma$ =
1) prior to running the models.

### Model evaluation

We first evaluated the goodness-of-fit of each model containing covariates
(climate, soil, disturbances and full) against the baseline model using
likelihood ratio tests [@jackson_multi-state_2011], which evaluate if the
addition of one or more new parameters significantly increases the likelihood of
the model. We also compared and ranked the models using the Akaike information
criterion [AIC; @burnham_model_2002]. The model with the lowest AIC was
considered to be the best model and thus used in further analyses.


### Model baseline and hazard ratios

We first evaluated the trends in recent forest transition dynamics. We used the
baseline hazards ($q_{rs.0}$) estimated by our best model as indicators of the
underlying forest response. For each pair of states, the baseline hazard
describes the risk to make a transition for a mean forest plot (when all
covariates are set to 0). We then investigated how  environmental covariates
influenced the transition dynamics (question 1) by comparing the estimated
hazard ratios derived from our best model ($exp(\beta_{rs})$).

### Transient dynamics and equilibrium

We further investigated how disturbances modify the long-term equilibrium
(question 2) and the forest transient dynamics (question 3). We computed
different properties on the Markov transition matrix along the latitudinal
temperature gradient and compared them among five disturbance scenarios
defined by disturbance type and severity: (1) **no or minor** disturbances, when the
covariates logging and natural were both fixed at 0; (2) moderate natural, with
the covariate natural fixed at 1 and logging fixed at 0; and (3) vice versa for
moderate logging; (4) major natural, with the covariate natural fixed at 2 and
logging fixed at 0; and (5) vice versa for major logging. The temperature
covariate was also allowed to vary from its lower 10th to its upper 90th
percentile, whereas all other covariates were fixed at the average conditions
found in the ecotone, the balsam fir-yellow birch domain (Fig. 1), to focus
solely on the effect of disturbances along the temperature gradient.

An extensive literature describes the multiple properties of discrete-time
Markov transition matrices [@caswell_matrix_2008; @hill_markov_2004] which can
be adapted to continuous-time models. We chose four informative and
complementary properties that fully characterise both the short and long-time
scale dynamics of our modelled system: (1) the steady-state distribution, which
corresponds to the potential long-term proportion of forest states at
equilibrium; (2) the half-life to equilibrium, which evaluates the time of
convergence to the steady-state and the length of the transient period; (3) the
turnover time, which measures the rate of transient successional changes; and
(4) the entropy, which captures the uncertainty regarding the next transitions.
While their absolute values should be interpreted with caution, their comparison
under various disturbance scenarios can highlight essential features of the
dynamics.

First, to measure the potential direction of forest dynamics under a given
scenario, we estimated the steady-state distribution, $\pi$. For a regular
Markov process, any initial state distribution converges to the same equilibrium
as time approaches infinity. The vector of equilibrium $\pi$ can be obtained by
taking the left eigenvector of the intensity matrix $Q$, which has an
eigenvalue of 0, normalised to sum to 1, or the left eigenvector of the
transition probability matrix $A$, which has an eigenvalue of 1, normalised to
sum to 1 [@norris_1997].

Then, the convergence rate to the equilibrium distribution can be measured
using the damping ratio [@hill_markov_2004]:

$$\rho = \lambda_{A1} / \lambda_{A2} = exp(\lambda_{Q1} - \lambda_{Q2})\text{,}$$

where $\lambda_{A1}$ and $\lambda_{A2}$ are the largest and second-largest
eigenvalues of $A$ ($\lambda_{A1}$ = 1 for stochastic $A$), whereas
$\lambda_{Q1}$ and $\lambda_{Q2}$ are the largest and second-largest eigenvalues
of $Q$ ($\lambda_{Q1}$ = 0 for stochastic $Q$). The convergence time was
approximated using the half-life to equilibrium:

$$t_{1/2} = log(2)/log(\rho)\text{.}$$

We also measured the turnover time in each forest state, also called the
sojourn time in multi-state models, which corresponds to the time spent in one
state before transitioning to a different state. The turnover time can be estimated by
$Turnover_{r} = −1/q_{rr}$, where $q_{rr}$ is the $r^{th}$ entry on the
diagonal of the estimated $Q$ matrix. The turnover of the whole system is
given by the average of each state turnover time over the steady-state
distribution:

$$Turnover = \displaystyle -\sum_r{\pi_r \times Turnover_{r}}\text{.}$$

Finally, @hill_markov_2004 proposed to use the entropy of a discrete-time
transition matrix as an index of the incertitude of successional changes. It
measures how uncertain we are about the next new state of a site knowing its
current state. For a continuous-time process, the entropy can be measured using
the jump matrix [@spencer_continuous-time_2005], which contains the
probabilities that the next state after state $r$ is state $s$:

$$j_{rs} = −q_{rs}/q_{rr}\text{.}$$

The entropy of state $s$ is then:

$$H(j_{.s}) = \displaystyle -\sum_r{j_{rs} \times log(j_{rs})}\text{.}$$

The normalised entropy of the whole system is the average of the entropies
over the steady state, divided by $H_{max} = log(n_{state} = 4)$:

$$\text{Entropy} = \displaystyle \frac{-\sum_r{\pi_r \times H(j_{.s})}}{H_{max}}\text{.}$$

Values of entropy closer to zero indicate more deterministic transition
dynamics whereas values closer to one indicate more random dynamics.


### Demographic processes

We finally decomposed the transition dynamics into its underlying demographic
components (question 4) for the most abundant species (i.e., three temperate,
*Acer rubrum*, *Acer saccharum* and *Betula alleghaniensis*; two boreal, *Abies
balsamea* and *Picea mariana*; two pioneer, *Betula papyrifera* and *Populus
tremuloides*). The transitions between states can result from various
combinations of increases in basal area through tree recruitment and growth and
decreases in basal area through mortality and logging. We measured recruitment
as the increase in basal area from new trees that had reached or exceeded the
threshold diameter of 9.1 cm. Growth was measured as the increase in tree basal
area between consecutive surveys. During the surveys, tree vitality was
characterised. We used this information to separate mortality as either due to
harvesting or to any other causes and measured the loss in basal area that
resulted from each of these two mortality processes.

Next, we used an indicator value analysis to quantify the contribution of each
demographic process and species to each of the sixteen forest transitions
[@dufrene_species_1997]. The indicator value ($IV_{jk}$), which measures the
exclusiveness of a process $j$ to a transition $k$, is given by the product of
the relative abundance (specificity; $RA_{jk}$) and the relative frequency
(fidelity; $RF_{jk}$):

$IV_{jk} = 100 \times RA_{jk} \times RF_{jk}\text{.}$


All analyses were performed using the R programming language version 3.6.1
[@r_core_team_r_2019]. The list of R packages that were used to carry out the
analyses is provided in the Supporting Information (Table S5). All data used
in the study, in addition to R scripts that reproduced the analyses and
figures, will be made available online on GitHub upon manuscript acceptance.


# Results

During the study period, a large fraction of Mixed forests transitioned to
Temperate forests (20.5%) but few did the opposite (6.3%). There were many
transitions from Boreal to Pioneer (13.0%), and more from Pioneer to Boreal
(19.3%). Temperate and Boreal forests were generally more persistent (90.3 and
84.9%, respectively, did not transition during the study period) than Mixed and
Pioneer forests (69.2 and 72.8%, respectively; Fig. 2a).

Overall, the full model, which includes climate, soil and disturbance variables,
had the best fit and predictive performances for the observed data (Table 2;
Fig. S5). The second-best model was the disturbance model, but it was far behind
with a difference in AIC of almost 1500 units from the full model (Table 2).
All variable subsets improved significantly the likelihood of the model (all
likelihood ratio tests were highly significant, *p* << 0.001; Table 2). Model
performance was also evaluated using a 10-fold cross-validation (see
Supplementary Methods); it revealed that including climate and disturbances
improved overall model predictive performance, while soil variables had a
negligible effect (Fig. S5). Thereafter, all inferences about transition
probabilities were derived from the full model.


## Baseline transition intensities

The baseline transition intensities of the full model provide insights about the
background rate of forest changes (Fig. 3). Forest dynamics over the whole study
area was largely dominated by transitions from Pioneer to Boreal (q~PB~ =
0.0270) and from Mixed to Temperate (q~MT~ = 0.0229; Fig. 3). Mixed
forests were 1.6 times (q~MT~ / q~TM~) more likely to transition to Temperate
than the reverse, indicating that temperate species had been successfully
colonising mixedwoods, outcompeting boreal species, during the study period. For
Boreal forests, regeneration from Pioneer to Boreal was 3.9 times (q~PB~ /
q~BP~) more likely than transition from Boreal to Pioneer.


## Effect of covariates on transition probabilities

The full multi-state model indicates that forest state transitions are
contingent upon environmental covariates. All transitions to Pioneer were highly
influenced by disturbances (Fig. 4, Table S3). As could be expected, major
disturbances exert stronger effects than moderate disturbances (for both natural
and logging), but, for each level of severity, logging had stronger effects
than natural disturbances. For example, the risk of transition from Boreal to
Pioneer has surged up to 213 times higher for plots that suffered major logging
(logging 2) and 37 times higher for plots that suffered major natural
disturbances (natural 2) compared to undisturbed plots (minor). Disturbances
of all types and severities favoured transitions from Mixed to Temperate
forests. Moderate disturbances (natural and logging) doubled the risk of this
type of transition, whereas major disturbances increased it by ca. 5 times
(Hazard Ratio (HR) = 5.76 and 5.32, for natural and logging, respectively).
Although the effect of major disturbances on the instantaneous risk of
transition from Mixed to Temperate was stronger than for moderate disturbances,
the probability of this event decreased with time (Fig. S6). Moderate
disturbances also favoured transitions from Boreal to Mixed (HR = 2.76 and 3.45,
respectively), while major disturbances had no significant effect on this type
of transition. Overall, the effects of disturbances are well reflected by the
radical change of structure of the 10-year transition probability matrix (Fig.
S7).

Climate variables also had a significant influence on most transitions (Fig. 4).
Warmer summer temperature (higher temperature) and higher humidity (higher CMI)
favoured transitions from Boreal to Mixed as well as from Pioneer to Mixed and
Pioneer to Temperate. Interestingly, warmer temperature did not significantly
influence the risk of transition from Mixed to Temperate and higher CMI had a
negative effect.

State transitions were also influenced by soil variables (Fig. 4, Table S3).
Holding the other covariates constant, the instantaneous risk of transition from
Boreal to Mixed and from Pioneer to Temperate decreased by 27% and 23%,
respectively, on poorer drainage (more humid), but the risk of transition from
Temperate to Mixed increased by 30% (HR = 0.73, 0.77 and 1.30, respectively).
Higher pH (acidic soil) had a negative effect on the transitions from Temperate
to Mixed (HR = 0.73). These changes in risk ratios associated to soil variables
appear almost irrelevant compared to the effect of disturbances, but a slight
increase in drainage can dampen the positive effect of disturbances. For
instance, under moderate natural disturbances, the instantaneous risk of
transition from Boreal to Mixed is 0.007 at moderate drainage but decreases to
0.003 when increasing drainage by 1 point.


## Effect of disturbances on long-term equilibrium

The potential state proportion at equilibrium was strongly influenced by
disturbances (Fig. 5a). For the undisturbed scenario (minor), the predicted
equilibrium at the ecotone was relatively close to the initial observed
proportions, with signs of regeneration from Pioneer to Boreal states and slight
increases in Mixed and Temperate states. The steady-state proportion of
Temperate almost doubled with moderate disturbances (minor: 33%; moderate
natural: 56%; moderate logging: 60%), while the boreal state was more than
halved. At major disturbances, Pioneer forests dominated the equilibrium
landscape, while the other states collapsed.  

The steady-state proportion also changed as expected along the temperature
gradient (Fig. 5b,c). The Boreal state dominates at low temperature (high
latitude) and the Temperate state dominates at high temperature (low latitude),
highlighting the position of the boundary between these two biomes at a growing
season temperature of about 12.9°C, which is found in the actual ecotone.
Moderate disturbances (both natural and logging) displaced the temperate-boreal
boundary at lower temperatures (ca. 12.2°C), hence further north of the current
ecotone (Fig. 5b,c). Because of the dominance of the Pioneer state, the boundary
modestly moved north with major natural disturbances (12.7°C), while it
retreated to the south with major logging (13.4°C).


## Effect of disturbances on transient dynamics

Disturbances affected forest transient dynamics with greater impact for higher
disturbance severity (Fig. 6). In the minor disturbance scenario, turnover
time was generally longer at low temperature, indicating slower transition
dynamics in northern forests (Fig. 6a,b). The turnover time then rapidly
declined to reach a minimum at ca. 13.25°C, at the southern limit of the
ecotone, and went back up after this point. This trough, where transition
dynamics is the fastest, is located just a little south of the boundary between
Boreal and Temperate dominances found in Figure 5. Major disturbances
accelerated transition dynamics all along the temperature gradient, while
moderate disturbances also decreased turnover time but more strongly in the
northern boreal region (Fig. 6a,b). The effect on turnover time was similar for
both disturbance types, except that the effect of major logging was much
stronger in northern boreal forests than natural disturbances (Fig. 6a,b). These
spatial patterns reflect the turnover time of the dominant state at each point
along the temperature gradient (Fig. S8).

At minor disturbances, the entropy of the system generally increased from north
to south and peaked at ca. 12.6°C, at the northern end of the ecotone (Fig.
6c,d). This peak illustrates where the transition dynamics is most uncertain
(transition to all states are possible at this point), while it is very
predictable in northern boreal forests (Boreal stays Boreal until it transitions
to Pioneer later on). The peak can be mainly attributed to the entropy of the
Boreal state at the ecotone, and the generally high values at low latitudes can
be principally attributed to the Temperate state (Fig. S9). This latitudinal
pattern of entropy is modified by disturbances. Moderate natural disturbances
decreased the entropy throughout the gradient, but especially where the peak is
found (Fig. 6c). With moderate logging, the peak disappeared, and entropy
increased monotonically from north to south (Fig. 6d). The peak of entropy was
displaced to the south when major disturbances were included, whether natural or
logging (Fig. 6c,d), where it was dominated by the entropy of the Pioneer state
(Fig. S9).

Half-life to equilibrium was the longest at ca. 11.8°C, north of the ecotone,
in the balsam fir-white birch domain, while it was the shortest in the
southernmost latitudes (Fig. 6e,f). Moderate disturbances flattened and shifted
this peak to the north and the effect of moderate logging (Fig. 6f) was stronger
than natural disturbances (Fig. 6e). In the balsam fir-white birch, the
half-life to reach equilibrium distribution was reduced almost by half by
moderate logging. With major disturbances, forests all along the temperature
gradient can reach very quickly their steady-state distribution (maximum of
about 8 years for major logging and 25 years for major natural disturbances).

## Contribution of demographic processes

Only the demographic processes of a few species contributed substantially to the
observed transition dynamics (Fig. 7). The importance of some processes were
expected. For example, transitions from Boreal to Pioneer were dominated by
mortality and logging of *Picea mariana*, while the transitions from Pioneer to
Boreal were characterised by recruitment and growth of *Picea mariana* and
*Abies balsamea*. Most interestingly, the transitions from Mixed to Temperate
were determined by the mortality of *Abies balsamea* and the growth of temperate
species, mainly *Acer rubrum* and *Betula alleghaniensis*, and to a lesser
extent *Acer saccharum*. The recruitment of temperate species was not indicator
of the Mixed to Temperate transitions, but rather of the transitions from
Pioneer to Temperate.


# Discussion

Our study reveals that **forest transition dynamics in the temperate-boreal
ecotone was predominantly controlled by natural and anthropogenic disturbances
and secondarily by climate, whereas local soil conditions exerted relatively
minor constraints. While major disturbances only promoted transitions to the
pioneer state, moderate disturbances increased the probability of transition
from mixed to temperate states.** Our analysis of the equilibrium further
highlights that the long-term forest dynamics under moderate disturbances
favours an increased proportion of temperate forests and thereby a northward
shift of the temperate-boreal ecotone. Disturbances also modify the forest
transient dynamics, accelerating both the turnover and convergence time and
making the dynamics more predictable. **Contrary to our expectation, transitions
from mixed to temperate forests were not driven by recruitment but mostly
mortality and growth.** In accordance with the hypothesis formulated in previous
studies [@johnstone_changes_2010; @johnstone_changing_2016;
@vissault_biogeographie_2016; @brice_disturbances_2019], our findings show that
moderate disturbances catalyse transitions to the alternate, temperate-dominated
forest state and could therefore promote regime shifts. **Moreover, our results
emphasise that forest dynamics are affected by multiple factors operating across
different spatial and temporal scales and that predicting range shifts during
climate change will thus require approaches that integrate multi-scale patterns
and processes [@allen_interactions_2007].**


## Trends in recent forest transition dynamics in Québec

Forest dynamics in Québec during the last 48 years was dominated by
transitions from pioneer to boreal and from mixed to temperate stands. The
important regeneration of boreal forests could be attributed to past natural
disturbances, notably the last spruce budworm outbreak. Indeed, the last
outbreak, which occurred during the 1970s, has caused major mortality in
coniferous species followed by important recruitment pulses and growth releases
[@bouchard_tree_2006].

Although we did not directly evaluate the impact of climate change, our results
suggest that recent climate warming may contribute to the forest transition
dynamics. The high baseline transition rate from mixed to temperate is
consistent with the expectation of a northward range shift of temperate trees
into the mixed and boreal forests. In our study, these transitions were caused
by the concomitant high mortality of an abundant boreal species, *Abies
balsamea*, and the increased growth of temperate species. Accordingly, the
warming trend of the last decades (Fig. S3) has been shown to increase growth
and reproductive rates of temperate species at their northern limit
[@reich_geographic_2015; @fisichelli_temperate_2014;
@boisvertmarsh_divergent_2019; @goldblum_tree_2005; @bolte_climate_2010], thus
providing a competitive advantage to temperate over boreal species.  

The increased transition rate to temperate forests is likely also a response to
historical disturbances and climate change. Comparisons of pre-settlement and
present-day forested landscapes of North America have highlighted an important
deciduous encroachment in response to historical human activities
[@danneyrolles_stronger_2019; @terrail_reorganization_2019;
@boucher_logging-induced_2006]. Moreover, evidence suggests that the ongoing
northern range expansion of some tree species is the result of delayed postglacial
migration [@svenning_disequilibrium_2013]. Historical legacies and recent
climate change are presumably mutually non-exclusive explanations. Indeed,
simulations by @boulanger_climate_2019 showed that the future climate-induced
expansion in temperate species to the detriment of boreal species would amplify
the already ongoing trend since preindustrial times.


## Disturbances catalyse forest state transition

Our study highlighted that moderate disturbances favour Mixed to Temperate
transitions following climate warming, whereas major disturbances merely promote pioneer states.
Disturbances directly remove trees, which leads to immediate and substantial
changes in forest composition [@brice_disturbances_2019]. Forests are expected
to be resilient to normally experienced disturbances and should thus return to
their preceding states [@gunderson_ecological_2000]. However, climate change
alters the conditions that initially supported the persistence of a given state,
making forests susceptible to transition to other states
[@johnstone_changing_2016].

**Following a disturbance, three mechanisms can contribute to the observed change in tree
cover: (1) the loss of a dominant species; (2) the growth release of advanced
regeneration of co-occurring species; and (3) the pulse establishment of new
species. Our results show that the first two mechanisms may operate
simultaneously, while the third had a limited influence.**

In the study area, both natural and anthropogenic disturbances
disproportionately affected *Abies balsamea*, which has suffered significant
mortality due to spruce budworm outbreaks and was also intensively harvested
[@duchesne_population_2008]. The canopy gaps created by the loss of this
ubiquitous and abundant boreal species probably allowed for the growth release
of co-occurring temperate species. These findings are in line with a study in
the temperate-boreal ecotone of Scandinavia where a boreal tree, *Picea abies*,
was particularly affected by a drought and an insect outbreak which then favour
the growth of a temperate species, *Fagus sylvatica* [@bolte_understory_2014].
**Combined effects of selective disturbances and climate warming may thus initiate
a shift in the competitive balance between boreal and temperate species
[@grundmann_impact_2011; @bolte_understory_2014].** We only found a weak
contribution of temperate tree recruitment to the Mixed to Temperate
transitions, likely because our analyses were based on tree basal area. However,
other studies analysing abundance data suggest that moderate disturbances may
also facilitate colonisation and establishment by opportunistic temperate
species under warmer conditions [@brice_disturbances_2019;
@leithead_northward_2010; @landhausser_disturbance_2010]. Moreover, it is
possible that, in the long run, the increased proportion of temperate species in
forest communities could alter soil properties and ultimately facilitate the
recruitment of even more temperate species.

In contrast to moderate disturbances, severe disturbances, primarily
clearcutting but also large fires in the study area (Fig. S1), may result in
large forest dieback and create openings of very large extent. These newly
opened landscapes can be colonised swiftly by early-successional species that
benefit from a long-distance seed dispersal and a fast growth, such as *Populus*
and *Betula* [@landhausser_disturbance_2010]. In contrast, temperate species
may be slower to come back following major disturbances because they dispersed
over shorter distances [maximum of ca. 200m for *Acer* compared to 5000m for
*Populus*; @boulanger_climate_2017]. Due to the increase in large-scale logging
during the last century, the proportion of young recently disturbed forests have
been found to have increased in North America [@boucher_logging-induced_2006;
@danneyrolles_anthropogenic_2018; @thompson_four_2013]. The expected increase in
frequency and severity of climate-induced disturbances in combination with
clearcuts may further promote the expansion of young pioneer forests in the
future.

Compared to the catalysing effect of disturbances, **local** soil
characteristics do not appear to represent a large impediment to state
transitions, but transitions may be slower on some soil types. Poor drainage
constrained climate-related transitions from Boreal to Mixed states, but not
from Mixed to Temperate. This indicates that temperate species can readily
colonise soils found in mixedwoods but may have more difficulty in colonising
hydric boreal soils. **Thus, local soils may be important to explain the low
transition rate from Boreal to Mixed.** Very poor drainage, often associated
with peatland and thick organic layer, is usually thought to be improper for the
regeneration of temperate species [@lafleur_response_2010]. Several studies
found that *Acer saccharum* regenerates well across the ecotone because of its
large tolerance to various soil conditions [@barras_supply_1998;
@goldblum_age_2002; @kellman_sugar_2004; @fisichelli_temperate_2014;
@collin_can_2018]. At their northern range limit, *A. saccharum* and *A.
rubrum*, the species contributing most to compositional changes in Québec
[@brice_disturbances_2019], are hypothesised to be mostly limited by cold soil
temperature [@barras_supply_1998; @goldblum_age_2002].


Moreover, disturbances may counteract any effect of soil properties. Indeed,
disturbances, such as logging and fire, often remove the surface organic layers
and expose mineral soil. They can, consequently, provide an appropriate seedbed
for temperate species recruitment [@archambault_fifty_2006;
@landhausser_disturbance_2010]. In combination with climate warming,
disturbances may also facilitate temperate migration by increasing understory
air and soil temperatures [@stevens_forest_2015; @de_frenne_microclimate_2013].


## Changes in potential long-term equilibrium and biome boundary

Our model highlights the potential role of disturbances in influencing the
position of the temperate-boreal boundary as well as the proportion of temperate
and boreal biomes at equilibrium. As a result of the increased replacement of
Mixed by Temperate states and a decline of Boreal to Pioneer states, the
equilibrium temperate-boreal boundary shifts northward with moderate
disturbances. While our results should not be interpreted as projections for the
future, they are useful to highlight the direction of forest dynamics under
different disturbance scenarios **and underscore that short-term changes in the
transition probabilities can impact long-term regional forest patterns**. Our
findings also support the simulations of @boulanger_climate_2019 where
harvesting under future climate warming was projected to promote further
invasions of pioneer species, such as *Populus*, and temperate species, such as
*Acer* and *Fagus*, in mixedwoods of Québec.


Based on their simulations, @liang_how_2018 and @vanderwel_how_2014 concluded
that logging would primarily accelerate the expansion of pioneer forests but
would have little or no effect on extensive biome shifts over the next century
in eastern United States. In contrast to their results, we found a clear range
shift of the Temperate state under moderate disturbances, whereas the Pioneer state
would have the advantage and become dominant at equilibrium only under major
disturbances. We hypothesise that the northern shift of the Temperate state
induced by moderate disturbances was mainly the result of the increased
dominance of temperate species in areas where they are already present. Indeed,
the current disturbance regime in our study area contributed to the decline of
one boreal species in particular, *Abies balsamea*, which in turn likely
benefited the growth of co-occurring temperate species. Moreover, because of its
positive response to past [@danneyrolles_stronger_2019], recent
[@brice_disturbances_2019] and future [@boulanger_climate_2019] disturbances in
Québec, *Acer rubrum* is likely to play a disproportionate role in the temperate
biome shift. However, the low probability of transition from Boreal to Mixed
suggests, like other studies, that migration of temperate trees into pure boreal forest will be
a much slower process [@vissault_biogeographie_2016; @solarik_priority_2019].

## Disturbances accelerate the transient dynamics

Beyond their impacts on the equilibrium, our results suggest that disturbances
may have a substantial influence on forest transient dynamics. **Disturbances
generally increased the rate of tree species replacement (reduced turnover time)
and induced a convergence of the dynamics (reduced entropy), thereby
accelerating transition dynamics toward a new equilibrium (reduced half-life; Fig. 6).
While disturbances are known to accelerate stand scale forest succession
[@abrams_disturbance-mediated_1989; @bolte_understory_2014], here we provided
evidence that their effects could translate to an acceleration of broad-scale
biome shifts.**

In the continuous boreal zone (spruce-moss domain), forests dominated by *Picea
mariana* are usually characterised by dynamics of stand self-replacement with
minimal compositional changes across disturbance cycles
[@goldblum_deciduous_2010]. Consistent with this dynamics, the turnover time of
undisturbed northern boreal forests was very long and the entropy very low in
our results. The turnover was shortened by disturbances, but the entropy
remained low, indicating that the dynamics was still very predictable (back and
forth transitions between Boreal and Pioneer states) and that there was no
directional shift associated with disturbances. Hence, boreal forests lose their
persistence when moderately disturbed but remain resilient as they return to
their previous boreal state. Under major disturbances, boreal forests collapsed
to Pioneer state and reached this new equilibrium swiftly (short half-life).
This observation is consistent with previous studies suggesting that boreal
forests can easily shift into an alternative treeless state in response to
severe or repeated disturbances [@payette_shift_2003;
@sanchez-pinillos_resistance_2019].

In contrast, the ecotone is characterised by a rapid turnover and high entropy
indicating abrupt compositional shift which can go in any direction. Compared
to northern boreal forests, the short turnover time implies a low persistence of
the forest states in this region even under minor disturbances. This result
corroborates the predictions made by @vissault_biogeographie_2016, where mixed
forests would undergo a swift conversion to temperate forests in the next
decades, whereas boreal forests would present a large inertia presumably because of
dispersal limitation. The dynamics of the ecotone appears unstable because it
is caught between two stable states, i.e. Boreal to the north and Temperate to
the south. Under moderate disturbances, the probability of transitioning to
Temperate increases to the detriment of the other possible states, hence the
entropy is decreased, and the dynamics becomes more predictable. Such a clear
directional shift strongly indicates non-equilibrium dynamics in this region.
Although turnover is fast, half-life to equilibrium is long because a forest may
not move towards equilibrium and may undergo multiple transitions.


## Ecological and management implications

**A common assumption is that factors determining species distribution are
hierarchical, such that climate would govern the distribution at the regional
scale while soil conditions would be more important at the local scale
[@pearson_predicting_2003]. However, our study provides empirical evidence that,
through their effect on demography, landscape disturbances, and to a lesser
extent, local soil factors, may interact with global warming to influence
regional shifts in forest types. Specifically, natural and anthropogenic
disturbances cause widespread mortality of a dominant species, while climate
warming likely increased the growth of co-occurring temperate species, thus
altering post-disturbance successional trajectory
[@anderson-teixeira_altered_2013; @duveneck_measuring_2016].**

A shift in dominant forest cover from conifer to deciduous broadleaf species
entails large changes in tree species diversity and composition
[@berteaux_cc-bio_2010] that **can accumulate through time and space** and
induce a complete transformation of regional forest dynamics and functions
[@peters_crossscale_2007]. In the long term, this regime shift could increase
carbon sequestration [@thurner_carbon_2014], modify disturbance regimes [reduced
flammability of broadleaf species @terrier_potential_2013; and reduced
sensitivity to current outbreak-prone pest @mffp_insectes_2018], alter soil
microbial activity [@laganiere_how_2010] and affect wildlife distribution
[@mizel_rapidly_2016].

Such regime shifts will impact strongly on forest management strategies in area
where silvicultural practices are tailored to the regional disturbance regimes
and rely on natural regeneration. In Québec, ecosystem-based forest management
seeks to maintain the composition and structure of a reference state, defined as
the preindustrial forest conditions [@pinna_amenagement_2009]. Yet,
@boulanger_climate_2019 showed that such management would fail to restore
historical forest conditions under future climate change, and that disturbances
would only exacerbate the gap. While trying to maintain a historical state is
likely impractical, our results emphasise that forest management should not only
consider the present system state, but also its most likely trajectory. Our
study also reveals the potential of moderate disturbances to catalyse
regional forest transitions. This suggests that partial cutting could be
used to reduce migration lags. Other studies additionally recommend planting temperate trees farther
north outside their current range, i.e. assisted migration, to facilitate
**range expansion** [Vieira et al. in prep.; @duveneck_measuring_2016]. However,
before implementing such silvicultural strategies, key questions need to be
answered. For instance, will multiple interacting disturbances exacerbate tree
mortality? And how will these rapid transitions impact the ecosystem processes
and functions? The rate of recent climate change already outpaces tree migration
capacity [@sittaro_tree_2017], and even more so the scientific capacity to
understand and mitigate its consequences [@jackson_balancing_2010]. Therefore,
in order to insure long-term forest health in the temperate-boreal ecotone, we
need to simultaneously limit global warming through drastic reduction of
greenhouse gas emissions and intensify research effort to develop effective
adaptation strategies for forest management.



## Acknowledgements

The authors declare no conflict of interest.

We are grateful to Kevin Cazelles for providing helpful suggestions and comments
that improved our analyses and manuscript. We also thank Guillaume Guénard for
useful advices on the model. Our thanks to the anonymous reviewers for their
constructive comments on an previous version of the manuscript. We gratefully
acknowledge the staff of the Ministère des Forêts, de la Faune et des Parcs du
Québec (MFFP) for their work on forest inventories. This research was supported
by Natural Sciences and Engineering Research Council of Canada (NSERC) research
grant no. 7738 to P. L. and no. 5134 to M.‐J. F.

\pagebreak

# Tables

**Table 1**: Description of the explanatory variables used in the multi-state models.


|**Variable name**|**Variable description**                                   |
|:----------------|:----------------------------------------------------------|
|**Climate**      |                                                           |
|Temperature      |Mean temperature during growing season, 10-year average prior to plot measurement (°C). |
|CMI              |Mean Climate Moisture Index from May to September, 10-year average prior to plot measurement (cm). |
|**Soil**         |                                                           |
|pH               |pH of the surface horizon                                  |
|Drainage         |6 classes of soil drainage, which range from excessive to very poor, that were treated as numeric.|
|**Disturbances** |                                                           |
|Logging          |Tree harvesting, including clearcutting, selection cutting, shelterwood cutting, seed-tree cutting, etc. None or minor (0), moderate (1) or major (2). |
|Natural          |Natural disturbances, including forest fires, insect outbreaks, windfall, etc. No or minor (0), moderate (1) or major (2). |


\pagebreak

**Table 2**: Comparisons of the five candidate multi-state models. The number of
parameters used in each model corresponds to the number of modelled transitions
(10) $\times$ the number of covariates - 1. The $\Delta$AIC is the difference
between the Akaike information criterion of each model (AIC~m~) and the
minimum of AIC among all models (AIC~min~): $\Delta$AIC = AIC~m~ – AIC~min~.
Models are presented in decreasing order of their $\Delta$AIC. Each model containing
covariates was compared to the baseline model using a Likelihood Ratio (LR) test.
The best model is the one in bold with $\Delta$AIC = 0.

&nbsp;

\begin{tabular}{llrrrl}
\toprule
\textbf{ } & \textbf{Covariates} & \textbf{Number of parameters} & \textbf{-2 Log-likelihood} & \textbf{Delta AIC} & \textbf{LR test}\\
\midrule
Baseline & Intercept & 10 & 37874.4 & 8298.4 & ---\\
Soil & Drainage, pH & 24 & 37713.7 & 8165.7 & < 0.001\\
Climate & Temperature, CMI & 24 & 36288.8 & 6740.8 & < 0.001\\
Disturbances & Natural, Logging & 50 & 30993.5 & 1497.5 & < 0.001\\
\textbf{Full} & \textbf{All} & \textbf{78} & \textbf{29440.0} & \textbf{0.0} & \textbf{< 0.001}\\
\bottomrule
\end{tabular}

\pagebreak

# Figures

![Locations of the 11,058 forest inventory plots in meridional Québec, Canada. Colours delimit the six bioclimatic domains. The two southernmost domains (red) are here combined. The number of  plots in each domain is shown in parentheses. The balsam fir-yellow birch domain (in bold) is the ecotone between the hardwood and boreal forests.](res/fig1_region.pdf)

\pagebreak

![Multi-state transition diagram (a), intensity matrix Q (b) and equations of our
full model (c). Directional arrows in the diagram (a) depict the allowed transitions between
states. The numbers represent the percentage of observed transitions between
states (nb~rs~/nb~r.~ $\times$ 100). Instantaneous transition from Boreal to
Temperate and vice versa are considered impossible in the model (hence the
absence of arrows in the diagram and the zeros in the Q matrix), however rare
transitions from Boreal to Temperate and from Temperate to Boreal were observed
in the data (less than 0.2%). The Q matrix (b) contains the instantaneous risk
to move from one state (row) to another (column), here: (B)oreal, (M)ixed,
(P)ioneer and (T)emperate, in that order. Transitions from any other state to
Pioneer were modelled as only dependent on disturbances
(c).](res/fig2_trans_diagram.pdf)

\pagebreak


![Baseline transition intensities as estimated from the best multi-state transition model. Arrows depict the direction of transitions between states. The numbers represent the estimated baseline hazards ($q_{rs.0}$), i.e., the instantaneous risk of moving from one state to another when all covariates are set to 0 (i.e., the means of standardised covariates and disturbance level 0).](res/fig3_baseline.pdf)

\pagebreak

![Hazard ratios (HR) and 95% confidence intervals as estimated from the best multi-state transition model. Each plot represents the estimated HR for transitions from row to column state, e.g., the plot on the first row, second column shows the HR for the Boreal to Mixed transition. The ordinate is in log scale. The HR of predictors are interpretable as multiplicative effects on the hazard, where values above 1 (in blue) indicate that the predictor is associated with a greater risk of state transition, whereas values below 1 (in red) indicate a lower risk of transition. Predictors statistically different from 1 are shown in dark blue or red. Numbers following disturbance predictors indicate their levels of intensity: 1 = moderate and 2 = major.](res/fig4_HR.pdf)


\pagebreak

![Changes in forest state proportions at equilibrium for different disturbance types (natural or logging) and intensity (no or minor, moderate, major). The barplot (a) compares the observed state proportion in the ecotone to the potential state proportion at equilibrium for different disturbance scenarios with all other covariates fixed at the average conditions found in the ecotone. The curved lines (b,c) show the proportions of Boreal (blue) and Temperate forests (red) at equilibrium along the temperature (latitudinal) gradient for no or minor (solid), moderate (dashed) and major (dotted) disturbances, with all other covariates fixed at the average conditions found in the ecotone. The light (no or minor), medium (moderate) and dark (major) grey circles indicate the positions of the boundary between dominance of Boreal and Temperate forests (i.e. the advancing front) while the corresponding arrows show how moderate and major disturbances move the boundary. The polygon approximates the position of the ecotone along the temperature gradient.](res/fig5_steady.pdf)

\pagebreak

![Changes in the characteristics of the forest transient dynamics along the temperature (latitudinal) gradient for different disturbance scenarios: no or minor (solid), moderate (dashed) and major (dotted) disturbances for both natural (a,c,e) and logging (b,d,f). All other covariates are fixed at the average conditions found in the ecotone to focus solely on the effect of disturbances along the temperature gradient. The turnover of the whole system (i.e. whole transition matrix) (a,b) corresponds to the time spent in a state before transitioning to the next and is given by the average of each state turnover time over the steady-state distribution. The entropy of the whole system (c,d) corresponds to the uncertainty of the next transition and is given by the average of each state entropy over the steady-state distribution. The half-life to equilibrium (e,f) is the time taken to reach 50% of the steady-state distribution, i.e. when the first eigenvalue becomes twice as large as the contribution of the second eigenvalue. The polygon approximates the positions of the ecotone along the temperature gradient.](res/fig6_transients.pdf)

\pagebreak

![Species and demographic process contribution to all observed state transitions across the study area. Letters on the *x* axis correspond to the four forest states: (B)oreal, (M)ixed, (P)ioneer and (T)emperate. Each pair of letters denotes a transition from one state (first letter) to the next (second letter), e.g., BB is Boreal to Boreal. The darker colours indicate higher indicator value.](res/fig7_demo_trans.png)

\pagebreak

# References
