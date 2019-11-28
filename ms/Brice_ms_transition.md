---
title: Moderate disturbances accelerate forest transition dynamics under climate change at the temperate-boreal ecotone
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
factors, such as disturbances and unsuitable soil conditions.

Using decades of forest inventory data, we model the state transition dynamics
of Québec's forests to identify the environmental conditions that promote or
prevent forest transitions. We further investigate how different disturbance
types and intensities impact forests' short-term transient dynamics and
long-term equilibrium.

We analysed over 10,000 forest inventory plots, sampled from 1970 to 2018 in
meridional Québec, Canada. We used a continuous-time multi-state Markov model to
quantify the transition probabilities between forest states (temperate, boreal,
mixed, pioneer) in relation to climate, soil conditions and disturbances. We
described transient dynamics and equilibrium under different disturbance
scenarios, using properties of Markov transition matrices.

Whereas the larger portion of the forest plots remained in the same state, most
transitions were conversions from mixed to temperate stands, or regeneration
from pioneer to boreal forests. Transition probabilities were mainly driven by
natural and anthropogenic disturbances and secondarily by climate, whereas soils
exerted only minor constraints. Moderate disturbances not only increased the
probability of transition to pioneer, but also from mixed to temperate states.
At low disturbance levels, boreal forests were characterised by great inertia
and predictable stand self-replacement dynamics, while mixed forests presented
rapid and uncertain transitions. Disturbances reduced turnover and convergence
times for all transitions. Yet, while major disturbances only promoted Pioneer
states, moderate disturbances induced a clear directional shift toward temperate
forests. In the long term, these changes in the transition dynamics promoted a
northward shift of the temperate forests at equilibrium. Hence, moderate
disturbances could catalyse rapid forest transitions and accelerate broad-scale
biome shifts.


## Keywords

Climate change,
Natural disturbances and logging,
Continuous-time multi-state Markov model,
Québec,
Temperate-boreal ecotone,
Transition probabilities,
Equilibrium,
Transient dynamics


\pagebreak

# Introduction

Global climate warming has led several temperate deciduous tree species to
slowly migrate northward, colonising conifer-dominated forests
[@fisichelli_temperate_2014; @evans_borealtemperate_2017;
@boisvert-marsh_shifting_2014; @sittaro_tree_2017]. As climate warms up and tips
the balance in favour of temperate over boreal species, forests at the ecotone
are expected to transition from coniferous to mixedwood and from mixedwood to
temperate deciduous [@boulanger_climate_2019; @price_anticipating_2013]. While
boreal forest dynamics is characterised by broad-scale disturbances (fires and
insect outbreaks) and slow decomposition of its acidic and nutrient-poor litter,
temperate forest dynamics is characterised by small-scale canopy gaps and rapid
decomposition of a rich litter [@goldblum_deciduous_2010]. As ecological
processes strongly differ among these biomes, climate-induced range shifts not
only impact species distributions, but also alter the structure of communities,
microclimates, biogeochemical cycles and ecosystem functioning and might trigger
a "regime shift" [@scheffer_catastrophic_2001].

Climatic niches for tree species are expected to shift northward by several
hundred kilometres by the end of the century [@mckenney_potential_2007].
However, many studies indicate that tree migration may not keep pace with global
warming [e.g. @zhu_failure_2012; @sittaro_tree_2017; @renwick_temporal_2015;
@talluto_extinction_2017; @vissault_biogeographie_2016]. Indeed, warmer climate
may improve recruitment, survival and growth of temperate species at their
northern range limit [@goldblum_tree_2005; @boisvertmarsh_divergent_2019;
@graignic_geographical_2014], while reducing growth and increasing mortality of
boreal species at their southern range limit [@goldblum_tree_2005;
@peng_drought-induced_2011]. These demographic changes alter colonisation and
extinction rates, which are driving recent species range shifts
[@talluto_extinction_2017] and should continue until forests reach a new
equilibrium with climate, i.e., the steady-state or long-term proportion of
forest states. However, because trees are long-lived species that disperse over
very short distances, colonisation and extinction events in response to
environmental changes are often delayed, and forests rarely reach their
equilibrium [@talluto_extinction_2017]. If forests are undisturbed, transition
rates between forest types following natural succession pathways will be mainly
limited by the persistence and turnover of resident species as well as the
dispersal and establishment rates of migrating species
[@neilson_transient_1993], resulting in transient dynamics that may last a very
long time [@hastings_transient_2018; @jackson_balancing_2010;
@talluto_extinction_2017]. While estimating the equilibrium proportion of forest
states provides insight into the potential long-term attractor of forest
dynamics under given environmental conditions, the transient dynamics depicts
the trajectory to reach it [@hastings_transient_2018]. The transient dynamics is
therefore more relevant to the changes likely to unfold during the 21st century
and, consequently, to the design of management strategies. Knowledge of
equilibrium and transient dynamics can thus expand our understanding of forest
responses at different time scales.


Given that forests are increasingly subjected to pervasive climatic stresses and
direct human disturbances, abrupt transitions are likely to play a key role in
driving the climate shift in biomes. Indeed, as climate change slowly modifies
the competitive balance among species, disturbances destroy the resident
community in whole or in part, thus providing establishment opportunities for
migrating species and making resources available for a fast growth.
Consequently, forest composition may shift to species that are better suited to
current conditions and fail to return to its previous state following a
disturbance [@johnstone_changing_2016; @renwick_temporal_2015;
@turner_disturbance_2010]. For example, canopy gaps have been shown to locally
facilitate establishment of temperate species in mixed forests of Ontario
[@leithead_northward_2010]. In Alaska, white spruce (*Picea glauca*) is invading
black spruce (*Picea mariana*) stands following fire and permafrost degradation
[@wirth_white_2008]. Similarly, moderate disturbances favoured the increase of
warm-adapted species and have led to a broad-scale community thermophilization of
forests in Québec [@brice_disturbances_2019].

While these examples suggest that disturbances have the potential to catalyse
shifts to an increasingly deciduous-dominated landscape, other simulation
studies have concluded that they are unlikely to drive extensive biome shifts in
the coming decades [@vanderwel_how_2014; @liang_how_2018]. In some cases,
disturbances have been shown to promote colonisation by early successional
species, which can then displace long-lived shade-tolerant species. For instance,
clearcutting has been found to favour the expansion of a pioneer species, the
trembling aspen (*Populus tremuloides*), in mixed and boreal stands of Québec
[@laquerre_augmentation_2009; @grondin_have_2018] and Alberta
[@landhausser_disturbance_2010]. These divergent results suggest that the effect
of disturbances on forest dynamics may depend on their intensity and type
(natural or anthropogenic). Therefore, more empirical evidence is essential to
disentangle the role of various intensities and types of disturbances in
broad-scale biome shifts.

The northward migration of temperate species may nevertheless be contingent on
their capacity to colonise different types of soil [@lafleur_response_2010;
@bennett_plant-soil_2017; @brown_non-climatic_2014]. Indeed, soils of cold
boreal forests generally have lower pH, lower microbial activity and slower
decomposition rates of organic matter than warmer southern temperate forest
soils [@goldblum_deciduous_2010]. These local and regional variations in soil
properties are expected to slow down or inhibit the establishment of temperate
trees into the boreal forest. For instance, transplant experimental studies have
shown that seedlings of sugar maple (*Acer saccharum*) in conifer-dominated
stands were negatively affected by seed predators and fungal pathogens
[@brown_non-climatic_2014] as well as by soil acidity through reduced foliar
nutrition [@collin_conifer_2017]. Yet, @kellman_sugar_2004 found that, after
high initial mortality due to seed predation, survival of *Acer saccharum*
seedlings in boreal stands was high, even superior to survival in temperate
stands, potentially because of increased light availability. Hence, it has been
suggested that soil properties in boreal forests may not be a major impediment
to the migration of temperate species showing broad ecological tolerance
[@lafleur_response_2010; @barras_supply_1998; @kellman_sugar_2004]. Nonetheless,
suboptimal soil conditions under a boreal canopy could delay forest transitions
under climate change [@solarik_priority_2019]. While experimental studies
provide valuable insights on the potential role of soils at local scales, we
need to test the generality of such constraints, or the lack thereof, on
long-term forest dynamics at regional scale and across species to better
anticipate future biome transitions.


One approach to investigating biome shifts in response to climate change is to
model transitions of forest plots among states as a stochastic process
influenced by their current state, as well as their current environmental
characteristics. Given the unequivocal distinction between temperate and boreal
forests, the dynamics of tree communities at the temperate-boreal ecotone can be
adequately characterised using discrete functional and successional states,
namely boreal, mixed, temperate and pioneer [@vissault_biogeographie_2016] and
thus can be formalised as a multi-state Markov model
[@jackson_multi-state_2018]. Markov models provide a useful framework for
modelling changes of state over time using longitudinal data. In epidemiology,
for example, these models are often used to describe the progression of diseases
[@van_den_hout_multi-state_2016]. In ecology, they have been used to study
processes such as forest succession [@runkle_gap_1981;
@waggoner_transition_1970; @lienard_data-intensive_2015], metapopulation
dynamics [@hanski_metapopulation_2003; @moilanen_patch_1999], landcover changes
[@yang_land_2014; @muller_markov_1994], or stage class transitions
[@caswell_matrix_2008].


Representation of forest dynamics with Markov chains allows the exploration of
ecological mechanisms [@wootton_prediction_2001] underlying biome shifts. For
example, transitions to pioneer reflect disturbance, transitions from pioneer
reflect colonisation, dispersal and recruitment limitation, and transitions
between the other states reflect competitive exclusion. In addition, multi-state
models can be used to investigate biome shifts from the perspective of both
transient dynamics and long-term equilibrium. Markov transition matrices can be
estimated from the model output and their well-established matrix properties can
then be compared under different scenarios [@hill_markov_2004;
@boulangeat_transient_2018]. For instance, the steady-state distribution can be
derived from a transition matrix and used to infer the potential long-term
forest composition under some environmental conditions [i.e., the attractor,
@scheffer_catastrophic_2001], providing insights about the direction of the
current forest dynamic [@waggoner_transition_1970; @hill_markov_2004]. Moreover,
transient periods can be described using the time of convergence to reach the
steady-state distribution, which measures the length of the transient period;
the turnover time indicates how fast the transitions occur and informs about the
persistence of forest states; and the entropy reveals the uncertainty about the
next transition. Contrasting empirically derived transition matrices and their
properties among disturbance scenarios can shed new light on forest dynamics
under climate change and may even provide hints about management measures.


Here, we investigate how forest dynamics is influenced by disturbances and soil
conditions under recent climate warming. In particular, we ask the following
questions: (1) What are the trends in recent forest transition dynamics? (2) How
do disturbances and soil characteristics influence the transition probabilities
among forest states? (3) Do different disturbance types and intensities impact
the potential long-term equilibrium distribution of forest states? And (4) how
do different disturbance types and intensities influence the short-term
transient dynamics under climate change? We answer those question by estimating
the influence of environmental covariates on transition probabilities among four
forest states (boreal, mixed, temperate and pioneer) using a continuous-time
Markov multi-state model. Using results from our model, we then examine the
impact of disturbances on forest equilibrium and transient dynamics by comparing
different complementary matrix properties.

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
result in a northward range shift of temperate forests at equilibrium.

# Methods

## Study area and forest inventory data

We used forest inventory plots in Québec, Canada, to investigate large-scale
transition dynamics in forest communities. Permanent plots have been sampled
approximately every ten years from 1970 to 2018 (and ongoing) by the *Ministère
des forêts, de la Faune et des Parcs* [@mffp_placettes-echantillons_2016] in
order to monitor changes in forest productivity and growth. The study area
extends from approximately 45° to 52° North latitude (ca. 795 000 km^2^). It
covers six bioclimatic domains (Fig. 1) and three different forest subzones; the
mixed forest, which corresponds to the balsam fir-yellow birch domain (from 47°N
to 48°N; hereafter, the ecotone), marks the transition between the hardwood
forest to the south, dominated by *Acer saccharum*, and the boreal forest to the
north, dominated by *Abies balsamea*  and *Picea mariana*. Plots were randomly
positioned across these three subzones with a decreasing sampling intensity
northward [@mffp_reseaux_2014].

We first selected all inventory plots that had been sampled at least twice. We
then disregarded plots that were subjected to active reforestation (i.e.,
plantation) during the study period because we were interested in transition
dynamics resulting from natural recolonisation processes. Finally, we kept plots
for which soil covariates were available. This yielded a total of 11,058 plots
analysed (Fig. 1). The time intervals between plot surveys varied from 3 to 39
years, with a mean interval of 11 years (SD = 3.45).


## Forest states

We classified the forest inventory plots into four forest states using species
basal area and composition at each sampling date. We first assigned each studied
species to a state: boreal, temperate or pioneer according to their functional
traits [Table S1; see @brice_disturbances_2019 for details]. For each plot, we
computed the total basal area of each species group and then classified the plot
to one of the four states similar to the @mffp_placettes-echantillons_2016
definitions; Boreal (boreal species representing >75% of the plot basal area),
Temperate (temperate species representing >75% of the plot basal area), Mixed
(when temperate and boreal species both occupy between >25% and <75% of the plot
basal area), and Pioneer (when the basal area of pioneer species is superior to
that of boreal and temperate species or when plot total basal area <5m^2^/ha).
We analysed state transitions between consecutive plot surveys. Based on this
classification, for the 42,633 observations (plots $\times$ number of years
measured), we recorded 31,690 state transitions (Fig. 2, Table S2).

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
plots, while CMI have shown no trends (Fig. S1).

We also collected information pertaining to natural and anthropogenic
disturbances that have affected the forest plots during the study period (Table
1; Figs. S2). At each plot, the type of disturbances (21 types) and their level
of intensity (moderate or major) were recorded during field surveys [Fig. S2;
@mffp_placettes-echantillons_2016]. The MFFP defined major disturbances as
events that have eliminated more than 75% of the total tree basal area, whereas
moderate disturbances have eliminated between 25% and 75%. For our multi-state
model, we differentiated two main types of disturbances: natural disturbances
and logging, with three levels of intensity each (minor, moderate or major).

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

Note that we did not include all potential ecological processes affecting the
forest dynamics. Rather, we focused on a set of variables that allowed us to
determine how disturbances and soils influence transition dynamics under recent
climate change. For instance, we decided not to include an index of propagule
availability, even though it is known to affect tree range shifts
[@pearson_climate_2006], as the neighbourhood composition is already very
strongly correlated with our climate covariates. Our model is therefore
well-suited for our research goals, however it is not designed to make future
range shift projections.



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
governed by a 4 $\times$ 4 transition intensity matrix, $Q$, where rows
are the current states and columns are the future states (Fig. 2b). For each state
$r,s \in {B, M, P, T}$, the transition intensity ($q_{rs}$) represents the
instantaneous risk that a plot transitions from state $r$ to state $s$. Because
the states were defined based on stand basal area, an instantaneous transition
from Boreal to Temperate and from Temperate to Boreal was considered impossible.
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
on the variable value at $t$, but can change between the intervals. Their
inclusion in the model allows one to fit a non-homogeneous Markov process.
Estimation of model parameters were obtained by maximising the log-likelihood
(see Supplementary Methods for details).

We built five models: one baseline model that solely includes the $q_{rs.0}$,
one model for each category of covariates independently (climate, soil and
disturbances), and one full model, which combines all covariates (Table 1).
Because multiple state transitions are estimated in a single model (all $q_{rs}$
in Fig. 2b), the number of parameters increases rapidly with the number of
covariates (number of modelled transitions (here 10) $\times$ (number of
covariates + 1)). Thus, to reduce the number of parameters, we hypothesised that
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
criterion [AIC; @burnham_model_2002].


### Model baseline and hazard ratios

We first evaluated the trends in recent forest transition dynamics (question
1). We used the baseline hazards ($q_{rs.0}$) estimated by our best model as
indicators of the underlying forest response. For each pair of states, the
baseline hazard describes the risk to make a transition for a mean forest plot
(when all covariates are set to 0).


We investigated the influence of environmental covariates on transition dynamics
(question 2). We compared the estimated hazard ratios derived from our best
model ($exp(\beta_{rs})$). We also computed the predicted 10-year transition
probabilities of forest plots under different disturbance scenarios, while
keeping all other covariates at the average found in the ecotonal zone (i.e.,
the balsam fir-yellow birch domain, Fig. 1), in order to facilitate visual
interpretation of the impacts of disturbances on the transition matrix
structure.

### Transient dynamics and equilibrium

We further investigated how disturbances modify the long-term equilibrium
(question 3) and the forest transient dynamics (question 4). We computed
different properties on the Markov transition matrix and compared them among
levels and types of disturbances along the latitudinal temperature gradient. An
extensive literature describes the multiple properties of discrete-time Markov
transition matrices [@caswell_matrix_2008; @hill_markov_2004] which can be
adapted to continuous-time models. We chose four informative and complementary
properties that fully characterise both the short and long-time scale dynamics
of our modelled system: (1) the steady-state distribution, which corresponds to
the potential long-term proportion of forest states at equilibrium; (2) the
half-life to equilibrium, which evaluates the time of convergence to the
steady-state and the length of the transient period; (3) the turnover time,
which measures the rate of transient successional changes; and (4) the entropy,
which captures the uncertainty regarding the next transitions. While their
absolute values should be interpreted with caution, their comparison under
various disturbance scenarios can highlight essential features of the dynamics.

First, to measure the potential direction of forest dynamics under a given
scenario, we estimated the steady-state distribution, $\pi$. For a regular
Markov process, any initial state distribution converges to the same equilibrium
as time approaches infinity. The vector of equilibrium $\pi$ can be obtained
by taking the left eigenvector of $Q$, which has an eigenvalue of 0, normalised
to sum to 1, or the left eigenvector of $P$, which has an eigenvalue of 1,
normalised to sum to 1 [@norris_1997].

Then, the convergence rate to the equilibrium distribution can be measured
using the damping ratio [@hill_markov_2004]:

$$\rho = \lambda_{P1} / \lambda_{P2} = exp(\lambda_{Q1} - \lambda_{Q2})\text{,}$$

where $\lambda_{P1}$ and $\lambda_{P2}$ are the largest and second-largest
eigenvalues of $P$ ($\lambda_{P1}$ = 1 for stochastic $P$), whereas
$\lambda_{Q1}$ and $\lambda_{Q2}$ are the largest and second-largest eigenvalues
of $Q$ ($\lambda_{Q1}$ = 0 for stochastic $Q$). The convergence time was
approximated using the half-life to equilibrium:

$$t_{1/2} = log(2)/log(\rho)\text{.}$$

We also measured the turnover time in each forest state, also called the
sojourn time in multi-state models, which corresponds to the time spent in one
state before transitioning to the next. The turnover time can be estimated by
$Turnover_{r} =−1/q_{rr}$, where $q_{rr}$ is the $r^{th}$ entry on the
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

## Model evaluation

Overall, the full model, which includes climate, soil and disturbance variables,
had the best fit and predictive performances for the observed data (Table 2;
Fig. S3). The second-best model was the disturbance model, but it was far behind
with a difference in AIC of almost 1500 units from the full model (Table 2).
All variable subsets improved significantly the likelihood of the model (all
likelihood ratio tests were highly significant, *p* << 0.001; Table 2). Model
performance was also evaluated using a 10-fold cross-validation (see
Supplementary Methods); it revealed that including climate and disturbances
improved overall model predictive performance, while soil variables had a
negligible effect (Fig. S3). Thereafter, all inferences about transition
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
q~BP~) more likely than a succession resetting caused by disturbances.



## Effect of covariates on transition probabilities

The full multi-state model indicates that forest state transitions are
contingent upon environmental covariates. All transitions to Pioneer were highly
influenced by disturbances (Fig. 4, Table S3). As could be expected, major
disturbances exert stronger effects than moderate disturbances (for both natural
and logging), but, for each level of intensity, logging had stronger effects
than natural disturbances. For example, the risk of transition from
Boreal to Pioneer has surged up to 213 times higher for plots that suffered
major logging (logging 2) compared to undisturbed plots (minor). Disturbances of
all types and intensities favoured transitions from Mixed to Temperate forests.
Moderate disturbances (natural and logging) doubled the risk of this type of
transition, whereas major disturbances increased it by ca. 5 times (Hazard Ratio
(HR) = 5.76 and 5.32, for natural and logging, respectively). Although the
effect of major disturbances on the instantaneous risk of transition from Mixed
to Temperate was stronger than for moderate disturbances, the probability of
this event decreased with time (Fig. S4). Moderate disturbances also favoured
transitions from Boreal to Mixed (HR = 2.76 and 3.45, respectively), while major
disturbances had no significant effect on this type of transition.

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


Disturbances completely altered the structure of the 10-year transition
probability matrix (Fig. 5). The largest values across most matrices were
generally associated with self-transitions (matrix diagonal), meaning that the
vast majority of forest plots (characterised by the average environmental
conditions found in the ecotone) remained in the same state after 10 years. For
undisturbed forest plots (minor), the self-transitions were very strong but
transitions from Pioneer to Boreal and from Mixed to Temperate were not trivial.
At moderate disturbances, probabilities of self-transitions decreased, while
transitions from Boreal to Pioneer, and from Mixed to Temperate increased the
most. Transitions from Mixed and Temperate to Pioneer did not increase much at
moderate disturbances, likely because such disturbances were less frequent and
less severe than in Boreal forests. The difference between natural disturbances
and logging was conspicuous only for major disturbances. For both types of major
disturbances, the transition probabilities to Pioneer showed a great increase
compared to moderate disturbances, but these values exploded in the major
logging transition matrix, exceeding self-transitions. Interestingly, the
estimated transition probability from Mixed to Temperate remained quite high after
major disturbances.



## Effect of disturbances on long-term equilibrium

The potential state proportion at equilibrium was strongly influenced by
disturbances (Fig. 6a). For the undisturbed scenario (minor), the predicted
equilibrium at the ecotone was relatively close to the initial observed
proportions, with signs of regeneration from Pioneer to Boreal states and slight
increases in Mixed and Temperate states. The steady-state proportion of
Temperate almost doubled with moderate disturbances (minor: 33%; moderate
natural: 56%; moderate logging: 60%), while the boreal state was more than
halved. At major disturbances, Pioneer forests dominated the equilibrium
landscape, while the other states collapsed.  

The steady-state proportion also changed as expected along the temperature
gradient (Fig. 6b,c). The Boreal state dominates at low temperature (high
latitude) and the Temperate state dominates at high temperature (low latitude),
highlighting the position of the boundary between these two biomes at a growing
season temperature of about 12.9°C, which is found in the actual ecotone. Moderate
disturbances (both natural and logging) displaced the temperate-boreal boundary
at lower temperatures (ca. 12.2°C), hence further north of the current ecotone
(Fig. 6b,c). Because of the dominance of the Pioneer state, the boundary
modestly moved north with major natural disturbances (12.7°C), while it
retreated to the south with major logging (13.4°C).


## Effect of disturbances on transient dynamics

Disturbances affected forest transient dynamics with greater impact for
higher disturbance intensity (Fig. 7). In the minor disturbance scenario,
turnover time was generally longer at low temperature, indicating slower
transition dynamics in northern forests (Fig. 7a,b). The turnover
time then rapidly declined to reach a minimum at ca. 13.25°C, at the southern limit
of the ecotone, and went back up after this point. This trough, where transition
dynamics is the fastest, is located just a little south of the boundary between
Boreal and Temperate dominances found in Figure 6. Major disturbances
accelerated transition dynamics all along the temperature gradient, while
moderate disturbances also decreased turnover time but more strongly in the
northern boreal region (Fig. 7a,b). The effect on turnover time was similar for
both disturbance types, except that the effect of major logging was much
stronger in northern boreal forests than natural disturbances (Fig. 7a,b). These
spatial patterns reflect the turnover time of the dominant state at each point
along the temperature gradient (Fig. S5).

At minor disturbances, the entropy of the system generally increased from north
to south and peaked at ca. 12.6°C, at the northern end of the ecotone (Fig. 7c,d).
This peak illustrates where the transition dynamics is most uncertain
(transition to all states are possible at this point), while it is very
predictable in northern boreal forests (Boreal stays Boreal until it transitions
to Pioneer later on). The peak can be mainly attributed to the entropy of the
Boreal state at the ecotone, and the generally high values at low latitudes can
be principally attributed to the Temperate state (Fig. S6). This latitudinal
pattern of entropy is modified by disturbances. Moderate natural disturbances
decreased the entropy throughout the gradient, but especially where the peak is
found (Fig. 7c). With moderate logging, the peak disappeared, and entropy
increased monotonically from north to south (Fig. 7d). The peak of entropy was
displaced to the south when major disturbances were included, whether natural or
logging (Fig. 7c,d), where it was dominated by the entropy of the Pioneer state
(Fig. S6).

Half-life to equilibrium was the longest at ca. 11.8°C, north of the ecotone,
in the balsam fir-white birch domain, while it was the shortest in the
southernmost latitudes (Fig. 7e,f). Moderate disturbances flattened and shifted
this peak to the north and the effect of moderate logging (Fig. 7f) was stronger
than natural disturbances (Fig. 7e). In the balsam fir-white birch, the
half-life to reach equilibrium distribution was reduced almost by half by
moderate logging. With major disturbances, forests all along the temperature
gradient can reach very quickly their steady-state distribution (maximum of
about 8 years for major logging and 25 years for major natural disturbances).


# Discussion

Our study reveals that disturbances are likely to accelerate forest responses to
climate change by promoting transitions from mixed to temperate forests. Our
analysis of the equilibrium further highlights that the long-term forest
dynamics under moderate disturbances favours an increased proportion of
temperate forests and a northward shift of the temperate-boreal ecotone.
Disturbances also modify the forest transient dynamics, accelerating both the
turnover and convergence time and making the dynamics more predictable. In
accordance with the hypothesis formulated in previous studies
[@johnstone_changes_2010; @johnstone_changing_2016;
@vissault_biogeographie_2016; @brice_disturbances_2019], our findings show that
moderate disturbances catalyse transitions to the alternate, temperate-dominated
forest state and could therefore promote regime shifts.



## Trends in recent forest transition dynamics

Forest dynamics in Québec during the last decades was dominated by transitions
from pioneer to boreal and from mixed to temperate stands. The important
regeneration of boreal forests could be attributed to past natural disturbances,
notably the last spruce budworm outbreak. Indeed, the last outbreak, which
occurred during the 1970s, has caused major mortality in coniferous species
followed by important recruitment pulses and growth releases
[@bouchard_tree_2006].

Although we did not directly evaluate the impact of climate change, our results
suggest that recent climate warming may contribute to the forest transition
dynamics. The high baseline transition rate from mixed to temperate is
consistent with the expectation of a northward range shift of temperate trees
into the mixed and boreal forests. These transitions were caused by a
concomitant increase in temperate species and decrease in boreal species (Fig.
S7). Accordingly, the warming trend of the last decades (Fig. S1) has been shown
to increase growth and reproductive rates of temperate species and reduce growth
of boreal species [@reich_geographic_2015; @fisichelli_temperate_2014;
@boisvertmarsh_divergent_2019; @goldblum_tree_2005], thus providing a
competitive advantage to temperate over boreal species.  

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

Our study highlighted that moderate disturbances favour climate-related
transitions, whereas major disturbances merely promote pioneer states.
Disturbances directly remove trees, which leads to immediate and substantial
changes in forest composition [@brice_disturbances_2019]. Without climate
change, forests are expected to be resilient to normally experienced
disturbances and should thus return to their preceding states. However, climate
change alters the conditions that initially supported the persistence of a given
state, making forests susceptible to transition to other states
[@johnstone_changing_2016]. Hence, moderate disturbances likely facilitate
colonisation and establishment by opportunistic temperate species under warmer
conditions [@brice_disturbances_2019; @leithead_northward_2010;
@landhausser_disturbance_2010]. In contrast, severe disturbances in the study
area, primarily clearcutting but also large fires (Fig. S2), create openings of
very large extent which are likely detrimental to temperate species and favour
early-successional species that can disperse seed over long distances, such as
*Populus sp* and *Betula sp* [@landhausser_disturbance_2010].


Compared to the catalysing effect of disturbances, soil characteristics do not
appear to represent a large impediment to state transitions, but transitions may
be slower on some soil types. Poor drainage constrained climate-related
transitions from Boreal to Mixed states, but not from Mixed to Temperate. This
indicates that temperate species can readily colonise soils found in mixedwoods
but may have more difficulty in colonising hydric boreal soils. Very poor
drainage, often associated with peatland and thick organic layer, is usually
thought to be improper for the regeneration of temperate species
[@lafleur_response_2010]. Several studies found that *Acer saccharum*
regenerates well across the ecotone because of its large tolerance to various
soil conditions [@barras_supply_1998; @goldblum_age_2002; @kellman_sugar_2004;
@fisichelli_temperate_2014; @collin_can_2018]. At their northern range limit,
*A. saccharum* and *A. rubrum*, the species contributing most to compositional
changes [@brice_disturbances_2019], are hypothesised to be mostly limited by
cold soil temperature [@barras_supply_1998; @goldblum_age_2002].


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
different disturbance scenarios. Our findings support the simulations of
@boulanger_climate_2019 where harvesting under future climate warming was
projected to promote further invasions of pioneer species, such as *Populus*,
and temperate species, such as *Acer* and *Fagus*, in mixedwoods of Québec. In
contrast, based on their simulations, @liang_how_2018 and @vanderwel_how_2014
concluded that logging would primarily accelerate the expansion of pioneer
forests, but would have little or no effect on extensive biome shifts over the
next century in eastern United States. Divergence in responses to disturbances
among tree species could explain these apparently conflicting results.
Disturbances may facilitate the range expansion of some species but hinder that
of others depending on their functional traits [@matthews_disturbance_2013;
@aubin_traits_2016]. For instance, because of its positive response to past
[@danneyrolles_stronger_2019], recent [@brice_disturbances_2019] and future
[@boulanger_climate_2019] disturbances in Québec, *Acer rubrum* is likely to
play a disproportionate role in the temperate biome shift.


## Disturbances accelerate the transient dynamics

Beyond their impacts on the equilibrium, our results suggest that disturbances
may have a substantial influence on forest transient dynamics. In the continuous
boreal zone (spruce-moss domain), forests dominated by *Picea mariana* are
usually characterised by dynamics of stand self-replacement with minimal
compositional changes across disturbance cycles [@goldblum_deciduous_2010].
Consistent with this dynamics, the turnover time of undisturbed northern boreal
forests was very long and the entropy very low in our results. The turnover was
shortened by disturbances, but the entropy remained low, indicating that the
dynamics was still very predictable (back and forth transitions between boreal
and pioneer states) and that there was no directional shift associated with
disturbances. Hence, boreal forests lose their persistence when moderately
disturbed but remain resilient as they return to their previous boreal state.
Under major disturbances, boreal forests collapsed to pioneer state and reached
this new equilibrium swiftly (short half-life). This observation is consistent
with previous studies suggesting that boreal forests can easily shift into an
alternative treeless state in response to severe or repeated disturbances
[@payette_shift_2003; @sanchez-pinillos_resistance_2019].

In contrast, the ecotone is characterised by a rapid turnover and high entropy
indicating abrupt compositional shift which can go in any direction. Compared
to northern boreal forests, the short turnover time implies a low persistence of
the forest states in this region even under minor disturbances. This result
corroborates the predictions made by @vissault_biogeographie_2016, where mixed
forests would undergo a swift conversion to temperate forests in the next
decades, whereas boreal forests would present a large inertia presumably because of
dispersal limitation. The dynamics of the ecotone appears unstable because it
is caught between two stable states, i.e. boreal to the north and temperate to
the south. Under moderate disturbances, the probability of transitioning to
Temperate increases to the detriment of the other possible states, hence the
entropy is decreased, and the dynamics becomes more predictable. Such a clear
directional shift strongly indicates non-equilibrium dynamics in this region.
Although turnover is fast, half-life to equilibrium is long because a forest may
not move towards equilibrium and may undergo multiple transitions.


## Ecological and management implications

Our study provides a strong empirical basis for predicting the types of changes
in forest dynamics that are likely to unfold in the 21st century. A shift in
dominant forest cover from conifer to deciduous broadleaf species not only
entails changes in tree species diversity and composition, but a complete
transformation of forest dynamics and functions. In the long term, this regime
shift could locally increase tree diversity [@berteaux_cc-bio_2010] and carbon
sequestration [@thurner_carbon_2014], modify disturbance regimes [reduced
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
climate-related transitions. This suggests that thoughtful logging practices
could be used to reduce extinction debt and colonisation credit, and thus reduce
migration lags. Other studies also recommend to plant temperate trees farther
north outside their current range to facilitate their migration
[Vieira et al. in prep.; @duveneck_measuring_2016]. However, before implementing
such silvicultural strategies, key questions need to be answered. For instance,
will multiple interacting disturbances exacerbate tree mortality? Which species
will be able to benefit from the opportunities created by canopy openings? And
how will these rapid transitions impact the ecosystem processes and functions?
The rate of recent climate change already outpaces tree migration capacity
[@sittaro_tree_2017], and even more so the scientific capacity to understand and
mitigate its consequences [@jackson_balancing_2010]. Therefore, in order to
insure long-term forest health at the temperate-boreal ecotone, we need to
simultaneously limit global warming through drastic reduction of greenhouse gas
emissions and intensify research effort to develop effective adaptation
strategies.



## Acknowledgements

The authors declare no conflict of interest.

We are grateful to Kevin Cazelles for providing helpful suggestions and
comments that improved our analyses and manuscript. We also thank Guillaume Guénard for
useful advices on the model. This research was supported by Natural Sciences and
Engineering Research Council of Canada (NSERC) research grant no. 7738 to P. L.
and no. 5134 to M.‐J. F.

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
| **Disturbances**|                                                           |
|Logging          |Tree harvesting, including clearcutting, selection cutting, shelterwood cutting, seed-tree cutting, etc. None or minor (0), moderate (1) or major (2). |
|Natural          |Natural disturbances, including forest fires, insect outbreaks, windfall, etc. None or minor (0), moderate (1) or major (2). |


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

![Hazard ratios (HR) and 95% confidence intervals as estimated from the best multi-state transition model. Each plot represents the estimated HR for transitions from row to column state, e.g., the plot on the first row, second column shows the HR for the Boreal to Mixed transition. The ordinate is in log scale. The HR of predictors are interpretable as multiplicative effects on the hazard, where values above 1 (in blue) indicate that the predictor is associated with a greater risk of state transition, whereas values below 1 (in red) indicate a lower risk of transition. Predictors statistically different from 1 are shonwed in dark blue or red. Numbers following disturbance predictors indicate their levels of intensity: 1 = moderate and 2 = major.](res/fig4_HR.pdf)

\pagebreak

![Predicted change in 10-year transition probabilities for different disturbance types and levels. All other covariates are fixed at the average conditions found in the ecotone. Letters correspond to the four forest states: (B)oreal, (M)ixed, (P)ioneer and (T)emperate. Numbers are the modelled transition probabilities from rows to columns and darker colour indicates higher transition probability.](res/fig5_pmatrix.pdf)

\pagebreak

![Changes in forest state proportions at equilibrium for different disturbance types (natural or logging) and intensity (minor, moderate, major). The barplot (a) compares the observed state proportion in the ecotone to the potential state proportion at equilibrium for different disturbance scenarios with all other covariates fixed at the average conditions found in the ecotone. The curved lines (b,c) show the proportions of Boreal (blue) and Temperate forests (red) at equilibrium along the temperature (latitudinal) gradient for minor (solid), moderate (dashed) and major (dotted) disturbances, with all other covariates fixed at the average conditions found in the ecotone. The light (minor), medium (moderate) and dark (major) grey circles indicate the positions of the boundary between dominance of Boreal and Temperate forests (i.e. the advancing front) while the corresponding arrows show how moderate and major disturbances move the boundary. The polygon approximates the position of the ecotone along the temperature gradient.](res/fig6_steady.pdf)

\pagebreak

![Changes in the characteristics of the forest transient dynamics along the temperature (latitudinal) gradient for different disturbance scenarios: minor (solid), moderate (dashed) and major (dotted) disturbances for both natural (a,c,e) and logging (b,d,f). All other covariates are fixed at the average conditions found in the ecotone to focus solely on the effect of disturbances along the temperature gradient. The turnover of the whole system (i.e. whole transition matrix) (a,b) corresponds to the time spent in a state before transitioning to the next and is given by the average of each state turnover time over the steady-state distribution. The entropy of the whole system (c,d) corresponds to the uncertainty of the next transition and is given by the average of each state entropy over the steady-state distribution. The half-life to equilibrium (e,f) is the time taken to reach 50% of the steady-state distribution, i.e. when the first eigenvalue becomes twice as large as the contribution of the second eigenvalue. The polygon approximates the positions of the ecotone along the temperature gradient.](res/fig7_transients.pdf)

\pagebreak

# References
