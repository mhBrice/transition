---
title: Supplementary Information
geometry: margin=1in
header-includes:
    - \usepackage{setspace}
    - \setstretch{1.5}
    - \usepackage{lscape}
    - \newcommand{\blandscape}{\begin{landscape}}
    - \newcommand{\elandscape}{\end{landscape}}
    - \usepackage{caption}
    - \usepackage{makecell}
    - \usepackage{colortbl}
    - \usepackage{xcolor}
    - \usepackage{float}
    - \floatplacement{figure}{H}
    - \renewcommand{\figurename}{\bfseries Figure}
    - \renewcommand\thefigure{S\arabic{figure}}
    - \renewcommand\thetable{S\arabic{table}}
    - \usepackage{setspace}
---

\begin{center}
{\Large Moderate disturbances accelerate forest transition dynamics under climate change in the temperate-boreal ecotone of eastern North America}
\end{center}

\begin{center}
{\large Marie-Hélène Brice, Steve Vissault, Willian Vieira, Dominique Gravel, Pierre Legendre, Marie-Josée Fortin}
\end{center}

# Supplementary Methods

## Multi-state models

For states $r,s \in {B, M, P, T}$ and time $t, \Delta{t} >= 0$, transition probabilities
($a_{rs}$) are defined as the probability that a plot in state $r$ at time $t$
is in state $s$ at time $t + \Delta{t}$ and can be written as:

$$A_{r,s}(t + \Delta{t}) = A(S_{t + \Delta{t}} = s | S_{t} = r)\text{.}$$

In a four-state transition model, the transition probability matrix $P(t +
\Delta{t})$, hereafter simplified to $A(t)$, is a 4 $\times$ 4 matrix, where the
rows are the current states and the columns the future states, containing the
transition probabilities $a_{rs}(t)$ for a specified time interval. For a
time-homogeneous model, $A(t)$ is solved by taking the matrix exponential of
the intensity matrix $Q$ scaled by the time interval:

$$A(t) = e^{tQ}\text{.}$$

The intensity matrix $Q$ contains transition intensities $q_{r,s}$ which
represent the instantaneous risk of moving from state $r$ to state $s$:

$$q_{r,s} = \lim_{\Delta \to 0} \frac{P(Y_{t+\Delta} = s | Y_t = r)}{\Delta}\text{, on
off-diagonal elements,}$$

$$q_{r,r} = − \sum_{s \neq r}{q_{rs}}\text{, on diagonal elements.}$$


We can define transition-specific hazard regression models for those states $r,s
\in {B, M, P, T}$ between which a direct transition is possible according to the
specified multi-state process (Fig. 2 in main text). The intensities $q_{r,s}$ can be
modelled as a product of a baseline hazard $q_{rs.0}$ and a log-linear effect of
the explanatory variables $x(t)$ and their coefficients $\beta_{rs}$:

$$q_{rs}(t) = q_{rs}(t|x(t)) = q_{rs.0}(t)exp(\beta_{rs}'x(t))\text{.}$$

In this model, $q_{rs.0}(t)$ is a baseline hazard function that describes the
risk for a reference plot $i$ with environment $x_i(t) = 0$, and
$exp(\beta_{rs}'x(t))$ is the relative increase or decrease in risk associated
with the set of characteristics $x_i(t)$. This model allows one to include the
effect of time-dependent covariates on transition intensities and therefore to
relax the time homogeneity assumption of Markov models. Time-dependent
covariates, such as climate and disturbances, are assumed to be
piecewise-constant, i.e., the hazard is constant within a specified time interval
$[t, t + \Delta{t}]$ and depends on the covariate value at $t$, but is allowed
to change between the intervals.

Estimation of model parameters can be obtained by maximising the log-likelihood
function using the transition probability matrix. The contribution of plot
$i$ at time $j$ to the likelihood is given by:

$$LL_i(\theta | s,x) = \prod\limits_{j=1}^{J} A(S_j=s_j|S_{j-1}=s_{j-1},\theta,x)\text{,}$$

where $\theta$ is the vector with all model parameters, $x$ denotes the
vector with the covariate values, and $s$ denotes the observed state trajectory
$s_1,...,s_J$ at times $t_1,...,t_J$. The full likelihood function is the
product of the contributions of all $N$ plots:

$$LL(\theta) = \prod\limits_{i=1}^{N} LL_i(\theta | s,x)\text{.}$$


\pagebreak

## Performance of candidate models

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

$$\hat{A}(r,s) = [\hat{A}(r|s)+ \hat{A}(s|r)]/2\text{.}$$

Averaging the pairwise AUC gives the overall multi-class AUC (hereafter mAUC) of
the model:

$$mAUC = \frac{2}{c(c-1)} \sum_{r<s} \hat{A}(r,s)\text{.}$$

The LS was proposed by Good (1952) and is often used in weather forecasts
[@gneiting_strictly_2007]. While AUC is a function of different classification
thresholds, LS measures the degree to which predicted probabilities are close to
the observed outcomes. We computed a global score for each model:

$$LS = \frac{1}{N} \sum_{i=1}^N -log(P(S_i = s_i))\text{,}$$

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



\pagebreak

# Supplementary tables

**Table S1**. List of the 46 species included in the analyses, their frequency and their corresponding group. The frequency corresponds to the number of forest plots in which they were observed. The species groups were defined using their trait values and knowledge of species ecology [see @brice_disturbances_2019 for details].

\begin{longtable}{>{\em}llr}
\toprule
\textbf{Species name} & \textbf{Vernacular name} & \textbf{Frequency}\\
\midrule
\endfirsthead
\multicolumn{3}{@{}l}{\textit{(continued)}}\\
\toprule
\textbf{Species name} & \textbf{Vernacular name} & \textbf{Frequency}\\
\midrule
\endhead
\
\endfoot
\bottomrule
\endlastfoot
\addlinespace[0.3em]
\multicolumn{3}{l}{\textbf{Boreal}}\\
\hspace{1em}Abies balsamea & Balsam fir & 7870\\
\hspace{1em}Larix laricina & Tamarack & 551\\
\hspace{1em}Picea glauca & White spruce & 4181\\
\hspace{1em}Picea mariana & Black spruce & 6082\\
\hspace{1em}Pinus banksiana & Jack pine & 1339\\
\addlinespace[0.3em]
\multicolumn{3}{l}{\textbf{Pioneer}}\\
\hspace{1em}Betula papyrifera & White birch & 5863\\
\hspace{1em}Betula populifolia & Grey birch & 226\\
\hspace{1em}Populus balsamifera & Balsam poplar & 216\\
\hspace{1em}Populus deltoides & Cottonwood & 2\\
\hspace{1em}Populus grandidentata & Large tooth aspen & 582\\
\hspace{1em}Populus tremuloides & Trembling aspen & 2504\\
\hspace{1em}Prunus pensylvanica & Pin cherry & 1495\\
\hspace{1em}Salix sp. & Willow & 499\\
\hspace{1em}Sorbus sp. & Mountain-ash & 515\\
\addlinespace[0.3em]
\multicolumn{3}{l}{\textbf{Temperate}}\\
\hspace{1em}Acer negundo & Manitoba maple & 1\\
\hspace{1em}Acer nigrum & Black maple & 3\\
\hspace{1em}Acer pensylvanicum & Striped maple & 719\\
\hspace{1em}Acer rubrum & Red maple & 3273\\
\hspace{1em}Acer saccharinum & Silver maple & 23\\
\hspace{1em}Acer saccharum & Sugar maple & 2190\\
\hspace{1em}Acer spicatum & Mountain maple & 206\\
\hspace{1em}Amelanchier sp. & Serviceberry & 33\\
\hspace{1em}Betula alleghaniensis & Yellow birch & 2582\\
\hspace{1em}Carpinus caroliniana & Blue beech & 6\\
\hspace{1em}Carya cordiformis & Bitternut hickory & 10\\
\hspace{1em}Fagus grandifolia & American beech & 928\\
\hspace{1em}Fraxinus americana & White ash & 261\\
\hspace{1em}Fraxinus nigra & Black ash & 465\\
\hspace{1em}Fraxinus pennsylvanica & Red ash & 33\\
\hspace{1em}Juglans cinerea & Butternut & 16\\
\hspace{1em}Ostrya virginiana & Ironwood & 493\\
\hspace{1em}Picea rubens & Red spruce & 861\\
\hspace{1em}Pinus resinosa & Red pine & 131\\
\hspace{1em}Pinus rigida & Pitch pine & 2\\
\hspace{1em}Pinus strobus & Eastern white pine & 762\\
\hspace{1em}Prunus serotina & Black cherry & 243\\
\hspace{1em}Quercus alba & White oak & 9\\
\hspace{1em}Quercus bicolor & Swamp white oak & 8\\
\hspace{1em}Quercus macrocarpa & Bur oak & 14\\
\hspace{1em}Quercus rubra & Red oak & 453\\
\hspace{1em}Thuja occidentalis & White cedar & 1269\\
\hspace{1em}Tilia americana & Basswood & 395\\
\hspace{1em}Tsuga canadensis & Eastern hemlock & 472\\
\hspace{1em}Ulmus americana & American elm & 228\\
\hspace{1em}Ulmus rubra & Red elm & 9\\
\hspace{1em}Ulmus thomasii & Rock elm & 3\\*
\end{longtable}


\pagebreak

**Table S2**. Frequency of all observed transitions between the four forest states during the study period. Transitions are from rows to columns.

\begin{tabular}{>{\bfseries}lccccc}
\toprule
\textbf{ } & \textbf{Boreal} & \textbf{Mixed} & \textbf{Pioneer} & \textbf{Temperate} & \textbf{Sum}\\
\midrule
Boreal & 9760 & 203 & 1454 & 3 & 11420\\
Mixed & 82 & 4996 & 393 & 682 & 6153\\
Pioneer & 1539 & 678 & 6460 & 198 & 8875\\
Temperate & 2 & 488 & 168 & 4584 & 5242\\
Sum & 11383 & 6365 & 8475 & 5467 & 31690\\
\bottomrule
\end{tabular}


\pagebreak
\blandscape

**Table S3**. Table of baseline transition intensities ($q_{rs,0}$ in first column) and Hazard ratios (HR) and their 95% confidence intervals as estimated from the best multi-state transition model. The HR of covariates are interpretable as multiplicative effects on the baseline hazard, where values above 1 indicate that the predictor is associated with a greater risk of state transition, whereas values below 1 indicate a lower risk of transition. Covariates statistically different from 1 are coloured in grey.

\begin{table}[H]
\setstretch{1.4}
\centering\begingroup\fontsize{6.8}{8.8}\selectfont

\begin{tabular}{>{\bfseries}l>{\centering\arraybackslash}p{1.9cm}>{\centering\arraybackslash}p{1.8cm}>{\centering\arraybackslash}p{1.8cm}>{\centering\arraybackslash}p{1.8cm}>{\centering\arraybackslash}p{1.8cm}>{\centering\arraybackslash}p{1.8cm}>{\centering\arraybackslash}p{1.8cm}>{\centering\arraybackslash}p{1.8cm}>{\centering\arraybackslash}p{1.8cm}}
\toprule
\textbf{Transitions} & \textbf{Baseline} & \textbf{Temperature} & \textbf{CMI} & \textbf{Drainage} & \textbf{pH} & \textbf{Natural1} & \textbf{Natural2} & \textbf{Logging1} & \textbf{Logging2}\\
\midrule
Boreal - Boreal & -0.009\newline (-0.01, -0.008) & \cellcolor{white}{ } & \cellcolor{white}{\cellcolor{white}{ }} & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ }\\
Boreal - Mixed & 0.002\newline (0.001, 0.003) & \cellcolor[HTML]{cccccc}{7.36\newline (5.28, 10.28)} & \cellcolor[HTML]{cccccc}{\cellcolor[HTML]{cccccc}{1.43\newline (1.11, 1.85)}} & \cellcolor[HTML]{cccccc}{0.74\newline (0.65, 0.84)} & \cellcolor{white}{1.12\newline (0.98, 1.29)} & \cellcolor[HTML]{cccccc}{2.76\newline (2.03, 3.75)} & \cellcolor{white}{2.34\newline (0.94, 5.82)} & \cellcolor[HTML]{cccccc}{3.45\newline (2.2, 5.42)} & \cellcolor{white}{0.72\newline (0.02, 27.18)}\\
Boreal - Pioneer & 0.007\newline (0.006, 0.008) & \cellcolor{white}{1.000} & \cellcolor{white}{\cellcolor{white}{1.000}} & \cellcolor{white}{1.000} & \cellcolor{white}{1.000} & \cellcolor[HTML]{cccccc}{4.76\newline (3.88, 5.84)} & \cellcolor[HTML]{cccccc}{37.01\newline (30.67, 44.66)} & \cellcolor[HTML]{cccccc}{8.45\newline (5.89, 12.12)} & \cellcolor[HTML]{cccccc}{213.52\newline (102.4, 445.24)}\\
Mixed - Boreal & 0.005\newline (0.003, 0.008) & \cellcolor{white}{1.07\newline (0.65, 1.75)} & \cellcolor[HTML]{cccccc}{\cellcolor[HTML]{cccccc}{1.56\newline (1.18, 2.06)}} & \cellcolor{white}{1\newline (0.8, 1.25)} & \cellcolor{white}{0.82\newline (0.65, 1.04)} & \cellcolor{white}{0.98\newline (0.52, 1.85)} & \cellcolor{white}{2.5\newline (0.33, 18.8)} & \cellcolor{white}{1.64\newline (0.94, 2.86)} & \cellcolor{white}{1.14\newline (0.02, 72.72)}\\
Mixed - Mixed & -0.033\newline (-0.037, -0.029) & \cellcolor{white}{ } & \cellcolor{white}{\cellcolor{white}{ }} & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ }\\
Mixed - Pioneer & 0.005\newline (0.004, 0.006) & \cellcolor{white}{1.000} & \cellcolor{white}{\cellcolor{white}{1.000}} & \cellcolor{white}{1.000} & \cellcolor{white}{1.000} & \cellcolor{white}{1.91\newline (0.98, 3.73)} & \cellcolor{white}{4.94\newline (0.6, 40.6)} & \cellcolor[HTML]{cccccc}{2.45\newline (1.22, 4.91)} & \cellcolor[HTML]{cccccc}{31.22\newline (19.07, 51.11)}\\
Mixed - Temperate & 0.023\newline (0.02, 0.026) & \cellcolor{white}{0.95\newline (0.76, 1.18)} & \cellcolor[HTML]{cccccc}{\cellcolor[HTML]{cccccc}{0.81\newline (0.7, 0.92)}} & \cellcolor{white}{0.96\newline (0.88, 1.05)} & \cellcolor{white}{0.94\newline (0.87, 1.01)} & \cellcolor[HTML]{cccccc}{2.52\newline (2.04, 3.11)} & \cellcolor[HTML]{cccccc}{5.76\newline (3.02, 10.98)} & \cellcolor[HTML]{cccccc}{3.39\newline (2.79, 4.11)} & \cellcolor[HTML]{cccccc}{5.32\newline (3.46, 8.16)}\\
Pioneer - Boreal & 0.027\newline (0.025, 0.029) & \cellcolor[HTML]{cccccc}{1.28\newline (1.17, 1.41)} & \cellcolor[HTML]{cccccc}{\cellcolor[HTML]{cccccc}{1.8\newline (1.65, 1.96)}} & \cellcolor{white}{0.98\newline (0.93, 1.02)} & \cellcolor{white}{0.96\newline (0.9, 1.01)} & \cellcolor{white}{0.94\newline (0.77, 1.15)} & \cellcolor[HTML]{cccccc}{0.38\newline (0.29, 0.49)} & \cellcolor[HTML]{cccccc}{1.9\newline (1.55, 2.33)} & \cellcolor{white}{1.84\newline (0.87, 3.89)}\\
Pioneer - Mixed & 0.004\newline (0.003, 0.004) & \cellcolor[HTML]{cccccc}{4.86\newline (3.7, 6.39)} & \cellcolor[HTML]{cccccc}{\cellcolor[HTML]{cccccc}{1.55\newline (1.25, 1.91)}} & \cellcolor{white}{0.9\newline (0.8, 1.01)} & \cellcolor{white}{0.99\newline (0.89, 1.1)} & \cellcolor{white}{1.28\newline (0.83, 1.96)} & \cellcolor{white}{0.52\newline (0.22, 1.22)} & \cellcolor[HTML]{cccccc}{2.54\newline (1.72, 3.76)} & \cellcolor{white}{0.97\newline (0.63, 1.48)}\\
Pioneer - Pioneer & -0.032\newline (-0.034, -0.029) & \cellcolor{white}{ } & \cellcolor{white}{\cellcolor{white}{ }} & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ }\\
Pioneer - Temperate & 0.001\newline (0.001, 0.001) & \cellcolor[HTML]{cccccc}{18.57\newline (12.5, 27.6)} & \cellcolor[HTML]{cccccc}{\cellcolor[HTML]{cccccc}{2.46\newline (1.88, 3.23)}} & \cellcolor[HTML]{cccccc}{0.77\newline (0.65, 0.91)} & \cellcolor{white}{1\newline (0.88, 1.13)} & \cellcolor{white}{0.2\newline (0.02, 2.05)} & \cellcolor{white}{0.54\newline (0.13, 2.3)} & \cellcolor{white}{1.25\newline (0.64, 2.42)} & \cellcolor{white}{0.31\newline (0.08, 1.15)}\\
Temperate - Mixed & 0.014\newline (0.012, 0.017) & \cellcolor[HTML]{cccccc}{0.54\newline (0.41, 0.7)} & \cellcolor{white}{\cellcolor{white}{0.96\newline (0.83, 1.12)}} & \cellcolor[HTML]{cccccc}{1.3\newline (1.15, 1.47)} & \cellcolor[HTML]{cccccc}{0.73\newline (0.65, 0.82)} & \cellcolor{white}{1.24\newline (0.78, 1.97)} & \cellcolor{white}{2.08\newline (0.57, 7.65)} & \cellcolor{white}{1.1\newline (0.84, 1.44)} & \cellcolor[HTML]{cccccc}{2.41\newline (1.12, 5.21)}\\
Temperate - Pioneer & 0.001\newline (0.001, 0.002) & \cellcolor{white}{1.000} & \cellcolor{white}{\cellcolor{white}{1.000}} & \cellcolor{white}{1.000} & \cellcolor{white}{1.000} & \cellcolor{white}{1.13\newline (0.12, 10.61)} & \cellcolor[HTML]{cccccc}{51.54\newline (22.41, 118.54)} & \cellcolor[HTML]{cccccc}{3.8\newline (2.13, 6.8)} & \cellcolor[HTML]{cccccc}{146.94\newline (94.03, 229.62)}\\
Temperate - Temperate & -0.016\newline (-0.019, -0.013) & \cellcolor{white}{ } & \cellcolor{white}{\cellcolor{white}{ }} & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ }\\
\bottomrule
\end{tabular}
\endgroup{}
\end{table}
\elandscape

\pagebreak
\blandscape

**Table S4**. Table of baseline transition intensities and hazard ratios as estimated from a full multi-state model ran on forest states defined with a different threshold than in the main manuscript. Here, plots are assigned to Boreal or Temperate states using a threshold of >85% (instead of >75%) of species dominance of the plot basal area. See Table S5 for details about the interpretation of the table.

\begin{table}[H]
\setstretch{1.4}
\centering\begingroup\fontsize{6.8}{8.8}\selectfont

\begin{tabular}{>{\bfseries}l>{\centering\arraybackslash}p{1.9cm}>{\centering\arraybackslash}p{1.8cm}>{\centering\arraybackslash}p{1.8cm}>{\centering\arraybackslash}p{1.8cm}>{\centering\arraybackslash}p{1.8cm}>{\centering\arraybackslash}p{1.8cm}>{\centering\arraybackslash}p{1.8cm}>{\centering\arraybackslash}p{1.8cm}>{\centering\arraybackslash}p{1.8cm}}
\toprule
\textbf{Transitions} & \textbf{Baseline} & \textbf{Temperature} & \textbf{CMI} & \textbf{Drainage} & \textbf{pH} & \textbf{Natural1} & \textbf{Natural2} & \textbf{Logging1} & \textbf{Logging2}\\
\midrule
Boreal - Boreal & -0.009\newline (-0.01, -0.008) & \cellcolor{white}{ } & \cellcolor{white}{\cellcolor{white}{ }} & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ }\\
Boreal - Mixed & 0.002\newline (0.001, 0.003) & \cellcolor[HTML]{cccccc}{5.29\newline (3.58, 7.83)} & \cellcolor{white}{\cellcolor{white}{1.01\newline (0.73, 1.39)}} & \cellcolor[HTML]{cccccc}{0.72\newline (0.62, 0.83)} & \cellcolor{white}{1.11\newline (0.94, 1.32)} & \cellcolor[HTML]{cccccc}{2.69\newline (1.92, 3.77)} & \cellcolor{white}{0.7\newline (0.05, 9.91)} & \cellcolor{white}{1.89\newline (0.87, 4.13)} & \cellcolor{white}{0.39\newline (0.01, 27.18)}\\
Boreal - Pioneer & 0.007\newline (0.006, 0.008) & \cellcolor{white}{1.000} & \cellcolor{white}{\cellcolor{white}{1.000}} & \cellcolor{white}{1.000} & \cellcolor{white}{1.000} & \cellcolor[HTML]{cccccc}{4.61\newline (3.73, 5.68)} & \cellcolor[HTML]{cccccc}{37.5\newline (31, 45.38)} & \cellcolor[HTML]{cccccc}{8.81\newline (6.08, 12.77)} & \cellcolor[HTML]{cccccc}{197.37\newline (110.34, 353.04)}\\
Mixed - Boreal & 0.002\newline (0.001, 0.004) & \cellcolor{white}{0.8\newline (0.41, 1.56)} & \cellcolor[HTML]{cccccc}{\cellcolor[HTML]{cccccc}{1.61\newline (1.1, 2.35)}} & \cellcolor{white}{1.08\newline (0.79, 1.48)} & \cellcolor{white}{0.69\newline (0.45, 1.07)} & \cellcolor{white}{0.88\newline (0.33, 2.33)} & \cellcolor{white}{1.19\newline (0.04, 38.52)} & \cellcolor{white}{1.05\newline (0.41, 2.68)} & \cellcolor{white}{0.76\newline (0, 201.91)}\\
Mixed - Mixed & -0.021\newline (-0.024, -0.019) & \cellcolor{white}{ } & \cellcolor{white}{\cellcolor{white}{ }} & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ }\\
Mixed - Pioneer & 0.005\newline (0.004, 0.006) & \cellcolor{white}{1.000} & \cellcolor{white}{\cellcolor{white}{1.000}} & \cellcolor{white}{1.000} & \cellcolor{white}{1.000} & \cellcolor[HTML]{cccccc}{3.2\newline (2.04, 5.04)} & \cellcolor[HTML]{cccccc}{12.6\newline (6.22, 25.5)} & \cellcolor[HTML]{cccccc}{2.7\newline (1.61, 4.54)} & \cellcolor[HTML]{cccccc}{46.45\newline (33.21, 64.97)}\\
Mixed - Temperate & 0.015\newline (0.013, 0.017) & \cellcolor{white}{0.89\newline (0.71, 1.11)} & \cellcolor[HTML]{cccccc}{\cellcolor[HTML]{cccccc}{0.76\newline (0.66, 0.88)}} & \cellcolor{white}{1.02\newline (0.94, 1.12)} & \cellcolor{white}{0.95\newline (0.88, 1.02)} & \cellcolor[HTML]{cccccc}{2.38\newline (1.92, 2.96)} & \cellcolor[HTML]{cccccc}{3.07\newline (1.6, 5.88)} & \cellcolor[HTML]{cccccc}{3.14\newline (2.61, 3.78)} & \cellcolor[HTML]{cccccc}{4.71\newline (3.02, 7.34)}\\
Pioneer - Boreal & 0.025\newline (0.023, 0.027) & \cellcolor[HTML]{cccccc}{1.21\newline (1.1, 1.33)} & \cellcolor[HTML]{cccccc}{\cellcolor[HTML]{cccccc}{1.8\newline (1.64, 1.96)}} & \cellcolor{white}{0.98\newline (0.94, 1.03)} & \cellcolor{white}{0.95\newline (0.9, 1.01)} & \cellcolor{white}{0.96\newline (0.78, 1.17)} & \cellcolor[HTML]{cccccc}{0.37\newline (0.28, 0.49)} & \cellcolor[HTML]{cccccc}{1.85\newline (1.5, 2.28)} & \cellcolor{white}{1.7\newline (0.94, 3.08)}\\
Pioneer - Mixed & 0.005\newline (0.005, 0.006) & \cellcolor[HTML]{cccccc}{5.31\newline (4.29, 6.58)} & \cellcolor[HTML]{cccccc}{\cellcolor[HTML]{cccccc}{1.82\newline (1.55, 2.14)}} & \cellcolor[HTML]{cccccc}{0.89\newline (0.81, 0.98)} & \cellcolor{white}{0.99\newline (0.91, 1.08)} & \cellcolor{white}{1.13\newline (0.8, 1.59)} & \cellcolor[HTML]{cccccc}{0.48\newline (0.25, 0.9)} & \cellcolor[HTML]{cccccc}{2.42\newline (1.82, 3.21)} & \cellcolor{white}{0.8\newline (0.57, 1.12)}\\
Pioneer - Pioneer & -0.031\newline (-0.033, -0.029) & \cellcolor{white}{ } & \cellcolor{white}{\cellcolor{white}{ }} & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ }\\
Pioneer - Temperate & 0\newline (0, 0.001) & \cellcolor[HTML]{cccccc}{22.53\newline (14.08, 36.06)} & \cellcolor[HTML]{cccccc}{\cellcolor[HTML]{cccccc}{2.48\newline (1.82, 3.38)}} & \cellcolor[HTML]{cccccc}{0.7\newline (0.57, 0.85)} & \cellcolor{white}{1.04\newline (0.91, 1.2)} & \cellcolor{white}{0.1\newline (0, 5.96)} & \cellcolor{white}{0.57\newline (0.14, 2.31)} & \cellcolor{white}{1.05\newline (0.46, 2.38)} & \cellcolor{white}{0.34\newline (0.09, 1.28)}\\
Temperate - Mixed & 0.017\newline (0.014, 0.021) & \cellcolor[HTML]{cccccc}{0.55\newline (0.42, 0.71)} & \cellcolor{white}{\cellcolor{white}{0.88\newline (0.76, 1.02)}} & \cellcolor[HTML]{cccccc}{1.35\newline (1.2, 1.52)} & \cellcolor[HTML]{cccccc}{0.79\newline (0.71, 0.87)} & \cellcolor[HTML]{cccccc}{1.78\newline (1.16, 2.74)} & \cellcolor{white}{1.61\newline (0.53, 4.91)} & \cellcolor[HTML]{cccccc}{1.32\newline (1.04, 1.67)} & \cellcolor[HTML]{cccccc}{4.35\newline (2.38, 7.93)}\\
Temperate - Pioneer & 0.001\newline (0.001, 0.002) & \cellcolor{white}{1.000} & \cellcolor{white}{\cellcolor{white}{1.000}} & \cellcolor{white}{1.000} & \cellcolor{white}{1.000} & \cellcolor{white}{0.54\newline (0.01, 46.28)} & \cellcolor[HTML]{cccccc}{37.81\newline (13.25, 107.92)} & \cellcolor[HTML]{cccccc}{3.41\newline (1.76, 6.62)} & \cellcolor[HTML]{cccccc}{142.39\newline (85.95, 235.89)}\\
Temperate - Temperate & -0.018\newline (-0.022, -0.015) & \cellcolor{white}{ } & \cellcolor{white}{\cellcolor{white}{ }} & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ } & \cellcolor{white}{ }\\
\bottomrule
\end{tabular}
\endgroup{}
\end{table}

\elandscape

\pagebreak


**Table S5**. List of R packages used.

|**Packages** |**Main functions** | **Uses**                     | **References**           |
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
| labdsv   | indval         | Calculates the indicator value       | @roberts_labdsv_2019 |


\pagebreak

# Supplementary figures

![Waffle charts representing the frequency of forest plots by disturbance type (natural disturbances and logging), level of severity (minor, moderate, major) and vegetation zone (boreal, ecotone and temperate). One square is one occurrence of a disturbance in a forest plot (a forest plot can be disturbed more than once). In each chart (except for the no or minor disturbances), the colours represent the 21 original disturbance types recorded in the field surveys. ](res/figS1_waffle.png)

\pagebreak

![Histograms of time intervals between surveys by bioclimatic domains.](res/figS2_hist_intervals.pdf)

\pagebreak

![Temporal trends in growing season temperatures (top) and annual climate moisture index (bottom). Grey lines represent averaged climate values across the 10,388 studied forest plots. Straight black lines show the fitted least-squared linear regression lines.](res/figS3_clim_trend.pdf)


\pagebreak


![Predicted proportion of temperate species in forest plots (temperate / (boreal + temperate)) measured in basal area using a negative binomial GLM with polynomial terms of two climate variables, temperature and climate moisture index during the growing season (R^2^ = 50%).](res/figS4_predBT.pdf)


\pagebreak



![Performance comparisons of the five candidate multi-state models using multi-class area under the curve (mAUC; a) and logarithmic skill score (LS; b) obtained through 10-fold cross-validation. Higher values of mean pairwise AUC (a) indicate a better capacity to discriminate between pairs of the four forest states: (B)oreal, (M)ixed, (P)ioneer and (T)emperate. The overall mAUC of each model is given next to the legend. Lower values of mean state-specific logarithmic scores (b) indicate better prediction accuracy for each of the four forest states. The overall logarithmic score of each model is given next to the legend.](res/figS5_cv.pdf)

Model evaluation using 10-fold cross-validation revealed that including climate
and disturbances improved overall model predictive performances, while soil
variables had a negligible effect. All models were good at
distinguishing Boreal from Temperate (high pairwise AUC). Soil variables
slightly help to predict Mixed and Temperate states. Including climate variables
help to distinguish Mixed from the other states, while including disturbances
help to distinguish Pioneer from the other states, especially Boreal.


\pagebreak

![Probability of transition between forest states through time for natural disturbances (left) and logging (right) and for different levels (line type) as predicted by the best multi-state model. All other covariates are fixed at the average conditions found in the ecotone, i.e. the balsam fir-yellow birch domain.](res/figS6_proba.pdf)


\pagebreak


![Predicted change in 10-year transition probabilities for different disturbance types and levels. All other covariates are fixed at the average conditions found in the ecotone. Letters correspond to the four forest states: (B)oreal, (M)ixed, (P)ioneer and (T)emperate. Numbers are the modelled transition probabilities from rows to columns and darker colour indicates higher transition probability.](res/figS7_pmatrix.pdf)

The largest values across most matrices were generally associated with
self-transitions (matrix diagonal), meaning that the vast majority of forest
plots remained in the same state after 10 years. At minor disturbances, the
self-transitions were very strong but transitions from Pioneer to Boreal and
from Mixed to Temperate were also important. At moderate disturbances,
probabilities of self-transitions decreased, while transitions from Boreal to
Pioneer, and from Mixed to Temperate increased the most. Transitions from Mixed
and Temperate to Pioneer did not increase much at moderate disturbances, likely
because such disturbances were less frequent and less severe than in Boreal
forests. The difference between natural disturbances and logging was conspicuous
only for major disturbances. For both types of major disturbances, the
transition probabilities to Pioneer showed a great increase compared to moderate
disturbances, but these values exploded in the major logging transition matrix,
exceeding self-transitions.

\pagebreak

![State contribution to forest turnover (see Fig. 7a,b in main text) along the temperature (latitudinal) gradient for different disturbance scenarios: minor (solid), moderate (dashed) and major (dotted) disturbances for both natural (a,c,e) and logging (b,d,f). All other covariates are fixed at the average conditions found in the ecotone, i.e. the balsam fir-yellow birch domain, to focus solely on the effect of disturbances along the temperature gradient. The turnover time of a state (or sojourn time) measures the time spent in this state before transitioning to the next. Long turnover time can translate to large resistance. Here, at any point along the gradient, state turnover time is scaled by the steady state distribution and the sum of all scaled state turnover gives the the turnover time of the transition matrix. The colors at the top of the plots approximate the position of the bioclimatic domains.](res/figS8_contrib2turnover.pdf)

\pagebreak

![State contribution to forest entropy (see Fig. 7c,d in main text) along the temperature (latitudinal) gradient for different disturbance scenarios: minor (solid), moderate (dashed) and major (dotted) disturbances for both natural (a,c,e) and logging (b,d,f). All other covariates are fixed at the average conditions found in the ecotone, i.e. the balsam fir-yellow birch domain, to focus solely on the effect of disturbances along the temperature gradient. The entropy of a state measures the incertitude of its next transition. Here, at any point along the gradient, state entropy is scaled by the steady state distribution and the sum of all scaled state entropy gives the the entropy of the transition matrix. The colors at the top of the plots approximate the position of the bioclimatic domains.](res/figS9_contrib2entropy.pdf)


\pagebreak




# References
