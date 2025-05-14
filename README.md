# Scripts and figures for Eastern screech-owl dynamic occupancy analysis in Fort Collins, CO

This repository contains the scripts to reproduce analyses and figures for the corresponding scientific publication in the [Wilson Journal of Ornithology (2025)](https://doi.org/10.1080/15594491.2025.2497148). The complete repository for this publication, including the raw / processed data and scripts to compile the data, fit the model, and summarize / plot model outputs, is on Zenodo (DOI: [10.5281/zenodo.15272794]( https://doi.org/10.5281/zenodo.15272794)).

We fit a dynamic occupancy model in the Bayesian framework to 9 years of occupancy data (2013–2021) from a participatory science monitoring program in Fort Collins, Colorado, USA to explore environmental factors associated with Eastern screech-owl (*Megascops asio maxwelliae*) breeding season occupancy dynamics in an understudied portion of its range. We found that Eastern screech-owl persistence probability decreased while colonization probability increased with aggregation of tree cover and average breeding season temperature. We also found colonization probability increased with cumulative winter precipitation, and persistence probability decreased with average winter minimum temperature. Overall, our results suggest that this population of Eastern screech-owls responds to climate before and during the breeding season in complex ways, and that the species appears to experience high turnover in this study area along the western edge of its continental range.

## Sharing/access information

The Eastern screech-owl monitoring program is ongoing and **potential users of the data are encouraged to contact [Rob Sparks](mailto:Rob.Sparks@birdconservancy.org?subject=inquiry%20about%20Eastern%20screech-owl%20project) for additional information or potential for collaboration.** Research permits to conduct surveys for Eastern screech-owl were secured with the City of Fort Collins Natural Areas Department.

## Workflow notes

**JAGS must be installed in order to fit the dynamic occupancy model.** Download JAGS at https://sourceforge.net/projects/mcmc-jags/.

**Open the easo-dyn-occ.Rproj file first to ensure that all folders within the project are properly loaded into the working directory.**

Relatedly, and in addition, **this project uses the [`renv` package](https://docs.posit.co/ide/user/ide/guide/environments/r/renv.html) to keep track of which R packages and versions were used** so that others can easily set up the same environment using `renv.lock` and `renv/`. If you open the `.Rproj` file first, the `renv` package will ask to install itself, and if that prompt is accepted, it will restore the required packages — i.e., you won't need to manually interact with `renv.lock` or `renv/` yourself.

## Description of the folder structure for scripts / plots

Descriptions of the raw and processed data and additional details about script outputs (e.g., model results) are in the associated Zenodo repository (DOI: [10.5281/zenodo.15272794]( https://doi.org/10.5281/zenodo.15272794)).

### Scripts
Within the subfolders of **scripts** &#x1F4C1;, scripts with the same numeric prefix can be run in any order.

* **data_prep** &#x1F4C1;
  + `1_data_prep_*.R` &#x1F4C4;. Code to calculate the climate, vegetation, and detection covariates.
  + `2_model_prep.R` &#x1F4C4;. Code to format the climate and vegetation covariates and generate variables needed for model fitting.
* **model_fitting** &#x1F4C1;
  + `3_model_fit_dynocc.R` &#x1F4C4;. Code to fit the dynamic occupancy model in the Bayesian framework.
  + `4_summarize_model_diagnostics.R` &#x1F4C4;. Code to summarize model diagnostics — Bayesian *p* value for the log likelihood and area under the receiver operating characteristic curve (AUC) — from the dynamic occupancy model.
  + `model_code_dynocc.txt` &#x1F4C4;. Text file of the dynamic occupancy model in JAGS. This file is called within the function to fit the model in the `3_model_fit_dynocc.R` script.
* **plots** &#x1F4C1;
  + `4_fig*.R` &#x1F4C4;. Code to create figures 2—4 in the publication. Outputs of these scripts are in **plots** &#x1F4C1; and include `fig*.png` (described in more detail below).
* **utils** &#x1F4C1;
  + `utils_*.R` &#x1F4C4;. Helper functions for calculating climate variables (`utils_calculate_daymet_climate_variables.R`), predicting data from the dynamic occupancy model (`utils_predict_data_from_model.R`), and generating figures (`utils_plot_occ_prob_by_year.R` and `utils_plot_predicted_data.R`).

### Plots
Figures from the manuscript are contained in the following folder:

* **plots** &#x1F4C1;. Plots for Figures 2—4 in the manuscript, which visualize parameter estimates for model coefficients (`fig2_coeff_param_ests.png`), model-predicted probability of persistence, colonization, and initial occupancy as a function of covariates (`fig3_coeff_param_preds.png`), and model-estimated probability of persistence, colonization, occupancy, and detection by year (`fig4_occ_proc_prob.png`).


## Code/Software

All scripts were run using R Core Team (2024). Loaded packages / software are indicated in the scripts and cited below:

- Bocinsky, R.K. 2024. _FedData: Download Geospatial Data Available from Several
  Federated Data Sources_. R package version 4.2.0,
  <https://CRAN.R-project.org/package=FedData>.
- Denwood, M.J. 2016. runjags: An R package providing interface
  utilities, model templates, parallel computing methods and additional
  distributions for MCMC models in JAGS. *Journal of Statistical Software* 71(9): 1-25. <https://doi.org/10.18637/jss.v071.i09>.
- Grolemund, G., and H. Wickham. 2011. Dates and times made easy with
  lubridate. *Journal of Statistical Software* 40(3): 1-25. <https://www.jstatsoft.org/v40/i03/>.
- Hesselbarth, M.H.K., M. Sciaini, K.A. With, K. Wiegand, and J. Nowosad. 2019. landscapemetrics: an open-source R tool to calculate landscape metrics. *Ecography* 42: 1648-1657 (ver. 2.2). <https://doi.org/10.1111/ecog.04617>.
- Hijmans, R. 2024. _terra: Spatial Data Analysis_. R package version 1.7-78,
  <https://CRAN.R-project.org/package=terra>.
- Lüdecke, D. 2018. ggeffects: Tidy data frames of marginal effects from regression
  models. *Journal of Open Source Software* 3(26): 772. <https://doi.org/10.21105/joss.00772>.
- Müller, K. 2020. _here: A Simpler Way to Find Your Files_. R package version 1.0.1,
  <https://CRAN.R-project.org/package=here>.
- Pebesma, E., & R. Bivand. 2023. *Spatial Data Science: With Applications in
  R*. Chapman and Hall/CRC. <https://doi.org/10.1201/9780429459016>
- Pebesma, E., 2018. Simple features for R: Standardized support for spatial
  vector data. *The R Journal* 10(1): 439-446. <https://doi.org/10.32614/RJ-2018-009>.
- Pierce, D. 2024. _ncdf4: Interface to Unidata netCDF (Version 4 or Earlier)
  Format Data Files_. R package version 1.23,
  <https://CRAN.R-project.org/package=ncdf4>.
- Plummer, M. 2003. JAGS: A Program for Analysis of Bayesian Graphical Models Using Gibbs Sampling. *Proceedings of the 3rd International Workshop on Distributed Statistical Computing* (DSC 2003), Vienna, 20-22 March 2003, 1-10.
- R Core Team. 2024. R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>.
- Robin X., N. Turck, A. Hainard, N. Tiberti, F. Lisacek, J.-C. Sanchez, and M. Müller. 2011. pROC: an open-source package for R and S+ to analyze and compare ROC curves. *BMC Bioinformatics* 12: 77. <https://doi.org/10.1186/1471-2105-12-77>.
- Thieurmel, B., and A. Elmarhraoui. 2022. _suncalc: Compute Sun Position, Sunlight
  Phases, Moon Position and Lunar Phase_. R package version 0.5.1,
  <https://CRAN.R-project.org/package=suncalc>.
- Ushey, K., and H. Wickham. 2025. *renv: Project Environments*. R package version
  1.1.2, <https://CRAN.R-project.org/package=renv>.
- Wickham, H. et al. 2019. Welcome to the tidyverse. *Journal of Open Source Software* 4(43): 1686. <https://doi.org/10.21105/joss.01686>.
- Wilke, C. 2024. _cowplot: Streamlined Plot Theme and Plot Annotations for
  'ggplot2'_. R package version 1.1.3,
  <https://CRAN.R-project.org/package=cowplot>.
