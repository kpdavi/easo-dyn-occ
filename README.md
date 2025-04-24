# Code for Eastern Screech-Owl dynamic occupancy analyses for Fort Collins, CO

This repository includes presence / absence survey data from 2013–2021 from Fort Collins, Colorado for Eastern screech-owl (*Megascops asio maxwelliae*) and vegetation and climate variables associated with survey locations, and scripts for analyzing and plotting these data in R. We fit a dynamic occupancy model in the Bayesian framework to explore environmental factors associated with the species' breeding season occupancy dynamics in an understudied portion of its range and to inform local land management efforts. We found that Eastern Screech-Owl persistence probability decreased while colonization probability increased with aggregation of tree cover and average breeding season temperature. We also found colonization probability increased with cumulative winter precipitation, and persistence probability decreased with average winter minimum temperature. Overall, our results suggest that this population of Eastern Screech-Owls responds to climate before and during the breeding season in complex ways, and that this species appears to experience high turnover in this study area along the western edge of its continental range. The code in this repository reproduces the statistics and figures of the corresponding paper.

## Workflow notes

**JAGS must be installed in order to fit the dynamic occupancy model.** Download JAGS at https://sourceforge.net/projects/mcmc-jags/.

**Open the easo-dyn-occ.Rproj file first to ensure that all folders within the project are properly loaded into the working directory.**


## Description of the data and folder structure

### Data
All data used in the analysis are contained in the folders described below. Please see the publication for more details about sources of raw data and how variables were calculated.

* **data** &#x1F4C1;. This folder contains the raw (**raw** &#x1F4C1;) and processed (**processed** &#x1F4C1;) data for the Eastern screech-owl surveys and environmental covariates.
  + **raw** &#x1F4C1;
    - `RAP_tree_cover_2013-2021.tif`. An image of tree cover rasters from the [Rangeland Analysis Platform](https://rangelands.app/rap/?biomass_t=herbaceous&ll=39.0000,-98.0000&z=5) from 2013–2021. This image is comprised of 9 remotely-sensed layers of continuous tree cover (one layer for each year) at a 30-m<sup>2</sup> spatial resolution (Allred et al. 2021).
    - `sites.csv`. The site IDs (*Transect*) and randomly-generated coordinates (*longitude* and *latitude*). The locations of the monitoring sites are not public, so coordinates were randomly generated so that the script would run. However, because the coordinates are not the actual coordinates for the sites, the values in the output file created in this script will not align with the data values used in the manuscript. Data with the real covariate values for tree cover are in `covs_tree.csv` in  **processed** &#x1F4C1;.
    - `surveys.csv`. The Eastern screech-owl survey data. Columns include site ID (*Transect*), year (*year*; 2013–2021), survey date (*date*), survey start time (*start.time*), presence / absence of the species (*presence*; 0 [absence] or 1 [presence]).
  + **processed** &#x1F4C1;
    - `coeff_names`. A dataframe that includes the column name of a covariate (*names_cov*), the real name of the covariate (*names_real*), the associated parameter in the model for the covariate (*param*), and the process the covariate is fit on (*proc*) for the initial occupancy, persistence, and colonization subprocesses of the dynamic occupancy model.
    - `covs_climate.csv`. The climate covariate data. Columns include year (*year*; 2013–2021), average breeding season minimum temperature (&deg;C) in the current year (*avg.breeding.tmin*) and the breeding season prior (*lag.avg.breeding.tmin*), average winter minimum temperature (&deg;C) preceding the breeding season (*avg.winter.tmin*) and the winter a year prior (*lag.avg.winter.tmin*), and cumulative winter precipitation (mm) preceding the breeding season (*sum.winter.prcp*).
    - `covs_lidar.csv`. The data for standard deviation in vegetation height in 2013. Columns include site ID (*Transect*), year (*year*; 2013), and standard deviation in vegetation height within 250-m-radii buffers surrounding the centroid of each site (meters; *sd.vegheight.m.2013*).
    - `covs_tree.csv`. The data for proportion and cohesion of tree cover. Columns include site ID (*Transect*), year (*year*; 2013–2021), proportion of tree cover (was multiplied by 100 so values are percentages; *tree.prop*), and cohesion of tree cover (%; *tree.coh*).
    - `model_data.RData`. Dataframes, arrays, and vectors needed to fit and / or plot results from the dynamic occupancy model.
      - `easo_obs`. A dataframe that includes the real (*Transect*) and a numeric form (*site_id*) values for site ID, the real (*year*) and numeric form (*year_id*) of year, presence / absence of the species (*presence*; 0 [absence] or 1 [presence]), an intercept column for the detection process (*intercept*; 1), and standardized values for moon phase (*moon.phase*) and ordinal date (*ord.date*).
      - `X_det`. A matrix of presence / absence of the species by survey (rows) and scaled detection covariates (columns). Column order -- intercept[1], ordinal date[2], and moon phase[3].
      - `cov_det_num`. The number of covariates included in the observation process (detection) in the dynamic occupancy model.
      - `cov_gamma_num`. The number of covariates included in the colonization subprocess in the dynamic occupancy model.
      - `cov_phi_num`. The number of covariates included in the persistence subprocess in the dynamic occupancy model.
      - `cov_psi1_num`. The number of covariates included in the initial occupancy subprocess in the dynamic occupancy model.
      - `position_lidar`. The number of the layer for standard deviation in vegetation height in `X_psi1`. Used for estimating hyperpriors in the model to account for missing data in that covariate.
      - `X_[gam, phi, psi1]`. Array of scaled covariates for the colonization, persisitence, or initial occupancy subprocesses in the dynamic occupancy model. Rows are sites, columns are years, and the layers are covariates.
    - `surveys_covs_dat.csv`. The Eastern screech-owl survey data and data for detection covariates. Columns include site ID (*Transect*), year (*year*; 2013–2021), survey date (*date*), survey start time (*start.time*), presence / absence of the species (*presence*; 0 [absence] or 1 [presence]), ordinal date (*date.ord*), and moon phase (%; *moon.phase*).


### Scripts
Within the subfolders of **scripts** &#x1F4C1;, scripts with the same numeric prefix can be run in any order.

* **data_prep** &#x1F4C1;
  + `1_data_prep_*.R` &#x1F4C4;. Code to calculate the climate, vegetation, and detection covariates. The output of these scripts are in **data/processed** &#x1F4C1; and include `cov_*.csv` for the climate and vegetation covariates and `surveys_covs_det.csv` for the detection covariates.
  + `2_model_prep.R` &#x1F4C4;. Code to format the climate and vegetation covariates and generate variables needed for model fitting. The output of the script is `model_data.RData` in **data/processed** &#x1F4C1;.
* **model_fitting** &#x1F4C1;
  + `3_model_fit_dynocc.R` &#x1F4C4;. Code to fit the dynamic occupancy model in the Bayesian framework. Outputs of the script are in **output** &#x1F4C1; and include parameter estimates for coefficients / scalars (`mod_out_coeffs.rds`), occupancy process parameters (`mod_out_occ.rds`), or parameters for calculating area under the receiver operating characteristic curve (AUC; `mod_out_auc.rds`).
  + `4_summarize_model_output_diagnostics.R` &#x1F4C4;. Code to summarize results — posterior means, 95% credible intervals, and probability of direction [*pd*] — for focal coefficient parameters, and model diagnostics — Bayesian *p* value for the log likelihood and area under the receiver operating characteristic curve [AUC] — from the dynamic occupancy model.
  + `model_code_dynocc.txt` &#x1F4C4;. Text file of the dynamic occupancy model in JAGS. This file is called within the function to fit the model in the `3_model_fit_dynocc.R` script.
* **plots** &#x1F4C1;
  + `4_fig*.R` &#x1F4C4;. Code to create figures 2--4 in the publication. Outputs of these scripts are in **plots** &#x1F4C1; and include `fig*.png` (described in more detail below).
* **utils** &#x1F4C1;
  + `utils_*.R` &#x1F4C4;. Helper functions for calculating climate variables (`utils_calculate_daymet_climate_variables.R`), predicting data from the dynamic occupancy model (`utils_predict_data_from_model.R`), and generating figures (`utils_plot_occ_prob_by_year.R` and `utils_plot_predicted_data.R`).

### Outputs
All model outputs created in the analyses are contained in the following folder:

* **output** &#x1F4C1;. Results from the dynamic occupancy model (`mod_*.rds`). Specifically, parameter estimates for coefficients / scalars (`mod_out_coeffs.rds`), occupancy process parameters (`mod_out_occ.rds`), and parameters for calculating area under the receiver operating characteristic curve (AUC; `mod_out_auc.rds`).

### Plots
Figures that were created from data in the manuscript are contained in the following folder:

* **plots** &#x1F4C1;. Plots for Figures 2--4 in the manuscript, which visualize parameter estimates for model coefficients (`fig2_coeff_param_ests.png`), model-predicted probability of persistence, colonization, and initial occupancy as a function of covariates (`fig3_coeff_param_preds.png`), and model-predicted probability of persistence, colonization, occupancy, and detection by year (`fig4_occ_proc_prob.png`).

## Sharing/Access information

Research permits to conduct surveys for Eastern screech-owl were secured with the City of Fort Collins Natural Areas Department.

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
- Hesselbarth, M.H.K., M. Sciaini, K.A. With, K. Wiegand, and J. Nowosad. 2019.     landscapemetrics: an open-source R tool to calculate landscape metrics. *Ecography* 42: 1648-1657 (v2.2).
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
  <http://www.biomedcentral.com/1471-2105/12/77/>
- Thieurmel, B., and A. Elmarhraoui. 2022. _suncalc: Compute Sun Position, Sunlight
  Phases, Moon Position and Lunar Phase_. R package version 0.5.1,
  <https://CRAN.R-project.org/package=suncalc>.
- Ushey, K., and H. Wickham. 2025. *renv: Project Environments*. R package version
  1.1.2, <https://CRAN.R-project.org/package=renv>.
- Wickham, H. et al. 2019. Welcome to the tidyverse. *Journal of Open Source Software* 4(43): 1686. <https://doi.org/10.21105/joss.01686>.
- Wilke, C. 2024. _cowplot: Streamlined Plot Theme and Plot Annotations for
  'ggplot2'_. R package version 1.1.3,
  <https://CRAN.R-project.org/package=cowplot>.
